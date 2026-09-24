"""
sigFitCond.py -- conditional signal fit: a double-sided Crystal Ball whose
width is set event by event from the per-event relative mass resolution.

Model
-----
    P(m, r) = P(m | r) * P(r)              r = sigma_m / m  (relative)

    P(m | r) = DSCB(m ; mean, sigma_CB(r), alphaL, nL, alphaR, nR)
    sigma_CB(r) = kappa * r * MH

    P(r)     = RooKeysPdf, from the MC distribution of r

CONDITIONAL-ONLY (Option A of PLAN_conditional_fit.md): the fit declares
`r` a conditional observable, so P(r) cancels from the likelihood and
never enters it. Consequences:

  * mean, alphaL, nL, alphaR, nR and kappa are determined by the
    conditional shape alone;
  * no background resolution model is needed, which is the whole reason
    to start here rather than with the full 2D fit;
  * P(r) is still built, because it is needed to PROJECT the model onto
    the mass axis for plotting, and later for toys.

Double-sided, not the single-tail RooCBShape: the 1D fits give
alpha_lo ~ 1.35 and alpha_hi ~ 1.65, so both tails carry weight and a
one-sided CB would fit worse than the model this is meant to improve on.

The number to look at
---------------------
    kappa ~ 1  means the per-event error is correctly calibrated in
    absolute terms. From the 1D fits (sigma = 1.495 GeV for ggHcat 2024)
    and the median r (~0.012), kappa = 1.495/(0.012*125) = 0.997 is the
    expectation. A kappa far from 1 says the propagated error is
    mis-scaled, and the conditional model inherits that.

The second number is chi2/ndf of the projected model, against 4.73 for
the 1D double-sided CB on the same sample. If the coherent wave through
121-135 GeV does not flatten, the single-width assumption was not causing
it and this does not fix it.

Usage
-----
    python sigFitCond.py -c ggHcat -y 2024 -b incl -s ggH
    python sigFitCond.py -c ggHcat -y 2024 -b incl -s ggH --max-events 200000
    python sigFitCond.py -c VBFcat -y Run3 -b incl -s qqH --pdf
"""

import argparse
import json
import os
from array import array
from datetime import datetime, timezone

import numpy as np
import ROOT

import sigFit as sf


ROOT.gROOT.SetBatch(True)

# ---------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------

MASSCOL = "HiggsCandCorrMass"
RELERRCOL = "HiggsCandCorrMassErr"      # relative, and paired with the above
WEIGHTCOL = "w_allSF"

# Range of r kept in the fit.
#
# Set WIDE by default, so the conditional fit uses the same events as the
# 1D fit and the two are comparable. The badly measured tail is not junk --
# it is exactly the population the conditional model handles better than a
# single width -- so cutting it would both bias the comparison and discard
# the method's advantage.
#
# The only hard requirement is r > 0, or sigma_CB = 0. The upper edge
# exists because a RooRealVar needs a range and RooKeysPdf is built over
# it; a very long sparse tail degrades the kernel bandwidth in the
# populated region. 0.10 is ~8x the median (~0.012), well clear of the
# bulk. Tighten it only for the kernel estimate, never for the fit, and
# check with --report-cut what any change actually costs.
R_LO, R_HI = 1e-4, 0.10

# narrower window used only for DRAWING the resolution distribution
R_PLOT_HI = 0.05
R_NBINS = 100

# Unbinned conditional fits are O(N) per likelihood call, and RooKeysPdf is
# O(N) per EVALUATION. Both need a cap.
# The binned 1D reference runs on the SAME subsample, so the comparison is
# fair at any size and neither fit needs full statistics. 200k keeps the
# unbinned fit to minutes while leaving plenty of power to see the
# percent-level shape structure the comparison is about.
MAX_EVENTS = 200000        # events loaded into the RooDataSet
KEYS_MAX = 30000           # events used to build the kernel estimate
KEYS_RHO = 1.0             # RooKeysPdf bandwidth scale

# kappa: the scale between the propagated per-event error and the fitted
# width. Expected near 1; the range is deliberately wide enough to show a
# gross mis-scaling rather than hide it on a bound.
KAPPA_START, KAPPA_LO, KAPPA_HI = 1.0, 0.3, 3.0

# Optional intrinsic floor added in quadrature:
#     sigma_CB = sqrt((kappa * r * MH)^2 + c^2)
# The plain model forces the width proportional to r through the origin.
# If the true width has a constant contribution that MinvRelErr does not
# propagate -- angular resolution, FSR, the beamspot correlation -- that
# rigidity makes the conditional fit WORSE than a free single width, even
# though it has the same number of parameters. A significant c says so.
OFFSET_START, OFFSET_LO, OFFSET_HI = 0.3, 0.0, 2.0

MODES = sf.MC_LIST

# sigFit sets the ratio panel range inline rather than as a constant, so
# keep a local one (and pick sigFit's up if it ever gains it).
RATIO_Y_RANGE = getattr(sf, "RATIO_Y_RANGE", (0.0, 2.0))

# BDT-bin cuts. This mirrors the selMVA* dicts inside
# prepareFits.getHisto, which are local to that function and so cannot be
# imported -- the same duplication resolutionScan.py carries. Hoisting
# them to module scope in prepareFits and importing from there would
# remove both copies.
SEL_MVA = {
    "VBFcat":  {"bdt2": "discrMVA0>=0.92", "bdt1": "discrMVA0>=0.64 && discrMVA0<0.92", "bdt0": "discrMVA0<0.64"},
    "ggHcat":  {"bdt2": "discrMVA0>=0.76", "bdt1": "discrMVA0>=0.5 && discrMVA0<0.76",  "bdt0": "discrMVA0<0.5"},
    "VLcat":   {"bdt2": "discrMVA0>=0.78", "bdt1": "discrMVA0>=0.32 && discrMVA0<0.78", "bdt0": "discrMVA0<0.32"},
    "VHcat":   {"bdt2": "discrMVA0>=0.94", "bdt1": "discrMVA0>=0.86 && discrMVA0<0.94", "bdt0": "discrMVA0<0.86"},
    "Zinvcat": {"bdt2": "discrMVA0>=0.98", "bdt1": "discrMVA0>=0.78 && discrMVA0<0.98", "bdt0": "discrMVA0<0.78"},
    "TTHcat":  {"bdt2": "discrMVA0>=0.98", "bdt1": "discrMVA0>=0.8 && discrMVA0<0.98",  "bdt0": "discrMVA0<0.8"},
    "TTLcat":  {"bdt1": "discrMVA0>=0.54", "bdt0": "discrMVA0<0.54"},
}


def bin_cut(category, binMVA):
    """RDataFrame cut for a BDT bin, or None for the inclusive bin."""
    if binMVA in ("incl", "", None):
        return None
    try:
        return SEL_MVA[category][binMVA]
    except KeyError:
        raise SystemExit(f"no cut defined for {category} bin {binMVA}; "
                         f"known: {sorted(SEL_MVA.get(category, {}))}")


# ---------------------------------------------------------------------------
# input
# ---------------------------------------------------------------------------

def load_dataset(x, r, w, category, binMVA, year, sig, max_events, seed=0):
    """Unbinned RooDataSet in (m, r), weighted.

    The conditional fit needs the per-event r, so sigFit's binned path does
    not carry over. Values come through numpy rather than a RooFit import
    loop, which is far faster for ~1e6 rows.
    """
    from prepareFits import getHisto           # noqa: F401  (import check)
    import prepareFits as pf

    # reuse prepareFits' file lookup if it exposes one, else glob as
    # resolutionScan does
    if hasattr(pf, "get_signal_files"):
        files = pf.get_signal_files(sig, category, year,
                                    getattr(pf, "ROOTFILES_DIR", ""))
    else:
        import glob
        base = getattr(pf, "ROOTFILES_DIR",
                       "/work/submit/kbai/HmumuRun3/ROOTFILES")
        tags = {"ggH": ["11"], "qqH": ["10"],
                "VH": ["12", "13", "14"], "ttH": ["15"]}[sig]
        files = []
        for t in tags:
            pat = (f"{base}/{category}/snapshot_mc_{t}_*_{category}.root"
                   if year == "Run3"
                   else f"{base}/{category}/snapshot_mc_{t}_{year}_{category}.root")
            files.extend(sorted(glob.glob(pat)))
    if not files:
        raise SystemExit(f"no snapshots for {category} {sig} {year}")
    print(f"{len(files)} file(s):")
    for f in files:
        print(f"    {f}")

    df = ROOT.RDataFrame("events", files)
    cols = set(str(c) for c in df.GetColumnNames())
    for need in (MASSCOL, RELERRCOL):
        if need not in cols:
            raise SystemExit(f"input has no {need}")
    weighted = WEIGHTCOL in cols
    if not weighted:
        print(f"warning: no {WEIGHTCOL}, fitting unweighted")

    cut = bin_cut(category, binMVA)
    # exactly the 1D fit's selection ...
    df = (df.Filter("mc >= 10 && mc <= 15", "signal MC")
            .Filter(f"!std::isnan({MASSCOL})", "valid mass")
            .Filter(f"{MASSCOL} >= {sf.XLOW} && {MASSCOL} < {sf.XHIGH}",
                    "fit window"))
    if cut:
        df = df.Filter(cut, "bdt bin")
    df_mass_only = df
    # ... plus the r window, which only this fit needs
    df = df.Filter(f"{RELERRCOL} > {R_LO} && {RELERRCOL} < {R_HI}",
                   "resolution window")

    take = [MASSCOL, RELERRCOL] + ([WEIGHTCOL] if weighted else [])
    arr = df.AsNumpy(columns=take)
    m_all = np.asarray(arr[MASSCOL], dtype=float)
    r_all = np.asarray(arr[RELERRCOL], dtype=float)
    w_all = (np.asarray(arr[WEIGHTCOL], dtype=float) if weighted
             else np.ones_like(m_all))

    n_tot = len(m_all)
    if n_tot == 0:
        raise SystemExit("no events survive the selection")

    # What the r window costs, relative to the mass window alone. The 1D
    # fit applies no r requirement, so this is the size of the only
    # selection difference between the two.
    n_mass = (df_mass_only.Count().GetValue()
              if df_mass_only is not None else None)
    if n_mass:
        lost = n_mass - n_tot
        print(f"r window ({R_LO:g}, {R_HI:g}) keeps {n_tot} of {n_mass} "
              f"events in the mass window "
              f"({100.0 * lost / n_mass:.3f} % lost)")
        if lost and 100.0 * lost / n_mass > 0.1:
            print("  NOTE: this is the only selection difference from the "
                  "1D fit; widen R_HI to make them identical")

    # Subsample if needed. Random, not the first N: the snapshots are
    # ordered by era, so a head slice would not represent the mixture.
    if n_tot > max_events:
        rng = np.random.default_rng(seed)
        idx = rng.choice(n_tot, size=max_events, replace=False)
        m_all, r_all, w_all = m_all[idx], r_all[idx], w_all[idx]
        # keep the total weight, so `norm` still means the expected yield
        w_all = w_all * (n_tot / float(max_events))
        print(f"subsampled {max_events} of {n_tot} events "
              f"(weights rescaled to preserve the yield)")

    data = build_dataset(x, r, w, m_all, r_all, w_all, category)

    print(f"loaded {data.numEntries()} entries, "
          f"sum of weights {data.sumEntries():.2f}")
    return data, m_all, r_all, w_all, files


def build_dataset(x, r, w, m_all, r_all, w_all, category):
    """RooDataSet from numpy arrays.

    RooDataSet.from_numpy (ROOT >= 6.30) fills in one call; the row-by-row
    loop below is the fallback and is minutes slower at 200k rows -- often
    slower than the fit itself.
    """
    name = f"data_{category}"
    if hasattr(ROOT.RooDataSet, "from_numpy"):
        try:
            return ROOT.RooDataSet.from_numpy(
                {x.GetName(): m_all, r.GetName(): r_all, w.GetName(): w_all},
                ROOT.RooArgSet(x, r, w), name=name, weight_name=w.GetName())
        except Exception as exc:
            print(f"  from_numpy unavailable ({exc}); filling row by row")

    data = ROOT.RooDataSet(name, "data", ROOT.RooArgSet(x, r, w),
                           ROOT.RooFit.WeightVar(w))
    for mi, ri, wi in zip(m_all, r_all, w_all):
        x.setVal(float(mi))
        r.setVal(float(ri))
        data.add(ROOT.RooArgSet(x, r), float(wi))
    return data


def effective_entries(w_values):
    """N_eff = (sum w)^2 / sum w^2.

    How many unweighted events would give the same statistical power. Not
    the same as the raw count when the weights vary, and not the same as
    the yield at all -- for these samples the weights are ~6e-4, so the
    yield is a few hundred while N_eff is in the hundreds of thousands.
    """
    sw = float(np.sum(w_values))
    sw2 = float(np.sum(np.square(w_values)))
    return (sw * sw / sw2) if sw2 > 0 else 0.0


def cache_path(wsdir, tag):
    return os.path.join(wsdir, f"cache_{tag}.npz")


def save_cache(wsdir, tag, m_all, r_all, w_all, files):
    os.makedirs(wsdir, exist_ok=True)
    p = cache_path(wsdir, tag)
    np.savez_compressed(p, m=m_all, r=r_all, w=w_all,
                        files=np.array(files, dtype=object))
    print(f"cached {len(m_all)} events -> {p}")


def load_cache(wsdir, tag):
    p = cache_path(wsdir, tag)
    if not os.path.isfile(p):
        return None
    d = np.load(p, allow_pickle=True)
    print(f"loaded {len(d['m'])} cached events from {p}")
    return d["m"], d["r"], d["w"], list(d["files"])


# ---------------------------------------------------------------------------
# model
# ---------------------------------------------------------------------------

def make_cond_pdf(x, r, tag, scale_unc=None, res_unc=None,
                  correlated=True, offset=False):
    """DSCB with a per-event width.

        mean        = (cb_mu + MH - MH_REF) * (1 + scale_unc * CMS_scale_m)
        sigma_CB(r) = kappa * r * MH * (1 + res_unc * CMS_res_m)

    The width is now a FUNCTION OF THE CONDITIONAL OBSERVABLE, not a fitted
    constant. Note MH and not the observable m: if the width depended on m
    the conditional PDF's normalisation would itself depend on m and the fit
    would be silently wrong.
    """
    scale_unc = sf.SCALE_UNC if scale_unc is None else scale_unc
    res_unc = sf.RES_UNC if res_unc is None else res_unc

    nom = {
        "mu": ROOT.RooRealVar(f"cb_mu_{tag}", "cb_mu",
                              *sf.BOUNDS["mu"]),
        "kappa": ROOT.RooRealVar(f"cb_kappa_{tag}", "kappa",
                                 KAPPA_START, KAPPA_LO, KAPPA_HI),
        "alphaL": ROOT.RooRealVar(f"cb_alphaL_{tag}", "cb_alphaL",
                                  *sf.BOUNDS["alphaL"]),
        "nL": ROOT.RooRealVar(f"cb_nL_{tag}", "cb_nL", *sf.BOUNDS["nL"]),
        "alphaR": ROOT.RooRealVar(f"cb_alphaR_{tag}", "cb_alphaR",
                                  *sf.BOUNDS["alphaR"]),
        "nR": ROOT.RooRealVar(f"cb_nR_{tag}", "cb_nR", *sf.BOUNDS["nR"]),
    }

    mh = ROOT.RooRealVar("MH", "Higgs mass hypothesis", sf.MH_REF,
                         120.0, 130.0)
    mh.setConstant(True)

    nsuf = "" if correlated else f"_{tag}"
    nuis = {
        "scale": ROOT.RooRealVar(f"CMS_scale_m{nsuf}", "muon scale",
                                 0.0, -5.0, 5.0),
        "res": ROOT.RooRealVar(f"CMS_res_m{nsuf}", "muon resolution",
                               0.0, -5.0, 5.0),
    }

    mean = ROOT.RooFormulaVar(
        f"cb_mean_{tag}", "mean",
        f"(@0 + @1 - {sf.MH_REF})*(1 + {scale_unc}*@2)",
        ROOT.RooArgList(nom["mu"], mh, nuis["scale"]))

    if offset:
        nom["c0"] = ROOT.RooRealVar(f"cb_c0_{tag}", "c0",
                                    OFFSET_START, OFFSET_LO, OFFSET_HI)
        # sqrt((kappa*r*MH)^2 + c0^2), then the resolution nuisance
        sigma = ROOT.RooFormulaVar(
            f"cb_sigmaCond_{tag}", "sigmaCond",
            f"sqrt(@0*@0*@1*@1*@2*@2 + @3*@3)*(1 + {res_unc}*@4)",
            ROOT.RooArgList(nom["kappa"], r, mh, nom["c0"], nuis["res"]))
    else:
        # kappa * r * MH * (1 + res_unc * res)
        sigma = ROOT.RooFormulaVar(
            f"cb_sigmaCond_{tag}", "sigmaCond",
            f"@0*@1*@2*(1 + {res_unc}*@3)",
            ROOT.RooArgList(nom["kappa"], r, mh, nuis["res"]))

    pdf = ROOT.RooCrystalBall(
        f"crystal_ball_{tag}", "crystal_ball_conditional",
        x, mean, sigma, nom["alphaL"], nom["nL"], nom["alphaR"], nom["nR"])

    return pdf, {"nom": nom, "nuis": nuis, "mh": mh,
                 "func": {"mean": mean, "sigma": sigma},
                 "cfg": {"scale_unc": scale_unc, "res_unc": res_unc,
                         "correlated": correlated, "mh_ref": sf.MH_REF,
                         "offset": offset,
                         "model": ("conditional DSCB, sigma = "
                                   + ("sqrt((kappa*r*MH)^2 + c0^2)"
                                      if offset else "kappa*r*MH"))}}


def make_keys_pdf(r, r_values, w_values, tag, n_max=KEYS_MAX, rho=KEYS_RHO,
                  seed=0):
    """RooKeysPdf for the relative resolution, from a subsample.

    Two practical points:
      * RooKeysPdf sums a kernel over EVERY event on every evaluation, so
        it is capped at n_max. It does not enter the conditional
        likelihood, but it is evaluated for plotting and toys.
      * MirrorLeft: r is bounded below and piles up near its mode, so
        without mirroring the kernel leaks past the boundary and biases the
        low-r edge downwards.
    """
    n = len(r_values)
    if n > n_max:
        rng = np.random.default_rng(seed)
        # weighted draw, so the subsample follows the weighted distribution
        p = np.clip(w_values, 0, None)
        p = p / p.sum() if p.sum() > 0 else None
        sel = rng.choice(n, size=n_max, replace=False, p=p)
        vals = r_values[sel]
        print(f"kernel estimate from {n_max} of {n} events "
              f"(weighted draw)")
    else:
        vals = r_values

    ds = ROOT.RooDataSet(f"rset_{tag}", "r", ROOT.RooArgSet(r))
    for v in vals:
        r.setVal(float(v))
        ds.add(ROOT.RooArgSet(r))

    keys = ROOT.RooKeysPdf(f"keys_r_{tag}", "keys_r", r, ds,
                           ROOT.RooKeysPdf.MirrorLeft, rho)
    return keys, ds


def mass_histo(m_values, w_values, name, nbins=None):
    """TH1 of the mass from the SAME arrays the unbinned fit uses.

    This is what makes the comparison honest: the binned and unbinned fits
    then see identical events with identical weights, so chi2/ndf is
    comparable without either having to run at full statistics.
    """
    nbins = nbins or int((sf.XHIGH - sf.XLOW) * sf.BINS_PER_GEV)
    h = ROOT.TH1D(name, "", nbins, sf.XLOW, sf.XHIGH)
    h.SetDirectory(0)
    h.Sumw2()
    for mi, wi in zip(m_values, w_values):
        h.Fill(float(mi), float(wi))
    return h


def fit_1d_reference(x1d, m_values, w_values, tag, label, plotdir, year,
                     save_pdf):
    """sigFit's ordinary binned 1D fit, on the same subsample.

    Uses sigFit.fit_one unchanged, so the reference really is the 1D model
    and not a reimplementation of it. A separate observable and a distinct
    tag keep the RooFit names from clashing with the conditional model's.
    """
    h = mass_histo(m_values, w_values, f"h1d_{tag}")
    out = sf.fit_one(x1d, h, f"{tag}_1d", label, plotdir,
                     freeze=False, year=year, save_pdf=save_pdf)
    if out is None:
        return None
    pdf1d, bundle, norm_var, record = out
    xs, ys, es, fs = binned_expectation(h, pdf1d, x1d, norm_var.getVal())
    return {"bins": (xs, ys, es, fs),
            "sigma": record["parameters"]["sigma"]["value"],
            "sigmaErr": record["parameters"]["sigma"]["error"],
            "mu": record["parameters"]["mu"]["value"],
            "muErr": record["parameters"]["mu"]["error"],
            "chi2_ndf": record["fit_quality"]["chi2_ndf"],
            "ok": record["fit_quality"]["ok"],
            "n_eff": record["normalization"]["effective_entries"],
            "record": record, "_keep": (bundle, pdf1d, h)}


def mismatch_metrics(data_y, data_e, fit_y):
    """Quantify the leftover mismatch, three ways.

    chi2_ndf    mean squared pull. Scales with N, so only comparable
                between fits on the SAME events.
    pull_width  sqrt(<p^2>); 1 for a correct model. Same N dependence.
    frac_rms    sqrt(<(d/f - 1)^2> - <e^2/f^2>), the statistics-subtracted
                RMS of the Template/Fit deviation. Reads as "the fit is
                wrong by X% RMS" and does NOT depend on N -- the number to
                quote when comparing models, or the same model at
                different statistics.
    runs_z      runs test on the pull signs. Near 0 means the deviations
                alternate randomly; strongly negative means they come in
                long same-sign stretches, i.e. coherent structure. This is
                what "did the wave flatten" actually means.
    acf1        lag-1 autocorrelation of the pulls, same question.
    """
    import math
    p, ratio_dev, stat_var = [], [], []
    for d, e, f in zip(data_y, data_e, fit_y):
        if e <= 0 or d <= 0 or f <= 0:
            continue
        p.append((d - f) / e)
        ratio_dev.append((d / f - 1.0) ** 2)
        stat_var.append((e / f) ** 2)
    n = len(p)
    if n < 5:
        return None

    chi2 = sum(v * v for v in p) / n
    frac = max(sum(ratio_dev) / n - sum(stat_var) / n, 0.0) ** 0.5

    signs = [1 if v > 0 else -1 for v in p]
    npos = sum(1 for sgn in signs if sgn > 0)
    nneg = n - npos
    runs = 1 + sum(1 for a, b in zip(signs, signs[1:]) if a != b)
    if npos and nneg:
        exp_runs = 1.0 + 2.0 * npos * nneg / n
        var = (2.0 * npos * nneg * (2.0 * npos * nneg - n)
               / (n * n * (n - 1.0)))
        z = (runs - exp_runs) / math.sqrt(var) if var > 0 else float("nan")
    else:
        z = float("nan")

    mean = sum(p) / n
    num = sum((p[i] - mean) * (p[i + 1] - mean) for i in range(n - 1))
    den = sum((v - mean) ** 2 for v in p)
    return {"n_bins": n, "chi2_ndf": chi2, "pull_width": chi2 ** 0.5,
            "frac_rms": frac, "runs": runs, "runs_z": z,
            "acf1": num / den if den > 0 else 0.0}


def binned_expectation(hist, pdf, x, norm):
    """Expected counts per bin for an ordinary 1D pdf, for the metrics."""
    b_lo, b_hi = sf.range_bins(hist)
    nset = ROOT.RooArgSet(x)
    saved = x.getVal()
    xs, ys, es, fs = [], [], [], []
    for i in range(b_lo, b_hi + 1):
        x.setVal(hist.GetBinCenter(i))
        xs.append(hist.GetBinCenter(i))
        ys.append(hist.GetBinContent(i))
        es.append(hist.GetBinError(i))
        fs.append(norm * pdf.getVal(nset) * hist.GetBinWidth(i))
    x.setVal(saved)
    return xs, ys, es, fs


def plot_ratio_overlay(cond, ref1d, label, outbase, year, save_pdf):
    """Both Template/Fit curves on one panel.

    The single most direct answer to the question: if the conditional
    model has flattened the wave, its points scatter about 1 where the 1D
    points still swing coherently.
    """
    keep = []
    frame = ROOT.TH1F("ovr_frame", ";m_{#mu#mu} [GeV];Template/Fit",
                      1, sf.XLOW, sf.XHIGH)
    frame.SetDirectory(0)
    frame.SetStats(0)
    frame.SetMinimum(0.80)
    frame.SetMaximum(1.20)
    sf._style_ratio_axis(frame.GetXaxis(), 1.15)
    sf._style_ratio_axis(frame.GetYaxis(), 1.60)

    c = _single_canvas("c_ratio_overlay")
    frame.Draw()

    leg = ROOT.TLegend(sf.PAD_LEFT_MARGIN + 0.04, 0.72,
                       sf.PAD_LEFT_MARGIN + 0.46, 0.86)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.028)

    for (xs, ys, es, fs), col, name in (
            (cond, sf.fit_color(), "conditional"),
            (ref1d, ROOT.kAzure + 2, "binned 1D")):
        g = ROOT.TGraphErrors()
        k = 0
        for xv, d, e, f in zip(xs, ys, es, fs):
            if e <= 0 or d <= 0 or f <= 0:
                continue
            g.SetPoint(k, xv, d / f)
            g.SetPointError(k, 0.0, e / f)
            k += 1
        g.SetMarkerStyle(sf.DATA_MARKER_STYLE)
        g.SetMarkerSize(sf.DATA_MARKER_SIZE)
        g.SetMarkerColor(col)
        g.SetLineColor(col)
        g.Draw("P SAME")
        leg.AddEntry(g, name, "p")
        keep.append(g)

    one = ROOT.TLine(sf.XLOW, 1.0, sf.XHIGH, 1.0)
    one.SetLineColor(ROOT.kBlack)
    one.SetLineWidth(2)
    one.Draw("same")
    leg.Draw()
    keep += [frame, leg, one, sf.cms_label(c, year)]

    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextFont(42)
    tex.SetTextSize(sf.ANN_TEXT_SIZE)
    dy = sf.LINE_SPACING * sf.ANN_TEXT_SIZE
    for i, ln in enumerate(sf.header_lines(label)):
        tex.DrawLatex(sf.PAD_LEFT_MARGIN + 0.04, 0.88 - i * dy, ln)
    keep.append(tex)

    c.SaveAs(f"{outbase}_ratio_overlay.png")
    if save_pdf:
        c.SaveAs(f"{outbase}_ratio_overlay.pdf")
    del c


# ---------------------------------------------------------------------------
# plots
# ---------------------------------------------------------------------------

def _single_canvas(name, right=None):
    c = ROOT.TCanvas(name, "", sf.CANVAS_W, sf.CANVAS_H)
    c.SetTopMargin(sf.PAD_TOP_MARGIN)
    c.SetBottomMargin(sf.PAD_BOTTOM_MARGIN)
    c.SetLeftMargin(sf.PAD_LEFT_MARGIN)
    c.SetRightMargin(sf.PAD_RIGHT_MARGIN if right is None else right)
    c.SetTickx()
    c.SetTicky()
    return c


def plot_mass_projection(x, r, data, pdf, nom, label, outbase, year,
                         save_pdf, n_eff=0.0, nbins=None):
    """Mass projection of the conditional model, with ratio and pull.

    The curve is the model AVERAGED OVER THE OBSERVED r VALUES
    (ProjWData), which is what a conditional model predicts for the
    marginal mass distribution -- not the model at any single r.

    The ratio and pull are read off the RooPlot's own histogram and curve,
    using RooCurve::average over each bin. That avoids re-deriving the
    normalisation by hand, and avoids the bin/curve-point misalignment that
    indexing into the curve would risk.
    """
    nbins = nbins or int((sf.XHIGH - sf.XLOW) * sf.BINS_PER_GEV)
    x.setBins(nbins)

    frame = x.frame(ROOT.RooFit.Title(""))
    data.plotOn(frame, ROOT.RooFit.Name("dat"),
                ROOT.RooFit.MarkerStyle(sf.DATA_MARKER_STYLE),
                ROOT.RooFit.MarkerSize(sf.DATA_MARKER_SIZE),
                ROOT.RooFit.MarkerColor(ROOT.kBlack),
                ROOT.RooFit.LineColor(ROOT.kBlack),
                ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
    pdf.plotOn(frame, ROOT.RooFit.Name("fit"),
               ROOT.RooFit.ProjWData(ROOT.RooArgSet(r), data),
               ROOT.RooFit.LineColor(sf.fit_color()),
               ROOT.RooFit.LineWidth(2))

    n_float = sum(0 if v.isConstant() else 1 for v in nom.values())
    chi2_ndf = frame.chiSquare("fit", "dat", n_float)

    h = frame.getHist("dat")
    curve = frame.getCurve("fit")
    bw = (sf.XHIGH - sf.XLOW) / nbins

    ratio = ROOT.TH1D(f"ratio_{outbase.split('/')[-1]}", "",
                      nbins, sf.XLOW, sf.XHIGH)
    ratio.SetDirectory(0)
    # GetPointX/GetPointY rather than GetPoint(i, x&, y&): the
    # pass-by-reference form needs ctypes.c_double in PyROOT and is easy
    # to get wrong.
    bin_x, bin_y, bin_e, bin_f = [], [], [], []
    for i in range(h.GetN()):
        xv = h.GetPointX(i)
        yv = h.GetPointY(i)
        ey = h.GetErrorY(i)
        f = curve.average(xv - 0.5 * bw, xv + 0.5 * bw)
        bin_x.append(xv)
        bin_y.append(yv)
        bin_e.append(ey)
        bin_f.append(f)
        b = ratio.FindBin(xv)
        if yv <= 0 or f <= 0:
            ratio.SetBinContent(b, -999.0)
            ratio.SetBinError(b, 0.0)
            continue
        ratio.SetBinContent(b, yv / f)
        ratio.SetBinError(b, ey / f)

    ratio.SetMarkerStyle(sf.DATA_MARKER_STYLE)
    ratio.SetMarkerSize(sf.DATA_MARKER_SIZE)
    ratio.SetMarkerColor(ROOT.kBlack)
    ratio.SetLineColor(ROOT.kBlack)
    ratio.GetYaxis().SetRangeUser(*RATIO_Y_RANGE)
    ratio.SetTitle(";m_{#mu#mu} [GeV];Template/Fit")
    sf._style_ratio_axis(ratio.GetXaxis(), sf.RATIO_X_TITLE_OFFSET)
    sf._style_ratio_axis(ratio.GetYaxis(), sf.RATIO_Y_TITLE_OFFSET)

    pull_frame = x.frame(ROOT.RooFit.Title(""))
    ph = frame.pullHist("dat", "fit")
    for i in range(ph.GetN() - 1, -1, -1):
        if ph.GetErrorYhigh(i) <= 0 and ph.GetErrorYlow(i) <= 0:
            ph.RemovePoint(i)
    ph.SetMarkerStyle(sf.DATA_MARKER_STYLE)
    ph.SetMarkerSize(sf.DATA_MARKER_SIZE)
    ph.SetMarkerColor(ROOT.kBlack)
    ph.SetLineColor(ROOT.kBlack)
    pull_frame.addPlotable(ph, "P")
    pull_frame.SetMinimum(-6.0)
    pull_frame.SetMaximum(6.0)
    # exactly sigFit's label, so the two pull panels are comparable at a
    # glance rather than needing the axis read each time
    pull_frame.GetYaxis().SetTitle("#frac{Template#minusFit}{#deltaTemplate}")
    pull_frame.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
    sf._style_ratio_axis(pull_frame.GetXaxis(), sf.RATIO_X_TITLE_OFFSET)
    sf._style_ratio_axis(pull_frame.GetYaxis(), sf.RATIO_Y_TITLE_OFFSET)

    frame.SetTitle(f";;Events / {bw:.2f} GeV")
    sf._style_ratio_axis(frame.GetYaxis(), sf.RATIO_Y_TITLE_OFFSET)
    frame.GetXaxis().SetLabelSize(0)
    frame.GetXaxis().SetTitleSize(0)
    frame.SetMaximum(1.1 * frame.GetMaximum()) # headroom

    os.makedirs(os.path.dirname(outbase), exist_ok=True)
    for kind, lower, central, guides, draw in (
            ("ratio", ratio, 1.0, sf.RATIO_GUIDES, "pe"),
            ("pull", pull_frame, 0.0, sf.PULL_GUIDES, "")):
        canvas, outer, pad1, pad2 = sf.make_canvas_pads(f"c_cond_{kind}")
        pad1.cd()
        frame.Draw()
        keep_txt = _annotate(nom, chi2_ndf, data.sumEntries(),
                             data.numEntries(), n_eff, label)
        keep_cms = sf.cms_label(pad1, year)
        keep_leg = sf.fit_legend(frame,
                                 model="DSCB model (#sigma = #kappa r M_{H})")

        pad2.cd()
        lower.Draw(draw)
        l0 = ROOT.TLine(sf.XLOW, central, sf.XHIGH, central)
        l0.SetLineColor(sf.fit_color())
        l0.SetLineWidth(2)
        l0.Draw("same")
        lines = [l0]
        for gy in guides:
            ln = ROOT.TLine(sf.XLOW, gy, sf.XHIGH, gy)
            ln.SetLineColor(11)
            ln.SetLineStyle(ROOT.kDashed)
            ln.Draw("same")
            lines.append(ln)

        canvas.Update()
        canvas.SaveAs(f"{outbase}_{kind}.png")
        if save_pdf:
            canvas.SaveAs(f"{outbase}_{kind}.pdf")
        del keep_txt, keep_cms, keep_leg, lines, pad1, pad2, outer, canvas

    return chi2_ndf, (bin_x, bin_y, bin_e, bin_f)


def _annotate(nom, chi2_ndf, yield_, n_mc, n_eff, label):
    """Parameter block. kappa first: it is the number this fit is for."""
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(sf.ANN_TEXT_SIZE)
    left = sf.PAD_LEFT_MARGIN + 0.05
    dy = sf.LINE_SPACING * sf.ANN_TEXT_SIZE
    for i, line in enumerate(sf.header_lines(label)):
        latex.DrawLatex(left, 0.85 - i * dy, line)

    labels = {"kappa": ("#kappa", ""), "c0": ("c_{0}", "GeV"),
              "mu": ("#mu", "GeV"),
              "alphaL": ("#alpha_{lo}", ""), "nL": ("n_{lo}", ""),
              "alphaR": ("#alpha_{high}", ""), "nR": ("n_{high}", "")}
    latex.SetTextSize(sf.PARAM_TEXT_SIZE)
    latex.SetTextAlign(12)
    right = 1.0 - sf.PAD_RIGHT_MARGIN - 0.36
    # three different things, and they are easy to confuse:
    #   yield  sum of weights = expected signal events (NOT fitted: this
    #          fit is not extended)
    #   N_MC   raw simulated events that entered the fit
    #   N_eff  (sum w)^2 / sum w^2, the statistical power
    lines = [f"yield = {yield_:.1f}",
             f"N_{{MC}} = {n_mc:.0f}",
             f"N_{{eff}} = {n_eff:.0f}"]
    order = ["kappa"] + (["c0"] if "c0" in nom else []) + \
            ["mu", "alphaL", "nL", "alphaR", "nR"]
    for k in order:
        v = nom[k]
        lab, unit = labels[k]
        txt = sf.fmt_with_unc(lab, v.getVal(), v.getError(), unit)
        if sf.near_bound(v):
            txt += " #color[2]{(bound)}"
        lines.append(txt)
    lines.append(f"#chi^{{2}}/ndf = {chi2_ndf:.2f}")
    pdy = sf.LINE_SPACING * sf.PARAM_TEXT_SIZE
    for i, line in enumerate(lines):
        latex.DrawLatex(right, 0.86 - i * pdy, line)
    return latex


def plot_resolution(r, r_values, w_values, keys, label, outbase, year,
                    save_pdf, n_curve=400):
    """Kernel estimate over the MC distribution of r -- a validation plot.

    Black points: the weighted MC distribution of sigma_m/m, normalised to
    a probability density. Red line: the RooKeysPdf built from a subsample
    of the same events. They should lie on top of each other; where they
    do not, the kernel is not describing the distribution and anything
    later that uses it (toys, a full 2D fit) inherits the discrepancy.

    Two things have to match for the comparison to be meaningful, and both
    were wrong in the first version of this function:

      * The NORMALISATION RANGE. RooKeysPdf is normalised over the whole
        of r's range, so the histogram must be too -- normalising it over
        a narrower display window would scale it up and the curve would
        sit below the points for no reason. Filled over the full range,
        then the x axis is zoomed for display.

      * The DRAWING. A RooPlot cannot be overlaid on a TH1 with "same":
        it brings its own frame and its own axis scaling. The kernel is
        evaluated into a TGraph instead, which is drawn in the pad's own
        coordinates like any other graph.
    """
    nb = 2 * R_NBINS
    h = ROOT.TH1D("h_r", "", nb, R_LO, R_HI)     # FULL range, as the pdf
    h.SetDirectory(0)
    for v, w in zip(r_values, w_values):
        h.Fill(float(v), float(w))
    if h.Integral() > 0:
        h.Scale(1.0 / h.Integral() / h.GetBinWidth(1))   # -> density

    h.SetMarkerStyle(sf.DATA_MARKER_STYLE)
    h.SetMarkerSize(sf.DATA_MARKER_SIZE)
    h.SetMarkerColor(ROOT.kBlack)
    h.SetLineColor(ROOT.kBlack)
    h.SetTitle(";#sigma_{m}/m;probability density")
    sf._style_ratio_axis(h.GetXaxis(), 1.15)
    sf._style_ratio_axis(h.GetYaxis(), 1.60)
    h.GetXaxis().SetRangeUser(R_LO, R_PLOT_HI)   # zoom, do not renormalise

    # the kernel, evaluated point by point
    nset = ROOT.RooArgSet(r)
    saved = r.getVal()
    g = ROOT.TGraph(n_curve)
    step = (R_PLOT_HI - R_LO) / (n_curve - 1)
    ymax = h.GetMaximum()
    for i in range(n_curve):
        v = R_LO + i * step
        r.setVal(v)
        y = keys.getVal(nset)
        g.SetPoint(i, v, y)
        ymax = max(ymax, y)
    r.setVal(saved)
    g.SetLineColor(sf.fit_color())
    g.SetLineWidth(3)

    h.SetMaximum(1.35 * ymax)

    c = _single_canvas("c_res")
    h.Draw("pe")
    g.Draw("L SAME")

    leg = ROOT.TLegend(1.0 - sf.PAD_RIGHT_MARGIN - 0.42, 0.68,
                       1.0 - sf.PAD_RIGHT_MARGIN - 0.03, 0.80)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.028)
    leg.AddEntry(h, "MC (weighted)", "pe")
    leg.AddEntry(g, "RooKeysPdf (MirrorLeft)", "l")
    leg.Draw()

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(sf.ANN_TEXT_SIZE)
    dy = sf.LINE_SPACING * sf.ANN_TEXT_SIZE
    med = float(np.median(r_values))
    for i, line in enumerate(sf.header_lines(label) +
                             [f"median = {med:.5f}"]):
        latex.DrawLatex(sf.PAD_LEFT_MARGIN + 0.05, 0.85 - i * dy, line)

    os.makedirs(os.path.dirname(outbase), exist_ok=True)
    c.SaveAs(f"{outbase}_resolution.png")
    if save_pdf:
        c.SaveAs(f"{outbase}_resolution.pdf")

    # how well does the kernel actually describe the points?
    dev = []
    for i in range(1, h.GetNbinsX() + 1):
        xv = h.GetBinCenter(i)
        if xv > R_PLOT_HI or h.GetBinContent(i) <= 0:
            continue
        r.setVal(xv)
        y = keys.getVal(nset)
        if y > 0:
            dev.append(h.GetBinContent(i) / y - 1.0)
    r.setVal(saved)
    if dev:
        rms = (sum(d * d for d in dev) / len(dev)) ** 0.5
        print(f"  kernel vs MC: RMS of (MC/kernel - 1) = {rms:.2%} "
              f"over {len(dev)} bins")

    del c, leg, g
    return h


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--category", default="ggHcat",
                    choices=sorted(MODES))
    ap.add_argument("-b", "--bin", default="incl")
    ap.add_argument("-s", "--sig", default=None,
                    help="production mode; default = the first for the "
                         "category")
    ap.add_argument("-y", "--year", default="Run3")
    ap.add_argument("--max-events", type=int, default=MAX_EVENTS,
                    help=f"cap on the unbinned dataset (default "
                         f"{MAX_EVENTS}). An unbinned conditional fit is "
                         f"O(N) per likelihood call, so full statistics is "
                         f"impractical -- but the 1D reference runs on the "
                         f"same subsample, so the comparison stays fair")
    ap.add_argument("--keys-max", type=int, default=KEYS_MAX,
                    help=f"cap on the kernel estimate (default {KEYS_MAX}); "
                         f"RooKeysPdf is O(N) per evaluation")
    ap.add_argument("--keys-rho", type=float, default=KEYS_RHO,
                    help="RooKeysPdf bandwidth scale (default 1.0)")
    ap.add_argument("--offset", action="store_true",
                   help="add an intrinsic floor in quadrature, "
                        "sigma = sqrt((kappa*r*MH)^2 + c0^2). Tests "
                        "whether the width really is proportional to the "
                        "per-event error through the origin.")
    ap.add_argument("--resume", action="store_true",
                   help="reuse the cached events and the fitted parameters "
                        "from a previous run of the same tag, and skip "
                        "both the RDataFrame pass and the fit. For "
                        "iterating on the plots.")
    ap.add_argument("--refit", action="store_true",
                   help="with --resume, reuse the cached events but fit "
                        "again")
    ap.add_argument("--no-cache", action="store_true",
                   help="do not write the event cache")
    ap.add_argument("--no-1d", action="store_true",
                    help="skip the binned 1D reference fit. It runs on the "
                         "SAME subsample, which is the point: it makes "
                         "chi2/ndf comparable without either fit needing "
                         "full statistics")
    ap.add_argument("--no-keys", action="store_true",
                    help="skip the kernel estimate; the conditional fit "
                         "does not need it, only the plots and toys do")
    ap.add_argument("--free-nuis", action="store_true",
                    help="leave the scale/resolution nuisances floating "
                         "(they are degenerate with mu and kappa -- for "
                         "debugging only)")
    ap.add_argument("--wsdir", default="WS_COND")
    ap.add_argument("--plotdir",
                    default=os.path.expanduser(
                        "~/public_html/HmumuFits/signal_fits_cond"))
    ap.add_argument("--pdf", action="store_true")
    ap.add_argument("--seed", type=int, default=0)
    args = ap.parse_args()

    sf.setup_style()
    sig = args.sig or MODES[args.category][0]
    tag = f"{args.category}_{args.bin}_{args.year}_{sig}"
    label = f"{args.category} {args.bin} {sig}"

    print(f"category : {args.category}"
          f"\nbin      : {args.bin}"
          f"\nyear     : {args.year}"
          f"\nsample   : {sig}"
          f"\nmass     : {MASSCOL}"
          f"\nrel. err : {RELERRCOL}  in ({R_LO}, {R_HI})"
          f"\nmodel    : DSCB, sigma_CB = kappa * r * MH")

    # ---- observables --------------------------------------------------
    x = ROOT.RooRealVar(f"mh{args.category}", "m_{#mu#mu}",
                        sf.XLOW, sf.XHIGH)
    x.setRange("full", sf.XLOW, sf.XHIGH)
    r = ROOT.RooRealVar("massRelErr", "#sigma_{m}/m", R_LO, R_HI)
    w = ROOT.RooRealVar("wgt", "weight", 1.0, -1e9, 1e9)

    cached = load_cache(args.wsdir, tag) if args.resume else None
    if cached is not None:
        m_all, r_all, w_all, files = cached
        data = build_dataset(x, r, w, m_all, r_all, w_all, args.category)
        print(f"rebuilt dataset: {data.numEntries()} entries, "
              f"sum of weights {data.sumEntries():.2f}")
    else:
        data, m_all, r_all, w_all, files = load_dataset(
            x, r, w, args.category, args.bin, args.year, sig,
            args.max_events, args.seed)
        if not args.no_cache:
            save_cache(args.wsdir, tag, m_all, r_all, w_all, files)

    n_eff = effective_entries(w_all)
    print(f"  N_MC  (raw events fitted) = {len(m_all)}")
    print(f"  N_eff (statistical power) = {n_eff:.0f}")
    print(f"  yield (sum of weights)    = {data.sumEntries():.2f}"
          "   <- expected signal events, not a fitted quantity")

    # ---- model ---------------------------------------------------------
    pdf, bundle = make_cond_pdf(x, r, tag, offset=args.offset)
    nom, nuis = bundle["nom"], bundle["nuis"]

    # The nuisances multiply exactly what mu and kappa already scale, so
    # they are pinned for the MC fit and released afterwards -- same
    # reasoning as sigFit.
    for v in nuis.values():
        v.setVal(0.0)
        v.setConstant(not args.free_nuis)

    fitopts = [
        ROOT.RooFit.Minimizer("Minuit2"),
        ROOT.RooFit.Strategy(2),
        ROOT.RooFit.Save(True),
        ROOT.RooFit.Range("full"),
        # r is CONDITIONAL: P(r) cancels and never enters the likelihood
        ROOT.RooFit.ConditionalObservables(ROOT.RooArgSet(r)),
        # correct errors for a weighted UNBINNED fit (unlike the binned
        # case, where SumW2Error is the right choice)
        ROOT.RooFit.AsymptoticError(True),
        ROOT.RooFit.PrintLevel(-1),
    ]
    jf = os.path.join(args.wsdir, f"SignalCond_{tag}_params.json")
    prev = None
    if args.resume and not args.refit and os.path.isfile(jf):
        with open(jf) as fh:
            prev = json.load(fh)

    # A cached fit is only valid for the SAME model. Restoring a
    # no-offset fit into an --offset model would leave c0 at its start
    # value and silently produce a shape nobody fitted.
    if prev:
        stored_offset = prev.get("shape_model", {}).get("offset", False)
        stored_pars = set(prev.get("parameters", {}))
        wanted_pars = set(nom)
        if stored_offset != bool(args.offset) or stored_pars != wanted_pars:
            print(f"\ncached fit in {jf} is for a different model "
                  f"(offset={stored_offset}, parameters={sorted(stored_pars)});"
                  f"\nrefitting instead of restoring it")
            prev = None

    if prev:
        print(f"\nrestoring the fit from {jf} (no refit)")
        for k, v in nom.items():
            if k in prev["parameters"]:
                v.setVal(prev["parameters"][k]["value"])
                v.setError(prev["parameters"][k]["error"])
        q = prev["fit_quality"]
        status, cov, edm = q["status"], q["cov_qual"], q["edm"]
        parked, ok = q["params_at_bound"], q["ok"]
    else:
        print("\nfitting (conditional on r)...")
        pdf.fitTo(data, *fitopts)
        result = pdf.fitTo(data, *fitopts)
        status, cov, edm = result.status(), result.covQual(), result.edm()
        parked = [k for k, v in nom.items()
                  if not v.isConstant() and sf.near_bound(v)]
        ok = status == 0 and cov == 3 and not parked

    print(f"\n  kappa  = {nom['kappa'].getVal():7.4f} "
          f"+/- {nom['kappa'].getError():.4f}"
          "      <-- expect ~1 if the per-event error is calibrated")
    print(f"  mu     = {nom['mu'].getVal():7.3f} "
          f"+/- {nom['mu'].getError():.3f}")
    if "c0" in nom:
        c0, c0e = nom["c0"].getVal(), nom["c0"].getError()
        sig = f"{c0 / c0e:.1f} sigma from zero" if c0e > 0 else "no error"
        print(f"  c0     = {c0:7.4f} +/- {c0e:.4f} GeV      <-- {sig}")
    for k in ("alphaL", "nL", "alphaR", "nR"):
        print(f"  {k:6s} = {nom[k].getVal():7.3f} "
              f"+/- {nom[k].getError():.3f}")
    print(f"  status={status} covQual={cov} edm={edm:.2e}"
          + ("" if ok else "   <-- CHECK THIS FIT"))
    if parked:
        print(f"  !! at bound: {', '.join(parked)}")

    # ---- save the fit NOW, before anything that can fail ---------------
    # A plotting error must not cost the fit. chi2_ndf is not known yet;
    # it is filled in after the projection and the file rewritten.
    os.makedirs(args.wsdir, exist_ok=True)
    record = {
        "provenance": {
            "written_utc": datetime.now(timezone.utc)
                           .isoformat(timespec="seconds"),
            "git_commit": sf.git_commit(),
            "category": args.category, "bin": args.bin,
            "year": args.year, "sig": sig, "files": files,
            "mass_column": MASSCOL, "relerr_column": RELERRCOL,
            "r_range": [R_LO, R_HI],
            "n_entries": data.numEntries(),
            "max_events": args.max_events,
            "keys_max": args.keys_max, "keys_rho": args.keys_rho,
            "conditional": True,
        },
        "shape_model": bundle["cfg"],
        "parameters": {k: {"value": v.getVal(), "error": v.getError(),
                           "min": v.getMin(), "max": v.getMax(),
                           "at_bound": sf.near_bound(v)}
                       for k, v in nom.items()},
        "normalization": {"sum_weights": data.sumEntries(),
                          "n_mc": data.numEntries(),
                          "n_eff": n_eff,
                          "note": ("sum_weights is the expected signal "
                                   "yield; this fit is not extended, so "
                                   "it is not a fitted quantity")},
        "fit_quality": {"chi2_ndf": None, "status": status,
                        "cov_qual": cov, "edm": edm,
                        "params_at_bound": parked, "ok": ok},
        "reference_1d": None,
    }
    with open(jf, "w") as fh:
        json.dump(record, fh, indent=2, sort_keys=True)
    print(f"saved the fit to {jf} (chi2 and the 1D reference follow)")

    # ---- resolution kernel --------------------------------------------
    keys = rset = None
    if not args.no_keys:
        print()
        keys, rset = make_keys_pdf(r, r_all, w_all, tag,
                                   args.keys_max, args.keys_rho, args.seed)

    # ---- plots ---------------------------------------------------------
    plotdir = os.path.join(args.plotdir, args.category)
    outbase = os.path.join(plotdir, f"cond_{tag}")
    chi2_ndf, cond_bins = plot_mass_projection(
        x, r, data, pdf, nom, label, outbase, args.year, args.pdf,
        n_eff=n_eff)
    print(f"\n  chi2/ndf (projected) = {chi2_ndf:.3f}"
          "     <-- compare with the 1D DSCB on the same sample")
    if keys is not None:
        plot_resolution(r, r_all, w_all, keys, label, outbase, args.year,
                        args.pdf)

    # ---- the binned 1D reference, on the same events -------------------
    ref = None
    if not args.no_1d:
        print("\n--- binned 1D reference (same subsample) ---")
        x1d = ROOT.RooRealVar(f"mh{args.category}_1d", "m_{#mu#mu}",
                              sf.XLOW, sf.XHIGH)
        x1d.setRange("full", sf.XLOW, sf.XHIGH)
        ref = fit_1d_reference(x1d, m_all, w_all, tag, label, plotdir,
                               args.year, args.pdf)

    if ref:
        # sigma_eff of the conditional model: the width it predicts for the
        # typical event. Not the same object as a fitted constant width, but
        # the closest comparable number.
        kap = nom["kappa"].getVal()
        r_med = float(np.median(r_all))
        sig_eff = kap * r_med * sf.MH_REF
        if "c0" in nom:
            sig_eff = (sig_eff ** 2 + nom["c0"].getVal() ** 2) ** 0.5
        print("\n" + "=" * 62)
        print(f"{'':<22}{'conditional':>18}{'binned 1D':>18}")
        print("-" * 62)
        print(f"{'N_MC (raw)':<22}{data.numEntries():>18d}"
              f"{data.numEntries():>18d}")
        print(f"{'N_eff':<22}{n_eff:>18.0f}{ref['n_eff']:>18.0f}")
        print(f"{'yield (sum w)':<22}{data.sumEntries():>18.2f}"
              f"{data.sumEntries():>18.2f}")
        print(f"{'mu [GeV]':<22}{nom['mu'].getVal():>18.4f}"
              f"{ref['mu']:>18.4f}")
        print(f"{'width [GeV]':<22}{sig_eff:>18.4f}{ref['sigma']:>18.4f}")
        print(f"{'chi2/ndf':<22}{chi2_ndf:>18.3f}"
              f"{ref['chi2_ndf']:>18.3f}")
        print("=" * 62)
        print(f"conditional width = kappa * median(r) * MH "
              f"= {kap:.4f} * {r_med:.5f} * {sf.MH_REF:g}")
        # ---- the metrics that answer the actual question -------------
        mc_ = mismatch_metrics(cond_bins[1], cond_bins[2], cond_bins[3])
        m1_ = mismatch_metrics(ref["bins"][1], ref["bins"][2],
                               ref["bins"][3])
        if mc_ and m1_:
            print()
            print(f"{'bins used':<22}{mc_['n_bins']:>18d}"
                  f"{m1_['n_bins']:>18d}")
            print(f"{'pull width':<22}{mc_['pull_width']:>18.3f}"
                  f"{m1_['pull_width']:>18.3f}")
            print(f"{'frac RMS of T/F':<22}{mc_['frac_rms']:>17.2%}"
                  f"{m1_['frac_rms']:>18.2%}")
            print(f"{'runs test z':<22}{mc_['runs_z']:>+18.1f}"
                  f"{m1_['runs_z']:>+18.1f}")
            print(f"{'lag-1 acf':<22}{mc_['acf1']:>+18.3f}"
                  f"{m1_['acf1']:>+18.3f}")
            print("=" * 62)
            print("frac RMS is statistics-independent: it is the honest "
                  "'how wrong is\nthe shape' number. runs z and the "
                  "autocorrelation say whether what\nis left is "
                  "STRUCTURED (both far from 0) or just noise (both ~0).")

            better_frac = mc_["frac_rms"] < m1_["frac_rms"]
            flatter = abs(mc_["runs_z"]) < abs(m1_["runs_z"])
            print(f"\n  shape mismatch     : "
                  f"{'smaller' if better_frac else 'LARGER'} "
                  f"({mc_['frac_rms']:.2%} vs {m1_['frac_rms']:.2%})")
            print(f"  coherent structure : "
                  f"{'weaker' if flatter else 'NOT weaker'} "
                  f"(|z| {abs(mc_['runs_z']):.1f} vs "
                  f"{abs(m1_['runs_z']):.1f})")

        plot_ratio_overlay(cond_bins, ref["bins"], label, outbase,
                           args.year, args.pdf)

    # ---- freeze and write ----------------------------------------------
    norm = data.sumEntries()
    for v in nom.values():
        v.setConstant(True)
    for v in nuis.values():
        v.setConstant(False)
    norm_var = ROOT.RooRealVar(f"{pdf.GetName()}_norm",
                               f"{pdf.GetName()}_norm", norm)
    norm_var.setConstant(True)

    os.makedirs(args.wsdir, exist_ok=True)
    ws = ROOT.RooWorkspace("w", "workspace")
    getattr(ws, "import")(pdf)
    getattr(ws, "import")(norm_var)
    if keys is not None:
        getattr(ws, "import")(keys)
    wsfile = os.path.join(args.wsdir, f"SignalCond_{tag}_workspace.root")
    ws.writeToFile(wsfile)

    record["fit_quality"]["chi2_ndf"] = chi2_ndf
    record["normalization"]["sum_weights"] = norm
    record["metrics"] = {
        "conditional": mismatch_metrics(cond_bins[1], cond_bins[2],
                                        cond_bins[3]),
        "binned_1d": (mismatch_metrics(ref["bins"][1], ref["bins"][2],
                                       ref["bins"][3]) if ref else None),
    }
    record["reference_1d"] = ({"mu": ref["mu"], "sigma": ref["sigma"],
                               "chi2_ndf": ref["chi2_ndf"],
                               "ok": ref["ok"], "n_eff": ref["n_eff"]}
                              if ref else None)
    with open(jf, "w") as fh:
        json.dump(record, fh, indent=2, sort_keys=True)

    print(f"\nwrote {wsfile}")
    print(f"wrote {jf}")
    print(f"plots {outbase}_{{ratio,pull,resolution}}.png")
    observables = {x.GetName(), r.GetName()}
    floating = [v.GetName() for v in ws.allVars()
                if not v.isConstant() and v.GetName() not in observables]
    print(f"floating in workspace: {floating or '(none)'}")
    if not ok:
        print("WARNING: this fit is flagged, check the plots")


if __name__ == "__main__":
    main()
