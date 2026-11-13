"""
bkgFit.py -- blinded background fit, with an F-test-selected envelope.

Blinding
--------
    Higgs SR   120 - 130 GeV        blinded
    sideband   110-120 + 130-150    what the fit sees

The normalisation is extrapolated from the sideband.

The full-range data is saved in the workspace as `datahist_<tag>` for the observed results.

Model selection
---------------

  1. F-test within each nested family: take the lowest order where
     adding the next parameter is not justified at p < F_THRESHOLD. Each
     order is fitted from its default start and from the fit of the order
     below (seed_values); the better one is kept, so a higher order does
     not lose to a local minimum. A failed fit (status != 0) ends its
     family: the order below it is kept.
     The families are those of the Run 2 analysis note, Sec. 5.3
     (RUN2_FAMILIES in pdfDefinitions.py): BWZ x Bernstein, BWZGamma,
     BWZRedux, FEWZ x Bernstein, Bernstein, exponential and power-law
     sums, exponential / power law + constant, Landau x Bernstein.
     --extra-models adds EXTRA_FAMILIES: the DY-core families dyexp and
     dypow and -- when mcCoreFit.py has written MCcore_<tag>.json -- the
     MC-tuned-core families mcexp, mcpow and mcdy.
  2. Goodness-of-fit across families: keep every family's choice whose
     chi2 p-value exceeds GOF_THRESHOLD.
  3. Order the survivors by penalised NLL, NLL + 0.5 * npar, as Combine
     compares them: the lowest is pdfindex 0, the default member.

Everything that survives goes into the RooMultiPdf envelope, which covers
the ignorance about the shape.

Outputs
-------
This script performs fits, but it does not write the Combine workspace.
It runs in conda (ROOT >= 6.28) and writes everything the workspace needs:

    Bkg_<tag>_params.json   envelope, fitted parameters of every
                            candidate, extrapolated norm
    Bkg_<tag>_data.root     TH1D data_<tag>, the FULL-range data

makeWorkspaces.py then rebuilds the PDFs from these inside the Combine
container and writes the RooMultiPdf workspace bwsHrare.py reads, with
tag = "<cat>_<bin>_<year>":

    multipdf_<tag>_bkg        RooMultiPdf
    pdfindex_<tag>            RooCategory
    multipdf_<tag>_bkg_norm   RooRealVar, floating
    datahist_<tag>            RooDataHist, FULL range

Note the tag is built exactly as bwsHrare does, so an empty --bin gives
the double underscore ("VBFcat__Run3") that the existing datacards have.

Usage
-----
    # fits, plots, JSON -- conda (mcCoreFit.py first, for the mc* families)
    python bkgFit.py -c ggHcat -y Run3 -b incl

    # ...then the cards and workspace -- inside the Combine container
    python buildWorkspace.py -c ggHcat -y Run3 -b incl

    python bkgFit.py -c ggHcat -y Run3 -b ""    # VBFcat__Run3 style
"""

import argparse
import json
import os
from datetime import datetime, timezone

import ROOT

from fitlib import config, mcCore, style
from fitlib.gof import (CHI2_MIN_BIN_SIGNIFICANCE, GOF_THRESHOLD,
                        merge_chi2_bins)
from fitlib import pdfDefinitions as pdefs


ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
for _i in range(ROOT.RooMsgService.instance().numStreams()):
    ROOT.RooMsgService.instance().getStream(_i).removeTopic(ROOT.RooFit.Plotting)

# ---------------------------------------------------------------------------
# blinding
# ---------------------------------------------------------------------------
# Blinded window for the BACKGROUND FIT (see the module docstring for why
# it is narrower than `isH`). Everything downstream -- the named ranges,
# the chi2 bin merging, the shaded box on the plot, f_sb and so the
# extrapolated norm -- follows from these two numbers.
SR_LO, SR_HI = 120.0, 130.0        # blinded Higgs signal region
R_SB_LOW = "sbLow"                 # 110 - 120
R_SB_HIGH = "sbHigh"               # 130 - 150
R_FULL = "full"
R_SR = "sr"

# ---------------------------------------------------------------------------
# model selection
# ---------------------------------------------------------------------------
# F-test: add a parameter only if 2*dNLL against chi2_(dnpar) gives
# p < F_THRESHOLD. 0.05 is the usual choice.
F_THRESHOLD = 0.05
# Penalty per free parameter on the sideband NLL, as RooMultiPdf applies it
# when Combine profiles pdfindex (its correction factor, 0.5). The envelope
# member with the lowest penalised NLL becomes pdfindex 0, the default.
NLL_PENALTY_PER_PAR = 0.5
# GOF_THRESHOLD and CHI2_MIN_BIN_SIGNIFICANCE: fitlib/gof.py, shared with
# mcCoreFit.py

# normalisation is extrapolated from the sideband, so give it room
NORM_RANGE = (0.2, 5.0)            # x the extrapolated value

# Data gets a FILLED circle. sigFit uses an open one (24) for simulation,
# so the two are distinguishable at a glance: open = MC, filled = data.
# The SIZE still comes from sigFit, so the points match in scale.
DATA_MARKER_STYLE = 20

# ---------------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------------
# y-title offset for the TALL upper pad; the ratio-pad offset is too small
# there and puts the title on top of the tick labels
MAIN_Y_TITLE_OFFSET = 1.45


# ---------------------------------------------------------------------------
# input
# ---------------------------------------------------------------------------

def load_data(category, binMVA, year, nbins):
    """Full-range data histogram, from prepareFits."""
    from fitlib.prepareFits import getHisto
    h = getHisto(nbins, config.XLOW, config.XHIGH, False, category, year,
                 False, binMVA or "incl")
    h.SetDirectory(0)
    return h


def setup_ranges(x):
    """Name the sideband pieces, the SR and the full range; return the
    sideband union as a range list."""
    x.setRange(R_FULL, config.XLOW, config.XHIGH)
    x.setRange(R_SB_LOW, config.XLOW, SR_LO)
    x.setRange(R_SB_HIGH, SR_HI, config.XHIGH)
    x.setRange(R_SR, SR_LO, SR_HI)
    # RooFit takes a comma-separated list of named ranges, so the union is
    # expressed at use time rather than as its own range
    return f"{R_SB_LOW},{R_SB_HIGH}"


def sideband_bins(h):
    """Bin indices inside the sideband, for counting and for chi2."""
    out = []
    for i in range(1, h.GetNbinsX() + 1):
        c = h.GetBinCenter(i)
        if config.XLOW <= c < SR_LO or SR_HI <= c < config.XHIGH:
            out.append(i)
    return out


def full_yield(n_sb, fit):
    """Full-window yield extrapolated from the sideband count.

    Every PDF is evaluated normalised over the FULL window (the fit's
    Range() does not persist on x), so this one number is what the chi2,
    the plot and the workspace norm all use.
    """
    f_sb = fit["sb_fraction"]
    return n_sb / f_sb if f_sb > 0 else n_sb


# ---------------------------------------------------------------------------
# fitting
# ---------------------------------------------------------------------------

def fit_candidate(x, data, name, tag, sb_range, mc_core=None, seed=None):
    """Fit one candidate on the sideband only.

    seed: {parameter name: value}, a second start (seed_values: the fit of
    the order below, carried over). The fit is run from the default start
    and, if given, from the seed; the better result is kept (status 0
    first, then the lower NLL). `start` records which one it was.

    Returns a dict with the masked NLL (for the F-test), the sideband
    chi2 and its p-value (for the GOF cut), and the fraction of the PDF
    that lies in the sideband (for the normalisation extrapolation).
    """
    pdf, keep = pdefs.create_bkg_pdf(x, tag, name, mc_core)
    npar = pdefs.BKG_NPAR[name]
    pars = keep["params"]

    opts = [ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.Strategy(2),
            ROOT.RooFit.Save(True),
            # the masked likelihood: only the sideband contributes, and the
            # PDF is normalised over the sideband
            ROOT.RooFit.Range(sb_range),
            ROOT.RooFit.SumW2Error(True),
            ROOT.RooFit.PrintLevel(-1),
            ROOT.RooFit.PrintEvalErrors(-1)]
    def run():
        pdf.fitTo(data, *opts)
        return pdf.fitTo(data, *opts)

    def state():
        return {p.GetName(): (p.getVal(), p.getError()) for p in pars}

    def load(vals):
        for p in pars:
            if p.GetName() in vals:
                v, e = vals[p.GetName()]
                p.setVal(min(max(v, p.getMin()), p.getMax()))
                if e is not None:
                    p.setError(e)

    def rank(r):
        return (r.status() != 0, r.minNll())

    init = state()
    res, start = run(), "default"
    if seed:
        kept = state()
        load(init)
        load({k: (v, None) for k, v in seed.items()})
        res_seed = run()
        if rank(res_seed) < rank(res):
            res, start = res_seed, "seeded"
        else:
            load(kept)

    # sideband fraction: what part of the PDF lies in the sideband once
    # normalised over the full range. Needed to turn the observed sideband
    # count into a full-range yield.
    nset = ROOT.RooArgSet(x)
    itg_sb = pdf.createIntegral(nset, ROOT.RooFit.NormSet(nset),
                                ROOT.RooFit.Range(sb_range)).getVal()

    # fitted values by parameter name, for makeWorkspaces.py to rebuild from
    params = {p.GetName(): {"value": p.getVal(), "error": p.getError(),
                            "min": p.getMin(), "max": p.getMax()}
              for p in keep["params"]}

    return {"name": name, "pdf": pdf, "keep": keep, "npar": npar,
            "params": params, "start": start,
            "nll": res.minNll(), "status": res.status(),
            "cov_qual": res.covQual(), "edm": res.edm(),
            "sb_fraction": itg_sb, "result": res}


def pdf_histogram(x, h_data, pdf, norm_full, name):
    """Expectation per bin on h_data's binning, as a TH1 for the merge.

    Same convention as sideband_chi2: pdf.getVal(nset) is normalised over
    the FULL window, so the expectation is norm_full * density * binwidth.
    """
    h = h_data.Clone(name)
    h.SetDirectory(0)
    h.Reset()
    nset = ROOT.RooArgSet(x)
    saved = x.getVal()
    for i in range(1, h_data.GetNbinsX() + 1):
        x.setVal(h_data.GetBinCenter(i))
        h.SetBinContent(i, norm_full * pdf.getVal(nset)
                        * h_data.GetBinWidth(i))
    x.setVal(saved)
    return h


def sideband_chi2(x, h, pdf, sb_idx, norm_full, npar, threshold=None):
    """chi2 over SIDEBAND BINS ONLY, with low-significance bins merged.

    A goodness-of-fit that included bins the fit never saw would not be a
    goodness-of-fit. pdf.getVal(nset) is normalised over the FULL window,
    so the expectation in a bin is norm_full * density * binwidth, with
    norm_full = full_yield(n_sb, fit) = n_sb / f_sb.

    The merging is what makes the p-value trustworthy in the sparse tail;
    set threshold = 0 to recover the unmerged, bin-by-bin chi2.
    """
    threshold = (CHI2_MIN_BIN_SIGNIFICANCE if threshold is None
                 else threshold)
    h_pdf = pdf_histogram(x, h, pdf, norm_full, f"hpdf_{pdf.GetName()}")
    chi2, nbins = merge_chi2_bins(h, h_pdf, SR_LO, SR_HI, threshold)
    ndf = nbins - npar
    if ndf <= 0:
        print(f"  warning: {pdf.GetName()} left {nbins} merged chi2 bins "
              f"for {npar} parameters -- no degrees of freedom")
    p = ROOT.TMath.Prob(chi2, ndf) if ndf > 0 else 0.0
    return {"chi2": chi2, "ndf": ndf,
            "chi2_ndf": chi2 / ndf if ndf > 0 else -1.0,
            "pvalue": p, "n_bins": nbins,
            "n_bins_raw": len(sb_idx),
            "chi2_threshold": threshold}


# Bernstein-type families: the prefix of their coefficient names. c0 is
# fixed to 1 and has no parameter; c1..cn are <prefix><k>_<tag>.
BERN_COEF = {"bernstein": "bern_c", "bwzbern": "bwzbern_c",
             "fewzbern": "fewzbern_c", "landbern": "landbern_c"}
# sums of exponentials / power laws: the fraction that keeps the components
# of the order below (<name><tag>, no underscore), per member
SUM_FRAC = {"exp2": "exp_c1", "exp3": "exp_c2",
            "pow2": "frac1", "pow3": "frac2"}
SEED_FRAC = 0.99


def lower_order():
    """{member: the order below it in its family}, from BKG_FAMILIES."""
    return {hi: (fam, lo) for fam, mem in pdefs.BKG_FAMILIES.items()
            for lo, hi in zip(mem, mem[1:])}


def seed_values(name, family, prev, tag):
    """A start for `name` at the fit of `prev`, the order below it.

    Parameters the two share keep prev's fitted values. Then, per family:
      Bernstein types   degree elevation, c'_k = k/(n+1) c_(k-1)
                        + (1 - k/(n+1)) c_k: the same polynomial, so the
                        start reproduces prev exactly;
      exp / pow sums    the components of prev keep SEED_FRAC of the
                        weight, the new one starts at its default: close to
                        prev, but not at the stationary point an exact copy
                        (new exponent = old) would sit on;
      dy* / mc*         the new coefficient at its default 0 is exact.
    """
    vals = {n: p["value"] for n, p in prev["params"].items()}
    if family in BERN_COEF:
        pre = BERN_COEF[family]
        n = 0
        while f"{pre}{n + 1}_{tag}" in vals:
            n += 1
        b = [1.0] + [vals[f"{pre}{k}_{tag}"] for k in range(1, n + 1)]
        for k in range(1, n + 2):
            r = k / (n + 1)
            vals[f"{pre}{k}_{tag}"] = (r * b[k - 1]
                                       + (1 - r) * (b[k] if k <= n else 0.0))
    elif name in SUM_FRAC:
        vals[f"{SUM_FRAC[name]}{tag}"] = SEED_FRAC
    return vals


def ftest(fits, family, failed=None):
    """Lowest order in `family` whose successor is not justified.

    2*(NLL_low - NLL_high) is compared with chi2 at the difference in
    parameter count. Both NLLs come from the same masked likelihood, so
    the comparison is valid; mixing a blinded and an unblinded fit here
    would not be.

    failed: the members whose fit failed. The chain starts at the lowest
    order that fitted and stops at the first failed order after it: the
    order below is kept, rather than comparing across the gap with an
    order whose predecessor was never tested.
    """
    failed = failed or {}
    names = [n for n in family if n in fits or n in failed]
    while names and names[0] not in fits:
        names.pop(0)
    if not names:
        return None, []
    steps = []
    chosen = names[0]
    for lo, hi in zip(names, names[1:]):
        if hi in failed:
            steps.append({"from": lo, "to": hi, "failed": True,
                          "status": failed[hi]["status"]})
            break
        dnll = fits[lo]["nll"] - fits[hi]["nll"]
        dnpar = fits[hi]["npar"] - fits[lo]["npar"]
        if dnpar <= 0:
            continue
        # a higher order cannot fit worse; guard against a failed fit
        chi2 = max(2.0 * dnll, 0.0)
        p = ROOT.TMath.Prob(chi2, dnpar)
        steps.append({"from": lo, "to": hi, "2dNLL": chi2,
                      "dnpar": dnpar, "pvalue": p,
                      "justified": p < F_THRESHOLD})
        if p < F_THRESHOLD:
            chosen = hi
        else:
            break
    return chosen, steps


# ---------------------------------------------------------------------------
# plots
# ---------------------------------------------------------------------------

def plot_envelope(x, data, fits, envelope, best, label, outbase, year,
                  sb_range, save_pdf, nbins, n_sb, labels=None):
    """Data with the SR suppressed, every envelope member drawn across the
    full range, and a ratio panel against the best-fit model.

    Two RooFit details matter here:
      * data.plotOn(..., CutRange(sideband)) ZEROES the SR bins but still
        draws them, so the SR points are removed from the RooHist
        afterwards. That removal is the blinding on the plot.
      * Each curve gets its yield explicitly, Normalization(n_sb / f_sb,
        NumEvent), AND NormRange(full). Without an explicit NormRange,
        RooFit sees the "fitrange" attribute fitTo(Range=sideband) left on
        x and normalises the curve to the sideband instead, putting it
        1/f_sb too high. The explicit yield is the same extrapolation the
        chi2 and the workspace norm use, so all three agree.

    The legend names each member as the Run 2 note does, with its order
    (pdfDefinitions.model_label); `labels` overrides that per member.
    """
    labels = labels or {}
    frame = x.frame(ROOT.RooFit.Title(""))
    data.plotOn(frame, ROOT.RooFit.Name("dat"),
                ROOT.RooFit.CutRange(sb_range),
                ROOT.RooFit.MarkerStyle(DATA_MARKER_STYLE),
                ROOT.RooFit.MarkerSize(style.DATA_MARKER_SIZE),
                ROOT.RooFit.MarkerColor(ROOT.kBlack),
                ROOT.RooFit.LineColor(ROOT.kBlack),
                ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))

    # blinding: drop the zeroed SR points CutRange leaves behind
    hdat = frame.getHist("dat")
    for i in reversed(range(hdat.GetN())):
        if SR_LO <= hdat.GetPointX(i) < SR_HI:
            hdat.RemovePoint(i)

    cols = [style.fit_color(), ROOT.kAzure + 2, ROOT.kGreen + 2,
            ROOT.kOrange + 7, ROOT.kViolet + 1, ROOT.kTeal + 3]
    for i, nm in enumerate(envelope):
        fits[nm]["pdf"].plotOn(
            frame, ROOT.RooFit.Name(nm),
            ROOT.RooFit.Range(R_FULL),
            # explicit, or RooFit silently normalises to the last fit
            # range (the sideband) and the curve comes out 1/f_sb too high
            ROOT.RooFit.NormRange(R_FULL),
            ROOT.RooFit.Normalization(full_yield(n_sb, fits[nm]),
                                      ROOT.RooAbsReal.NumEvent),
            ROOT.RooFit.LineColor(cols[i % len(cols)]),
            ROOT.RooFit.LineStyle(ROOT.kSolid if nm == best
                                  else ROOT.kDashed),
            ROOT.RooFit.LineWidth(2))

    bw = (config.XHIGH - config.XLOW) / nbins
    frame.SetTitle(f";;Events / {bw:.2f} GeV")
    style._style_ratio_axis(frame.GetYaxis(), style.RATIO_Y_TITLE_OFFSET)
    frame.GetYaxis().SetTitleOffset(MAIN_Y_TITLE_OFFSET)
    frame.GetXaxis().SetLabelSize(0)
    frame.GetXaxis().SetTitleSize(0)
    frame.SetMinimum(0.0)
    frame.SetMaximum(1.35 * frame.GetMaximum())

    # ratio of data and of each member to the best-fit model
    cbest = frame.getCurve(best)
    rgraphs = []
    g = ROOT.TGraphErrors()
    k = 0
    for i in range(hdat.GetN()):
        xv, yv = hdat.GetPointX(i), hdat.GetPointY(i)
        ey = hdat.GetErrorY(i)
        f = cbest.average(xv - 0.5 * bw, xv + 0.5 * bw)
        if yv <= 0 or f <= 0:
            continue
        g.SetPoint(k, xv, yv / f)
        g.SetPointError(k, 0.0, ey / f)
        k += 1
    g.SetMarkerStyle(DATA_MARKER_STYLE)
    g.SetMarkerSize(style.DATA_MARKER_SIZE)
    g.SetMarkerColor(ROOT.kBlack)
    g.SetLineColor(ROOT.kBlack)
    rgraphs.append(g)

    for i, nm in enumerate(envelope):
        if nm == best:
            continue
        c = frame.getCurve(nm)
        gr = ROOT.TGraph()
        npt = 0
        xv = config.XLOW + 0.5 * bw
        while xv < config.XHIGH:
            fb = cbest.average(xv - 0.5 * bw, xv + 0.5 * bw)
            fn = c.average(xv - 0.5 * bw, xv + 0.5 * bw)
            if fb > 0:
                gr.SetPoint(npt, xv, fn / fb)
                npt += 1
            xv += bw
        gr.SetLineColor(cols[(i) % len(cols)])
        gr.SetLineStyle(ROOT.kDashed)
        gr.SetLineWidth(2)
        rgraphs.append(gr)

    rframe = ROOT.TH1F("rframe", "", 1, config.XLOW, config.XHIGH)
    rframe.SetDirectory(0)
    rframe.SetStats(0)
    rframe.SetMinimum(0.80)
    rframe.SetMaximum(1.20)
    rframe.SetTitle(";m_{#mu#mu} [GeV];Data/Template")
    style._style_ratio_axis(rframe.GetXaxis(), style.RATIO_X_TITLE_OFFSET)
    style._style_ratio_axis(rframe.GetYaxis(), style.RATIO_Y_TITLE_OFFSET)

    canvas, outer, pad1, pad2 = style.make_canvas_pads("c_bkg")
    pad1.cd()
    frame.Draw()

    # shade the blinded region so it cannot be mistaken for empty data
    box = ROOT.TBox(SR_LO, 0.0, SR_HI, frame.GetMaximum())
    box.SetFillColorAlpha(ROOT.kGray, 0.35)
    box.SetLineWidth(0)
    box.Draw("same")
    frame.Draw("same")

    keep = [box]
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(style.ANN_TEXT_SIZE)
    dy = style.LINE_SPACING * style.ANN_TEXT_SIZE
    for i, ln in enumerate(style.header_lines(label) +
                           [f"blinded {SR_LO:g}-{SR_HI:g} GeV"]):
        latex.DrawLatex(style.PAD_LEFT_MARGIN + 0.05, 0.85 - i * dy, ln)
    keep.append(latex)

    # wide enough for the note-style names: from just right of the text block
    # on the left to the frame's edge
    leg = ROOT.TLegend(style.PAD_LEFT_MARGIN + 0.31,
                       1 - style.PAD_TOP_MARGIN - 0.04
                       - 0.048 * (len(envelope) + 1),
                       1 - style.PAD_RIGHT_MARGIN - 0.01,
                       1 - style.PAD_TOP_MARGIN - 0.04)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.026)
    leg.SetMargin(0.12)
    leg.AddEntry(frame.getHist("dat"), "Data (sideband)", "pe")
    for nm in envelope:
        leg.AddEntry(frame.getCurve(nm),
                     labels.get(nm, pdefs.model_label(nm))
                     + ("  (best)" if nm == best else ""), "l")
    leg.Draw()
    keep.append(leg)
    keep.append(style.cms_label(pad1, year, full_year=True))

    pad2.cd()
    rframe.Draw()
    for gr in rgraphs[1:]:
        gr.Draw("L SAME")
    rgraphs[0].Draw("P SAME")
    one = ROOT.TLine(config.XLOW, 1.0, config.XHIGH, 1.0)
    one.SetLineColor(style.fit_color())
    one.SetLineWidth(2)
    one.Draw("same")
    rbox = ROOT.TBox(SR_LO, 0.80, SR_HI, 1.20)
    rbox.SetFillColorAlpha(ROOT.kGray, 0.35)
    rbox.SetLineWidth(0)
    rbox.Draw("same")
    rframe.Draw("axis same")
    keep += [one, rbox] + rgraphs

    os.makedirs(os.path.dirname(outbase), exist_ok=True)
    canvas.Update()
    canvas.SaveAs(f"{outbase}_envelope.png")
    if save_pdf:
        canvas.SaveAs(f"{outbase}_envelope.pdf")
    del keep, pad1, pad2, outer, canvas


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--category", default="ggHcat",
                    choices=sorted(config.MC_LIST))
    ap.add_argument("-b", "--bin", default="incl",
                    help="BDT bin label; pass \"\" to reproduce the "
                         "double-underscore tag the existing datacards use")
    ap.add_argument("-y", "--year", default="Run3")
    ap.add_argument("--extra-models", action="store_true",
                    help="add the families beyond the Run 2 note "
                         "(EXTRA_FAMILIES in pdfDefinitions): dyexp, "
                         "dypow, and mcexp/mcpow/mcdy if mcCoreFit.py "
                         "has run")
    ap.add_argument("--candidates", nargs="+", default=None,
                    help="fit exactly these members instead (overrides "
                         "--extra-models)")
    ap.add_argument("--wsdir", default=None,
                    help="fit JSONs and workspaces (default WS_<year>)")
    ap.add_argument("--mc-core", default=None,
                    help="frozen MC cores from mcCoreFit.py (default "
                         "<wsdir>/MCcore_<tag>.json, if it exists)")
    ap.add_argument("--no-mc-core", action="store_true",
                    help="with --extra-models: leave out the "
                         "mcexp/mcpow/mcdy families")
    ap.add_argument("--plotdir", default=None,
                    help="default ~/public_html/HmumuFits/bkg_fits")
    ap.add_argument("--pdf", action="store_true")
    ap.add_argument("--bins-per-gev", type=float, default=None)
    args = ap.parse_args()
    args.wsdir = args.wsdir or config.default_wsdir(args.year)
    args.plotdir = args.plotdir or os.path.expanduser(
        "~/public_html/HmumuFits/bkg_fits")
    if args.candidates:
        model_set = "candidates"
    elif args.extra_models:
        model_set = "run2+extra"
    else:
        model_set = "run2"
    # the MC cores only matter when mc* members can be asked for
    use_mc = model_set != "run2" and not args.no_mc_core

    style.setup_style()
    if args.bins_per_gev:
        config.BINS_PER_GEV = args.bins_per_gev
    nbins = int((config.XHIGH - config.XLOW) * config.BINS_PER_GEV)

    # exactly bwsHrare's tag, double underscore and all
    tag = f"{args.category}_{args.bin}_{args.year}"
    label = f"{args.category} {args.bin or 'incl'} data"

    print(f"category : {args.category}"
          f"\nbin      : {args.bin or '(empty)'}"
          f"\nyear     : {args.year}"
          f"\ntag      : {tag}"
          f"\nwindow   : {config.XLOW:g}-{config.XHIGH:g} GeV, {nbins} bins"
          f"\nblinded  : {SR_LO:g}-{SR_HI:g} GeV"
          f"\nsideband : {config.XLOW:g}-{SR_LO:g} + {SR_HI:g}-{config.XHIGH:g} GeV")

    # ---- data ----------------------------------------------------------
    h = load_data(args.category, args.bin, args.year, nbins)
    x = ROOT.RooRealVar(f"mh{args.category}", "m_{#mu#mu}",
                        config.XLOW, config.XHIGH)
    sb_range = setup_ranges(x)
    data = ROOT.RooDataHist(f"datahist_{tag}", "data",
                            ROOT.RooArgList(x), h)

    sb_idx = sideband_bins(h)
    n_sb = sum(h.GetBinContent(i) for i in sb_idx)
    n_full = h.Integral(1, h.GetNbinsX())
    print(f"\nevents   : {n_full:.0f} in the window, {n_sb:.0f} in the "
          f"sideband ({100 * n_sb / n_full:.1f} %)")
    print(f"sideband bins used for chi2: {len(sb_idx)} of {nbins}")

    # ---- MC-tuned cores (mcCoreFit.py) ---------------------------------
    mc_file = args.mc_core or os.path.join(args.wsdir, f"MCcore_{tag}.json")
    mc_core = None
    if model_set == "run2":
        print("\nMC cores : not used (Run 2 families; --extra-models "
              "adds the MC-tuned cores)")
    elif not use_mc:
        print("\nMC cores : off (--no-mc-core)")
    elif os.path.isfile(mc_file):
        with open(mc_file) as fh:
            mc_core = json.load(fh)["cores"] or None
        print(f"\nMC cores : {', '.join(mc_core or []) or 'none'} "
              f"from {mc_file}")
    elif args.mc_core:
        raise SystemExit(f"--mc-core {args.mc_core}: no such file")
    else:
        print(f"\nMC cores : none ({mc_file} not found; run mcCoreFit.py "
              f"for the mcexp/mcpow/mcdy families)")

    # ---- fit every candidate -------------------------------------------
    def available(nm):
        return not mcCore.is_member(nm) or \
            nm.rstrip("0123456789") in (mc_core or {})
    if args.candidates:
        names = args.candidates
    elif args.extra_models:
        names = pdefs.family_members(pdefs.RUN2_FAMILIES
                                     + pdefs.EXTRA_FAMILIES)
    else:
        names = pdefs.family_members(pdefs.RUN2_FAMILIES)
    absent = [nm for nm in names if not available(nm)]
    if args.candidates and absent:
        raise SystemExit(f"no MC core for {absent}; run mcCoreFit.py")
    names = [nm for nm in names if available(nm)]
    print(f"\nfitting {len(names)} candidates on the sideband")
    print(f"{'model':<11}{'npar':>5}{'NLL':>14}{'chi2/ndf':>10}"
          f"{'p(chi2)':>9}{'f_sb':>8}{'status':>8}")
    fits, failed = {}, {}
    below = lower_order()
    for nm in names:
        # a second start at the fit of the order below, if that one fitted
        fam, prev = below.get(nm, (None, None))
        seed = (seed_values(nm, fam, fits[prev], tag) if prev in fits
                else None)
        try:
            f = fit_candidate(x, data, nm, tag, sb_range, mc_core, seed)
        except Exception as exc:
            print(f"{nm:<11}  failed: {exc}")
            failed[nm] = {"status": -1, "error": str(exc)}
            continue
        f.update(sideband_chi2(x, h, f["pdf"], sb_idx,
                               full_yield(n_sb, f), f["npar"]))
        flag = "" if f["status"] == 0 and f["cov_qual"] == 3 else "  <-- check"
        if f["start"] == "seeded":
            flag += "  (from the order below)"
        print(f"{nm:<11}{f['npar']:>5}{f['nll']:>14.2f}"
              f"{f['chi2_ndf']:>10.3f}{f['pvalue']:>9.3f}"
              f"{f['sb_fraction']:>8.3f}{f['status']:>8}{flag}")
        # a failed fit (status != 0) never enters the F-test or the envelope
        if f["status"] != 0:
            failed[nm] = f
        else:
            fits[nm] = f

    if failed:
        print("\nWARNING: failed fits, excluded from the envelope: "
              + ", ".join(f"{nm} (status {f['status']})"
                          for nm, f in failed.items()))
    if not fits:
        raise SystemExit("no candidate fitted")

    # ---- F-test within families ----------------------------------------
    print(f"\nF-test within families (p < {F_THRESHOLD} justifies the "
          f"extra parameter)")
    chosen, ftests = {}, {}
    for fam, members in pdefs.BKG_FAMILIES.items():
        pick, steps = ftest(fits, members, failed)
        if pick is None:
            continue
        chosen[fam] = pick
        ftests[fam] = steps
        for st in steps:
            if st.get("failed"):
                print(f"  {fam:<12} {st['from']:>6} -> {st['to']:<6} "
                      f"fit failed (status {st['status']}): stop here")
                continue
            print(f"  {fam:<12} {st['from']:>6} -> {st['to']:<6} "
                  f"2dNLL={st['2dNLL']:7.2f} dnpar={st['dnpar']} "
                  f"p={st['pvalue']:.4f}  "
                  f"{'justified' if st['justified'] else 'NOT justified'}")
        print(f"  {fam:<12} -> {pick}")

    # ---- GOF across families -------------------------------------------
    envelope = [nm for nm in chosen.values()
                if fits[nm]["pvalue"] > GOF_THRESHOLD]
    rejected = [nm for nm in chosen.values() if nm not in envelope]
    print(f"\nGOF cut (p > {GOF_THRESHOLD})")
    for nm in sorted(set(chosen.values())):
        print(f"  {nm:<11} p={fits[nm]['pvalue']:.4f}  "
              f"{'IN' if nm in envelope else 'rejected'}")
    if not envelope:
        raise SystemExit("no candidate passed the goodness-of-fit cut; "
                         "inspect the fits before lowering GOF_THRESHOLD")

    # ---- default member: lowest penalised NLL --------------------------
    for f in fits.values():
        f["nll_pen"] = f["nll"] + NLL_PENALTY_PER_PAR * f["npar"]
    envelope = sorted(set(envelope), key=lambda n: fits[n]["nll_pen"])
    print(f"\npenalised NLL = NLL + {NLL_PENALTY_PER_PAR:g} * npar")
    for nm in envelope:
        print(f"  {nm:<11} {fits[nm]['nll_pen']:14.2f}  "
              f"(+{fits[nm]['nll_pen'] - fits[envelope[0]]['nll_pen']:.2f})")
    # the member Combine would pick on the sideband: pdfindex 0, the source
    # of the Asimov and of the extrapolated norm, and the reference for the
    # ratio panel. It is NOT 'the' background model: the whole envelope goes
    # into the multipdf and Combine profiles over pdfindex.
    best = envelope[0]
    print(f"\nenvelope : {', '.join(envelope)}   (best fit: {best})")

    # ---- normalisation, extrapolated -----------------------------------
    # The fit only saw the sideband, so the full-range yield is the
    # observed sideband count divided by the fraction of the PDF that lies
    # in the sideband. NOT the observed total, which would defeat the
    # blinding, and not the sideband count, which would be too low.
    f_sb = fits[best]["sb_fraction"]
    norm_val = full_yield(n_sb, fits[best])
    spread = [full_yield(n_sb, fits[nm]) for nm in envelope
              if fits[nm]["sb_fraction"] > 0]
    print(f"\nnorm     : {norm_val:.1f} from {best} "
          f"(sideband {n_sb:.0f} / f_sb {f_sb:.4f})")
    if len(spread) > 1:
        print(f"           envelope spread {min(spread):.1f} - "
              f"{max(spread):.1f}  "
              f"({100 * (max(spread) - min(spread)) / norm_val:.1f} %)")
        print("           this spread is the extrapolation uncertainty; "
              "Combine\n           profiles it through the floating norm "
              "and pdfindex")

    # ---- plots ---------------------------------------------------------
    plotdir = os.path.join(args.plotdir, args.category)
    outbase = os.path.join(plotdir, f"bkg_{tag}")
    plot_envelope(x, data, fits, envelope, best, label, outbase,
                  args.year, sb_range, args.pdf, nbins, n_sb)

    # ---- data, for makeWorkspaces.py -----------------------------------
    # A plain TH1D, which every ROOT version reads; a RooDataHist written
    # here would not be guaranteed to read back in the Combine container.
    os.makedirs(args.wsdir, exist_ok=True)
    datafile = os.path.join(args.wsdir, f"Bkg_{tag}_data.root")
    fout = ROOT.TFile(datafile, "RECREATE")
    h.Write(f"data_{tag}")
    fout.Close()
    print(f"\nwrote {datafile}")

    # ---- json ----------------------------------------------------------
    rec = {
        "provenance": {
            "written_utc": datetime.now(timezone.utc)
                           .isoformat(timespec="seconds"),
            "git_commit": config.git_commit(),
            "category": args.category, "bin": args.bin,
            "year": args.year, "tag": tag,
            "window": [config.XLOW, config.XHIGH], "n_bins": nbins,
            "blind_region": [SR_LO, SR_HI],
            "f_threshold": F_THRESHOLD, "gof_threshold": GOF_THRESHOLD,
            "nll_penalty_per_par": NLL_PENALTY_PER_PAR,
            "model_set": model_set,
            "observable": x.GetName(), "data_file": os.path.basename(datafile),
        },
        "data": {"n_full": n_full, "n_sideband": n_sb,
                 "n_sideband_bins": len(sb_idx)},
        "candidates": {nm: {k: v for k, v in f.items()
                            if k not in ("pdf", "keep", "result")}
                       for nm, f in fits.items()},
        "failed_fits": {nm: f["status"] for nm, f in failed.items()},
        "ftest": ftests,
        "family_choice": chosen,
        "envelope": envelope,
        "rejected_by_gof": rejected,
        "best": best,
        # the frozen cores the mc* candidates were built on, so that
        # makeWorkspaces.py rebuilds them without the MC
        "mc_core": {"file": mc_file, "cores": mc_core} if mc_core else None,
        "normalization": {"value": norm_val, "sideband_count": n_sb,
                          "sb_fraction": f_sb, "norm_range": NORM_RANGE,
                          "envelope_spread": [min(spread), max(spread)]
                          if len(spread) > 1 else None},
    }
    jf = os.path.join(args.wsdir, f"Bkg_{tag}_params.json")
    with open(jf, "w") as fh:
        json.dump(rec, fh, indent=2, sort_keys=True)
    print(f"wrote {jf}")
    print(f"plots {outbase}_envelope.png")
    bin_arg = args.bin or '""'
    sig = os.path.join(args.wsdir, f"Signal_{tag}_params.json")
    if not os.path.isfile(sig):
        print(f"\nnote: {sig} is missing; buildWorkspace.py needs it "
              f"(sigFit.py --wsdir {args.wsdir})")
    print(f"\nnext, inside the Combine container:\n  python buildWorkspace.py "
          f"-c {args.category} -b {bin_arg} -y {args.year} "
          f"--wsdir {args.wsdir}")


if __name__ == "__main__":
    main()
