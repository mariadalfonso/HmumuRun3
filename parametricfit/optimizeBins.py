"""
optimizeBins.py -- the BDT bin boundaries of a category that maximise the
expected significance, as the Run 2 note chose its ggH categories.

    python optimizeBins.py -c ggHcat -y Run3
    python optimizeBins.py -c ggHcat -y Run3 --bkg data --min-bkg 2000

The Run 2 note (AN-2019/124, Sec. 6.2, "Category optimization"):
  * the figure of merit is the expected significance of a combined fit of
    the categories, signal (all modes) and background modelled over
    110-150 GeV, the background from simulation
  * the boundary is scanned in signal-efficiency quantiles, 1% steps
  * categories are added one at a time until the next gains < 1%
  * odd event numbers of the signal MC for the optimisation, even ones for
    the final signal model
It found 5 ggH categories, at 70, 35, 20 and 5% signal efficiency.

Here, per candidate bin (a range of discrMVA):
  * s_k: the signal MC in 0.25 GeV mass bins over 110-150 GeV, every mode
    (mc 10-15), with sigFit.py's per-era yield scaling
  * b_k: the background over 110-150 GeV, a smooth fit evaluated in the
    same bins. --bkg mc (default, as the note): the background MC of
    mcCoreFit.py (the same era scaling; --bkg-samples: the Drell-Yan ones,
    prepareFits.BKG_DY_TAGS, by default in ggHcat, where DY is the
    background; every sample, BKG_MC_TAGS, in the other categories), a
    Poisson fit of exp(cubic in m) to its 1 GeV histogram (the note used
    BWZRedux; any smooth fit serves here). --bkg data: the data SIDEBANDS
    only (110-120, 130-150 GeV), exp(quadratic in m), extrapolated under
    the peak; the 120-130 GeV data are never read.
  * Z^2 = sum_k 2((s_k+b_k) ln(1+s_k/b_k) - s_k), the Asimov significance
    of the binned shape fit with the shapes fixed: what the note's
    combined fit gives without shape uncertainties. The total is
    sqrt(sum over bins of Z^2).

The search is exact: dynamic programming over the grid finds the best split
for each number of bins (the note's greedy search adds one boundary at a
time and can only do as well). --grid quantile (default) puts the grid at
1% steps of signal efficiency, as the note; --grid score at --step in
discrMVA. The recommended number of bins is the last one before a bin that
gains less than --min-gain (default 1%, the note's rule).

Every bin needs --min-bkg background events over 110-150 GeV (the
background MC's weighted sum, or the sideband data count with --bkg data;
default 1000) and, with --bkg mc, --min-mc-raw unweighted MC events
(default 200), so its background fit is stable. The purest bin tends to end
on that limit, so compare a few values.

--signal-parity odd uses odd event numbers only (weights x2), as the note;
the final signal model should then use the even ones, which sigFit.py does
not do yet, so the default is all events.

A guide to the boundaries, not the result: no background-shape uncertainty
(the envelope or core-PDF profiling) and no systematics. Check the best
candidates with the full chain (sigFit, bkgFit, the cards, combine -M
Significance).

Adapted from analysis/MACRO/optimize_categories.py, which counted signal and
data in |m - 125| < 25 GeV (the blinded region included), with s^2/b, at
least 20 background events per bin, and kept only the overall best.
Runs in conda, from parametricfit/. Writes <wsdir>/BinScan_<cat>_<year>.json
and, in <plotdir>/<cat>/ (--plot-nbins picks the split, default the
recommended one; --no-plots for none):
    binscan_<cat>_<year>_score[_log].png   signal and background score at
        unit area, the split's boundaries (signal efficiency above each)
        and the current SELMVA ones
    binscan_<cat>_<year>_gain.png          Z against the number of bins
    binscan_<cat>_<year>_fits_<k>bins.png  per bin, the background histogram,
        the fit the significance uses (chi2/ndf) and the signal scaled up
"""

import argparse
import json
import math
import os
import re

import numpy as np
import ROOT

from fitlib import config, prepareFits as pf, style

ROOT.gROOT.SetBatch(True)

MASS = "HiggsCandCorrMass"
SCORE = "discrMVA"
SB_LO = (config.XLOW, 120.0)          # the data sidebands, --bkg data
SB_HI = (130.0, config.XHIGH)
FIT_BIN = 1.0                          # GeV, the background fit's histogram
# --bkg-samples per category: ggH's background is Drell-Yan
DEFAULT_BKG_SAMPLES = {"ggHcat": "dy"}
Z_BIN = 0.25                           # GeV, the significance's mass bins
N_FINE = 1000                          # score histogram bins over [0, 1]
POLY = {"mc": 3, "data": 2}            # degree of ln(b) in m


# ---------------------------------------------------------------------------
# input: one pass per file group, (score x mass) histograms
# ---------------------------------------------------------------------------

def hist2d_array(h):
    """TH2 contents as an array [score bin, mass bin], the score's under-
    and overflow folded into the first and last bin."""
    nx, ny = h.GetNbinsX(), h.GetNbinsY()
    a = np.zeros((nx, ny))
    for i in range(nx + 2):
        ii = min(max(i - 1, 0), nx - 1)
        for j in range(1, ny + 1):
            a[ii, j - 1] += h.GetBinContent(i, j)
    return a


def stale_files(files, branches=(SCORE, MASS)):
    """the files whose events tree lacks one of `branches`: snapshots made
    before the BDT score was added, say"""
    stale = []
    for path in files:
        f = ROOT.TFile.Open(path)
        t = f.Get("events") if f and not f.IsZombie() else None
        if not t or any(not t.GetBranch(b) for b in branches):
            stale.append(path)
        if f:
            f.Close()
    return stale


def require_branches(files, what):
    """stop with the list of files that cannot be read for `what`"""
    stale = stale_files(files)
    if stale:
        raise SystemExit(
            f"{len(stale)} of {len(files)} {what} snapshot(s) have no "
            f"{SCORE} (or {MASS}) branch; remake them with the current "
            "analysis code:\n  " + "\n  ".join(stale)
            + ("\n(--bkg-samples dy uses the Drell-Yan samples only)"
               if what.startswith("background") else ""))


def fill(files, filters, nmass, name, weight=True):
    df = ROOT.RDataFrame("events", files)
    for f in filters:
        df = df.Filter(f)
    df = df.Define("wgt", "w_allSF" if weight else "1.0")
    h = df.Histo2D((name, "", N_FINE, 0.0, 1.0, nmass, config.XLOW,
                    config.XHIGH), SCORE, MASS, "wgt")
    return hist2d_array(h.GetValue())


def parity_filter(parity):
    return {"all": [], "odd": ["event % 2 == 1"],
            "even": ["event % 2 == 0"]}[parity]


def load_signal(args):
    """[fine score bin, Z_BIN mass bin], summed over modes and MC eras."""
    nz = int(round((config.XHIGH - config.XLOW) / Z_BIN))
    tags = sorted({t for s in args.modes for t in pf.SIGNAL_TAGS[s]})
    scale_parity = 1.0 if args.signal_parity == "all" else 2.0
    total = np.zeros((N_FINE, nz))
    for era, scale in pf.signal_yield_scales(args.year).items():
        files = pf._collect_files([(era, tags)], args.category,
                                  args.rootfiles_dir)
        require_branches(files, f"signal {era}")
        total += scale * scale_parity * fill(
            files, [f"mc >= 10 && mc <= 15 && !isnan({MASS})"]
            + parity_filter(args.signal_parity), nz, f"s_{era}")
        print(f"  signal {era}: {len(files)} file(s), x {scale:.3f}"
              + ("" if scale_parity == 1 else
                 f", {args.signal_parity} events x 2"))
    return total


def load_background(args):
    """[fine score bin, FIT_BIN mass bin]: background MC over 110-150 GeV
    (weighted, and unweighted for the raw count), or the data sidebands."""
    nm = int(round((config.XHIGH - config.XLOW) / FIT_BIN))
    if args.bkg == "data":
        files = pf.get_data_files(args.year, args.category,
                                  args.rootfiles_dir)
        require_branches(files, "data")
        sb = (f"({MASS} >= {SB_LO[0]} && {MASS} < {SB_LO[1]}) || "
              f"({MASS} >= {SB_HI[0]} && {MASS} < {SB_HI[1]})")
        a = fill(files, [f"mc < 0 && mc > -100 && !isnan({MASS})", sb], nm,
                 "b_data", weight=False)
        print(f"  data: {len(files)} file(s), {a.sum():.0f} sideband events "
              "(120-130 GeV not read)")
        return a, a
    eras = pf.signal_yield_scales(args.year)
    dy_only = args.bkg_samples == "dy"
    per_era = {era: pf.get_bkg_mc_files(args.category, era,
                                        args.rootfiles_dir, dy_only=dy_only)
               for era in eras}
    # check every era first, so one run lists all the stale files
    require_branches([f for fs in per_era.values() for f in fs],
                     "background MC")
    w = np.zeros((N_FINE, nm))
    raw = np.zeros((N_FINE, nm))
    for era, scale in eras.items():
        files = per_era[era]
        cut = [f"mc > 15 && !isnan({MASS})"]
        w += scale * fill(files, cut, nm, f"b_{era}")
        raw += fill(files, cut, nm, f"braw_{era}", weight=False)
        print(f"  background MC {era}: {len(files)} file(s), x {scale:.3f}")
    print(f"  background MC: {w.sum():.0f} expected, {raw.sum():.0f} raw "
          "events over 110-150 GeV")
    return w, raw


# ---------------------------------------------------------------------------
# the significance of one candidate bin
# ---------------------------------------------------------------------------

def background_density(counts, mode):
    """Poisson fit of exp(poly(u)), u = m - 130, to the FIT_BIN histogram:
    over 110-150 (mc) or the sidebands only (data). Returns the fitted
    density per GeV as a function, or None."""
    centres = config.XLOW + FIT_BIN * (np.arange(len(counts)) + 0.5)
    use = np.ones(len(counts), bool) if mode == "mc" else \
        (centres < SB_LO[1]) | (centres >= SB_HI[0])
    u, n = centres[use] - 130.0, counts[use]
    deg = POLY[mode]
    if (n > 0).sum() < deg + 2:
        return None
    X = np.vstack([u ** p for p in range(deg + 1)]).T
    pos = n > 0
    w0 = np.sqrt(n[pos])
    beta, *_ = np.linalg.lstsq(X[pos] * w0[:, None], np.log(n[pos]) * w0,
                               rcond=None)
    # Poisson regression (IRLS); quasi-Poisson for the weighted MC
    for _ in range(30):
        mu = np.exp(np.clip(X @ beta, -50, 50))
        z = X @ beta + (n - mu) / mu
        sw = np.sqrt(mu)
        new, *_ = np.linalg.lstsq(X * sw[:, None], z * sw, rcond=None)
        done = np.max(np.abs(new - beta)) < 1e-10
        beta = new
        if done:
            break

    def dens(m):
        uu = np.asarray(m) - 130.0
        return np.exp(sum(beta[p] * uu ** p for p in range(deg + 1))) / FIT_BIN
    return dens


def bin_stats(sig, bkg, raw, args):
    """(s, b over 110-150, Z^2, s in 120-130, b in 120-130) of one bin, or
    None if it fails the minimum statistics."""
    nb = bkg.sum()
    if nb < args.min_bkg or (args.bkg == "mc" and raw.sum() < args.min_mc_raw):
        return None
    dens = background_density(bkg, args.bkg)
    if dens is None:
        return None
    edges = config.XLOW + Z_BIN * np.arange(len(sig) + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    b = dens(centres) * Z_BIN
    s = sig
    good = (b > 0) & (s > 0)
    z2 = float(np.sum(2.0 * ((s[good] + b[good]) * np.log1p(s[good] / b[good])
                             - s[good])))
    win = (centres >= 120.0) & (centres < 130.0)
    return {"s": float(s.sum()), "b_full": float(b.sum()), "Z2": z2,
            "s_window": float(s[win].sum()), "b_window": float(b[win].sum()),
            "n_bkg": float(nb), "n_raw": float(raw.sum())}


# ---------------------------------------------------------------------------
# the grid and the search
# ---------------------------------------------------------------------------

def make_grid(sig, args):
    """Cell edges as fine score-bin indices, from the top down: 1% steps in
    signal efficiency (quantile) or --step in score."""
    if args.grid == "score":
        step = int(round(args.step * N_FINE))
        return list(range(0, N_FINE + 1, step))
    per_bin = sig.sum(axis=1)
    above = np.cumsum(per_bin[::-1])[::-1]          # signal at or above bin i
    total = above[0]
    edges = {0, N_FINE}
    for q in np.arange(args.quantile_step, 1.0 - 1e-9, args.quantile_step):
        # the lowest fine edge with at most a fraction q of the signal above
        i = int(np.searchsorted(-above, -q * total, side="left"))
        if 0 < i < N_FINE:
            edges.add(i)
    return sorted(edges)


def best_partitions(grid, stats, max_bins):
    """Exact best split of the grid into k = 1..max_bins ranges maximising
    the sum of Z^2. stats(i, j) gives a range's statistics or None.
    Returns {k: (Z, [grid edge indices])}."""
    n = len(grid) - 1
    table = {}
    for i in range(n):
        for j in range(i + 1, n + 1):
            r = stats(grid[i], grid[j])
            if r is not None:
                table[i, j] = r
    best = [[-1.0] * (n + 1) for _ in range(max_bins + 1)]
    prev = [[None] * (n + 1) for _ in range(max_bins + 1)]
    best[0][0] = 0.0
    for k in range(1, max_bins + 1):
        for j in range(1, n + 1):
            for i in range(j):
                if best[k - 1][i] < 0 or (i, j) not in table:
                    continue
                v = best[k - 1][i] + table[i, j]["Z2"]
                if v > best[k][j]:
                    best[k][j], prev[k][j] = v, i
    out = {}
    for k in range(1, max_bins + 1):
        if best[k][n] < 0:
            continue
        idx, j = [n], n
        for kk in range(k, 0, -1):
            j = prev[kk][j]
            idx.append(j)
        out[k] = (math.sqrt(best[k][n]), [grid[i] for i in idx[::-1]])
    return out


def current_edges(category):
    """The cut values of fitlib/prepareFits.SELMVA for the category."""
    vals = set()
    for b, cut in pf.SELMVA.get(category, {}).items():
        if re.fullmatch(r"bdt\d+", b):
            vals.update(float(v) for v in re.findall(r"[0-9]*\.[0-9]+", cut))
    return sorted(vals)


def describe(fine_edges, stats, sig_total_above):
    rows = []
    for lo, hi in zip(fine_edges[:-1], fine_edges[1:]):
        r = dict(stats(lo, hi))
        r.update(lo=lo / N_FINE, hi=hi / N_FINE,
                 eff_lo=float(sig_total_above(hi)),
                 eff_hi=float(sig_total_above(lo)), Z=math.sqrt(r["Z2"]))
        rows.append(r)
    return rows


def print_rows(rows):
    print(f"    {'score range':<16}{'signal eff.':>14}{'s':>9}{'b (110-150)':>13}"
          f"{'s/sqrt(b) 120-130':>19}{'Z':>8}")
    for r in rows:
        swb = r["s_window"] / math.sqrt(r["b_window"]) if r["b_window"] else 0
        eff = f"{100 * r['eff_hi']:.0f}-{100 * r['eff_lo']:.0f}%"
        print(f"    [{r['lo']:.3f}, {r['hi']:.3f}){'':<2}{eff:>14}{r['s']:>9.1f}"
              f"{r['b_full']:>13.0f}{swb:>19.3f}{r['Z']:>8.3f}")


def selmva(edges):
    """Score edges (fractions of 1) as a prepareFits.SELMVA entry."""
    inner = edges[1:-1]
    cuts = {}
    for k in range(len(inner) + 1):
        parts = ([f"discrMVA>={inner[k - 1]:g}"] if k > 0 else []) + \
                ([f"discrMVA<{inner[k]:g}"] if k < len(inner) else [])
        cuts[f"bdt{k}"] = " && ".join(parts)
    return cuts


# ---------------------------------------------------------------------------
# validation plots
# ---------------------------------------------------------------------------

SCORE_PLOT_BINS = 50                   # score bins of the distribution plot
SIG_COLOR = ROOT.kRed + 1
BKG_FILL, BKG_LINE = ROOT.kAzure - 9, ROOT.kAzure + 2
EDGE_COLOR = ROOT.kBlack
CUR_COLOR = ROOT.kGray + 2
FIT_SIG_COLOR = ROOT.kAzure + 2        # signal in the fit panels (fit is red)


def _latex(size, align=11, font=42):
    t = ROOT.TLatex()
    t.SetNDC()
    t.SetTextFont(font)
    t.SetTextSize(size)
    t.SetTextAlign(align)
    return t


def _save(c, stem, args):
    os.makedirs(os.path.dirname(stem), exist_ok=True)
    for ext in ("png", "pdf") if args.pdf else ("png",):
        c.SaveAs(f"{stem}.{ext}")
    print(f"  wrote {stem}.png")


def _bkg_label(args):
    if args.bkg == "data":
        return "data sidebands"
    return "background MC" + (" (DY)" if args.bkg_samples == "dy" else "")


def plot_score(sig, bkg, edges, effs, cur, args, stem):
    """Signal and background score, unit area, in 110-150 GeV (the data
    sidebands only with --bkg data), with the recommended boundaries (solid,
    the signal efficiency above each) and the current ones (dashed)."""
    reb = N_FINE // SCORE_PLOT_BINS
    hs = ROOT.TH1D("score_sig", "", SCORE_PLOT_BINS, 0.0, 1.0)
    hb = ROOT.TH1D("score_bkg", "", SCORE_PLOT_BINS, 0.0, 1.0)
    for h, a in ((hs, sig.sum(axis=1)), (hb, bkg.sum(axis=1))):
        for i, v in enumerate(a.reshape(-1, reb).sum(axis=1)):
            h.SetBinContent(i + 1, v)
        h.Scale(1.0 / max(h.Integral(), 1e-30))
        h.SetDirectory(0)
    hb.SetFillColor(BKG_FILL)
    hb.SetLineColor(BKG_LINE)
    hb.SetLineWidth(2)
    hs.SetLineColor(SIG_COLOR)
    hs.SetLineWidth(3)
    keep = [hs, hb]
    for logy in (True, False):
        c = ROOT.TCanvas(f"c_score_{logy}", "", 800, 650)
        c.SetMargin(style.PAD_LEFT_MARGIN, style.PAD_RIGHT_MARGIN,
                    style.PAD_BOTTOM_MARGIN, style.PAD_TOP_MARGIN)
        c.SetLogy(logy)
        ymax = max(hs.GetMaximum(), hb.GetMaximum())
        pos = [h.GetBinContent(b) for h in (hs, hb)
               for b in range(1, SCORE_PLOT_BINS + 1) if h.GetBinContent(b) > 0]
        ylo = 0.3 * min(pos) if logy else 0.0
        ytop = ylo * (ymax / ylo) ** 1.6 if logy else 1.55 * ymax
        frame = hb.Clone(f"frame_{logy}")
        frame.Reset()
        frame.SetMinimum(ylo)
        frame.SetMaximum(ytop)
        frame.GetXaxis().SetTitle("BDT score (discrMVA)")
        frame.GetYaxis().SetTitle(f"Fraction / {1.0 / SCORE_PLOT_BINS:g}")
        frame.GetYaxis().SetTitleOffset(1.3)
        frame.Draw("AXIS")
        hb.Draw("HIST SAME")
        hs.Draw("HIST SAME")
        # boundaries up to the curves' height, so the legend stays clear
        yline = ymax * (3 if logy else 1.08)
        for x in cur:
            ln = ROOT.TLine(x, ylo, x, yline)
            ln.SetLineColor(CUR_COLOR)
            ln.SetLineStyle(2)
            ln.SetLineWidth(2)
            ln.Draw()
            keep.append(ln)
        lab = _latex(0.026, 21)
        last_xn, lift = -1.0, 0.0
        for x, e in zip(edges[1:-1], effs[1:-1]):
            ln = ROOT.TLine(x, ylo, x, yline)
            ln.SetLineColor(EDGE_COLOR)
            ln.SetLineWidth(2)
            ln.Draw()
            keep.append(ln)
            # NDC position of the label: just above the line
            xn = c.GetLeftMargin() + (1 - c.GetLeftMargin()
                                      - c.GetRightMargin()) * x
            fy = (math.log(yline / ylo) / math.log(ytop / ylo) if logy
                  else yline / ytop)
            yn = c.GetBottomMargin() + (1 - c.GetBottomMargin()
                                        - c.GetTopMargin()) * fy + 0.012
            # neighbouring labels closer than their width: stack the next
            # one above, so the cut values never overlap
            lift = lift + 0.065 if xn - last_xn < 0.075 else 0.0
            last_xn = xn
            yn += lift
            lab.DrawLatex(xn, yn, f"{x:.3f}")
            lab.DrawLatex(xn, yn + 0.03, f"{100 * e:.0f}%")
        leg = ROOT.TLegend(0.19, 0.70, 0.62, 0.88)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.030)
        leg.AddEntry(hs, "signal (all modes)", "l")
        leg.AddEntry(hb, _bkg_label(args), "f")
        dummy = ROOT.TLine()
        dummy.SetLineWidth(2)
        leg.AddEntry(dummy, f"best {len(edges) - 1} bins (cut, sig. eff. above)",
                     "l")
        if cur:
            dcur = ROOT.TLine()
            dcur.SetLineColor(CUR_COLOR)
            dcur.SetLineStyle(2)
            dcur.SetLineWidth(2)
            leg.AddEntry(dcur, "current SELMVA", "l")
            keep.append(dcur)
        leg.Draw()
        keep += [frame, leg, dummy]
        c.RedrawAxis()
        style.cms_label(c, args.year, full_year=True)
        _latex(0.032, 13).DrawLatex(0.65, 0.87,
                                    style.CAT_LABEL.get(args.category,
                                                        args.category))
        _latex(0.028, 13).DrawLatex(0.65, 0.82, "110 < m_{#mu#mu} < 150 GeV")
        _save(c, stem + ("_log" if logy else ""), args)
    return keep


def plot_gain(results, recommended, cur_z, args, stem):
    """Z of the best split against the number of bins, the gain of each
    added bin, the current split's Z as a dashed line."""
    ks = sorted(results)
    g = ROOT.TGraph(len(ks))
    for i, k in enumerate(ks):
        g.SetPoint(i, k, results[k][0])
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.3)
    g.SetLineWidth(2)
    c = ROOT.TCanvas("c_gain", "", 800, 650)
    c.SetMargin(style.PAD_LEFT_MARGIN, style.PAD_RIGHT_MARGIN,
                style.PAD_BOTTOM_MARGIN, style.PAD_TOP_MARGIN)
    zs = [results[k][0] for k in ks] + ([cur_z] if cur_z else [])
    lo, hi = min(zs), max(zs)
    frame = c.DrawFrame(0.5, lo - 0.15 * (hi - lo), ks[-1] + 0.5,
                        hi + 0.35 * (hi - lo))
    frame.GetXaxis().SetTitle("number of BDT bins")
    frame.GetYaxis().SetTitle("expected significance Z (shape, no syst.)")
    frame.GetXaxis().SetNdivisions(len(ks) + 1, 0, 0)
    frame.GetYaxis().SetTitleOffset(1.3)
    keep = [g, frame]
    if cur_z:
        ln = ROOT.TLine(0.5, cur_z, ks[-1] + 0.5, cur_z)
        ln.SetLineColor(CUR_COLOR)
        ln.SetLineStyle(2)
        ln.SetLineWidth(2)
        ln.Draw()
        keep.append(ln)
        t = _latex(0.028, 11)
        t.SetNDC(False)
        t.SetTextColor(CUR_COLOR)
        t.DrawLatex(0.6, cur_z + 0.015 * (hi - lo),
                    f"current SELMVA: {cur_z:.3f}")
    g.Draw("PL")
    if recommended in results:
        r = ROOT.TGraph(1)
        r.SetPoint(0, recommended, results[recommended][0])
        r.SetMarkerStyle(24)
        r.SetMarkerSize(2.6)
        r.SetMarkerColor(SIG_COLOR)
        r.SetLineColor(SIG_COLOR)
        r.Draw("P")
        keep.append(r)
    t = _latex(0.026, 21)
    t.SetNDC(False)
    for k0, k1 in zip(ks[:-1], ks[1:]):
        z0, z1 = results[k0][0], results[k1][0]
        t.DrawLatex(k1, z1 + 0.06 * (hi - lo), f"+{100 * (z1 / z0 - 1):.1f}%")
    style.cms_label(c, args.year, full_year=True)
    _latex(0.032, 13).DrawLatex(0.19, 0.87,
                                style.CAT_LABEL.get(args.category,
                                                    args.category))
    _latex(0.028, 13).DrawLatex(
        0.19, 0.82, f"circled: recommended ({recommended} bins, next "
                    f"< {100 * args.min_gain:g}%)")
    _save(c, stem, args)
    return keep


def _poisson_chi2(n, err, mu):
    good = err > 0
    return float(np.sum(((n[good] - mu[good]) / err[good]) ** 2)), int(good.sum())


def plot_bin_fits(fine_edges, sig_c, bkg_c, raw_c, args, stem):
    """Per bin of the split: the background histogram (1 GeV) the fit sees,
    the fitted exp(poly) the significance uses, and the signal scaled up;
    chi2/ndf of the fit in each pad. With --bkg data the 120-130 GeV bins
    are not read and not drawn."""
    nbin = len(fine_edges) - 1
    ncol = min(nbin, 3)
    nrow = (nbin + ncol - 1) // ncol
    c = ROOT.TCanvas("c_fits", "", 460 * ncol, 420 * nrow)
    c.Divide(ncol, nrow, 0.002, 0.002)
    keep = []
    nm = bkg_c.shape[1]
    lo_m = config.XLOW
    for k, (a, b) in enumerate(zip(fine_edges[:-1], fine_edges[1:])):
        pad = c.cd(k + 1)
        pad.SetMargin(0.17, 0.04, 0.13, 0.10)
        n = bkg_c[b] - bkg_c[a]
        raw = raw_c[b] - raw_c[a]
        s = sig_c[b] - sig_c[a]
        dens = background_density(n, args.bkg)
        centres = lo_m + FIT_BIN * (np.arange(nm) + 0.5)
        blind = ((centres >= SB_LO[1]) & (centres < SB_HI[0])
                 if args.bkg == "data" else np.zeros(nm, bool))
        err = np.where(raw > 0, n / np.sqrt(np.maximum(raw, 1e-30)), 0.0)
        err[blind] = 0.0
        hb = ROOT.TH1D(f"fit_b{k}", "", nm, lo_m, config.XHIGH)
        for i in range(nm):
            if not blind[i]:
                hb.SetBinContent(i + 1, n[i])
                hb.SetBinError(i + 1, err[i])
        hb.SetMarkerStyle(style.DATA_MARKER_STYLE)
        hb.SetMarkerSize(style.DATA_MARKER_SIZE)
        mu = dens(centres) * FIT_BIN if dens else np.zeros(nm)
        hf = ROOT.TH1D(f"fit_f{k}", "", 400, lo_m, config.XHIGH)
        fx = lo_m + (config.XHIGH - lo_m) * (np.arange(400) + 0.5) / 400
        if dens:
            for i, v in enumerate(dens(fx) * FIT_BIN):
                hf.SetBinContent(i + 1, v)
        hf.SetLineColor(style.fit_color())
        hf.SetLineWidth(2)
        # the signal per 1 GeV, scaled to a round factor that makes its
        # peak about a fifth of the background under it
        s1 = s.reshape(nm, -1).sum(axis=1)
        peak_b = float(np.interp(125.0, centres, mu)) or 1.0
        raw_k = 0.2 * peak_b / max(s1.max(), 1e-30)
        mag = 10 ** math.floor(math.log10(max(raw_k, 1e-30)))
        kfac = max(1, int(round(raw_k / mag)) * mag) if raw_k >= 1 else 1
        hs = ROOT.TH1D(f"fit_s{k}", "", nm, lo_m, config.XHIGH)
        for i, v in enumerate(s1 * kfac):
            hs.SetBinContent(i + 1, v)
        hs.SetLineColor(FIT_SIG_COLOR)
        hs.SetFillColorAlpha(FIT_SIG_COLOR, 0.25)
        hs.SetLineWidth(2)
        top = max(hb.GetMaximum(), hf.GetMaximum(), hs.GetMaximum())
        frame = pad.DrawFrame(lo_m, 0.0, config.XHIGH, 1.45 * top)
        frame.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
        frame.GetYaxis().SetTitle(f"Events / {FIT_BIN:g} GeV")
        for ax in (frame.GetXaxis(), frame.GetYaxis()):
            ax.SetTitleSize(0.050)
            ax.SetLabelSize(0.045)
        frame.GetYaxis().SetTitleOffset(1.6)
        frame.GetYaxis().SetMaxDigits(3)
        hs.Draw("HIST SAME")
        hf.Draw("HIST L SAME")
        hb.Draw("E1 P SAME")
        chi2, npts = _poisson_chi2(n[~blind], err[~blind], mu[~blind])
        ndf = npts - (POLY[args.bkg] + 1)
        t = _latex(0.050, 13)
        t.DrawLatex(0.21, 0.86, f"bin {k}: [{a / N_FINE:.3f}, {b / N_FINE:.3f})")
        t2 = _latex(0.042, 13)
        t2.DrawLatex(0.21, 0.79, f"#chi^{{2}}/ndf = {chi2:.0f}/{ndf}")
        t2.DrawLatex(0.21, 0.73, f"b = {n.sum():.4g}, s = {s.sum():.1f}")
        leg = ROOT.TLegend(0.60, 0.70, 0.95, 0.88)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.040)
        leg.AddEntry(hb, "MC" if args.bkg == "mc" else "data", "pe")
        leg.AddEntry(hf, f"exp(poly{POLY[args.bkg]}) fit", "l")
        leg.AddEntry(hs, f"signal #times {kfac:g}", "f")
        leg.Draw()
        pad.RedrawAxis()
        keep += [hb, hf, hs, frame, t, t2, leg]
    c.cd(0)
    _save(c, stem, args)
    return keep


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--category", default="ggHcat",
                    choices=sorted(pf.SELMVA))
    ap.add_argument("-y", "--year", default="Run3")
    ap.add_argument("--bkg", choices=("mc", "data"), default="mc",
                    help="background from the background MC (default, as "
                         "the note) or the data sidebands")
    ap.add_argument("--bkg-samples", choices=("all", "dy"), default=None,
                    help="--bkg mc: every background MC sample "
                         "(BKG_MC_TAGS) or the Drell-Yan ones only "
                         "(BKG_DY_TAGS); default "
                         + ", ".join(f"{k} for {c}" for c, k in
                                     DEFAULT_BKG_SAMPLES.items())
                         + ", all otherwise")
    ap.add_argument("--grid", choices=("quantile", "score"),
                    default="quantile",
                    help="boundaries at signal-efficiency steps (default) "
                         "or at --step in discrMVA")
    ap.add_argument("--quantile-step", type=float, default=0.01,
                    help="signal-efficiency step of --grid quantile "
                         "(default 0.01, as the note)")
    ap.add_argument("--step", type=float, default=0.01,
                    help="discrMVA step of --grid score (default 0.01)")
    ap.add_argument("--max-bins", type=int, default=6)
    ap.add_argument("--min-gain", type=float, default=0.01,
                    help="recommend the last number of bins before one that "
                         "gains less than this (default 0.01, the note)")
    ap.add_argument("--min-bkg", type=float, default=1000,
                    help="minimum background events per bin over 110-150 "
                         "GeV (MC weighted sum, or sideband data count; "
                         "default 1000)")
    ap.add_argument("--min-mc-raw", type=float, default=200,
                    help="--bkg mc: minimum unweighted MC events per bin "
                         "(default 200)")
    ap.add_argument("--signal-parity", choices=("all", "odd", "even"),
                    default="all",
                    help="signal MC events used (odd: the note's choice, "
                         "weights x2; default all)")
    ap.add_argument("--modes", nargs="+", default=sorted(pf.SIGNAL_TAGS),
                    choices=sorted(pf.SIGNAL_TAGS),
                    help="signal modes (default all)")
    ap.add_argument("--rootfiles-dir", default=pf.ROOTFILES_DIR)
    ap.add_argument("--wsdir", default=None,
                    help="where the JSON goes (default WS_<year>)")
    ap.add_argument("--plotdir", default=None,
                    help="validation plots go to <plotdir>/<cat>/ "
                         "(default ~/public_html/HmumuFits/bin_scan)")
    ap.add_argument("--plot-nbins", type=int, default=None,
                    help="the split the plots show (default the "
                         "recommended one)")
    ap.add_argument("--no-plots", action="store_true")
    ap.add_argument("--pdf", action="store_true",
                    help="also write the plots as PDF")
    args = ap.parse_args()
    args.bkg_samples = args.bkg_samples or DEFAULT_BKG_SAMPLES.get(
        args.category, "all")
    args.wsdir = args.wsdir or config.default_wsdir(args.year)
    args.plotdir = args.plotdir or os.path.expanduser(
        "~/public_html/HmumuFits/bin_scan")
    if args.grid == "score":
        step = args.step * N_FINE
        if abs(step - round(step)) > 1e-9 or N_FINE % int(round(step)):
            ap.error(f"--step must be a multiple of {1 / N_FINE:g} dividing 1")

    print(f"category : {args.category}, {args.year}\n"
          f"background: {_bkg_label(args)}"
          f", fitted over 110-150 GeV\n"
          f"grid     : "
          + (f"{100 * args.quantile_step:g}% signal efficiency"
             if args.grid == "quantile" else f"{args.step:g} in discrMVA")
          + f", at least {args.min_bkg:g} background events per bin, up to "
          f"{args.max_bins} bins")
    sig = load_signal(args)
    bkg, raw = load_background(args)
    print(f"  signal over 110-150 GeV: {sig.sum():.1f}")

    sig_c = np.vstack([np.zeros(sig.shape[1]), np.cumsum(sig, axis=0)])
    bkg_c = np.vstack([np.zeros(bkg.shape[1]), np.cumsum(bkg, axis=0)])
    raw_c = np.vstack([np.zeros(raw.shape[1]), np.cumsum(raw, axis=0)])
    cache = {}

    def stats(lo, hi):
        if (lo, hi) not in cache:
            cache[lo, hi] = bin_stats(sig_c[hi] - sig_c[lo],
                                      bkg_c[hi] - bkg_c[lo],
                                      raw_c[hi] - raw_c[lo], args)
        return cache[lo, hi]

    per_bin = sig.sum(axis=1)
    above = np.concatenate([np.cumsum(per_bin[::-1])[::-1], [0.0]])

    def eff_above(i):
        return above[i] / above[0]

    grid = make_grid(sig, args)
    results = best_partitions(grid, stats, args.max_bins)

    out = {"category": args.category, "year": args.year, "bkg": args.bkg, "bkg_samples": args.bkg_samples,
           "grid": args.grid, "min_bkg": args.min_bkg,
           "signal_parity": args.signal_parity, "modes": args.modes,
           "best": {}}
    cur = current_edges(args.category)
    cur_z = None
    if cur:
        fe = [0] + [int(round(v * N_FINE)) for v in cur] + [N_FINE]
        if all(stats(a, b) is not None for a, b in zip(fe[:-1], fe[1:])):
            rows = describe(fe, stats, eff_above)
            z = cur_z = math.sqrt(sum(r["Z2"] for r in rows))
            print(f"\ncurrent  : edges {', '.join(f'{v:g}' for v in cur)} "
                  f"(prepareFits.SELMVA)   Z = {z:.3f}")
            print_rows(rows)
            out["current"] = {"edges": [0.0] + cur + [1.0], "Z": z,
                              "bins": rows}
        else:
            print(f"\ncurrent  : edges {cur}: a bin fails the minimum "
                  "statistics")

    prevz, recommended = None, None
    for k, (z, fe) in sorted(results.items()):
        edges = [i / N_FINE for i in fe]
        rows = describe(fe, stats, eff_above)
        gain = None if prevz is None else z / prevz - 1
        if gain is not None and gain < args.min_gain and recommended is None:
            recommended = k - 1
        print(f"\n{k} bin(s): Z = {z:.3f}"
              + ("" if gain is None else f"   (+{100 * gain:.1f}%)")
              + "   edges " + ", ".join(f"{e:g}" for e in edges[1:-1])
              + "   signal eff. "
              + ", ".join(f"{100 * eff_above(i):.0f}%" for i in fe[1:-1]))
        print_rows(rows)
        out["best"][k] = {"edges": edges, "Z": z, "bins": rows,
                          "signal_eff_edges": [eff_above(i) for i in fe],
                          "selmva": selmva(edges)}
        prevz = z
    if results:
        if recommended is None:
            recommended = max(results)
            print(f"\nevery added bin gains >= {100 * args.min_gain:g}% up to "
                  f"{recommended}; raise --max-bins to look further")
        out["recommended"] = recommended
        print(f"\nrecommended: {recommended} bin(s) (the next gains < "
              f"{100 * args.min_gain:g}%), as a prepareFits.SELMVA entry:\n  "
              + json.dumps(out["best"][recommended]["selmva"]))
    if results and not args.no_plots:
        k = args.plot_nbins or out["recommended"]
        if k not in results:
            print(f"\nno {k}-bin split to plot (have {sorted(results)})")
        else:
            fe = results[k][1]
            stem = os.path.join(args.plotdir, args.category,
                                f"binscan_{args.category}_{args.year}")
            print(f"\nvalidation plots, the {k}-bin split:")
            style.setup_style()
            plot_score(sig, bkg, [i / N_FINE for i in fe],
                       [eff_above(i) for i in fe], cur, args, stem + "_score")
            plot_gain(results, out["recommended"], cur_z, args,
                      stem + "_gain")
            plot_bin_fits(fe, sig_c, bkg_c, raw_c, args,
                          stem + f"_fits_{k}bins")
    os.makedirs(args.wsdir, exist_ok=True)
    path = os.path.join(args.wsdir,
                        f"BinScan_{args.category}_{args.year}.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
