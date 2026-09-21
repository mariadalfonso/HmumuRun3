"""
fitMVABins.py -- fit the signal mass shape in bins of the MVA score.

Splits signal MC into bins of `discrMVA` (default 10 bins of width 0.1),
fits a double-sided Crystal Ball in each, and summarises how mu, sigma, the
yield and the effective statistics evolve with MVA score.

Reuses sigFit.py for the PDF, the plotting and the diagnostics, so the
per-bin plots are identical in form to the standard signal fits.

Why this is worth doing
-----------------------
The MVA is trained on kinematics that correlate with mass resolution (muon
eta, pT, FSR activity), so sigma is NOT expected to be constant across
bins. That is the point of the scan: if sigma varies, a single shape for
the whole category is wrong, and each bin needs its own -- which is exactly
how the analysis is categorised downstream.

The catch: the highest MVA bins hold few events. Six free CB parameters
need a few thousand effective entries to be constrained (a closure scan at
N = 200 gave sigma anywhere from 1.58 to 2.32 GeV). Watch N_eff and the
`at bound` flags, and use --fix-tails when the thin bins misbehave.

Usage
-----
    python fitMVABins.py -c ggHcat -y 2024 -s ggH
    python fitMVABins.py -c ggHcat -y 2024 -s ggH --fix-tails
    python fitMVABins.py -c VBFcat -y Run3 -s qqH --edges 0 0.3 0.6 0.8 1.0
    python fitMVABins.py -c ggHcat -y 2024 -s ggH --mvavar discrMVA0
"""

import argparse
import glob
import inspect
import json
import os

import ROOT

import sigFit as sf

# sigFit.plot_fit gained an optional `result` argument when the fit-error
# band was added, and that patch may or may not be in the copy on disk.
# Check once rather than assuming either way.
_PLOT_FIT_TAKES_RESULT = (
    "result" in inspect.signature(sf.plot_fit).parameters)


# ---------------------------------------------------------------------------
# multi-line annotation header, without touching sigFit.py
# ---------------------------------------------------------------------------
# sigFit._annotate renders its header via sigFit.header_lines(label), which
# accepts only a string -- so an MVA bin would end up as one over-long line.
# Rather than edit sigFit.py, extend the function here: a list is taken as
# the lines themselves, anything else falls through to the original.
#
# This has to be installed before any plot_fit call, which it is: the
# assignment runs at import.
_sf_header_lines = sf.header_lines


def _header_lines(label):
    if isinstance(label, (list, tuple)):
        return ["H #rightarrow #mu#mu"] + list(label)
    return _sf_header_lines(label)


sf.header_lines = _header_lines


# Signal MC codes, matching prepareFits.signal_files
SIGNAL_MC = {
    "ggH": ["11"],
    "qqH": ["10"],
    "VH":  ["12", "13", "14"],
    "ttH": ["15"],
}

DIR_DEFAULT = "/work/submit/kbai/HmumuRun3/ROOTFILES"

TAIL_PARAMS = ("alphaL", "nL", "alphaR", "nR")

# Gap from the left edge of "CMS" to "Internal", in units of the CMS text
# size ("CMS" in font 61 is ~2.2 wide). Reproduces sigFit's 0.10 NDC on its
# 600x420 upper pad, and gives 0.144 on the 600x600 canvases here.
CMS_EXTRA_OFFSET = 2.4

# Pull overlay: +/-6 clipped the excursions near the peak, and a legend
# inside the frame sat on top of ten noisy series.
PULL_OVERLAY_MAX = 8.0
# The pull overlay carries ~400 points per slice across 40 GeV, so it gets a
# wider canvas than the square ones. Axis text is in absolute pixels (font
# 43) and ROOT scales NDC text by the shorter side, so nothing stretches --
# only the mass axis gets more room. The right margin can then be a smaller
# FRACTION while still leaving the legend the same number of pixels.
PULL_CANVAS_W = 1000
PULL_RIGHT_MARGIN = 0.16


# ---------------------------------------------------------------------------
# input
# ---------------------------------------------------------------------------

def collect_files(category, sig, year, basedir):
    files = []
    for tag in SIGNAL_MC[sig]:
        pat = (f"{basedir}/{category}/snapshot_mc_{tag}_*_{category}.root"
               if year == "Run3"
               else f"{basedir}/{category}/snapshot_mc_{tag}_{year}_{category}.root")
        files.extend(sorted(glob.glob(pat)))
    if not files:
        raise SystemExit(f"no snapshots for {category} {sig} {year} under {basedir}")
    return files


def mva_histo(files, lo, hi, mvavar, nbins, name, last=False):
    """Weighted mass histogram for one MVA slice.

    Bins are half-open [lo, hi) so they partition the range without
    double-counting a score sitting exactly on a boundary -- except the
    LAST bin, which is closed [lo, hi]. Without that, a score of exactly
    1.0 lands in no bin at all, and XGBoost's binary:logistic output can
    saturate to 1.0 in float32.
    """
    upper = "<=" if last else "<"
    df = ROOT.RDataFrame("events", files)
    df = (df.Filter("mc >= 10 && mc <= 15", "signal MC")
            .Filter("!isnan(HiggsCandCorrMass)", "valid mass")
            .Filter(f"{mvavar} >= {lo} && {mvavar} {upper} {hi}", "mva bin")
            .Define("weight", "w_allSF"))
    h = df.Histo1D((name, "", nbins, sf.XLOW, sf.XHIGH),
                   "HiggsCandCorrMass", "weight").GetValue()
    h.SetDirectory(0)
    return h


# ---------------------------------------------------------------------------
# the fit
# ---------------------------------------------------------------------------

def fit_bin(x, hist, tag, label, plotdir, year, save_pdf, fixed_tails=None):
    """Fit one MVA slice. Mirrors sigFit.fit_one, with the option to hold the
    tail parameters at values measured on the inclusive sample.

    Returns a dict, or None if the slice is empty.
    """
    norm = sf.range_integral(hist)
    n_eff = sf.effective_entries(hist)
    if norm <= 0:
        print(f"  !! empty slice for {tag}, skipping")
        return None, None

    data = ROOT.RooDataHist(f"datahist_{tag}", "data", ROOT.RooArgList(x), hist)
    pdf, bundle = sf.make_signal_pdf(x, tag)
    nom, nuis = bundle["nom"], bundle["nuis"]

    # nuisances are degenerate with cb_mu / cb_sigma: pin them for the fit
    for v in nuis.values():
        v.setVal(0.0)
        v.setConstant(True)

    if fixed_tails:
        for k in TAIL_PARAMS:
            nom[k].setVal(fixed_tails[k])
            nom[k].setConstant(True)

    fitopts = [
        ROOT.RooFit.Minimizer("Minuit2"),
        ROOT.RooFit.Strategy(2),
        ROOT.RooFit.Save(True),
        ROOT.RooFit.Range("full"),
        ROOT.RooFit.SumW2Error(True),
        ROOT.RooFit.PrintLevel(-1),
    ]
    pdf.fitTo(data, *fitopts)
    result = pdf.fitTo(data, *fitopts)

    n_float = sum(0 if v.isConstant() else 1 for v in nom.values())
    plot_kw = {"year": year, "save_pdf": save_pdf}
    if _PLOT_FIT_TAKES_RESULT:
        plot_kw["result"] = result
    chi2_ndf = sf.plot_fit(x, data, pdf, nom, hist, label,
                           f"{plotdir}/signal_{tag}", n_float, norm, n_eff,
                           **plot_kw)
    if isinstance(chi2_ndf, tuple):        # some versions also return the
        chi2_ndf = chi2_ndf[0]             # band summary

    status, cov = result.status(), result.covQual()
    parked = [k for k, v in nom.items() if not v.isConstant()
              and sf.near_bound(v)]
    ok = (status == 0 and cov == 3 and not parked)

    flag = "" if ok else "   <-- CHECK"
    print(f"  mu = {nom['mu'].getVal():8.3f} +/- {nom['mu'].getError():.3f}   "
          f"sigma = {nom['sigma'].getVal():6.3f} +/- {nom['sigma'].getError():.3f}   "
          f"yield = {norm:8.2f}   N_eff = {n_eff:9.0f}   "
          f"chi2/ndf = {chi2_ndf:5.2f}{flag}")
    if parked:
        print(f"    !! at bound: {', '.join(parked)}")

    record = {
        "params": {k: {"value": v.getVal(), "error": v.getError(),
                       "fixed": v.isConstant(),
                       "at_bound": (not v.isConstant()) and sf.near_bound(v)}
                   for k, v in nom.items()},
        "yield": norm,
        "n_eff": n_eff,
        "chi2_ndf": chi2_ndf,
        "status": status,
        "cov_qual": cov,
        "n_float": n_float,
        "params_at_bound": parked,
        "ok": ok,
    }
    # ROOT objects are kept separately: `record` is JSON-serialised, and
    # PyROOT would garbage-collect the pdf the moment this frame exits.
    objects = {"pdf": pdf, "bundle": bundle, "hist": hist,
               "norm": norm, "data": data, "result": result}
    return record, objects


# ---------------------------------------------------------------------------
# summary plots
# ---------------------------------------------------------------------------

def cms_label(pad, year):
    """sigFit.cms_label, with the "Internal" offset scaled to this pad.

    sigFit puts the extra text a fixed 0.10 NDC after "CMS", which suits its
    600x420 upper pad. ROOT scales text by a pad's SHORTER side while NDC x
    is a fraction of its width, so on the square canvases used here "CMS" is
    wider in NDC and the two overlap. Text, sizes and luminosity still come
    from sigFit. Same fix as fitEtaBins.py.
    """
    pad.cd()
    t = pad.GetTopMargin()
    l = pad.GetLeftMargin()
    r = pad.GetRightMargin()
    y = 1 - t + 0.2 * t
    cms_size = sf.CMS_TEXT_SIZE_FRAC * t

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(ROOT.kBlack)

    latex.SetTextFont(61)
    latex.SetTextAlign(11)
    latex.SetTextSize(cms_size)
    latex.DrawLatex(l, y, sf.CMS_TEXT)

    if sf.CMS_EXTRA_TEXT:
        w_px = pad.GetWw() * pad.GetAbsWNDC()
        h_px = pad.GetWh() * pad.GetAbsHNDC()
        dx = (CMS_EXTRA_OFFSET * cms_size * min(w_px, h_px) / w_px
              if w_px > 0 else 0.10)
        latex.SetTextFont(52)
        latex.SetTextSize(sf.EXTRA_OVER_CMS_TEXT_SIZE * cms_size)
        latex.DrawLatex(l + dx, y, sf.CMS_EXTRA_TEXT)

    lumi = sf.get_lumi(year) if year else 0.0
    if lumi <= 0:
        print(f"warning: no luminosity found for year tag '{year}'")
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.SetTextSize(sf.LUMI_TEXT_SIZE_FRAC * t)
    if lumi > 0:
        txt = "%.1f fb^{-1} (%.1f TeV)" % (lumi, sf.SQRT_S_TEV)
    else:
        txt = "(%.1f TeV)" % sf.SQRT_S_TEV
    latex.DrawLatex(1 - r, y, txt)
    return latex


def _single_canvas(name, right_margin=None, width=None):
    """`right_margin` widens the space to the right of the frame, so a
    legend can sit OUTSIDE the plot rather than on top of the data.
    `width` overrides the canvas width in pixels."""
    c = ROOT.TCanvas(name, "", sf.CANVAS_W if width is None else width,
                     sf.CANVAS_H)
    c.SetTopMargin(sf.PAD_TOP_MARGIN)
    c.SetBottomMargin(sf.PAD_BOTTOM_MARGIN)
    c.SetLeftMargin(sf.PAD_LEFT_MARGIN)
    c.SetRightMargin(sf.PAD_RIGHT_MARGIN if right_margin is None
                     else right_margin)
    c.SetTickx()
    c.SetTicky()
    return c


def _trend_graph(centres, widths, vals, errs):
    from array import array
    n = len(centres)
    g = ROOT.TGraphErrors(n, array('d', centres), array('d', vals),
                          array('d', widths), array('d', errs))
    g.SetMarkerStyle(sf.DATA_MARKER_STYLE)
    g.SetMarkerSize(sf.DATA_MARKER_SIZE)
    g.SetMarkerColor(ROOT.kBlack)
    g.SetLineColor(ROOT.kBlack)
    g.SetLineWidth(sf.DATA_LINE_WIDTH)
    return g


def _bin_colors(n):
    """n distinct colours from a perceptually ordered palette, so the legend
    reads as a progression in MVA score rather than an arbitrary set."""
    ROOT.gStyle.SetPalette(ROOT.kRainBow)
    ncol = ROOT.TColor.GetNumberOfColors()
    return [ROOT.TColor.GetColorPalette(int(i * (ncol - 1) / max(n - 1, 1)))
            for i in range(n)]


def overlay_models(x, objs, results, edges, label, outbase, year, save_pdf):
    """All fitted CB shapes on one canvas, each normalised to unit area.

    Unit-normalised on purpose: the question is whether the SHAPE changes
    with MVA score, and the yields differ by orders of magnitude across
    bins, so plotting them at their fitted scale would show nothing but the
    yield trend.
    """
    keys = [f"{edges[i]:.2f}_{edges[i + 1]:.2f}" for i in range(len(edges) - 1)]
    keys = [k for k in keys if objs.get(k)]
    if not keys:
        print("no fits to overlay")
        return

    colors = _bin_colors(len(keys))
    keep = []

    x.setRange("zoom", 115.0, 135.0)
    frame = x.frame(ROOT.RooFit.Range("zoom"), ROOT.RooFit.Title(""))

    leg = ROOT.TLegend(0.62, 0.42, 0.94, 0.86)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.028)
    leg.SetHeader("MVA score:  #sigma [GeV]")

    for i, k in enumerate(keys):
        objs[k]["pdf"].plotOn(frame,
                              ROOT.RooFit.LineColor(colors[i]),
                              ROOT.RooFit.LineWidth(2),
                              ROOT.RooFit.Range("zoom"),
                              ROOT.RooFit.NormRange("zoom"),
                              ROOT.RooFit.Name(k))
        lo, hi = k.split("_")
        sig_v = results[k]["params"]["sigma"]["value"]
        leg.AddEntry(frame.findObject(k),
                     f"[{lo}, {hi}]:  {sig_v:.3f}", "l")

    frame.SetTitle(";m_{#mu#mu} [GeV];a.u. (unit area)")
    sf._style_ratio_axis(frame.GetXaxis(), 1.15)
    sf._style_ratio_axis(frame.GetYaxis(), 1.60)

    c = _single_canvas("c_overlay_models")
    frame.Draw()
    leg.Draw()
    keep += [frame, leg, cms_label(c, year)]

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(sf.ANN_TEXT_SIZE)
    dy = sf.LINE_SPACING * sf.ANN_TEXT_SIZE
    for i, ln in enumerate(["H #rightarrow #mu#mu", label,
                            "fitted models, unit area"]):
        latex.DrawLatex(sf.PAD_LEFT_MARGIN + 0.05, 0.86 - i * dy, ln)
    keep.append(latex)

    c.Update()
    c.SaveAs(f"{outbase}_overlay_models.png")
    if save_pdf:
        c.SaveAs(f"{outbase}_overlay_models.pdf")
    del c


def _native_pulls(x, objs, k):
    """(centres, pulls) at native binning; undefined bins skipped.

    An empty data bin has zero error, so its pull is undefined -- skipped
    rather than drawn at zero, which would read as perfect agreement.
    """
    hist, pdf, norm = objs[k]["hist"], objs[k]["pdf"], objs[k]["norm"]
    b_lo, b_hi = sf.range_bins(hist)
    nset = ROOT.RooArgSet(x)
    saved = x.getVal()
    cs, ps = [], []
    for b in range(b_lo, b_hi + 1):
        err = hist.GetBinError(b)
        if err <= 0 or hist.GetBinContent(b) <= 0:
            continue
        x.setVal(hist.GetBinCenter(b))
        exp = norm * pdf.getVal(nset) * hist.GetBinWidth(b)
        cs.append(hist.GetBinCenter(b))
        ps.append((hist.GetBinContent(b) - exp) / err)
    x.setVal(saved)
    return cs, ps


def overlay_pulls(x, objs, edges, label, outbase, year, save_pdf,
                  pull_max=None, draw_lines=False, width=None):
    """Raw pull versus mass for every MVA bin, at NATIVE binning.

    No rebinning and no smoothing: every point is one 0.1 GeV bin, so the
    points are statistically independent and nothing is attenuated. Markers
    only by default -- connecting 400 noisy points per series with lines is
    what made this unreadable.

    Pulls come from the PDF evaluated per bin, not read off a RooCurve.
    """
    keys = [f"{edges[i]:.2f}_{edges[i + 1]:.2f}" for i in range(len(edges) - 1)]
    keys = [k for k in keys if objs.get(k)]
    if not keys:
        return

    if pull_max is None:
        pull_max = PULL_OVERLAY_MAX

    colors = _bin_colors(len(keys))
    keep = []

    graphs = []
    for i, k in enumerate(keys):
        cs, ps = _native_pulls(x, objs, k)
        if not ps:
            continue
        g = ROOT.TGraph(len(cs))
        for j in range(len(cs)):
            g.SetPoint(j, cs[j], ps[j])
        g.SetMarkerColor(colors[i])
        g.SetLineColor(colors[i])
        g.SetMarkerStyle(20)
        g.SetMarkerSize(0.4)
        g.SetLineWidth(1)
        graphs.append((k, g))

    bw = objs[keys[0]]["hist"].GetBinWidth(1)
    frame = ROOT.TH2F("ovpull_frame", "", 10, sf.XLOW, sf.XHIGH,
                      10, -pull_max, pull_max)
    frame.SetStats(0)
    frame.SetTitle(";m_{#mu#mu} [GeV];Pull")
    sf._style_ratio_axis(frame.GetXaxis(), 1.15)
    sf._style_ratio_axis(frame.GetYaxis(), 1.60)

    # Legend outside the frame: with ten series there is no clear patch of
    # canvas left inside it.
    leg = ROOT.TLegend(1.0 - PULL_RIGHT_MARGIN + 0.01, 0.30, 0.995, 0.88)
    leg.SetFillStyle(0)
    leg.SetBorderSize(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.026)
    leg.SetHeader("MVA score")

    c = _single_canvas("c_overlay_pulls", right_margin=PULL_RIGHT_MARGIN,
                       width=PULL_CANVAS_W if width is None else width)
    frame.Draw()

    # A pull has unit error BY CONSTRUCTION: (data - fit)/sigma_data, so its
    # own uncertainty is sigma/sigma = 1. Every point's error bar would be
    # +/-1; this band is that statement drawn once instead of 4000 times.
    band = ROOT.TBox(sf.XLOW, -1.0, sf.XHIGH, 1.0)
    band.SetFillColorAlpha(ROOT.kGray + 1, 0.35)
    band.SetLineWidth(0)
    band.Draw("same")
    keep.append(band)

    opt = "lp same" if draw_lines else "p same"
    for k, g in graphs:
        g.Draw(opt)
        lo, hi = k.split("_")
        leg.AddEntry(g, f"[{lo}, {hi}]", "p")

    for y in (0.0, 2.0, -2.0):
        ln = ROOT.TLine(sf.XLOW, y, sf.XHIGH, y)
        ln.SetLineColor(ROOT.kBlack if y == 0 else 11)
        ln.SetLineStyle(ROOT.kSolid if y == 0 else ROOT.kDashed)
        ln.SetLineWidth(2 if y == 0 else 1)
        ln.Draw("same")
        keep.append(ln)
    leg.AddEntry(band, "#pm1 (expected)", "f")

    leg.Draw()
    keep += [frame, leg] + [g for _, g in graphs] + [cms_label(c, year)]

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(sf.ANN_TEXT_SIZE)
    dy = sf.LINE_SPACING * sf.ANN_TEXT_SIZE
    for i, ln in enumerate(["H #rightarrow #mu#mu", label,
                            f"raw pull, {bw:.2f} GeV bins"]):
        latex.DrawLatex(sf.PAD_LEFT_MARGIN + 0.04, 0.86 - i * dy, ln)
    keep.append(latex)

    c.Update()
    c.SaveAs(f"{outbase}_overlay_pulls.png")
    if save_pdf:
        c.SaveAs(f"{outbase}_overlay_pulls.pdf")
    del c


def pull_map(x, objs, edges, label, outbase, year, save_pdf,
             mrange=(115.0, 140.0), zmax=4.0):
    """2D map of the pull: mass on x at NATIVE binning, MVA bin on y.

    Ten overlaid 400-point series are unreadable; a map shows all of them at
    full resolution with no overlap. Read it for orientation:
      - vertical stripes spanning all rows  -> mismatch common to every bin
      - patches confined to a few rows      -> MVA-dependent mismatch

    Restricted to `mrange` because outside it the thin high-MVA slices have
    empty bins, and an empty bin has an undefined pull that would render as
    a colour meaning "agreement".
    """
    keys = [f"{edges[i]:.2f}_{edges[i + 1]:.2f}" for i in range(len(edges) - 1)]
    keys = [k for k in keys if objs.get(k)]
    if not keys:
        return

    ref = objs[keys[0]]["hist"]
    b_lo = ref.FindBin(mrange[0] + 1e-6)
    b_hi = ref.FindBin(mrange[1] - 1e-6)
    nx = b_hi - b_lo + 1

    h2 = ROOT.TH2D("pullmap", "", nx,
                   ref.GetBinLowEdge(b_lo),
                   ref.GetBinLowEdge(b_hi) + ref.GetBinWidth(b_hi),
                   len(keys), 0.0, float(len(keys)))
    h2.SetDirectory(0)

    nset = ROOT.RooArgSet(x)
    saved_x = x.getVal()
    for row, k in enumerate(keys, start=1):
        hist, pdf, norm = objs[k]["hist"], objs[k]["pdf"], objs[k]["norm"]
        for j, b in enumerate(range(b_lo, b_hi + 1), start=1):
            err = hist.GetBinError(b)
            if err <= 0 or hist.GetBinContent(b) <= 0:
                continue
            x.setVal(hist.GetBinCenter(b))
            exp = norm * pdf.getVal(nset) * hist.GetBinWidth(b)
            h2.SetBinContent(j, row, (hist.GetBinContent(b) - exp) / err)
        lo, hi = k.split("_")
        h2.GetYaxis().SetBinLabel(row, f"{lo}-{hi}")
    x.setVal(saved_x)

    h2.SetStats(0)
    h2.SetTitle(";m_{#mu#mu} [GeV];MVA score;pull")
    h2.SetMinimum(-zmax)
    h2.SetMaximum(zmax)
    sf._style_ratio_axis(h2.GetXaxis(), 1.15)
    sf._style_ratio_axis(h2.GetYaxis(), 1.60)
    h2.GetZaxis().SetLabelFont(43)      # absolute px, as _style_ratio_axis
    h2.GetZaxis().SetLabelSize(16)

    # symmetric diverging palette: zero must read as neutral
    ROOT.gStyle.SetPalette(ROOT.kTemperatureMap)

    c = _single_canvas("c_pullmap", right_margin=0.16)
    h2.Draw("COLZ")
    keep = [h2, cms_label(c, year)]

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(sf.ANN_TEXT_SIZE)
    dy = sf.LINE_SPACING * sf.ANN_TEXT_SIZE
    for i, ln in enumerate(["H #rightarrow #mu#mu", label,
                            f"pull, {ref.GetBinWidth(1):.2f} GeV bins"]):
        latex.DrawLatex(sf.PAD_LEFT_MARGIN + 0.04, 0.86 - i * dy, ln)
    keep.append(latex)

    c.Update()
    c.SaveAs(f"{outbase}_pullmap.png")
    if save_pdf:
        c.SaveAs(f"{outbase}_pullmap.pdf")
    del c
    ROOT.gStyle.SetPalette(ROOT.kBird)      # restore


def oscillation_size(x, objs, edges, core=(120.0, 132.0)):
    """Quantify the mismatch per MVA bin, at native binning.

    Returns {key: {...}} with three numbers, all over the `core` window
    where every slice is populated:

    rms_pull      sqrt(<p^2>). A correct model gives 1.
    excess_pull   sqrt(<p^2> - 1), the part not explained by statistics.
                  NOT comparable across bins: a fixed fractional mismatch
                  scales as sqrt(N_eff), so a fat bin looks worse.
    frac_rms      sqrt(<(r-1)^2> - <sigma^2/f^2>), the statistics-subtracted
                  RMS of the Template/Fit deviation. Reads as "the fit is
                  wrong by X% RMS" and IS comparable across bins. This is
                  the one to use for "which MVA bin is described worst".
    """
    out = {}
    nset = ROOT.RooArgSet(x)
    saved_x = x.getVal()

    for i in range(len(edges) - 1):
        k = f"{edges[i]:.2f}_{edges[i + 1]:.2f}"
        if not objs.get(k):
            continue
        hist, pdf, norm = objs[k]["hist"], objs[k]["pdf"], objs[k]["norm"]
        b_lo = hist.FindBin(core[0] + 1e-6)
        b_hi = hist.FindBin(core[1] - 1e-6)

        sp2 = sr2 = sn2 = 0.0
        n = 0
        for b in range(b_lo, b_hi + 1):
            err = hist.GetBinError(b)
            c = hist.GetBinContent(b)
            if err <= 0 or c <= 0:
                continue
            x.setVal(hist.GetBinCenter(b))
            f = norm * pdf.getVal(nset) * hist.GetBinWidth(b)
            if f <= 0:
                continue
            sp2 += ((c - f) / err) ** 2
            sr2 += (c / f - 1.0) ** 2
            sn2 += (err / f) ** 2
            n += 1

        if n == 0:
            continue
        rms_pull = (sp2 / n) ** 0.5
        excess = max(sp2 / n - 1.0, 0.0) ** 0.5
        frac = max(sr2 / n - sn2 / n, 0.0) ** 0.5
        out[k] = {"n_bins": n, "rms_pull": rms_pull,
                  "excess_pull": excess, "frac_rms": frac,
                  "core": list(core)}

    x.setVal(saved_x)
    return out


def plot_oscillation(osc, edges, label, outbase, year, save_pdf):
    """frac_rms and excess_pull versus MVA score."""
    keys = [f"{edges[i]:.2f}_{edges[i + 1]:.2f}" for i in range(len(edges) - 1)]
    keys = [k for k in keys if k in osc]
    if not keys:
        return
    centres, widths = [], []
    for k in keys:
        lo, hi = (float(v) for v in k.split("_"))
        centres.append(0.5 * (lo + hi))
        widths.append(0.5 * (hi - lo))

    for name, ytitle, vals in (
            ("fracrms", "RMS of (Template/Fit #minus 1)",
             [osc[k]["frac_rms"] for k in keys]),
            ("excesspull", "#sqrt{#LTp^{2}#GT #minus 1}",
             [osc[k]["excess_pull"] for k in keys])):
        c = _single_canvas(f"c_osc_{name}")
        g = _trend_graph(centres, widths, vals, [0.0] * len(keys))
        g.SetTitle(f";MVA score;{ytitle}")
        sf._style_ratio_axis(g.GetXaxis(), 1.15)
        sf._style_ratio_axis(g.GetYaxis(), 1.60)
        g.GetXaxis().SetLimits(edges[0], edges[-1])
        g.SetMinimum(0.0)
        g.Draw("ap")
        keep = [g, cms_label(c, year)]

        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextFont(42)
        latex.SetTextSize(sf.ANN_TEXT_SIZE)
        dy = sf.LINE_SPACING * sf.ANN_TEXT_SIZE
        note = ("statistics subtracted; comparable across bins"
                if name == "fracrms"
                else "grows as #sqrt{N_{eff}}; not comparable")
        for i, ln in enumerate(["H #rightarrow #mu#mu", label, note]):
            latex.DrawLatex(sf.PAD_LEFT_MARGIN + 0.05, 0.86 - i * dy, ln)
        keep.append(latex)

        c.Update()
        c.SaveAs(f"{outbase}_osc_{name}.png")
        if save_pdf:
            c.SaveAs(f"{outbase}_osc_{name}.pdf")
        del c


def plot_trends(results, edges, mvavar, label, outbase, year, save_pdf,
                inclusive=None):
    """sigma, mu, yield and N_eff versus MVA score."""
    keep = []
    centres, widths = [], []
    keys = []
    for i in range(len(edges) - 1):
        k = f"{edges[i]:.2f}_{edges[i + 1]:.2f}"
        if k not in results or results[k] is None:
            continue
        keys.append(k)
        centres.append(0.5 * (edges[i] + edges[i + 1]))
        widths.append(0.5 * (edges[i + 1] - edges[i]))

    if not keys:
        print("no successful fits to plot")
        return

    panels = [
        ("sigma", "#sigma [GeV]",
         [results[k]["params"]["sigma"]["value"] for k in keys],
         [results[k]["params"]["sigma"]["error"] for k in keys]),
        ("mu", "#mu [GeV]",
         [results[k]["params"]["mu"]["value"] for k in keys],
         [results[k]["params"]["mu"]["error"] for k in keys]),
        ("yield", "yield",
         [results[k]["yield"] for k in keys], [0.0] * len(keys)),
        ("neff", "N_{eff}",
         [results[k]["n_eff"] for k in keys], [0.0] * len(keys)),
        ("chi2", "#chi^{2}/ndf",
         [results[k]["chi2_ndf"] for k in keys], [0.0] * len(keys)),
    ]

    for name, ytitle, vals, errs in panels:
        c = _single_canvas(f"c_{name}")
        g = _trend_graph(centres, widths, vals, errs)
        g.SetTitle(f";{mvavar};{ytitle}")
        sf._style_ratio_axis(g.GetXaxis(), 1.15)
        sf._style_ratio_axis(g.GetYaxis(), 1.60)
        g.GetXaxis().SetLimits(edges[0], edges[-1])
        g.Draw("ap")
        keep.append(g)

        # inclusive reference, where it makes sense
        if inclusive and name in ("sigma", "mu"):
            ref = inclusive["params"][name]["value"]
            line = ROOT.TLine(edges[0], ref, edges[-1], ref)
            line.SetLineColor(sf.fit_color())
            line.SetLineStyle(ROOT.kDashed)
            line.SetLineWidth(2)
            line.Draw("same")
            keep.append(line)

        keep.append(cms_label(c, year))
        latex = ROOT.TLatex()
        latex.SetNDC()
        latex.SetTextFont(42)
        latex.SetTextSize(sf.ANN_TEXT_SIZE)
        dy = sf.LINE_SPACING * sf.ANN_TEXT_SIZE
        left = sf.PAD_LEFT_MARGIN + 0.05
        lines = ["H #rightarrow #mu#mu", label]
        if inclusive and name in ("sigma", "mu"):
            lines.append("#color[2]{dashed: inclusive fit}")
        for i, ln in enumerate(lines):
            latex.DrawLatex(left, 0.86 - i * dy, ln)
        keep.append(latex)

        c.Update()
        c.SaveAs(f"{outbase}_vs_mva_{name}.png")
        if save_pdf:
            c.SaveAs(f"{outbase}_vs_mva_{name}.pdf")
        del c


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--category", default="ggHcat",
                    choices=sorted(sf.MC_LIST))
    ap.add_argument("-s", "--sig", default=None,
                    help="production mode; default = the dominant one for "
                         "the category")
    ap.add_argument("-y", "--year", default="Run3")
    ap.add_argument("--mvavar", default="discrMVA",
                    help="MVA branch name (prepareFits.py uses discrMVA0 "
                         "for its own binning -- check which your snapshots "
                         "carry)")
    ap.add_argument("--edges", nargs="+", type=float, default=None,
                    help="bin edges; default 0, 0.1, ... 1.0")
    ap.add_argument("--fix-tails", action="store_true",
                    help="fit the inclusive sample first, then hold "
                         "alphaL/nL/alphaR/nR at those values in every bin, "
                         "floating only mu and sigma. Use when the thin "
                         "high-MVA bins will not converge.")
    ap.add_argument("--basedir", default=DIR_DEFAULT)
    ap.add_argument("--plotdir",
                    default=os.path.expanduser(
                        "~/public_html/HmumuFits/mva_bins"))
    ap.add_argument("--outdir", default="MVA_BINS",
                    help="JSON summary output directory")
    ap.add_argument("--pull-max", type=float, default=None,
                    help=f"y range of the pull overlay (default "
                         f"{PULL_OVERLAY_MAX:g})")
    ap.add_argument("--pull-width", type=int, default=PULL_CANVAS_W,
                    help=f"canvas width in px for the pull overlay "
                         f"(default {PULL_CANVAS_W}; the square plots stay "
                         f"at {{sf.CANVAS_W}})".replace(
                             "{sf.CANVAS_W}", str(sf.CANVAS_W)))
    ap.add_argument("--pull-lines", action="store_true",
                    help="connect the pull points with lines (off by "
                         "default: 400 noisy points per series is "
                         "unreadable joined up)")
    ap.add_argument("--map-range", nargs=2, type=float,
                    default=[115.0, 140.0],
                    help="mass window for the 2D pull map (outside it the "
                         "thin slices have empty, undefined-pull bins)")
    ap.add_argument("--core", nargs=2, type=float, default=[120.0, 132.0],
                    help="mass window used for the oscillation-size metrics")
    ap.add_argument("--pdf", action="store_true")
    ap.add_argument("--bins-per-gev", type=float, default=None)
    args = ap.parse_args()

    if args.bins_per_gev:
        sf.BINS_PER_GEV = args.bins_per_gev

    sig = args.sig or sf.MC_LIST[args.category][0]
    # round: 0.1*i gives 0.30000000000000004, which is ugly in the printed
    # edge list and in any downstream cut string
    edges = args.edges if args.edges else [round(0.1 * i, 10) for i in range(11)]
    edges = [round(e, 10) for e in edges]
    if len(edges) < 2:
        raise SystemExit("need at least two edges")

    nbins = int((sf.XHIGH - sf.XLOW) * sf.BINS_PER_GEV)
    files = collect_files(args.category, sig, args.year, args.basedir)
    print(f"{len(files)} file(s) for {args.category} {sig} {args.year}")
    print(f"MVA variable: {args.mvavar}")
    print(f"edges: {edges}")

    catdir = os.path.join(args.plotdir, args.category)
    os.makedirs(catdir, exist_ok=True)
    x = ROOT.RooRealVar(f"mh{args.category}", "m_{#mu#mu}", sf.XLOW, sf.XHIGH)
    x.setRange("full", sf.XLOW, sf.XHIGH)

    base_tag = f"{args.category}_{args.year}_{sig}"

    # ---- inclusive fit: reference, and the tail source for --fix-tails ----
    print(f"\n--- inclusive ({edges[0]} <= {args.mvavar} < {edges[-1]}) ---")
    h_incl = mva_histo(files, edges[0], edges[-1], args.mvavar, nbins,
                       f"h_{base_tag}_incl", last=True)
    cat_lab = sf.CAT_LABEL.get(args.category, args.category)
    sig_lab = sf.MODE_LABEL.get(sig, sig)

    objs = {}
    inclusive, obj_incl = fit_bin(x, h_incl, f"{base_tag}_mvaincl",
                        [f"{cat_lab} inclusive",
                         f"MVA score: all",
                         f"Signal: {sig_lab}"],
                        catdir, args.year, args.pdf)

    fixed_tails = None
    if args.fix_tails:
        if inclusive is None:
            raise SystemExit("--fix-tails needs a successful inclusive fit")
        fixed_tails = {k: inclusive["params"][k]["value"] for k in TAIL_PARAMS}
        print("\nholding tails at the inclusive values: " +
              "  ".join(f"{k}={v:.3f}" for k, v in fixed_tails.items()))

    # ---- per-bin fits ------------------------------------------------------
    results = {}
    for i in range(len(edges) - 1):
        lo, hi = edges[i], edges[i + 1]
        key = f"{lo:.2f}_{hi:.2f}"
        print(f"\n--- {args.mvavar} in [{lo:.2f}, {hi:.2f}) ---")
        h = mva_histo(files, lo, hi, args.mvavar, nbins,
                      f"h_{base_tag}_{key}", last=(i == len(edges) - 2))
        tag = f"{base_tag}_mva{lo:.2f}to{hi:.2f}".replace(".", "p")
        closing = "]" if i == len(edges) - 2 else ")"
        results[key], objs[key] = fit_bin(
            x, h, tag,
            [f"{cat_lab}",
             f"MVA score: [{lo:.2f}, {hi:.2f}{closing}",
             f"Signal: {sig_lab}"],
            catdir, args.year, args.pdf, fixed_tails)

    # ---- trends and summary ----------------------------------------------
    trend_label = f"{cat_lab}, Signal: {sig_lab} ({args.year})"
    overlay_models(x, objs, results, edges, trend_label,
                   os.path.join(catdir, base_tag), args.year, args.pdf)
    overlay_pulls(x, objs, edges, trend_label,
                  os.path.join(catdir, base_tag), args.year, args.pdf,
                  args.pull_max, args.pull_lines, args.pull_width)

    pull_map(x, objs, edges, trend_label,
             os.path.join(catdir, base_tag), args.year, args.pdf,
             tuple(args.map_range))
    osc = oscillation_size(x, objs, edges, tuple(args.core))
    plot_oscillation(osc, edges, trend_label,
                     os.path.join(catdir, base_tag), args.year, args.pdf)

    plot_trends(results, edges, args.mvavar,
                f"{args.category} {sig} ({args.year})",
                os.path.join(catdir, base_tag), args.year, args.pdf,
                inclusive)

    os.makedirs(args.outdir, exist_ok=True)
    outjson = os.path.join(args.outdir, f"mvabins_{base_tag}.json")
    with open(outjson, "w") as fh:
        json.dump({"category": args.category, "sig": sig, "year": args.year,
                   "mvavar": args.mvavar, "edges": edges,
                   "fix_tails": bool(args.fix_tails),
                   "fixed_tails": fixed_tails,
                   "inclusive": inclusive, "bins": results,
                   "oscillation": osc}, fh,
                  indent=2, sort_keys=True)

    # ---- table -----------------------------------------------------------
    print(f"\n{'bin':>14} {'yield':>10} {'N_eff':>10} {'mu':>9} "
          f"{'sigma':>8} {'chi2/ndf':>9} {'ok':>5}")
    tot = 0.0
    for i in range(len(edges) - 1):
        key = f"{edges[i]:.2f}_{edges[i + 1]:.2f}"
        r = results.get(key)
        if r is None:
            print(f"{key:>14} {'--':>10}")
            continue
        tot += r["yield"]
        print(f"{key:>14} {r['yield']:>10.2f} {r['n_eff']:>10.0f} "
              f"{r['params']['mu']['value']:>9.3f} "
              f"{r['params']['sigma']['value']:>8.3f} "
              f"{r['chi2_ndf']:>9.2f} {str(r['ok']):>5}")
    if inclusive:
        print(f"{'inclusive':>14} {inclusive['yield']:>10.2f} "
              f"{inclusive['n_eff']:>10.0f} "
              f"{inclusive['params']['mu']['value']:>9.3f} "
              f"{inclusive['params']['sigma']['value']:>8.3f} "
              f"{inclusive['chi2_ndf']:>9.2f} "
              f"{str(inclusive['ok']):>5}")
        print(f"\nsum of bin yields = {tot:.2f}, inclusive = "
              f"{inclusive['yield']:.2f}  "
              f"(closure {tot / inclusive['yield']:.4f})")

    if osc:
        print(f"\nmismatch size over m in [{args.core[0]:.0f}, "
              f"{args.core[1]:.0f}] GeV at native binning")
        print(f"{'bin':>14} {'RMS pull':>10} {'excess':>8} {'frac RMS':>10}")
        for i in range(len(edges) - 1):
            k = f"{edges[i]:.2f}_{edges[i + 1]:.2f}"
            if k not in osc:
                continue
            o = osc[k]
            print(f"{k:>14} {o['rms_pull']:>10.2f} {o['excess_pull']:>8.2f} "
                  f"{o['frac_rms']:>10.2%}")

    n_bad = sum(1 for r in results.values() if r and not r["ok"])
    if n_bad:
        print(f"\nWARNING: {n_bad} bin(s) flagged. If they are the high-MVA "
              f"thin ones, try --fix-tails.")
    print(f"\nwrote {outjson}")
    print(f"plots  {catdir}/{base_tag}_vs_mva_*.png")


if __name__ == "__main__":
    sf.setup_style()
    main()
