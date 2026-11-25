"""
dyCoreCompare.py -- compare Drell-Yan m_mumu shape predictions with the
FEWZ spline the fit uses now (fitlib/fewz.py, fewz_spline.root).

    key        label             source
    fewz_cur   FEWZ (current)    fitlib/fewz_spline.root, fewz_spline
    fewz_new   FEWZ (new)        fitlib/fewz_smooth.csv
    dyturbo    DYTurbo           fitlib/dyturbo_smooth.csv

All three sit in fitlib/ next to fewz.py (--csv-dir reads the two CSVs
from elsewhere). The CSVs (columns mass, sigma, 0.01 GeV apart) become
TSpline3s with knots every KNOT_STEP GeV, like the current file; that
reproduces the CSV to < 1e-5 in the fit window. The curves and fits use
only the SHAPE, each normalised over the window; DYcore_shapes.png also
labels each integral over the window relative to the current FEWZ (only
the current file states its units, pb/GeV) and the k of a fit of m^-k.

Plots, in --plotdir:

  DYcore_shapes.png          the three shapes, normalised to 1 over the
                             window, and their ratio to FEWZ (current);
                             needs no MC, so runs anywhere

and per BDT bin, with the DY MC (as mcCoreFit.py, unless --shapes-only):

  DYcore_raw_<tag>.png       MC with each shape normalised to it, not
                             fitted; below, shape / MC with the MC stat.
                             band. chi2 against the MC, ndf = bins - 1
  DYcore_bern<n>_<tag>.png   each shape x Bern_n fitted to the MC (as
                             fewzbern<n>, the order of --order), pulls
                             below, as the mc_core plots
  DYcore_<tag>.json          the chi2 / p values and fitted coefficients

The chi2 is fitlib/gof.py's, through mcCoreFit.gof_vs_mc: bins merged to
~9 effective MC entries each.

Usage, from parametricfit/ (conda, ROOT >= 6.28):

    python dyCoreCompare.py --shapes-only --plotdir ../../DY_core/plots
    python dyCoreCompare.py -c VBFcat -b bdt0 bdt1 bdt2 bdt3 -y Run3
    python dyCoreCompare.py -c ggHcat -b bdt0 bdt1 bdt2 -y Run3
"""

import argparse
import json
import os
import random
from datetime import datetime, timezone

import numpy as np
import ROOT

import mcCoreFit
from fitlib import config, fewz, style
from fitlib.bwz import BERN_INIT, BERN_RANGE, BERN_STEP, bernstein

ROOT.gROOT.SetBatch(True)

DEFAULT_CSV_DIR = os.path.dirname(fewz.FEWZ_FILE)     # fitlib/
CSV_FILES = {"fewz_new": "fewz_smooth.csv", "dyturbo": "dyturbo_smooth.csv"}
KNOT_STEP = 0.25               # GeV between spline knots, as fewz_spline

# key -> (legend label, colour); the order is the drawing order
PREDICTIONS = {
    "fewz_cur": ("FEWZ (current)", (0, 0, 0)),
    "fewz_new": ("FEWZ (new)", (42, 120, 214)),
    "dyturbo": ("DYTurbo", (235, 104, 52)),
}
REFERENCE = "fewz_cur"

SUBSTEPS = 20                  # midpoint rule per plot bin
SHAPE_STEP = 0.05              # GeV, for the drawn curves
RATIO_RANGE = (0.82, 1.12)
PRINT_MASSES = (110, 115, 120, 125, 130, 135, 140, 145, 150)


def colour(key):
    return ROOT.TColor.GetColor(*PREDICTIONS[key][1])


def label(key):
    return PREDICTIONS[key][0]


# --------------------------------------------------------------- the shapes

def csv_spline(path, name):
    """TSpline3 through the CSV rows on the KNOT_STEP grid."""
    d = np.loadtxt(path, delimiter=",", skiprows=1)
    m, y = d[:, 0], d[:, 1]
    on_grid = np.abs(m / KNOT_STEP - np.round(m / KNOT_STEP)) < 1e-6
    m, y = np.ascontiguousarray(m[on_grid]), np.ascontiguousarray(y[on_grid])
    if len(m) < 4 or np.any(y <= 0):
        raise RuntimeError(f"{path}: need >= 4 positive points on the "
                           f"{KNOT_STEP} GeV grid")
    return ROOT.TSpline3(name, m, y, len(m))


def load_predictions(csv_dir):
    """{key: TSpline3}, the current FEWZ first."""
    f = ROOT.TFile.Open(fewz.FEWZ_FILE)
    sp = f.Get(fewz.SPLINE).Clone(REFERENCE)
    f.Close()
    out = {REFERENCE: sp}
    for key, fname in CSV_FILES.items():
        out[key] = csv_spline(os.path.join(csv_dir, fname), key)
    for key, s in out.items():
        if not s.GetXmin() <= config.XLOW < config.XHIGH <= s.GetXmax():
            raise RuntimeError(f"{key} covers {s.GetXmin()}-{s.GetXmax()} "
                               f"GeV, not the window")
    return out


def density(sp, m):
    """sp normalised to 1 over the window, at the masses m."""
    step = (config.XHIGH - config.XLOW) / 8000
    grid = config.XLOW + (np.arange(8000) + 0.5) * step
    norm = step * sum(sp.Eval(v) for v in grid)
    return np.array([sp.Eval(v) for v in m]) / norm


def shape_stats(sp):
    """(integral over the window in the file's own units, k of a fit of
    m^-k, exponential slope in %/GeV); the fits are straight lines in
    ln(sigma) against ln(m) and m."""
    step = (config.XHIGH - config.XLOW) / 8000
    grid = config.XLOW + (np.arange(8000) + 0.5) * step
    y = np.array([sp.Eval(v) for v in grid])
    k = -np.polyfit(np.log(grid), np.log(y), 1)[0]
    slope = -100 * np.polyfit(grid, np.log(y), 1)[0]
    return step * y.sum(), k, slope


def shape_hist(sp, name, nbins, total):
    """sp integrated in the plot bins, normalised to `total`."""
    h = ROOT.TH1D(name, "", nbins, config.XLOW, config.XHIGH)
    h.SetDirectory(0)
    for i in range(1, nbins + 1):
        lo, w = h.GetBinLowEdge(i), h.GetBinWidth(i)
        step = w / SUBSTEPS
        h.SetBinContent(i, step * sum(sp.Eval(lo + (j + 0.5) * step)
                                      for j in range(SUBSTEPS)))
    h.Scale(total / h.Integral())
    return h


def header(pad, lumi_year):
    """CMS label; for the theory-only plot no luminosity."""
    if lumi_year:
        style.cms_label(pad, lumi_year, full_year=True)
        return
    t, l, r = pad.GetTopMargin(), pad.GetLeftMargin(), pad.GetRightMargin()
    y = 1 - t + 0.2 * t
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(61)
    lat.SetTextSize(style.CMS_TEXT_SIZE_FRAC * t)
    lat.DrawLatex(l, y, style.CMS_TEXT)
    lat.SetTextFont(52)
    lat.SetTextSize(style.EXTRA_OVER_CMS_TEXT_SIZE
                    * style.CMS_TEXT_SIZE_FRAC * t)
    lat.DrawLatex(l + 0.10, y, style.CMS_EXTRA_TEXT)
    lat.SetTextFont(42)
    lat.SetTextAlign(31)
    lat.SetTextSize(style.LUMI_TEXT_SIZE_FRAC * t)
    lat.DrawLatex(1 - r, y, "(%.1f TeV)" % style.SQRT_S_TEV)


def annotate(lines, size=style.ANN_TEXT_SIZE):
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextSize(size)
    dy = style.LINE_SPACING * size
    for i, line in enumerate(lines):
        lat.DrawLatex(style.PAD_LEFT_MARGIN + 0.05, 0.85 - i * dy, line)
    return lat


def ratio_guides(lo, hi, values=(1.0,)):
    out = []
    for v in values:
        ln = ROOT.TLine(lo, v, hi, v)
        ln.SetLineStyle(2)
        ln.SetLineColor(ROOT.kGray + 2)
        ln.Draw()
        out.append(ln)
    return out


def plot_shapes(splines, outbase, pdf_too):
    """The normalised shapes and their ratio to the current FEWZ."""
    m = np.arange(config.XLOW, config.XHIGH + 1e-9, SHAPE_STEP)
    dens = {k: density(sp, m) for k, sp in splines.items()}

    print(f"\nshape / {label(REFERENCE)}, each normalised over "
          f"{config.XLOW:g}-{config.XHIGH:g} GeV")
    print("  m [GeV]" + "".join(f"{label(k):>16}" for k in splines
                                if k != REFERENCE))
    for mass in PRINT_MASSES:
        i = int(round((mass - config.XLOW) / SHAPE_STEP))
        print(f"  {mass:7.1f}" + "".join(
            f"{dens[k][i] / dens[REFERENCE][i]:16.4f}"
            for k in splines if k != REFERENCE))

    stats = {k: shape_stats(sp) for k, sp in splines.items()}
    sig_ref = stats[REFERENCE][0]
    win = f"{config.XLOW:g}-{config.XHIGH:g}"
    print(f"\n{'':16}{'integral ' + win:>18}{'/ ' + label(REFERENCE):>18}"
          f"{'k (m^-k)':>10}{'slope %/GeV':>13}")
    for k, (sig, kpow, slope) in stats.items():
        print(f"{label(k):16}{sig:18.4g}{sig / sig_ref:18.4g}{kpow:10.2f}"
              f"{slope:13.2f}")

    c, outer, pad1, pad2 = style.make_canvas_pads("c_shapes")
    keep = []
    pad1.cd()
    mg = ROOT.TMultiGraph()
    leg = ROOT.TLegend(1 - style.PAD_RIGHT_MARGIN - 0.60, 0.42,
                       1 - style.PAD_RIGHT_MARGIN - 0.02, 0.70)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.030)
    for k in splines:
        g = ROOT.TGraph(len(m), m, np.ascontiguousarray(dens[k]))
        g.SetLineColor(colour(k))
        g.SetLineWidth(2)
        mg.Add(g, "L")
        sig, kpow, _ = stats[k]
        norm = (f"#sigma = {sig:.3g} pb" if k == REFERENCE else
                f"#sigma/#sigma_{{cur}} = {sig / sig_ref:.3g}")
        leg.AddEntry(g, f"{label(k)}:  {norm},  k = {kpow:.2f}", "l")
        keep.append(g)
    leg.AddEntry(ROOT.nullptr, f"#sigma: integral over {win} GeV", "")
    leg.AddEntry(ROOT.nullptr, "k: fit of m^{-k} over the same range", "")
    mg.SetTitle(";;1/#sigma d#sigma/dm_{#mu#mu} [GeV^{-1}]")
    mg.Draw("A")
    mg.GetXaxis().SetLimits(config.XLOW, config.XHIGH)
    mg.GetXaxis().SetLabelSize(0)
    mg.SetMinimum(0.0)
    mg.SetMaximum(1.3 * max(d.max() for d in dens.values()))
    style._style_ratio_axis(mg.GetYaxis(), style.RATIO_Y_TITLE_OFFSET)
    mg.GetYaxis().SetTitleOffset(1.45)
    leg.Draw()
    header(pad1, None)
    lat = annotate(["Drell-Yan m_{#mu#mu} shape",
                    f"normalised in {config.XLOW:g}-{config.XHIGH:g} GeV",
                    f"ratio to {label(REFERENCE)}"])

    pad2.cd()
    mg2 = ROOT.TMultiGraph()
    for k in splines:
        r = np.ascontiguousarray(dens[k] / dens[REFERENCE])
        g = ROOT.TGraph(len(m), m, r)
        g.SetLineColor(colour(k))
        g.SetLineWidth(2)
        mg2.Add(g, "L")
        keep.append(g)
    mg2.SetTitle(";m_{#mu#mu} [GeV];Ratio")
    mg2.Draw("A")
    mg2.GetXaxis().SetLimits(config.XLOW, config.XHIGH)
    mg2.SetMinimum(RATIO_RANGE[0])
    mg2.SetMaximum(RATIO_RANGE[1])
    mg2.GetYaxis().SetNdivisions(505)
    style._style_ratio_axis(mg2.GetXaxis(), style.RATIO_X_TITLE_OFFSET)
    style._style_ratio_axis(mg2.GetYaxis(), style.RATIO_Y_TITLE_OFFSET)
    guides = ratio_guides(config.XLOW, config.XHIGH, (0.9, 1.1))

    os.makedirs(os.path.dirname(outbase) or ".", exist_ok=True)
    c.SaveAs(outbase + ".png")
    if pdf_too:
        c.SaveAs(outbase + ".pdf")
    del keep, guides, lat


# ------------------------------------------------------------ against the MC

def make_pdfs(x, key, sp, order, tag):
    """(raw pdf, fitted pdf, its free coefficients) for one prediction.

    raw: the spline alone, no free parameter. fitted: spline x Bern_order,
    built as fitlib/fewz.py builds fewzbern<order>."""
    lo, hi = x.getMin(), x.getMax()
    core = fewz.spline_to_expr(sp, lo, hi)
    raw = fewz.precise_integration(ROOT.RooGenericPdf(
        f"{key}_raw_{tag}", label(key), core, ROOT.RooArgList(x)))
    pars = []
    for k in range(1, order + 1):
        v = ROOT.RooRealVar(f"{key}_c{k}_{tag}", f"c{k}", BERN_INIT,
                            *BERN_RANGE)
        v.setError(BERN_STEP)
        pars.append(v)
    args = ROOT.RooArgList(x)
    for v in pars:
        args.add(v)
    u = f"((@0-{lo!r})/{hi - lo!r})"
    coefs = ["1.0"] + [f"@{k}" for k in range(1, order + 1)]
    fitted = fewz.precise_integration(ROOT.RooGenericPdf(
        f"{key}_bern{order}_{tag}", f"{label(key)} x Bern{order}",
        f"{core}*{bernstein(u, coefs)}", args))
    return raw, fitted, pars


def fit_to_mc(pdf, pars, data, tries, rng):
    """Fit, retrying from random starts until status 0 (as mcCoreFit)."""
    opts = [ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.Save(True),
            ROOT.RooFit.SumW2Error(True), ROOT.RooFit.PrintLevel(-1),
            ROOT.RooFit.PrintEvalErrors(-1)]
    best = None
    for itry in range(tries):
        for v in pars:
            if itry == 0:
                v.setVal(BERN_INIT)
            else:
                lo, hi = BERN_RANGE
                v.setVal(min(max(BERN_INIT + mcCoreFit.RANDOM_SPREAD
                                 * (hi - lo) * rng.uniform(-1, 1), lo), hi))
        strat = ROOT.RooFit.Strategy(1 if itry == 0 else 2)
        pdf.fitTo(data, strat, *opts)
        res = pdf.fitTo(data, strat, *opts)
        cand = {"status": res.status(), "nll": res.minNll(), "try": itry,
                "params": {v.GetTitle(): {"value": v.getVal(),
                                          "error": v.getError()}
                           for v in pars}}
        if best is None or mcCoreFit.rank(cand) > mcCoreFit.rank(best):
            best = cand
        if cand["status"] == 0:
            break
    for v in pars:
        v.setVal(best["params"][v.GetTitle()]["value"])
    return best


def gof_text(name, g, status=0):
    if status > 1:
        return f"{name}: fit failed"
    if g["ndf"] <= 0:
        return f"{name}: low MC stats"
    p = "p < 0.001" if g["pvalue"] < 1e-3 else f"p = {g['pvalue']:.3f}"
    return f"{name}: #chi^{{2}}/ndf = {g['chi2_ndf']:.2f}, {p}"


def frame_axes(fr):
    fr.SetTitle(f";;Weighted events / "
                f"{1 / mcCoreFit.PLOT_BINS_PER_GEV:g} GeV")
    style._style_ratio_axis(fr.GetYaxis(), style.RATIO_Y_TITLE_OFFSET)
    fr.GetYaxis().SetTitleOffset(1.45)
    ROOT.TGaxis.SetMaxDigits(3)
    ROOT.TGaxis.SetExponentOffset(-0.10, 0.01, "y")
    fr.GetXaxis().SetLabelSize(0)
    fr.GetXaxis().SetTitleSize(0)
    fr.SetMinimum(0.0)
    fr.SetMaximum(1.35 * fr.GetMaximum())


def mc_legend(nrows):
    top = 0.88
    leg = ROOT.TLegend(1 - style.PAD_RIGHT_MARGIN - 0.52,
                       top - 0.055 * nrows,
                       1 - style.PAD_RIGHT_MARGIN - 0.03, top)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.027)
    return leg


def save(c, outbase, pdf_too):
    os.makedirs(os.path.dirname(outbase) or ".", exist_ok=True)
    c.SaveAs(outbase + ".png")
    if pdf_too:
        c.SaveAs(outbase + ".pdf")


def plot_raw(x, h, splines, raw_gof, cat_label, what, year, outbase,
             pdf_too):
    """MC with each shape normalised to it (not fitted); shape / MC below,
    with the MC stat. uncertainty as a band around 1."""
    nplot = int((config.XHIGH - config.XLOW) * mcCoreFit.PLOT_BINS_PER_GEV)
    tagname = os.path.basename(outbase)
    c, outer, pad1, pad2 = style.make_canvas_pads(f"c_{tagname}")
    fr = x.frame(ROOT.RooFit.Bins(nplot))
    mc_fill, mc_band = mcCoreFit.mc_hists(h.Clone(), nplot, tagname)
    fr.addTH1(mc_fill, "hist")
    fr.addTH1(mc_band, "E2")
    leg = mc_legend(2 + len(splines))
    leg.AddEntry(mc_fill, "MC", "f")
    leg.AddEntry(mc_band, "MC stat. unc.", "f")
    total = mc_fill.Integral()
    ratios = []
    for k, sp in splines.items():
        hs = shape_hist(sp, f"{k}_{tagname}", nplot, total)
        hs.SetLineColor(colour(k))
        hs.SetLineWidth(2)
        hs.SetLineStyle(2 if k == REFERENCE else 1)
        fr.addTH1(hs, "hist")
        leg.AddEntry(hs, gof_text(label(k), raw_gof[k]), "l")
        r = hs.Clone(f"ratio_{k}_{tagname}")
        r.SetDirectory(0)
        r.Divide(mc_fill.Clone())
        for i in range(1, r.GetNbinsX() + 1):
            r.SetBinError(i, 0)
        ratios.append(r)

    pad1.cd()
    frame_axes(fr)
    fr.Draw()
    leg.Draw()
    style.cms_label(pad1, year, full_year=True)
    lat = annotate(["H #rightarrow #mu#mu", cat_label, what,
                    "not fitted"])

    pad2.cd()
    band = mc_fill.Clone(f"relband_{tagname}")
    band.SetDirectory(0)
    for i in range(1, band.GetNbinsX() + 1):
        v = mc_fill.GetBinContent(i)
        band.SetBinContent(i, 1.0)
        band.SetBinError(i, mc_fill.GetBinError(i) / v if v > 0 else 0)
    band.SetFillColor(ROOT.kGray + 2)
    band.SetFillStyle(mcCoreFit.MC_BAND_STYLE)
    band.SetLineWidth(0)
    band.SetMarkerSize(0)
    band.SetTitle(";m_{#mu#mu} [GeV];Shape / MC")
    band.SetMinimum(RATIO_RANGE[0])
    band.SetMaximum(RATIO_RANGE[1])
    band.GetYaxis().SetNdivisions(505)
    style._style_ratio_axis(band.GetXaxis(), style.RATIO_X_TITLE_OFFSET)
    style._style_ratio_axis(band.GetYaxis(), style.RATIO_Y_TITLE_OFFSET)
    band.Draw("E2")
    guides = ratio_guides(config.XLOW, config.XHIGH)
    for r in ratios:
        r.Draw("hist same")
    save(c, outbase, pdf_too)
    del lat, guides


def plot_fitted(x, data, h, fits, order, cat_label, what, year, outbase,
                pdf_too):
    """Each shape x Bern_order fitted to the MC, pulls below (mc_core)."""
    nplot = int((config.XHIGH - config.XLOW) * mcCoreFit.PLOT_BINS_PER_GEV)
    tagname = os.path.basename(outbase)
    c, outer, pad1, pad2 = style.make_canvas_pads(f"c_{tagname}")
    fr = x.frame(ROOT.RooFit.Bins(nplot))
    mc_fill, mc_band = mcCoreFit.mc_hists(h.Clone(), nplot, tagname)
    fr.addTH1(mc_fill, "hist")
    fr.addTH1(mc_band, "E2")
    data.plotOn(fr, ROOT.RooFit.Name("mc"),
                ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2),
                ROOT.RooFit.Binning(nplot), ROOT.RooFit.Invisible())
    leg = mc_legend(2 + len(fits))
    leg.AddEntry(mc_fill, "MC", "f")
    leg.AddEntry(mc_band, "MC stat. unc.", "f")
    pulls = []
    for k, f in fits.items():
        col = colour(k)
        f["pdf"].plotOn(fr, ROOT.RooFit.Name(k),
                        ROOT.RooFit.LineColor(col),
                        ROOT.RooFit.LineWidth(2),
                        ROOT.RooFit.LineStyle(2 if k == REFERENCE else 1))
        leg.AddEntry(fr.findObject(k),
                     gof_text(f"{label(k)} #times Bern{order}", f["gof"],
                              f["status"]), "l")
        ph = fr.pullHist("mc", k)
        ph.SetLineColor(col)
        ph.SetMarkerColor(col)
        ph.SetMarkerSize(0.6)
        pulls.append(ph)

    pad1.cd()
    frame_axes(fr)
    fr.Draw()
    leg.Draw()
    style.cms_label(pad1, year, full_year=True)
    lat = annotate(["H #rightarrow #mu#mu", cat_label, what,
                    "shape fit only"])

    pad2.cd()
    fr2 = x.frame()
    for ph in pulls:
        fr2.addPlotable(ph, "P")
    fr2.SetMinimum(-5)
    fr2.SetMaximum(5)
    fr2.SetTitle(";m_{#mu#mu} [GeV];Pull")
    fr2.GetYaxis().SetNdivisions(505)
    style._style_ratio_axis(fr2.GetXaxis(), style.RATIO_X_TITLE_OFFSET)
    style._style_ratio_axis(fr2.GetYaxis(), style.RATIO_Y_TITLE_OFFSET)
    fr2.Draw()
    guides = ratio_guides(config.XLOW, config.XHIGH, (0.0,))
    save(c, outbase, pdf_too)
    del lat, guides


def load_mc(args, binMVA):
    """The DY (or total background) MC histogram, as mcCoreFit.py."""
    from fitlib.prepareFits import getHistoBkgMC
    nbins = int((config.XHIGH - config.XLOW) * config.BINS_PER_GEV)
    return getHistoBkgMC(nbins, config.XLOW, config.XHIGH, args.category,
                         args.year, binMVA or "incl", dy_only=args.dy_only)


def run_bin(args, binMVA, splines, rng):
    tag = f"{args.category}_{binMVA}_{args.year}"
    what = "DY MC" if args.dy_only else "background MC"
    cat_label = f"{args.category} {binMVA or 'incl'}"
    print(f"\n===== {tag}: {what} =====")
    h, _ = load_mc(args, binMVA)
    if h is None or h.Integral() <= 0:
        raise SystemExit(f"{tag}: no MC in the window")
    x = ROOT.RooRealVar(f"mh{args.category}", "m_{#mu#mu}",
                        config.XLOW, config.XHIGH)
    data = ROOT.RooDataHist(f"mc_{tag}", what, ROOT.RooArgList(x), h)
    print(f"sum of weights {h.Integral():.1f}, "
          f"effective entries {h.GetEffectiveEntries():.0f}")

    raw_gof, fits, keep = {}, {}, []
    print(f"{'shape':<16}{'raw chi2/ndf':>13}{'p':>9}"
          f"{f'xBern{args.order} chi2/ndf':>20}{'p':>9}{'status':>8}"
          f"   coefficients")
    for k, sp in splines.items():
        raw, fitted, pars = make_pdfs(x, k, sp, args.order, tag)
        keep += [raw, fitted, pars]
        raw_gof[k] = mcCoreFit.gof_vs_mc(x, h, raw, 0)
        f = fit_to_mc(fitted, pars, data, args.tries, rng)
        f["gof"] = mcCoreFit.gof_vs_mc(x, h, fitted, len(pars))
        f["pdf"] = fitted
        fits[k] = f
        cs = " ".join(f"{p['value']:+.3f}" for p in f["params"].values())
        print(f"{label(k):<16}{raw_gof[k]['chi2_ndf']:>13.3f}"
              f"{raw_gof[k]['pvalue']:>9.3f}{f['gof']['chi2_ndf']:>20.3f}"
              f"{f['gof']['pvalue']:>9.3f}{f['status']:>8}   {cs}")

    plot_raw(x, h, splines, raw_gof, cat_label, what, args.year,
             os.path.join(args.plotdir, f"DYcore_raw_{tag}"), args.pdf)
    plot_fitted(x, data, h, fits, args.order, cat_label, what, args.year,
                os.path.join(args.plotdir,
                             f"DYcore_bern{args.order}_{tag}"), args.pdf)

    rec = {
        "provenance": {
            "written_utc": datetime.now(timezone.utc)
                           .isoformat(timespec="seconds"),
            "git_commit": config.git_commit(), "tag": tag,
            "window": [config.XLOW, config.XHIGH],
            "dy_only": args.dy_only, "order": args.order,
            "sources": {REFERENCE: fewz.FEWZ_FILE,
                        **{k: os.path.join(args.csv_dir, v)
                           for k, v in CSV_FILES.items()}},
        },
        "mc": {"sum_weights": h.Integral(),
               "effective_entries": h.GetEffectiveEntries()},
        "raw": raw_gof,
        f"bern{args.order}": {k: {kk: v for kk, v in f.items()
                                  if kk != "pdf"}
                              for k, f in fits.items()},
    }
    jf = os.path.join(args.plotdir, f"DYcore_{tag}.json")
    with open(jf, "w") as fh:
        json.dump(rec, fh, indent=2, sort_keys=True)
    print(f"wrote {jf}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--category", default="ggHcat")
    ap.add_argument("-b", "--bin", nargs="+", default=["incl"],
                    help="BDT bin label(s), as given to mcCoreFit")
    ap.add_argument("-y", "--year", default="Run3")
    ap.add_argument("--csv-dir", default=DEFAULT_CSV_DIR,
                    help="folder with fewz_smooth.csv and dyturbo_smooth.csv "
                         "(default %(default)s)")
    ap.add_argument("--plotdir",
                    default=os.path.expanduser(
                        "~/public_html/HmumuFits/dy_core"))
    ap.add_argument("--shapes-only", action="store_true",
                    help="only DYcore_shapes.png, no MC")
    ap.add_argument("--order", type=int, default=3,
                    help="Bernstein order of the shape x Bern_n fit "
                         "(default %(default)s, as fewzbern3 in corePdf)")
    ap.add_argument("--all-bkg-mc", action="store_true",
                    help="compare with the total background MC instead of "
                         "the Drell-Yan MC alone")
    ap.add_argument("--tries", type=int, default=5)
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--pdf", action="store_true")
    args = ap.parse_args()
    args.dy_only = not args.all_bkg_mc

    ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
    fewz.quiet_bin_integrator()
    style.setup_style()
    splines = load_predictions(args.csv_dir)
    plot_shapes(splines, os.path.join(args.plotdir, "DYcore_shapes"),
                args.pdf)
    if args.shapes_only:
        return
    rng = random.Random(args.seed)
    for b in args.bin:
        run_bin(args, b, splines, rng)


if __name__ == "__main__":
    main()
