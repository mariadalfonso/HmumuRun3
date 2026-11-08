"""
signal.py -- the signal shape fit: the model, one fit of one histogram
(fit_one) and its ratio / pull plots (plot_fit).

The style and the fit window come from style.py and config.py, looked up as
module attributes at call time. Studies that swap in their own pieces
assign them here (sigFitDCBG: signal.make_signal_pdf, signal._annotate,
signal.PARAM_LABEL; fitMVABins: style.header_lines).
"""

import os

import ROOT

from . import config, style
from . import pdfDefinitions as pdefs

def make_signal_pdf(x, tag, scale_unc=pdefs.SCALE_UNC, res_unc=pdefs.RES_UNC,
                    correlated=True):
    """The DSCB with MH hypothesis + scale/resolution nuisance hooks, from
    pdfDefinitions.create_signal_pdf -- the same model makeWorkspaces.py
    rebuilds in the Combine container.

    Kept as a module-level name because fit_one calls it through this
    module: fitEtaBins/fitMVABins use it, and sigFitDCBG swaps in its own
    model by assigning signal.make_signal_pdf.

    Returns (pdf, bundle).
    """
    return pdefs.create_signal_pdf(x, tag, scale_unc=scale_unc,
                                   res_unc=res_unc, correlated=correlated)


def near_bound(var, tol=0.01):
    """True if a floating parameter ended up parked on a limit."""
    lo, hi = var.getMin(), var.getMax()
    if hi <= lo:
        return False
    frac = (var.getVal() - lo) / (hi - lo)
    return frac < tol or frac > (1.0 - tol)


PARAM_LABEL = {
    "mu":     ("#mu",            "GeV"),
    "sigma":  ("#sigma",         "GeV"),
    "alphaL": ("#alpha_{lo}",    ""),
    "nL":     ("n_{lo}",         ""),
    "alphaR": ("#alpha_{high}",  ""),
    "nR":     ("n_{high}",       ""),
}


def fmt_with_unc(label, val, unc, unit=""):
    """value +/- uncertainty with precision matched to the uncertainty.

    TODO: Check this later
    """
    if unc < 0.025:
        s = f"{label} = {val:.3f} #pm {unc:.3f}"
    elif unc < 0.250:
        s = f"{label} = {val:.2f} #pm {unc:.2f}"
    elif unc < 2.500:
        s = f"{label} = {val:.1f} #pm {unc:.1f}"
    else:
        s = f"{label} = {val:.0f} #pm {unc:.0f}"
    return s + (f" {unit}" if unit else "")


def _make_ratio(hist, pdf, x, norm):
    """Template/Fit per bin, inside the fit window.

    The expected count in a bin is norm * (normalised density at the bin
    centre) * (bin width).
    """
    b_lo, b_hi = range_bins(hist)
    ratio = ROOT.TH1D(f"ratio_{hist.GetName()}", "",
                      b_hi - b_lo + 1,
                      hist.GetBinLowEdge(b_lo),
                      hist.GetBinLowEdge(b_hi) + hist.GetBinWidth(b_hi))
    ratio.SetDirectory(0)
    nset = ROOT.RooArgSet(x)

    saved = x.getVal()
    for i in range(b_lo, b_hi + 1):
        j = i - b_lo + 1
        content = hist.GetBinContent(i)
        x.setVal(hist.GetBinCenter(i))
        expected = norm * pdf.getVal(nset) * hist.GetBinWidth(i)

        # Empty data bin set to off-scale so nothing is drawn.
        if content <= 0 or expected <= 0:
            ratio.SetBinContent(j, -999.0)
            ratio.SetBinError(j, 0.0)
            continue

        ratio.SetBinContent(j, content / expected)
        ratio.SetBinError(j, hist.GetBinError(i) / expected)

    x.setVal(saved)          # getVal() mutated x; restore it
    return ratio


def _annotate(nom, chi2_ndf, norm, n_eff, label):
    """Category / yield / fitted-parameter block on the upper pad.

    Returned so the caller keeps the TLatex alive until the canvas is written
    """
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(style.ANN_TEXT_SIZE)

    left = style.PAD_LEFT_MARGIN + 0.05
    y0, dy = 0.85, style.LINE_SPACING * style.ANN_TEXT_SIZE
    for i, line in enumerate(style.header_lines(label)):
        latex.DrawLatex(left, y0 - i * dy, line)

    # right-hand column: yield and fit parameters
    latex.SetTextSize(style.PARAM_TEXT_SIZE)
    latex.SetTextAlign(12)
    right = 1.0 - style.PAD_RIGHT_MARGIN - 0.36
    y = 0.85
    lines = [f"yield = {norm:.1f}", f"N_{{eff}} = {n_eff:.0f}"]
    for k in ("mu", "sigma", "alphaL", "nL", "alphaR", "nR"):
        v = nom[k]
        lab, unit = PARAM_LABEL[k]
        txt = fmt_with_unc(lab, v.getVal(), v.getError(), unit)
        if near_bound(v):
            txt += " #color[2]{(bound)}"
        lines.append(txt)
    lines.append(f"#chi^{{2}}/ndf = {chi2_ndf:.2f}")
    pdy = style.LINE_SPACING * style.PARAM_TEXT_SIZE
    for i, line in enumerate(lines):
        latex.DrawLatex(right, y - i * pdy, line)

    return latex

def fit_legend(frame, x1=None, y1=None, x2=None, y2=None,
               model="DSCB model", data="Simulation"):
    """Legend for a mass-fit panel: the points, then the curve.

    Takes the drawn objects off the RooPlot by the names given to
    plotOn ("dat", "fit"), so the marker and line styles in the legend are
    whatever was actually drawn rather than a second set that can drift.
    """
    x1 = style.PAD_LEFT_MARGIN + 0.03 if x1 is None else x1
    x2 = style.PAD_LEFT_MARGIN + 0.25 if x2 is None else x2
    y2 = 0.65 if y2 is None else y2
    y1 = 0.5 if y1 is None else y1

    leg = ROOT.TLegend(x1, y1, x2, y2)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.034)
    d = frame.findObject("dat")
    f = frame.findObject("fit")
    if d:
        leg.AddEntry(d, data, "pe")
    if f:
        leg.AddEntry(f, model, "l")
    leg.Draw()
    return leg


def plot_fit(x, data, pdf, nom, hist, label, outbase, n_float, norm, n_eff,
             year=None, save_pdf=False, scale=1.0):
    """Write two 2-panel canvases and return chi2/ndf.

    <outbase>_ratio.pdf/.png   fit + Template/Fit
    <outbase>_pull.pdf/.png    fit + pull

    scale: the full-year yield over the MC integral. When it is not 1 the
    template is drawn scaled by it and labelled with the data luminosity,
    so the plot shows the yield the workspace carries; the pulls, ratio and
    chi2 do not change under a constant scale.
    """
    if scale != 1.0:
        hist = hist.Clone(f"{hist.GetName()}_full_year")
        hist.SetDirectory(0)
        hist.Scale(scale)
        data = ROOT.RooDataHist(f"{data.GetName()}_full_year", "data",
                                ROOT.RooArgList(x), hist)
        norm *= scale
    frame = x.frame(ROOT.RooFit.Title(""))
    data.plotOn(frame, ROOT.RooFit.Name("dat"),
                ROOT.RooFit.MarkerStyle(style.DATA_MARKER_STYLE),
                ROOT.RooFit.MarkerSize(style.DATA_MARKER_SIZE),
                ROOT.RooFit.MarkerColor(ROOT.kBlack),
                ROOT.RooFit.LineColor(ROOT.kBlack),
                ROOT.RooFit.LineWidth(style.DATA_LINE_WIDTH),
                ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
    pdf.plotOn(frame, ROOT.RooFit.Name("fit"),
               ROOT.RooFit.LineColor(style.fit_color()),
               ROOT.RooFit.LineWidth(2))

    chi2_ndf = frame.chiSquare("fit", "dat", n_float)

    frame.SetTitle(f";;Events / {hist.GetBinWidth(1):.2f} GeV")
    # Absolute fonts on both pads, with the same y title offset.
    style._style_ratio_axis(frame.GetYaxis(), style.RATIO_Y_TITLE_OFFSET)
    frame.GetXaxis().SetLabelSize(0)      # shared axis: labels on lower pad
    frame.GetXaxis().SetTitleSize(0)
    frame.SetMaximum(1.1 * frame.GetMaximum())   # headroom

    ratio = _make_ratio(hist, pdf, x, norm)
    style._style_lower_marker(ratio)
    ratio.GetYaxis().SetTitle("Template/Fit")
    ratio.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
    ratio.GetYaxis().SetRangeUser(0.0, 2.0)
    ratio.GetYaxis().SetNdivisions(505)
    style._style_ratio_axis(ratio.GetXaxis(), style.RATIO_X_TITLE_OFFSET)
    style._style_ratio_axis(ratio.GetYaxis(), style.RATIO_Y_TITLE_OFFSET)

    pull_frame = x.frame(ROOT.RooFit.Title(""))
    pull_hist = frame.pullHist("dat", "fit")

    # An empty data bin has zero error, so its pull is undefined
    # RooFit sets these points to 0, removed from ther plot.
    for _i in range(pull_hist.GetN() - 1, -1, -1):
        if (pull_hist.GetErrorYhigh(_i) <= 0.0
                and pull_hist.GetErrorYlow(_i) <= 0.0):
            pull_hist.RemovePoint(_i)
    style._style_lower_marker(pull_hist)
    pull_frame.addPlotable(pull_hist, "P")
    pull_frame.SetMinimum(-6.0)
    pull_frame.SetMaximum(6.0)
    pull_frame.GetYaxis().SetTitle("#frac{Template#minusFit}{#deltaTemplate}")
    pull_frame.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
    pull_frame.GetYaxis().SetNdivisions(409)
    style._style_ratio_axis(pull_frame.GetXaxis(), style.RATIO_X_TITLE_OFFSET)
    style._style_ratio_axis(pull_frame.GetYaxis(), style.RATIO_Y_TITLE_OFFSET)

    os.makedirs(os.path.dirname(outbase), exist_ok=True)

    for kind, lower, central, guides, draw in (
            ("ratio", ratio,      1.0, style.RATIO_GUIDES, "pe"),
            ("pull",  pull_frame, 0.0, style.PULL_GUIDES,  "")):

        canvas, outer, pad1, pad2 = style.make_canvas_pads(f"c_{kind}")

        pad1.cd()
        frame.Draw()
        keep_txt = _annotate(nom, chi2_ndf, norm, n_eff, label)
        keep_cms = style.cms_label(pad1, year, full_year=(scale != 1.0))
        keep_leg = fit_legend(frame)

        pad2.cd()
        lower.Draw(draw)

        up, down = guides
        l0 = ROOT.TLine(config.XLOW, central, config.XHIGH, central)
        l0.SetLineColor(style.fit_color())
        l0.SetLineWidth(2)
        lu = ROOT.TLine(config.XLOW, up, config.XHIGH, up)
        ld = ROOT.TLine(config.XLOW, down, config.XHIGH, down)
        for ln in (lu, ld):
            ln.SetLineColor(11)
            ln.SetLineWidth(1)
            ln.SetLineStyle(ROOT.kDashed)
        l0.Draw("same")
        lu.Draw("same")
        ld.Draw("same")

        canvas.Update()
        canvas.SaveAs(f"{outbase}_{kind}.png")
        if save_pdf:
            canvas.SaveAs(f"{outbase}_{kind}.pdf")
        del keep_txt, keep_cms, keep_leg, pad1, pad2, outer, canvas

    return chi2_ndf


# ---------------------------------------------------------------------------
# the fit
# ---------------------------------------------------------------------------

def range_bins(hist, lo=config.XLOW, hi=config.XHIGH):
    """First and last bin strictly inside [lo, hi).

    The +/-1e-6 avoids the off-by-one at the edges: FindBin(hi) returns the
    bin *starting* at hi, which is outside the window.
    """
    return hist.FindBin(lo + 1e-6), hist.FindBin(hi - 1e-6)


def range_integral(hist, lo=config.XLOW, hi=config.XHIGH):
    """Integral strictly inside [lo, hi).

    Avoids GetSumOfWeights() silently including anything outside the fit
    window if the histogram axis is ever wider than the fit range.
    """
    b_lo, b_hi = range_bins(hist, lo, hi)
    return hist.Integral(b_lo, b_hi)


def effective_entries(hist, lo=config.XLOW, hi=config.XHIGH):
    """N_eff = (sum content)^2 / sum(error^2), inside the fit window.

    Computed from bin errors rather than TH1::GetEffectiveEntries(), which
    relies on ROOT's internal sum-of-weights bookkeeping. That bookkeeping
    is only filled by Fill(); a histogram built with SetBinContent (e.g.
    RooAbsData::createHistogram) reports a meaningless value.
    """
    b_lo, b_hi = range_bins(hist, lo, hi)
    s = sum(hist.GetBinContent(i) for i in range(b_lo, b_hi + 1))
    s2 = sum(hist.GetBinError(i) ** 2 for i in range(b_lo, b_hi + 1))
    return (s * s / s2) if s2 > 0 else 0.0


def fit_one(x, hist, tag, title, plotdir, freeze=True, year=None,
            save_pdf=False, yield_info=None):
    """Fit one histogram. Returns (pdf, bundle, norm_var, record) or None.

    `yield_info` (from sigFit.run_category) carries the signal yield extrapolated
    to the full year; it replaces the MC integral in the norm and the JSON,
    and the plot draws the template scaled to it.
    """
    norm = range_integral(hist)
    n_eff = effective_entries(hist)
    if norm <= 0:
        print(f"  !! empty histogram for {tag}, skipping")
        return None

    data = ROOT.RooDataHist(f"datahist_{tag}", "data", ROOT.RooArgList(x), hist)
    pdf, bundle = make_signal_pdf(x, tag)
    nom, nuis = bundle["nom"], bundle["nuis"]

    # Nuisances pinned at nominal: they are exactly degenerate with
    # cb_mu / cb_sigma, so letting them float here would make the fit
    # wander along a flat direction and report meaningless errors.
    for v in nuis.values():
        v.setVal(0.0)
        v.setConstant(True)

    fitopts = [
        ROOT.RooFit.Minimizer("Minuit2"),
        ROOT.RooFit.Strategy(2),
        ROOT.RooFit.Save(True),
        ROOT.RooFit.Range("full"),
        # Correct errors for weighted BINNED data. AsymptoticError is for
        # weighted *unbinned* fits and RooFit warns if used on a RooDataHist.
        ROOT.RooFit.SumW2Error(True),
        ROOT.RooFit.PrintLevel(-1),
    ]
    pdf.fitTo(data, *fitopts)                # first pass
    result = pdf.fitTo(data, *fitopts)       # refit from the result

    n_float = sum(0 if v.isConstant() else 1 for v in nom.values())
    sig_yield = norm if yield_info is None else yield_info["value"]
    chi2_ndf = plot_fit(x, data, pdf, nom, hist, title,
                        f"{plotdir}/signal_{tag}", n_float, norm, n_eff,
                        year=year, save_pdf=save_pdf, scale=sig_yield / norm)

    status, cov, edm = result.status(), result.covQual(), result.edm()
    parked = [k for k, v in nom.items() if near_bound(v)]
    ok = (status == 0 and cov == 3 and not parked)

    print(f"  mu     = {nom['mu'].getVal():7.3f} +/- {nom['mu'].getError():.3f}")
    print(f"  sigma  = {nom['sigma'].getVal():7.3f} +/- {nom['sigma'].getError():.3f}")
    print(f"  norm   = {norm:10.2f}   (N_eff = {n_eff:.1f})")
    if yield_info is not None:
        print(f"  yield  = {sig_yield:10.2f}   extrapolated to {year} "
              f"(x{sig_yield / norm:.3f})")
    print(f"  chi2/ndf = {chi2_ndf:6.3f}   status={status} covQual={cov} "
          f"edm={edm:.2e}")
    if parked:
        print(f"  !! parameters at bound: {', '.join(parked)}")
    if not ok:
        print("  <-- CHECK THIS FIT")

    try:
        corr_mu_sigma = result.correlation(nom["mu"], nom["sigma"])
    except Exception:
        corr_mu_sigma = None

    record = {
        "parameters": {
            k: {"value": v.getVal(), "error": v.getError(),
                "min": v.getMin(), "max": v.getMax(),
                "at_bound": near_bound(v)}
            for k, v in nom.items()
        },
        "shape_model": {
            "formula_mean":
                f"(cb_mu + MH - {pdefs.MH_REF})*(1 + {bundle['cfg']['scale_unc']}*CMS_scale_m)",
            "formula_sigma":
                f"cb_sigma*(1 + {bundle['cfg']['res_unc']}*CMS_res_m)",
            **bundle["cfg"],
        },
        "normalization": {
            "integral": sig_yield,          # the datacard rate
            "effective_entries": n_eff,
            **({"mc_integral": norm, "extrapolation": yield_info}
               if yield_info is not None else {}),
        },
        "fit_quality": {
            "chi2_ndf": chi2_ndf,
            "status": status,
            "cov_qual": cov,
            "edm": edm,
            "n_float": n_float,
            "params_at_bound": parked,
            "corr_mu_sigma": corr_mu_sigma,
            "ok": ok,
        },
    }

    norm_var = ROOT.RooRealVar(f"{pdf.GetName()}_norm",
                               f"{pdf.GetName()}_norm", sig_yield)
    norm_var.setConstant(True)

    if freeze:
        for v in nom.values():          # shape from MC is now fixed ...
            v.setConstant(True)
        for v in nuis.values():         # ... nuisances released for Combine
            v.setConstant(False)

    return pdf, bundle, norm_var, record
