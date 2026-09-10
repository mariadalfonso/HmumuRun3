"""
sigFit.py -- standalone H->mumu signal shape fits.

No Combine dependency for the fitting itself: uses ROOT's built-in
RooCrystalBall (root >= 6.26), so this runs in the plain conda env.

Signal shape model
------------------
Six CB parameters are fitted to MC per (category, BDT bin, production mode)
and then frozen. The PDF does not see cb_mu / cb_sigma directly; it sees:

    mean  = (cb_mu + MH - 125) * (1 + scale_unc * CMS_scale_m)
    sigma =  cb_sigma          * (1 + res_unc   * CMS_res_m)

  MH           Higgs mass hypothesis (GeV). Additive shift, so the same
               workspace can be re-evaluated at 124, 126, ... for a mass
               scan. Combine drives this with -m. NOT a systematic.
  CMS_scale_m  unit-Gaussian nuisance for muon momentum scale
  CMS_res_m    unit-Gaussian nuisance for muon momentum resolution

The nuisances are exactly degenerate with cb_mu / cb_sigma, so they are
pinned to zero during the MC fit and only released afterwards.

Usage
-----
    python sigFit.py --selftest
    python sigFit.py -c VBFcat -b bdt0 -s qqH
    python sigFit.py -c VBFcat -b incl bdt0 bdt1 bdt2
    python sigFit.py -c ggHcat VBFcat VLcat VHcat Zinvcat TTLcat TTHcat \
                     -b bdt0 bdt1 bdt2
"""

import argparse
import json
import os
import subprocess
from datetime import datetime, timezone

import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
# Empty tail bins have zero error, so their pull is undefined
# Supress these erros for EVERY stream
for _i in range(ROOT.RooMsgService.instance().numStreams()):
    ROOT.RooMsgService.instance().getStream(_i).removeTopic(ROOT.RooFit.Plotting)
ROOT.TH1.SetDefaultSumw2(True)

# ---------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------

XLOW, XHIGH = 110.0, 150.0
BINS_PER_GEV = 10         # 0.1 GeV bins, matching SIGfits.py and fitBkg.
                          # Checked with a closure-test scan (see README):
                          # unbiased at -0.11% +/- 0.25% on sigma. Anything
                          # from 0.1 to 0.5 GeV is fine; 1.0 GeV biases sigma
                          # by +2%. Override with --bins-per-gev.
MH_REF = 125.0            # reference mass the shape was fitted at

# TODO: Fractional uncertainties folded into the shape formulas.
# Replace with Run 3 muon POG numbers or from objScaleSmear().
SCALE_UNC = 0.002         # 0.2% muon momentum scale
RES_UNC   = 0.05          # 5% muon momentum resolution

# ---------------------------------------------------------------------------
# Plot style.
# ---------------------------------------------------------------------------
CMS_TEXT = "CMS"
CMS_EXTRA_TEXT = "Internal"      # drawn in font 52 (italic)
SQRT_S_TEV = 13.6

# Official CMS label size formulas: size = FRAC * pad top margin.
CMS_TEXT_SIZE_FRAC = 0.75
EXTRA_OVER_CMS_TEXT_SIZE = 0.76
LUMI_TEXT_SIZE_FRAC = 0.6

# Integrated luminosity.
#
# TODO: check the correct Run 3 lumi
RUN3_ERAS = ["_12022", "_22022", "_12023", "_22023", "_2024", "_2025", "_2026"]

LUMIS = {
    '_12016': 19.52,  # APV (B-F for 2016 pre)
    '_22016': 16.80,  # postVFP
    '_2016': 35.9,
    '_2017': 41.5,
    '_12017': 7.7,    # (F for 2017) for VBF
    '_2018': 59.70,
    '_12018': 39.54,
    '_all': 86.92,    # 19.52 + 7.7 + 59.70
    '_Run2': 138.,
    #
    '_12022': 7.99,   # C-D
    '_22022': 26.68,  # E, F, G
    '_12023': 17.96,  # C
    '_22023': 9.68,   # D
    '_2024': 109.82,  # C-I
    '_2025': 110.59,  # C-G
    '_2026': 25.31,   # C, B
    #
    '_2025C': 21.63,
    '_2025D': 25.52,
    '_2025E': 14.15,
    '_2025F': 26.89,
    '_2025G': 22.40,
}

# canvas / pad geometry (plot_style values)
CANVAS_W, CANVAS_H = 600, 600
RATIO_SPLIT = 0.30
# Dashed guide lines
RATIO_GUIDES = (1.2, 0.8)   # Template/Fit, +/-0.2 around 1
PULL_GUIDES = (2.0, -2.0)   # +/-2 sigma
PAD_TOP_MARGIN = 0.08
PAD_BOTTOM_MARGIN = 0.13
PAD_LEFT_MARGIN = 0.16
PAD_RIGHT_MARGIN = 0.05

# Lower-pad text in ABSOLUTE pixel sizes (font code ending in 3), so it
# renders consistently despite the pad's small height.
RATIO_X_TITLE_OFFSET = 1.0
RATIO_Y_TITLE_OFFSET = 1.4

# On-plot annotation. Sizes are fractions of the UPPER pad's height (font 42),
# so they scale with RATIO_SPLIT. Line spacing is derived from the text size
# rather than fixed, so raising the size does not make the lines overlap.
ANN_TEXT_SIZE = 0.040        # left block: "H -> mumu" and the category
PARAM_TEXT_SIZE = 0.040      # right block: yield, N_eff, fitted parameters
LINE_SPACING = 1.35          # dy = LINE_SPACING * text size
PARAM_BLOCK_WIDTH = 0.42     # left edge of the right block, from the frame's
                             # right edge. Widen if long lines get clipped.

DATA_MARKER_STYLE = 20
DATA_MARKER_SIZE = 1.2
DATA_LINE_WIDTH = 2

FIT_COLOR_RGB = (237, 41, 57)     # plot_style "redMed"

BOUNDS = {
    "mu":     (MH_REF, 120.0, 130.0),
    "sigma":  (2.0,      0.5,   6.0),
    "alphaL": (1.5,      0.1,   3.0),
    "nL":     (5.0,      0.1,  50.0),
    "alphaR": (1.5,      0.1,   3.0),
    "nR":     (5.0,      0.1,  50.0),
}
MC_LIST = {
    "ggHcat":  ["ggH", "qqH", "VH", "ttH"],
    "VBFcat":  ["ggH", "qqH", "VH", "ttH"],
    "Zinvcat": ["qqH", "VH", "ttH"],
    "VLcat":   ["VH", "ttH"],
    "VHcat":   ["VH", "ttH"],
    "TTLcat":  ["VH", "ttH"],
    "TTHcat":  ["VH", "ttH"],
}
CAT_LABEL = {
    "ggHcat":  "ggH cat.",
    "VBFcat":  "VBF cat.",
    "VHcat":   "VH cat.",
    "VLcat":   "V+lepton cat.",
    "Zinvcat": "Z(#nu#nu) cat.",
    "TTHcat":  "t#bar{t}H cat.",
    "TTLcat":  "t#bar{t}+lepton cat.",
}
BIN_LABEL = {
    "incl": "inclusive",
    "bdt0": "BDT bin 0",
    "bdt1": "BDT bin 1",
    "bdt2": "BDT bin 2",
    "":     "inclusive",
}
MODE_LABEL = {
    "ggH": "ggH",
    "qqH": "VBF",       # mc=10 sample, confirmed VBF
    "VH":  "VH",
    "ttH": "t#bar{t}H",
}


def header_lines(label):
    """Readable annotation header from a terse "{cat} {bin} {sig}" label."""
    parts = label.split()
    if (len(parts) == 3 and parts[0] in CAT_LABEL
            and parts[1] in BIN_LABEL and parts[2] in MODE_LABEL):
        cat, binMVA, sig = parts
        return [
            "H #rightarrow #mu#mu",
            f"{CAT_LABEL[cat]} {BIN_LABEL[binMVA]}",
            f"Signal: {MODE_LABEL[sig]}",
        ]
    return ["H #rightarrow #mu#mu", label]


# ---------------------------------------------------------------------------
# PDF
# ---------------------------------------------------------------------------

def make_signal_pdf(x, tag, scale_unc=SCALE_UNC, res_unc=RES_UNC,
                    correlated=True):
    """Build the CB with MH hypothesis + scale/resolution nuisance hooks.

    Returns (pdf, bundle). Everything is kept in the bundle because PyROOT
    will garbage-collect any RooFit object you do not hold a reference to.
    """
    # --- shape parameters measured from MC ---
    nom = {}
    for k, (start, lo, hi) in BOUNDS.items():
        nom[k] = ROOT.RooRealVar(f"cb_{k}_{tag}", f"cb_{k}", start, lo, hi)

    # --- mass hypothesis: shared across every category and bin ---
    # Named MH because that is what Combine's -m option looks for.
    mh = ROOT.RooRealVar("MH", "Higgs mass hypothesis", MH_REF, 120.0, 130.0)
    mh.setConstant(True)

    # --- nuisance parameters: unit Gaussians, constrained by Combine ---
    # correlated=True -> one shared parameter across all categories, which
    # is right for muon scale/resolution: it is the same detector effect.
    nsuf = "" if correlated else f"_{tag}"
    nuis = {
        "scale": ROOT.RooRealVar(f"CMS_scale_m{nsuf}", "muon scale",      0.0, -5.0, 5.0),
        "res":   ROOT.RooRealVar(f"CMS_res_m{nsuf}",   "muon resolution", 0.0, -5.0, 5.0),
    }

    # --- the shape the PDF actually uses ---
    mean = ROOT.RooFormulaVar(
        f"cb_mean_{tag}", "mean",
        f"(@0 + @1 - {MH_REF})*(1 + {scale_unc}*@2)",
        ROOT.RooArgList(nom["mu"], mh, nuis["scale"]),
    )
    sigma = ROOT.RooFormulaVar(
        f"cb_sigmaEff_{tag}", "sigmaEff",
        f"@0*(1 + {res_unc}*@1)",
        ROOT.RooArgList(nom["sigma"], nuis["res"]),
    )

    pdf = ROOT.RooCrystalBall(
        f"crystal_ball_{tag}", "crystal_ball",
        x, mean, sigma,
        nom["alphaL"], nom["nL"], nom["alphaR"], nom["nR"],
    )

    bundle = {
        "nom": nom,
        "nuis": nuis,
        "mh": mh,
        "func": {"mean": mean, "sigma": sigma},
        "cfg": {"scale_unc": scale_unc, "res_unc": res_unc,
                "correlated": correlated, "mh_ref": MH_REF},
    }
    return pdf, bundle


# Combined Run 3 Lumis, derived from RUN3_ERAS
# TODO: this sums to 308.03, but
#       prepareFits.lumis['_Run3'] hardcodes 312, and
#       limitPlot.py passes 286500 pb^-1 (= 286.5).
LUMIS['_Run3'] = sum(LUMIS[e] for e in RUN3_ERAS)

# ---------------------------------------------------------------------------
# diagnostics
# ---------------------------------------------------------------------------

def get_lumi(year):
    """Integrated luminosity in fb^-1 for a '_<year>' tag."""
    key = year if str(year).startswith("_") else f"_{year}"

    return LUMIS.get(key, 0.0)
    
def fit_color():
    """ROOT colour index for the fit curve."""
    return ROOT.TColor.GetColor(*FIT_COLOR_RGB)


def setup_style():
    """Minimal global style: no stat box, no title, ticks on all four sides."""
    ROOT.gStyle.SetOptStat(0)
    ROOT.gStyle.SetOptTitle(0)
    ROOT.gStyle.SetPadTickX(1)
    ROOT.gStyle.SetPadTickY(1)
    ROOT.gStyle.SetLabelFont(42, "xyz")
    ROOT.gStyle.SetTitleFont(42, "xyz")
    ROOT.gStyle.SetLegendBorderSize(0)
    ROOT.gStyle.SetHistLineColor(ROOT.kBlack)
    ROOT.gStyle.SetMarkerColor(ROOT.kBlack)
    ROOT.gStyle.SetLineColor(ROOT.kBlack)


def cms_label(pad, year):
    """CMS + extra text top-left, lumi + energy top-right, above the frame."""
    pad.cd()
    t = pad.GetTopMargin()
    l = pad.GetLeftMargin()
    r = pad.GetRightMargin()
    y = 1 - t + 0.2 * t

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(ROOT.kBlack)

    latex.SetTextFont(61)
    latex.SetTextAlign(11)
    latex.SetTextSize(CMS_TEXT_SIZE_FRAC * t)
    latex.DrawLatex(l, y, CMS_TEXT)

    if CMS_EXTRA_TEXT:
        latex.SetTextFont(52)
        latex.SetTextSize(EXTRA_OVER_CMS_TEXT_SIZE * CMS_TEXT_SIZE_FRAC * t)
        latex.DrawLatex(l + 0.10, y, CMS_EXTRA_TEXT)

    lumi = get_lumi(year) if year else 0.0
    if lumi <= 0:
        print(f"warning: no luminosity found for year tag '{year}'")
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.SetTextSize(LUMI_TEXT_SIZE_FRAC * t)
    if lumi > 0:
        txt = "%.1f fb^{-1} (%.1f TeV)" % (lumi, SQRT_S_TEV)
    else:
        txt = "(%.1f TeV)" % SQRT_S_TEV      # e.g. the closure test
    latex.DrawLatex(1 - r, y, txt)
    return latex


def make_canvas_pads(name):
    """Canvas with upper (fit) and lower (ratio/pull) pads glued together"""
    c = ROOT.TCanvas(name, "", CANVAS_W, CANVAS_H)
    outer = ROOT.TPad(f"outer_{name}", "", 0, 0, 1, 1)
    outer.SetTickx(False)
    outer.SetTicky(False)
    outer.Draw()
    outer.cd()

    pad1 = ROOT.TPad(f"upper_{name}", "", 0.0, RATIO_SPLIT, 1.0, 1.0)
    pad1.SetTopMargin(PAD_TOP_MARGIN)
    pad1.SetBottomMargin(0.0)
    pad1.SetLeftMargin(PAD_LEFT_MARGIN)
    pad1.SetRightMargin(PAD_RIGHT_MARGIN)

    pad2 = ROOT.TPad(f"lower_{name}", "", 0.0, 0.0, 1.0, RATIO_SPLIT)
    pad2.SetTopMargin(0.0)
    pad2.SetBottomMargin(PAD_BOTTOM_MARGIN / RATIO_SPLIT)
    pad2.SetLeftMargin(PAD_LEFT_MARGIN)
    pad2.SetRightMargin(PAD_RIGHT_MARGIN)

    pad1.Draw()
    pad2.Draw()
    return c, outer, pad1, pad2


def _style_ratio_axis(axis, offset):
    """Absolute-pixel lower-pad text."""
    axis.SetTitleFont(43)
    axis.SetLabelFont(43)
    axis.SetTitleSize(24)
    axis.SetLabelSize(20)
    axis.SetTitleOffset(offset)


def _style_lower_marker(obj):
    """Black points for a lower-pad plottable."""
    obj.SetMarkerStyle(DATA_MARKER_STYLE)
    obj.SetMarkerSize(DATA_MARKER_SIZE)
    obj.SetMarkerColor(ROOT.kBlack)
    obj.SetLineColor(ROOT.kBlack)


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
    latex.SetTextSize(ANN_TEXT_SIZE)

    left = PAD_LEFT_MARGIN + 0.05
    y0, dy = 0.85, LINE_SPACING * ANN_TEXT_SIZE
    for i, line in enumerate(header_lines(label)):
        latex.DrawLatex(left, y0 - i * dy, line)

    # right-hand column: yield and fit parameters
    latex.SetTextSize(PARAM_TEXT_SIZE)
    latex.SetTextAlign(12)
    right = 1.0 - PAD_RIGHT_MARGIN - 0.36
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
    pdy = LINE_SPACING * PARAM_TEXT_SIZE
    for i, line in enumerate(lines):
        latex.DrawLatex(right, y - i * pdy, line)

    return latex

def plot_fit(x, data, pdf, nom, hist, label, outbase, n_float, norm, n_eff,
             year=None, save_pdf=False):
    """Write two 2-panel canvases and return chi2/ndf.

    <outbase>_ratio.pdf/.png   fit + Template/Fit
    <outbase>_pull.pdf/.png    fit + pull
    """
    frame = x.frame(ROOT.RooFit.Title(""))
    data.plotOn(frame, ROOT.RooFit.Name("dat"),
                ROOT.RooFit.MarkerStyle(DATA_MARKER_STYLE),
                ROOT.RooFit.MarkerSize(DATA_MARKER_SIZE),
                ROOT.RooFit.MarkerColor(ROOT.kBlack),
                ROOT.RooFit.LineColor(ROOT.kBlack),
                ROOT.RooFit.LineWidth(DATA_LINE_WIDTH),
                ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))
    pdf.plotOn(frame, ROOT.RooFit.Name("fit"),
               ROOT.RooFit.LineColor(fit_color()),
               ROOT.RooFit.LineWidth(2))

    chi2_ndf = frame.chiSquare("fit", "dat", n_float)

    frame.SetTitle(f";;Events / {hist.GetBinWidth(1):.2f} GeV")
    # Absolute fonts on both pads, with the same y title offset.
    _style_ratio_axis(frame.GetYaxis(), RATIO_Y_TITLE_OFFSET)
    frame.GetXaxis().SetLabelSize(0)      # shared axis: labels on lower pad
    frame.GetXaxis().SetTitleSize(0)
    frame.SetMaximum(1.1 * frame.GetMaximum())   # headroom

    ratio = _make_ratio(hist, pdf, x, norm)
    _style_lower_marker(ratio)
    ratio.SetLineWidth(DATA_LINE_WIDTH)
    ratio.GetYaxis().SetTitle("Template/Fit")
    ratio.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
    ratio.GetYaxis().SetRangeUser(0.0, 2.0)
    ratio.GetYaxis().SetNdivisions(505)
    _style_ratio_axis(ratio.GetXaxis(), RATIO_X_TITLE_OFFSET)
    _style_ratio_axis(ratio.GetYaxis(), RATIO_Y_TITLE_OFFSET)

    pull_frame = x.frame(ROOT.RooFit.Title(""))
    pull_hist = frame.pullHist("dat", "fit")

    # An empty data bin has zero error, so its pull is undefined
    # RooFit sets these points to 0, removed from ther plot.
    for _i in range(pull_hist.GetN() - 1, -1, -1):
        if (pull_hist.GetErrorYhigh(_i) <= 0.0
                and pull_hist.GetErrorYlow(_i) <= 0.0):
            pull_hist.RemovePoint(_i)
    _style_lower_marker(pull_hist)
    pull_frame.addPlotable(pull_hist, "P")
    pull_frame.SetMinimum(-6.0)
    pull_frame.SetMaximum(6.0)
    pull_frame.GetYaxis().SetTitle("#frac{Template#minusFit}{#deltaTemplate}")
    pull_frame.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
    pull_frame.GetYaxis().SetNdivisions(409)
    _style_ratio_axis(pull_frame.GetXaxis(), RATIO_X_TITLE_OFFSET)
    _style_ratio_axis(pull_frame.GetYaxis(), RATIO_Y_TITLE_OFFSET)

    os.makedirs(os.path.dirname(outbase), exist_ok=True)

    for kind, lower, central, guides, draw in (
            ("ratio", ratio,      1.0, RATIO_GUIDES, "pe"),
            ("pull",  pull_frame, 0.0, PULL_GUIDES,  "")):

        canvas, outer, pad1, pad2 = make_canvas_pads(f"c_{kind}")

        pad1.cd()
        frame.Draw()
        keep_txt = _annotate(nom, chi2_ndf, norm, n_eff, label)
        keep_cms = cms_label(pad1, year)

        pad2.cd()
        lower.Draw(draw)

        up, down = guides
        l0 = ROOT.TLine(XLOW, central, XHIGH, central)
        l0.SetLineColor(fit_color())
        l0.SetLineWidth(2)
        lu = ROOT.TLine(XLOW, up, XHIGH, up)
        ld = ROOT.TLine(XLOW, down, XHIGH, down)
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
        del keep_txt, keep_cms, pad1, pad2, outer, canvas

    return chi2_ndf


def git_commit():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "unknown"


# ---------------------------------------------------------------------------
# the fit
# ---------------------------------------------------------------------------

def range_bins(hist, lo=XLOW, hi=XHIGH):
    """First and last bin strictly inside [lo, hi).

    The +/-1e-6 avoids the off-by-one at the edges: FindBin(hi) returns the
    bin *starting* at hi, which is outside the window.
    """
    return hist.FindBin(lo + 1e-6), hist.FindBin(hi - 1e-6)


def range_integral(hist, lo=XLOW, hi=XHIGH):
    """Integral strictly inside [lo, hi).

    Avoids GetSumOfWeights() silently including anything outside the fit
    window if the histogram axis is ever wider than the fit range.
    """
    b_lo, b_hi = range_bins(hist, lo, hi)
    return hist.Integral(b_lo, b_hi)


def effective_entries(hist, lo=XLOW, hi=XHIGH):
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
            save_pdf=False):
    """Fit one histogram. Returns (pdf, bundle, norm_var, record) or None."""
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
    chi2_ndf = plot_fit(x, data, pdf, nom, hist, title,
                        f"{plotdir}/signal_{tag}", n_float, norm, n_eff,
                        year=year, save_pdf=save_pdf)

    status, cov, edm = result.status(), result.covQual(), result.edm()
    parked = [k for k, v in nom.items() if near_bound(v)]
    ok = (status == 0 and cov == 3 and not parked)

    print(f"  mu     = {nom['mu'].getVal():7.3f} +/- {nom['mu'].getError():.3f}")
    print(f"  sigma  = {nom['sigma'].getVal():7.3f} +/- {nom['sigma'].getError():.3f}")
    print(f"  norm   = {norm:10.2f}   (N_eff = {n_eff:.1f})")
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
                f"(cb_mu + MH - {MH_REF})*(1 + {bundle['cfg']['scale_unc']}*CMS_scale_m)",
            "formula_sigma":
                f"cb_sigma*(1 + {bundle['cfg']['res_unc']}*CMS_res_m)",
            **bundle["cfg"],
        },
        "normalization": {
            "integral": norm,
            "effective_entries": n_eff,
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
                               f"{pdf.GetName()}_norm", norm)
    norm_var.setConstant(True)

    if freeze:
        for v in nom.values():          # shape from MC is now fixed ...
            v.setConstant(True)
        for v in nuis.values():         # ... nuisances released for Combine
            v.setConstant(False)

    return pdf, bundle, norm_var, record


def run_category(category, binMVA, year, modes, plotdir, wsdir, freeze,
                 save_pdf=False):
    from prepareFits import getHisto   # imported late: needs snapshot files

    tag_cat = f"{category}_{binMVA}_{year}"
    x = ROOT.RooRealVar(f"mh{category}", "m_{#mu#mu}", XLOW, XHIGH)
    x.setRange("full", XLOW, XHIGH)

    w = ROOT.RooWorkspace("w", "workspace")
    nbins = int((XHIGH - XLOW) * BINS_PER_GEV)

    # plots go to <plotdir>/<category>/ so a run over all categories does not
    # dump ~30 files into one directory
    catdir = os.path.join(plotdir, category)

    dump = {
        "provenance": {
            "written_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "git_commit": git_commit(),
            "category": category, "bin": binMVA, "year": year,
            "fit_range": [XLOW, XHIGH],
            "n_bins": nbins,
            "bins_per_gev": BINS_PER_GEV,
        },
        "fits": {},
    }

    keep = []   # hold references so PyROOT does not collect mid-loop
    for sig in modes:
        print(f"\n--- {category} {binMVA} {sig} ---")
        hist = getHisto(nbins, XLOW, XHIGH, False,
                        category, year, True, binMVA, sig)

        out = fit_one(x, hist, f"{tag_cat}_{sig}",
                      f"{category} {binMVA} {sig}", catdir, freeze,
                      year=year, save_pdf=save_pdf)
        if out is None:
            dump["fits"][sig] = {"skipped": "empty histogram"}
            continue
        pdf, bundle, norm_var, record = out
        keep.append((pdf, bundle, norm_var))

        getattr(w, "import")(pdf)
        getattr(w, "import")(norm_var)
        dump["fits"][sig] = record

    os.makedirs(wsdir, exist_ok=True)
    
    # Name by the modes being fit, so single-mode runs don't clobber each
    # other. When the full default set for the category is run, keep the
    # plain per-category name (that's the combined workspace Combine wants).
    if set(modes) == set(MC_LIST[category]):
        out_tag = tag_cat                       # e.g. ggHcat_incl_Run3
    else:
        out_tag = f"{tag_cat}_{'_'.join(modes)}" # e.g. ggHcat_incl_Run3_ggH

    wsfile = f"{wsdir}/Signal_{out_tag}_workspace.root"
    w.writeToFile(wsfile)

    jsonfile = f"{wsdir}/Signal_{out_tag}_params.json"
    with open(jsonfile, "w") as fh:
        json.dump(dump, fh, indent=2, sort_keys=True)

    print(f"\nwrote {wsfile}")
    print(f"wrote {jsonfile}")

    # Sanity check: the only floating variables should be the nuisances.
    floating = [v.GetName() for v in w.allVars()
                if not v.isConstant() and not v.GetName().startswith("mh")]
    print(f"floating in workspace: {floating or '(none)'}")

    n_bad = sum(1 for r in dump["fits"].values()
                if isinstance(r.get("fit_quality"), dict)
                and not r["fit_quality"]["ok"])
    if n_bad:
        print(f"WARNING: {n_bad} fit(s) flagged for review in {tag_cat}")


# ---------------------------------------------------------------------------
# selftest
# ---------------------------------------------------------------------------

TRUTH = {"mu": 125.10, "sigma": 2.00, "alphaL": 1.30,
         "nL": 5.00, "alphaR": 1.60, "nR": 8.00}


def closure(plotdir, n_gen=500000, tag="selftest", verbose=True):
    """Generate from a RooCrystalBall with known parameters and fit it back.

    Returns dict of {param: (fit, error, truth)} plus 'chi2_ndf', 'ok',
    or None if the fit failed.

    Generating from the same functional form is the point. A target that is
    not a CB (e.g. a sum of Gaussians) makes chi2/ndf meaningless and
    strongly binning-dependent, so it cannot serve as a pass criterion.
    """
    x = ROOT.RooRealVar(f"mh_toy_{tag}", "m_{#mu#mu}", XLOW, XHIGH)
    x.setRange("full", XLOW, XHIGH)
    x.setBins(int((XHIGH - XLOW) * BINS_PER_GEV))

    # Bounded and constant: RooRealVar(name, title, value) alone gets
    # [-inf, inf], and RooCrystalBall warns that sigma/alpha/n must be > 0.
    tv = {}
    for k, v in TRUTH.items():
        lo, hi = (v - 5.0, v + 5.0) if k == "mu" else (1e-3, 5.0 * v)
        tv[k] = ROOT.RooRealVar(f"true_{k}_{tag}", k, v, lo, hi)
        tv[k].setConstant(True)

    pdf_true = ROOT.RooCrystalBall(
        f"cb_true_{tag}", "cb_true", x,
        tv["mu"], tv["sigma"], tv["alphaL"], tv["nL"],
        tv["alphaR"], tv["nR"])

    ds = pdf_true.generateBinned(ROOT.RooArgSet(x), n_gen)
    h = ds.createHistogram(f"h_toy_{tag}", x)
    h.SetDirectory(0)

    out = fit_one(x, h, tag, f"closure {tag}",
                  os.path.join(plotdir, "closure"), freeze=False)
    if out is None:
        return None

    _, bundle, _, record = out
    nom = bundle["nom"]
    res = {k: (nom[k].getVal(), nom[k].getError(), v) for k, v in TRUTH.items()}
    res["chi2_ndf"] = record["fit_quality"]["chi2_ndf"]
    res["ok"] = record["fit_quality"]["ok"]

    if verbose:
        print("\n  parameter recovery (pull = (fit - truth)/error):")
        for k in TRUTH:
            v, e, t = res[k]
            pull = (v - t) / e if e > 0 else float("inf")
            print(f"    {k:7s} fit={v:8.3f} +/- {e:6.3f}   truth={t:7.3f}   "
                  f"pull={pull:+6.2f}")
    return res


def selftest(plotdir, n_gen=500000):
    print("=== selftest: closure test against a generated CB ===")
    print("  truth: " + "  ".join(f"{k}={v:.3f}" for k, v in TRUTH.items()))

    res = closure(plotdir, n_gen)
    if res is None:
        print("SELFTEST FAILED: fit returned nothing")
        return 1

    def pull(k):
        v, e, t = res[k]
        return abs((v - t) / e) if e > 0 else float("inf")

    # mu / sigma propagate into the limit; the tail parameters are weakly
    # constrained and correlated, so they get a looser tolerance.
    core_ok = all(pull(k) < 3.0 for k in ("mu", "sigma"))
    tail_ok = all(pull(k) < 5.0 for k in ("alphaL", "nL", "alphaR", "nR"))
    chi2_ok = res["chi2_ndf"] < 2.0

    good = core_ok and tail_ok and chi2_ok
    print(f"\n  core (mu, sigma) within 3 sigma : {core_ok}")
    print(f"  tails within 5 sigma            : {tail_ok}")
    print(f"  chi2/ndf < 2                    : {chi2_ok}  ({res['chi2_ndf']:.2f})")
    print(f"\nselftest {'PASSED' if good else 'FAILED'}")
    return 0 if good else 1


# ---------------------------------------------------------------------------

def main():
    global BINS_PER_GEV          # must precede any use of the name below
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--category", nargs="+", default=["VBFcat"],
                    choices=sorted(MC_LIST))
    ap.add_argument("-b", "--bins", nargs="+", default=["bdt0"],
                    help="BDT bin label(s): bdt0 bdt1 bdt2 incl")
    ap.add_argument("-s", "--sig", nargs="+", default=None,
                    help="production mode(s); default = all for the category")
    ap.add_argument("-y", "--year", default="Run3")
    ap.add_argument("--wsdir", default="WS_LOCAL")
    ap.add_argument("--plotdir",
                    default=os.path.expanduser(
                        "~/public_html/HmumuFits/signal_fits"),
                    help="plots land in <plotdir>/<category>/ (default "
                         "~/public_html/HmumuFits/signal_fits)")
    ap.add_argument("--pdf", action="store_true",
                    help="also write PDF alongside the PNG")
    ap.add_argument("--no-freeze", action="store_true",
                    help="leave shape parameters floating in the workspace")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--bins-per-gev", type=float, default=None,
                    help=f"override binning (default {BINS_PER_GEV})")
    ap.add_argument("--ngen", type=int, default=500000,
                    help="events generated for --selftest")
    args = ap.parse_args()

    if args.bins_per_gev:
        BINS_PER_GEV = args.bins_per_gev
        print(f"binning override: {BINS_PER_GEV} bins/GeV "
              f"({int((XHIGH - XLOW) * BINS_PER_GEV)} bins)")

    if args.selftest:
        raise SystemExit(selftest(args.plotdir, args.ngen))

    for cat in args.category:
        modes = args.sig if args.sig else MC_LIST[cat]
        for b in args.bins:
            run_category(cat, b, args.year, modes,
                         args.plotdir, args.wsdir, not args.no_freeze,
                         save_pdf=args.pdf)


if __name__ == "__main__":
    setup_style()
    main()
