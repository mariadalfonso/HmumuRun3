"""
style.py -- the CMS plot style shared by every script: label text and sizes,
canvas and pad geometry, markers, colours, the category / bin / mode labels,
and the helpers that draw them.

Importing it also sets ROOT up for batch plotting (batch mode, RooFit
messages, Sumw2), as importing sigFit.py used to.
"""

import ROOT

from . import config

ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
# Empty tail bins have zero error, so their pull is undefined
# Supress these erros for EVERY stream
for _i in range(ROOT.RooMsgService.instance().numStreams()):
    ROOT.RooMsgService.instance().getStream(_i).removeTopic(ROOT.RooFit.Plotting)
ROOT.TH1.SetDefaultSumw2(True)

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

# Integrated luminosity: LUMIS and the era groups live in prepareFits.py,
# so the label is always computed from the same eras the files come from.

# canvas / pad geometry (plot_style values)
CANVAS_W, CANVAS_H = 600, 600
RATIO_SPLIT = 0.30
# Dashed guide lines
RATIO_GUIDES = (1.2, 0.8)   # Template/Fit, +/-0.2 around 1
PULL_GUIDES = (2.0, -2.0)   # +/-2 sigma
PAD_TOP_MARGIN = 0.08
PAD_BOTTOM_MARGIN = 0.13
PAD_LEFT_MARGIN = 0.15
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

DATA_MARKER_STYLE = 24
DATA_MARKER_SIZE = 0.8
DATA_LINE_WIDTH = 1

FIT_COLOR_RGB = (237, 41, 57)     # plot_style "redMed"

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


def cms_label(pad, year, full_year=False):
    """CMS + extra text top-left, lumi + energy top-right, above the frame.
    full_year: label the data luminosity instead of the signal MC's."""
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

    lumi = config.get_lumi(year, full_year) if year else 0.0
    if lumi <= 0:
        print(f"warning: no luminosity found for year tag '{year}'")
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.SetTextSize(LUMI_TEXT_SIZE_FRAC * t)
    if lumi > 0:
        txt = "%.1f fb^{-1} (%.1f TeV)" % (lumi, SQRT_S_TEV)
    else:
        txt = "(%.1f TeV)" % SQRT_S_TEV      # e.g. the sigFit.closure test
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
    obj.SetLineWidth(DATA_LINE_WIDTH)
