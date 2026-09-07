"""
plot_style.py — central plotting style for the HmumuRun3 analysis.

Single source of truth for everything visual: global ROOT style, the
process→colour map, the luminosity table, the CMS label text, and canvas/pad
geometry. Both prepareHisto.py and SummaryPlots.py import from here instead of
defining style inline.

cmsstyle integration
--------------------
If the official `cmsstyle` package is available it is used for the global style
and the CMS label (setCMSStyle / SetExtraText / SetLumi / SetEnergy). If it is
not installed, everything falls back to the hand-rolled equivalents so the
scripts still run. Check availability with USE_CMSSTYLE.

To adopt more of cmsstyle incrementally (cmsDiCanvas, cmsLeg, cmsDraw, the
Petroff palettes p6/p8/p10), do it here — the rest of the code only sees the
helpers this module exposes.
"""

import ROOT

# --- optional cmsstyle -------------------------------------------------------
try:
    import cmsstyle as CMS
    USE_CMSSTYLE = True
except ImportError:
    CMS = None
    USE_CMSSTYLE = False


# ---------------------------------------------------------------------------
# Luminosity (fb^-1), keyed by the internal "_<year>" convention
# ---------------------------------------------------------------------------
lumis = {
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


# ---------------------------------------------------------------------------
# CMS label
# ---------------------------------------------------------------------------
CMS_TEXT = "CMS"
EXTRA_TEXT = "Internal"          # e.g. "Preliminary", "Simulation", "Internal"
SQRT_S_TEV = 13.6                # centre-of-mass energy

# Official CMS label size formulas (from the standard CMS_lumi macro used
# across CMS tutorials/tdrstyle). Sizes are `FRAC * t`, where t = pad top
# margin. These are the exact CMS-guideline multipliers — no extra scaling
# is applied; the label grows/shrinks only via PAD_TOP_MARGIN below.
CMS_TEXT_SIZE_FRAC = 0.75
EXTRA_OVER_CMS_TEXT_SIZE = 0.76      # extra text size = this * cmsTextSize
LUMI_TEXT_SIZE_FRAC = 0.6


# ---------------------------------------------------------------------------
# Colour palette (RGB). Kept as the analysis's current hand-tuned values so
# plots are unchanged. To switch to CMS Petroff colours, replace the values in
# PROCESS_COLORS below with e.g. CMS.p10.kBlue (requires USE_CMSSTYLE).
# ---------------------------------------------------------------------------
_RGB = {
    "darkBlue":   (0, 0, 255),
    "redDark":    (191, 34, 41),
    "redMed":     (237, 41, 57),
    "redLight":   (255, 82, 82),
    "orange":     (255, 204, 153),
    "orangeDark": (255, 150, 79),
    "gray":       (136, 139, 141),
    "azure":      (100, 192, 232),
    "azureDark":  (96, 147, 172),
    "green":      (144, 238, 144),
    "greenDark":  (98, 172, 141),
    "gold":       (212, 175, 55),
    "violet":     (181, 100, 227),
}

# process key -> palette name. Single source of truth for which process is
# which colour (previously an implicit zip() in prepareHisto.py).
PROCESS_COLORS = {
    "hDY":    "azure",
    "hTT2L":  "green",
    "hTop":   "greenDark",
    "hVV":    "orange",
    "hEWK":   "gold",
    "hVBFH":  "redMed",
    "hggH":   "redDark",
    "hWH":    "redLight",
    "hZH":    "redLight",
    "hTTH":   "redDark",
    "hZg":    "orangeDark",
}


def color(name):
    """ROOT colour index for a palette name (e.g. 'azure')."""
    return ROOT.TColor.GetColor(*_RGB[name])


def process_color(hname):
    """ROOT colour index for a process histogram name (e.g. 'hDY')."""
    return color(PROCESS_COLORS[hname])


# ---------------------------------------------------------------------------
# Canvas / pad geometry — official CMS values (from cmsstyle.C setCMSStyle)
# ---------------------------------------------------------------------------
CANVAS_W = 600
CANVAS_H = 600
RATIO_SPLIT = 0.30         # lower (ratio) pad occupies bottom 30% (was 0.2)

# Official CMS TStyle values (cmsstyle.C, setCMSStyle). Used by apply_cms_tstyle()
# and available to plotting code that styles axes directly.
CMS_FONT = 42
PAD_TOP_MARGIN = 0.08      # CMS official is 0.05; enlarged so the CMS/Internal/
                            # lumi label (sized as a fraction of this) reads
                            # bigger. This is the single knob for label size.
PAD_BOTTOM_MARGIN = 0.13
PAD_LEFT_MARGIN = 0.16
PAD_RIGHT_MARGIN = 0.05    # CMS official is 0.02; widened — the default left
                            # content (axis numbers, legend) too close to the
                            # edge.
AXIS_TITLE_SIZE = 0.06     # official (was 0.045 in the old hand-rolled version)
AXIS_LABEL_SIZE = 0.05     # official (was 0.04)
AXIS_LABEL_OFFSET = 0.012
X_TITLE_OFFSET = 1.1
Y_TITLE_OFFSET = 1.35      # official (was 1.1)
TICK_LENGTH = 0.03
N_DIVISIONS = 510

# ratio pad text: ABSOLUTE (pixel) sizes, so text renders consistently
# regardless of the pad's small height, instead of fragile pad-relative
# fractions that need retuning whenever pad geometry changes (font code
# ending in "3" = size given in pixels, per CMS/ATLAS convention).
RATIO_FONT_ABS = 43
RATIO_TITLE_SIZE_PX = 24
RATIO_LABEL_SIZE_PX = 20
RATIO_X_TITLE_OFFSET = 1.0    # pixel-mode offset: ~1.0 is the usual working
                               # value; the old 3.2 (tuned for fraction mode)
                               # pushed the title past the pad edge -> clipped
RATIO_Y_TITLE_OFFSET = 1.4

# data marker (official CMS marker style is 20; sizes kept from analysis)
DATA_MARKER_STYLE = 20
DATA_MARKER_SIZE = 1.2
DATA_LINE_WIDTH = 2

# legend
LEGEND_POS = (0.5, 0.60, 0.95, 0.88)   # widened: the old (0.45,0.60,0.95,0.88)
                                          # box was too narrow for the longest
                                          # entry, causing column overlap
LEGEND_NCOLUMNS = 2
LEGEND_TEXT_SIZE = 0.032                 # slightly smaller for extra margin
LEGEND_COLUMN_SEP = 0.02                 # gap between the two columns (ROOT
                                          # default is wider); smaller pulls
                                          # the left column closer to the right


# ---------------------------------------------------------------------------
# Global style setup — call once at program start
# ---------------------------------------------------------------------------
def apply_cms_tstyle():
    """Build and apply the official CMS TStyle, translated verbatim from
    cmsstyle.C (setCMSStyle). Used when the cmsstyle package is not importable
    so plots still match the official style. Returns the TStyle."""
    # Delete any pre-existing style of the same name first — creating two
    # TStyle objects named "cmsStyle" corrupts ROOT's registry (the C++
    # setCMSStyle does the same delete-then-create).
    existing = ROOT.gROOT.GetStyle("cmsStyle")
    if existing:
        ROOT.gROOT.GetListOfStyles().Remove(existing)

    st = ROOT.TStyle("cmsStyle", "Style for P-CMS")
    ROOT.gROOT.SetStyle(st.GetName())
    ROOT.gROOT.ForceStyle(False)

    # Canvas
    st.SetCanvasBorderMode(0)
    st.SetCanvasColor(ROOT.kWhite)
    st.SetCanvasDefH(CANVAS_H)
    st.SetCanvasDefW(CANVAS_W)
    st.SetCanvasDefX(0)
    st.SetCanvasDefY(0)
    st.SetPadBorderMode(0)
    st.SetPadColor(ROOT.kWhite)
    st.SetPadGridX(False)
    st.SetPadGridY(False)
    st.SetGridColor(0)
    st.SetGridStyle(3)
    st.SetGridWidth(1)

    # Frame
    st.SetFrameBorderMode(0)
    st.SetFrameBorderSize(1)
    st.SetFrameFillColor(0)
    st.SetFrameFillStyle(0)
    st.SetFrameLineColor(1)
    st.SetFrameLineStyle(1)
    st.SetFrameLineWidth(1)

    # Histogram
    st.SetHistLineColor(1)
    st.SetHistLineStyle(0)
    st.SetHistLineWidth(1)
    st.SetEndErrorSize(2)
    st.SetMarkerStyle(DATA_MARKER_STYLE)

    # Fit/function
    st.SetOptFit(1)
    st.SetFitFormat("5.4g")
    st.SetFuncColor(2)
    st.SetFuncStyle(1)
    st.SetFuncWidth(1)

    st.SetOptDate(0)

    # Statistics box
    st.SetOptFile(0)
    st.SetOptStat(0)
    st.SetStatColor(ROOT.kWhite)
    st.SetStatFont(CMS_FONT)
    st.SetStatFontSize(0.025)
    st.SetStatTextColor(1)
    st.SetStatFormat("6.4g")
    st.SetStatBorderSize(1)
    st.SetStatH(0.1)
    st.SetStatW(0.15)

    # Margins
    st.SetPadTopMargin(PAD_TOP_MARGIN)
    st.SetPadBottomMargin(PAD_BOTTOM_MARGIN)
    st.SetPadLeftMargin(PAD_LEFT_MARGIN)
    st.SetPadRightMargin(PAD_RIGHT_MARGIN)

    # Global title
    st.SetOptTitle(0)
    st.SetTitleFont(CMS_FONT)
    st.SetTitleColor(1)
    st.SetTitleTextColor(1)
    st.SetTitleFillColor(10)
    st.SetTitleFontSize(0.05)

    # Axis titles
    st.SetTitleColor(1, "XYZ")
    st.SetTitleFont(CMS_FONT, "XYZ")
    st.SetTitleSize(AXIS_TITLE_SIZE, "XYZ")
    st.SetTitleXOffset(X_TITLE_OFFSET)
    st.SetTitleYOffset(Y_TITLE_OFFSET)

    # Axis labels
    st.SetLabelColor(1, "XYZ")
    st.SetLabelFont(CMS_FONT, "XYZ")
    st.SetLabelOffset(AXIS_LABEL_OFFSET, "XYZ")
    st.SetLabelSize(AXIS_LABEL_SIZE, "XYZ")

    # Axis
    st.SetAxisColor(1, "XYZ")
    st.SetStripDecimals(True)
    st.SetTickLength(TICK_LENGTH, "XYZ")
    st.SetNdivisions(N_DIVISIONS, "XYZ")
    st.SetPadTickX(1)
    st.SetPadTickY(1)

    # Log plots
    st.SetOptLogx(0)
    st.SetOptLogy(0)
    st.SetOptLogz(0)

    # Postscript
    st.SetPaperSize(20.0, 20.0)
    st.SetHatchesLineWidth(2)
    st.SetHatchesSpacing(1.3)

    st.cd()
    return st


# module-level handle to keep the TStyle alive
_cms_tstyle = None
_style_done = False


def setup_style():
    """Apply the global CMS style once. Prefers the official cmsstyle package;
    if unavailable, applies the faithful translation in apply_cms_tstyle().
    Idempotent — safe (and cheap) to call from more than one module."""
    global _cms_tstyle, _style_done
    if _style_done:
        return
    ROOT.gROOT.SetBatch()
    if USE_CMSSTYLE:
        CMS.setCMSStyle()
        CMS.SetExtraText(EXTRA_TEXT)
        CMS.SetEnergy(SQRT_S_TEV)
    else:
        _cms_tstyle = apply_cms_tstyle()
    _style_done = True


def cms_label(pad, year):
    """Draw the CMS label + lumi/energy on `pad` for the given '_<year>' key.

    With cmsstyle: sets the lumi then calls CMS_lumi (guideline-compliant).
    Without: a faithful hand-rolled version using the official CMS fonts
    (61 for "CMS", 52 for the extra text) and the standard layout — "CMS" +
    extra text top-left, lumi + energy top-right, above the frame.
    Returns the TLatex kept alive (caller need not use it)."""
    lumi = lumis.get(year, 0.0)
    if USE_CMSSTYLE:
        CMS.SetLumi(lumi, unit="fb", run="", round_lumi=2)
        CMS.CMS_lumi(pad, iPosX=11)
        return None

    # fallback: official CMS fonts + layout (cmsTextFont=61, extraTextFont=52)
    pad.cd()
    t = pad.GetTopMargin()
    l = pad.GetLeftMargin()
    r = pad.GetRightMargin()
    y = 1 - t + 0.2 * t          # just above the frame

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(ROOT.kBlack)

    # "CMS" (font 61), left, above frame
    latex.SetTextFont(61)
    latex.SetTextAlign(11)
    latex.SetTextSize(CMS_TEXT_SIZE_FRAC * t)
    latex.DrawLatex(l, y, CMS_TEXT)

    # extra text (font 52) after "CMS"
    if EXTRA_TEXT:
        latex.SetTextFont(52)
        latex.SetTextSize(EXTRA_OVER_CMS_TEXT_SIZE * CMS_TEXT_SIZE_FRAC * t)
        latex.DrawLatex(l + 0.10, y, EXTRA_TEXT)

    # lumi + energy (font 42), right-aligned
    lumi_text = "%.2f fb^{-1} (%.1f TeV)" % (lumi, SQRT_S_TEV)
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.SetTextSize(LUMI_TEXT_SIZE_FRAC * t)
    latex.DrawLatex(1 - r, y, lumi_text)
    return latex


def make_canvas_pads(doLog):
    """Create the canvas + upper (main) and lower (ratio) pads, glued together
    with no vertical gap and a single x-axis shown only on the bottom pad.

    Structure mirrors the original createCanvasPads (outer wrapper pad, then
    pad1/pad2 as its children) to avoid a PyROOT object-ownership/double-free
    issue seen with a flatter pad hierarchy. Only the margins/fill style are
    changed from the original to glue pad1/pad2 together."""
    c = ROOT.TCanvas("c", "", CANVAS_W, CANVAS_H)
    pad = ROOT.TPad("outer_pad", "", 0, 0, 1, 1)
    if doLog:
        pad.SetLogy(1)
    pad.SetTickx(False)
    pad.SetTicky(False)
    pad.Draw()
    pad.cd()

    ydiv = RATIO_SPLIT
    pad1 = ROOT.TPad("upper_pad", "", 0., ydiv, 1., 1.)
    pad1.SetTopMargin(PAD_TOP_MARGIN)
    pad1.SetBottomMargin(0.0)          # glue: no gap at the shared edge
    pad1.SetLeftMargin(PAD_LEFT_MARGIN)
    pad1.SetRightMargin(PAD_RIGHT_MARGIN)
    if doLog:
        pad1.SetLogy(1)

    pad2 = ROOT.TPad("lower_pad", "", 0., 0., 1., ydiv)
    pad2.SetTopMargin(0.0)             # glue: touches the upper pad
    pad2.SetBottomMargin(PAD_BOTTOM_MARGIN / ydiv)
    pad2.SetLeftMargin(PAD_LEFT_MARGIN)
    pad2.SetRightMargin(PAD_RIGHT_MARGIN)

    pad1.Draw()
    pad2.Draw()

    return c, pad1, pad2
