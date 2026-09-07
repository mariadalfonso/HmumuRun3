"""
SummaryPlotsShapes.py -- shape-comparison plots from Hmm.py snapshots,
restricted to the SR sideband (m_mumu in [110,120] or [130,150] GeV).

Companion to SummaryPlots.py (which shows absolute stacked yields). This
script answers a different question -- "what does the SHAPE of each
dominant process look like, on top of each other" -- by normalizing each
selected histogram to unit integral and overlaying them as unfilled line
histograms (no stack, no ratio panel).

Every histogram (background, signal, AND data) is restricted to the SAME
SR-sideband window as SummaryPlots.py's SR_sideband region -- deliberately
excluding the innermost [120,130] GeV core, which is why data is always
shown unblinded here: the restriction itself is the safety mechanism, not
a separate blind cut on top of it (no --blind/--unblind flag).

There are 6 registered backgrounds and 5 signal processes, plus data --
overlaying all 12 curves on one plot is unreadable. So only a curated
subset is drawn: the top --n-bkg backgrounds and top --n-sig signals BY
INTEGRAL (within the SR-sideband window) for that specific variable/
category (default 2 + 2), rather than a hardcoded "DY is always dominant"
assumption, since dominance can shift between categories and variables.

Plots every variable SummaryPlots.py actively draws for the given category
(see get_active_vars() below) -- no variable list to type.

Usage:
    python SummaryPlotsShapes.py ggHcat 2024
    python SummaryPlotsShapes.py VBFcat 2024 --n-bkg 3 --n-sig 1
"""

import ROOT
import os
import sys
import argparse
import getpass

from LoadTree import loadTree
from prepareHisto import getHisto
import plot_style
from plot_style import lumis
import plot_vars

plot_style.setup_style()

# ---------------------------------------------------------------------------
# Configuration / CLI
# ---------------------------------------------------------------------------

CATEGORIES = ["VBFcat", "ggHcat", "VLcat", "TTLcat", "TTHcat",
              "VHcat", "Zinvcat"]

DEFAULT_INDIR  = f"/work/submit/{getpass.getuser()}/HmumuRun3/ROOTFILES/"
DEFAULT_OUTDIR = f"/home/submit/{getpass.getuser()}/public_html/HmumuRun3/"

# SR-sideband mass restriction: m_mumu in [110,120] or [130,150] GeV,
# deliberately excluding the innermost [120,130] core. Same definition as
# SummaryPlots.py's REGIONS["SR_sideband"]["region_filter"] -- kept in sync
# by hand since the two scripts don't currently share a config module.
SR_SIDEBAND_FILTER = ("((HiggsCandCorrMass>=110 && HiggsCandCorrMass<=120) || "
                       "(HiggsCandCorrMass>=130 && HiggsCandCorrMass<=150))")

# Display names for the mass-range label's category line. Same mapping as
# SummaryPlots.py's CATEGORY_LABELS -- kept in sync by hand.
CATEGORY_LABELS = {
    "VBFcat":  "VBF cat.",
    "ggHcat":  "ggF cat.",
    "VLcat":   "VH-lep cat.",
    "TTLcat":  "ttH-lep cat.",
    "TTHcat":  "ttH-had cat.",
    "VHcat":   "VH-had cat.",
    "Zinvcat": "Zinv cat.",
}

BKG_PROCS = ["hTT2L", "hTop", "hZg", "hVV", "hEWK", "hDY"]
SIG_PROCS = ["hWH", "hTTH", "hZH", "hVBFH", "hggH"]

BKG_LABELS = {
    "hDY":   "DY QCD",
    "hEWK":  "DY+jets (EWK)",
    "hVV":   "VV(V)",
    "hTT2L": "t#bar{t} 2l",
    "hTop":  "Top (1l, tW/tZq, ttV/4t)",
    "hZg":   "H#rightarrowZ#gamma + jets",
}
SIG_LABELS = {
    "hTTH":  "ttH",
    "hWH":   "WH",
    "hZH":   "ZH",
    "hVBFH": "VBF H",
    "hggH":  "ggH",
}

# Colour palette local to this script only -- SummaryPlots.py's stacked
# plots keep their original colours (plot_style.PROCESS_COLORS), unchanged.
#
# Colour is assigned by RANK (0 = most dominant, by integral, for THIS
# variable/category), not by fixed process identity. Rank 0 = the two lists'
# first entry, rank 1 = the second. If more than 2 backgrounds/signals are
# requested (--n-bkg/--n-sig > 2), extra ones reuse the last colour.
BKG_COLORS = [ROOT.kAzure + 7, ROOT.kAzure + 2]  # rank 0, rank 1
SIG_COLORS = [ROOT.kRed + 2, ROOT.kRed - 4]      # rank 0, rank 1

# x-axis title offset for THIS script's standalone canvas. Deliberately
# separate from plot_style.RATIO_X_TITLE_OFFSET (=1.0), which is tuned for
# SummaryPlots.py's ratio pad and shouldn't be touched from here -- a bigger
# value here to push the title further from the tick labels.
X_TITLE_OFFSET_SHAPES = 1.4
Y_TITLE_OFFSET_SHAPES = 2.0

# Legend position local to this script (deliberately not touching
# plot_style.LEGEND_POS, which is shared with SummaryPlots.py's stacked
# plots). Back to a compact 2-column box -- one row per process now (no more
# central/sideband pairing to keep adjacent), so it doesn't need to be tall.
MAIN_LEGEND_POS = (0.60, 0.68, 0.93, 0.90)


def _bkg_color(rank):
    return BKG_COLORS[min(rank, len(BKG_COLORS) - 1)]


def _sig_color(rank):
    return SIG_COLORS[min(rank, len(SIG_COLORS) - 1)]


def parse_args():
    p = argparse.ArgumentParser(
        description="Shape-comparison plots (unit-normalized), restricted to "
                     "the SR sideband (m_mumu in [110,120] or [130,150] GeV)."
    )
    p.add_argument("category", choices=CATEGORIES,
                    help="analysis category (snapshot subfolder)")
    p.add_argument("year", help="data-taking year, e.g. 2024, 2025, 12022 ...")
    p.add_argument("-i", "--indir", default=DEFAULT_INDIR,
                    help="input directory with snapshot ROOT files (default: %(default)s)")
    p.add_argument("-o", "--outdir", default=DEFAULT_OUTDIR,
                    help="output directory for PNG plots (default: %(default)s)")
    p.add_argument("--n-bkg", type=int, default=2,
                    help="how many dominant backgrounds to show, ranked by "
                         "integral (default: 2)")
    p.add_argument("--n-sig", type=int, default=2,
                    help="how many dominant signals to show, ranked by "
                         "integral (default: 2)")
    p.add_argument("--log", action="store_true",
                    help="use log scale on the y-axis for all plots (default: linear)")
    p.add_argument("--no-data", action="store_true",
                    help="don't overlay the data shape")
    return p.parse_args()


args = parse_args()

category = args.category
year = "_" + args.year
dirLOCAL_ = args.indir if args.indir.endswith("/") else args.indir + "/"
baseOutDir = args.outdir if args.outdir.endswith("/") else args.outdir + "/"
# Nest under <category><year>/SR_sideband/shapes/, as a sibling of
# SummaryPlots.py's SR_sideband/<group>/ trees under the same category+year
# folder, since this script now exclusively covers that region.
myOutDir = os.path.join(baseOutDir, f"{category}{year}", "shapes", "SR_sideband") + "/"
os.makedirs(myOutDir, exist_ok=True)

if year not in lumis:
    sys.exit(f"ERROR: year '{args.year}' not found in lumis table "
              f"(known: {[k.lstrip('_') for k in lumis if k.startswith('_')]})")

mytree = ROOT.TChain('events')
mytree = loadTree(mytree, dirLOCAL_, category, year)


def get_active_vars(category):
    """The variable names SummaryPlots.py actually plots for `category`.

    Mirrors its GROUP_FUNCS dispatch -- the "mass", "mva", "category", and
    "muons" groups are the only ones in ALL_GROUPS (run by default / via
    --plots all); the per-category reference functions (plotVBF, plotVHlep,
    etc.) are NOT part of the default dispatch and are mostly commented out,
    so they contribute no additional active variables today.

    NOTE: this list is maintained BY HAND to match SummaryPlots.py's active
    plot(...) calls -- there is currently no shared source of truth between
    the two scripts. If you add/remove an active plot(...) call in
    SummaryPlots.py's draw_mass/draw_mva/draw_category/plotMuons, update this
    function to match, or the shape script will silently fall out of sync.
    """
    names = ["dimu_mass", "mva"]

    if category == "VLcat":
        names.append("category_vlcat")
    elif category == "TTLcat":
        names.append("category_ttlcat")
    elif category == "TTHcat":
        names.append("category_tthcat")

    names += [
        "muon1_pt", "muon2_pt", "muon1_eta", "muon2_eta",
        "dimuon_pt", "dimuon_eta", "dimuon_rapidity",
        "muon1_norm_pt", "muon2_norm_pt", "costhetacs", "phistarcs",
    ]
    if category in ("VLcat", "TTLcat", "TTHcat"):
        names += ["muon1_sip3d", "muon2_sip3d"]
    if category in ("ggHcat", "TTHcat"):
        names.append("njets")
    names.append("deta_muons")

    return names


# ---------------------------------------------------------------------------
# Drawing
# ---------------------------------------------------------------------------

def _pick_dominant(hists, proc_list, n):
    """Return up to n (name, hist, integral) tuples from proc_list, ranked
    by integral descending, skipping empty/missing ones."""
    scored = []
    for name in proc_list:
        h = hists.get(name)
        if not h:
            continue
        integral = h.Integral()
        if integral > 0:
            scored.append((name, h, integral))
    scored.sort(key=lambda t: t[2], reverse=True)
    return scored[:n]


def plot_shape(varname):

    nbin, low, high = plot_vars.get_binning(varname)
    titleX = plot_vars.get_xlabel(varname)

    # ONE getHisto() call, restricted to the SR-sideband window, used for
    # backgrounds and data (both dominance ranking and the curves drawn).
    # blind=False unconditionally: the restriction itself is what makes
    # unblinded data safe here (it already excludes the innermost [120,130]
    # core), matching SummaryPlots.py's SR_sideband region.
    listHisto = getHisto(mytree, category, varname, year, nbin, low, high,
                          blind=False, region_filter=SR_SIDEBAND_FILTER)
    if not listHisto:
        print(f"   -> skipped '{varname}': input not available")
        return

    hists = {obj.name: obj.hOBJ.GetValue() for obj in listHisto}

    # Signal gets a SEPARATE, always-unrestricted call: a narrow resonance's
    # signal concentrates almost entirely in [120,130], the exact core the
    # SR-sideband window excludes, so region-restricting it here would leave
    # only a tiny, misleading sliver of the true predicted shape/yield.
    # Ranking (chosen_sig) is also based on this unrestricted signal, not the
    # SR-sideband-truncated version.
    sig_listHisto = getHisto(mytree, category, varname, year, nbin, low, high, blind=False)
    sig_hists = {obj.name: obj.hOBJ.GetValue() for obj in sig_listHisto} if sig_listHisto else {}

    chosen_bkg = _pick_dominant(hists, BKG_PROCS, args.n_bkg)
    chosen_sig = _pick_dominant(sig_hists, SIG_PROCS, args.n_sig)
    print(f"[{varname}] dominant backgrounds: {[n for n, _, i in chosen_bkg]}")
    print(f"[{varname}] dominant signals:     {[n for n, _, i in chosen_sig]}")

    show_data = (not args.no_data) and hists.get("hData") is not None \
                and hists["hData"].Integral() > 0

    c = ROOT.TCanvas("c", "", plot_style.CANVAS_W, plot_style.CANVAS_H)
    # cms_label() sizes the CMS/Internal/lumi text as a fraction of *this
    # pad's* height. SummaryPlots.py draws it on pad1, a sub-pad that's only
    # (1-RATIO_SPLIT) of the full canvas height; this script draws directly
    # on the full canvas, so reusing PAD_TOP_MARGIN as-is renders the label
    # ~1/(1-RATIO_SPLIT) times bigger than in the stacked plots -- crowding
    # the top margin band and overlapping the frame. Scale the top margin so
    # the label ends up the same ABSOLUTE size in both scripts.
    shapes_top_margin = plot_style.PAD_TOP_MARGIN * (1 - plot_style.RATIO_SPLIT)
    c.SetTopMargin(shapes_top_margin)
    c.SetBottomMargin(plot_style.PAD_BOTTOM_MARGIN)
    c.SetLeftMargin(plot_style.PAD_LEFT_MARGIN)
    c.SetRightMargin(plot_style.PAD_RIGHT_MARGIN)
    c.SetLogy(args.log)

    legend = ROOT.TLegend(*MAIN_LEGEND_POS)
    legend.SetNColumns(plot_style.LEGEND_NCOLUMNS)
    legend.SetColumnSeparation(plot_style.LEGEND_COLUMN_SEP)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(plot_style.LEGEND_TEXT_SIZE)
    legend.SetTextAlign(12)

    drawn = []
    ymax = 0.0

    def _normalize(name, h, label, color, is_data=False):
        nonlocal ymax
        hc = h.Clone(f"{name}_shape")
        integral = hc.Integral()
        if integral > 0:
            hc.Scale(1.0 / integral)
        hc.SetFillStyle(0)
        hc.SetLineWidth(4)
        hc.SetLineColor(color)
        if is_data:
            hc.SetMarkerStyle(20)
            hc.SetMarkerSize(1.0)
            hc.SetMarkerColor(color)
        ymax = max(ymax, hc.GetMaximum())
        drawn.append((hc, label, is_data))

    for rank, (name, h, integral) in enumerate(chosen_bkg):
        _normalize(name, h, BKG_LABELS.get(name, name), _bkg_color(rank))
    for rank, (name, h, integral) in enumerate(chosen_sig):
        _normalize(name, h, SIG_LABELS.get(name, name), _sig_color(rank))
    if show_data:
        _normalize("hData", hists["hData"], "Data", ROOT.kBlack, is_data=True)

    if not drawn:
        print(f"[{varname}] nothing to draw (no non-empty histograms) -- skipping")
        return

    for i, (hc, label, is_data) in enumerate(drawn):
        if args.log:
            hc.SetMaximum(ymax * 10)
            hc.SetMinimum(max(ymax * 1e-6, 1e-8))
        else:
            hc.SetMaximum(ymax * 1.3)
            hc.SetMinimum(0.0)
        hc.GetXaxis().SetTitle(titleX)
        # The inherited relative-mode default (title size as a fraction of
        # pad height, offset 1.1) was tuned for pad1 in the stacked plot,
        # which never actually shows its own x-axis (it's suppressed there;
        # only the ratio pad's x-axis is visible, using pixel-mode fonts).
        # Applied directly to this script's full-height standalone canvas,
        # that combination doesn't reliably fit inside the bottom margin and
        # pushes the title below the canvas edge. Reuse the same proven
        # pixel-mode fix already used for the ratio pad's x-axis instead.
        hc.GetXaxis().SetTitleFont(plot_style.RATIO_FONT_ABS)
        hc.GetXaxis().SetTitleSize(plot_style.RATIO_TITLE_SIZE_PX)
        hc.GetXaxis().SetTitleOffset(X_TITLE_OFFSET_SHAPES)
        hc.GetYaxis().SetTitle("Normalized Entries")
        # Same root cause as the x-axis fix above: the inherited relative-mode
        # title size (a fraction of pad HEIGHT) was tuned for pad1's reduced
        # sub-pad height in the stacked plot, not this script's full-height
        # standalone canvas -- so it renders oversized here regardless of the
        # offset. Switch to the same proven pixel-mode fix instead of keeping
        # this a fraction.
        hc.GetYaxis().SetTitleFont(plot_style.RATIO_FONT_ABS)
        hc.GetYaxis().SetTitleSize(plot_style.RATIO_TITLE_SIZE_PX)
        hc.GetYaxis().SetTitleOffset(Y_TITLE_OFFSET_SHAPES)
        draw_opt = "PE" if is_data else "HIST"
        if i > 0:
            draw_opt += " SAME"
        hc.Draw(draw_opt)

    for hc, label, is_data in drawn:
        legend.AddEntry(hc, label, "pe" if is_data else "l")
    legend.Draw()

    plot_style.cms_label(c, year)

    # Mass-range label, matching SummaryPlots.py's SR_sideband region label
    # exactly: unconditional for every variable, since the SR-sideband
    # restriction (not a per-variable blind check) is what makes unblinded
    # data safe here.
    cat_label = CATEGORY_LABELS.get(category, category)
    mass_label_lines = [
        "H #rightarrow #mu#mu",
        f"{cat_label} ({year.lstrip('_')})",
        "m_{#mu#mu}#in[110, 120]#cup[130, 150] GeV",
    ]
    _masslabel = ROOT.TLatex()
    _masslabel.SetNDC()
    _masslabel.SetTextFont(42)
    _masslabel.SetTextSize(0.035)
    y0, dy = 0.85, 0.045
    for i, line in enumerate(mass_label_lines):
        _masslabel.DrawLatex(plot_style.PAD_LEFT_MARGIN + 0.05, y0 - i * dy, line)

    outpath = f"{myOutDir}{varname}_{category}{year}_Shape.png"
    c.SaveAs(outpath)
    print(f"   -> {outpath}")


def main():
    active_vars = get_active_vars(category)
    print(f"[shapes] category={category} year={year.lstrip('_')} "
          f"n_bkg={args.n_bkg} n_sig={args.n_sig} "
          f"scale={'log' if args.log else 'linear'} "
          f"data={'off' if args.no_data else 'on'}")
    print(f"[shapes] region: SR_sideband -- m_mumu in [110,120] or [130,150] GeV "
          f"(unconditional, data always unblinded)")
    print(f"[shapes] output dir: {myOutDir}")
    print(f"[shapes] plotting {len(active_vars)} variables: {active_vars}")
    for v in active_vars:
        plot_shape(v)


if __name__ == "__main__":
    main()
