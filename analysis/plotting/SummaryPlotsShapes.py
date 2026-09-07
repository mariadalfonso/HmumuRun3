"""
SummaryPlotsShapes.py -- shape-comparison plots from Hmm.py snapshots.

Companion to SummaryPlots.py (which shows absolute stacked yields). This
script answers a different question -- "what does the SHAPE of each
dominant process look like, on top of each other" -- by normalizing each
selected histogram to unit integral and overlaying them as unfilled line
histograms (no stack, no ratio panel).

There are 6 registered backgrounds and 5 signal processes, plus data --
overlaying all 12 curves on one plot is unreadable. So only a curated
subset is drawn: the top --n-bkg backgrounds and top --n-sig signals BY
INTEGRAL for that specific variable/category (default 2 + 2), rather than
a hardcoded "DY is always dominant" assumption, since dominance can shift
between categories and variables.

Data is included by default, EXCEPT for a blinded variable (e.g. dimu_mass)
where a partially-blinded shape would be misleading -- pass --unblind to
include it there.

Plots every variable SummaryPlots.py actively draws for the given category
(see get_active_vars() below) -- no variable list to type.

Usage:
    python SummaryPlotsShapes.py ggHcat 2024
    python SummaryPlotsShapes.py ggHcat 2024 --unblind
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
DEFAULT_OUTDIR = f"/home/submit/{getpass.getuser()}/public_html/HmumuRun3/shapes/"

BKG_PROCS = ["hTT2L", "hTop", "hZg", "hVV", "hEWK", "hDY"]
SIG_PROCS = ["hWH", "hTTH", "hZH", "hVBFH", "hggH"]

BKG_LABELS = {
    "hDY":   "DY+jets (QCD)",
    "hEWK":  "DY+jets (EWK)",
    "hVV":   "VV + VVV",
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
BKG_COLORS = [ROOT.kAzure + 1, ROOT.kAzure + 3]  # rank 0, rank 1
SIG_COLORS = [ROOT.kRed + 2, ROOT.kRed - 4]      # rank 0, rank 1

# x-axis title offset for THIS script's standalone canvas. Deliberately
# separate from plot_style.RATIO_X_TITLE_OFFSET (=1.0), which is tuned for
# SummaryPlots.py's ratio pad and shouldn't be touched from here -- a bigger
# value here to push the title further from the tick labels.
X_TITLE_OFFSET_SHAPES = 1.4

# Legend positions local to this script (deliberately not touching
# plot_style.LEGEND_POS, which is shared with SummaryPlots.py's stacked
# plots). Two SEPARATE boxes: MAIN_LEGEND_POS for the per-process entries,
# STYLE_KEY_LEGEND_POS for the solid/dashed line-style key, positioned to
# its left -- so the style key reads as a distinct, general note rather
# than being mixed into the process list. Both raised higher than the
# shared default so they don't crowd the curves near the top of the frame.
MAIN_LEGEND_POS = (0.40, 0.68, 0.97, 0.90)
STYLE_KEY_LEGEND_POS = (0.16, 0.68, 0.38, 0.90)


def _bkg_color(rank):
    return BKG_COLORS[min(rank, len(BKG_COLORS) - 1)]


def _sig_color(rank):
    return SIG_COLORS[min(rank, len(SIG_COLORS) - 1)]


def parse_args():
    p = argparse.ArgumentParser(
        description="Shape-comparison plots (unit-normalized) from Hmm.py snapshots."
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
    p.add_argument("--unblind", action="store_true",
                    help="use unblinded data (needed to see the data shape "
                         "for a blinded variable, e.g. dimu_mass)")
    return p.parse_args()


args = parse_args()

category = args.category
year = "_" + args.year
dirLOCAL_ = args.indir if args.indir.endswith("/") else args.indir + "/"
myOutDir  = args.outdir if args.outdir.endswith("/") else args.outdir + "/"
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

    # Data: if the variable is blinded and --unblind wasn't passed, data is
    # drawn in its blinded (sideband-only) form, labeled clearly, rather than
    # dropped from the plot.
    blind_range = plot_vars.get_blind_range(varname)
    use_blind = blind_range is not None and not args.unblind
    is_sideband = use_blind   # True whenever data is actually blinded here

    listHisto = getHisto(mytree, category, varname, year, nbin, low, high, blind=use_blind)
    if not listHisto:
        print(f"   -> skipped '{varname}': input not available")
        return

    hists = {obj.name: obj.hOBJ.GetValue() for obj in listHisto}

    chosen_bkg = _pick_dominant(hists, BKG_PROCS, args.n_bkg)
    chosen_sig = _pick_dominant(hists, SIG_PROCS, args.n_sig)
    print(f"[{varname}] dominant backgrounds: {[n for n, _, i in chosen_bkg]}")
    print(f"[{varname}] dominant signals:     {[n for n, _, i in chosen_sig]}")

    show_data = (not args.no_data) and hists.get("hData") is not None \
                and hists["hData"].Integral() > 0
    if show_data and is_sideband:
        print(f"[{varname}] data is blinded -- showing sideband-only shape, "
              f"labeled 'Data (sideband)' (pass --unblind for full data)")
    data_label = "Data (sideband)" if is_sideband else "Data"

    # MC (backgrounds + signals): split into Higgs-mass sideband vs mH window,
    # independent of data blinding -- MC can be evaluated in both regions
    # freely, since blinding is a data-only policy. Reuses the SAME mass
    # window bounds as data blinding (plot_vars.get_blind_range("dimu_mass"))
    # as the single source of truth for what "the mH window" means.
    mass_lo, mass_hi = plot_vars.get_blind_range("dimu_mass")
    sideband_filter = f"(HiggsCandCorrMass<{mass_lo} || HiggsCandCorrMass>{mass_hi})"
    window_filter = f"(HiggsCandCorrMass>={mass_lo} && HiggsCandCorrMass<={mass_hi})"

    sideband_list = getHisto(mytree, category, varname, year, nbin, low, high,
                              blind=False, region_filter=sideband_filter)
    window_list = getHisto(mytree, category, varname, year, nbin, low, high,
                            blind=False, region_filter=window_filter)
    hists_sideband = {obj.name: obj.hOBJ.GetValue() for obj in sideband_list} if sideband_list else {}
    hists_window = {obj.name: obj.hOBJ.GetValue() for obj in window_list} if window_list else {}

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

    def _normalize(name, h, label, color, is_data=False, dashed=False, add_legend=True, region=""):
        nonlocal ymax
        tag = f"{name}_shape_{region}" if region else f"{name}_shape"
        hc = h.Clone(tag)
        integral = hc.Integral()
        if integral > 0:
            hc.Scale(1.0 / integral)
        hc.SetFillStyle(0)
        hc.SetLineWidth(4)
        hc.SetLineColor(color)
        if dashed:
            hc.SetLineStyle(2)   # dashed = mH window (backgrounds only)
        if is_data:
            hc.SetMarkerStyle(20)
            hc.SetMarkerSize(1.0)
            hc.SetMarkerColor(color)
        ymax = max(ymax, hc.GetMaximum())
        drawn.append((hc, label, is_data, add_legend))

    # Backgrounds: split into sideband (solid) / mH window (dashed), same
    # colour -- a real shape-comparison question (does the background shape
    # change between the sideband and the signal region?).
    for rank, (name, h, integral) in enumerate(chosen_bkg):
        color = _bkg_color(rank)
        label = BKG_LABELS.get(name, name)
        h_side = hists_sideband.get(name)
        h_win = hists_window.get(name)
        if h_side is not None and h_side.Integral() > 0:
            _normalize(name, h_side, label, color, dashed=True,
                       add_legend=False, region="sideband")
        if h_win is not None and h_win.Integral() > 0:
            _normalize(name, h_win, label, color,
                       add_legend=True, region="window")

    # Signal: a real resonance, essentially all inside the mH window by
    # construction -- there's no meaningful "sideband shape" to compare
    # against, so draw one solid curve, restricted to the window, per signal.
    for rank, (name, h, integral) in enumerate(chosen_sig):
        color = _sig_color(rank)
        label = SIG_LABELS.get(name, name)
        h_win = hists_window.get(name)
        if h_win is not None and h_win.Integral() > 0:
            # solid to match the "solid = mH window/SR" convention below --
            # signal is always in that region by construction, so it should
            # never appear dashed (dashed is reserved for sideband, which
            # only exists for background).
            _normalize(name, h_win, label, color, add_legend=True, region="window")

    if show_data:
        _normalize("hData", hists["hData"], data_label, ROOT.kBlack, is_data=True)

    if not drawn:
        print(f"[{varname}] nothing to draw (no non-empty histograms) -- skipping")
        return

    for i, (hc, label, is_data, add_legend) in enumerate(drawn):
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
        # Normalized-entries tick labels (small decimals, e.g. 0.0500) can be
        # wider than typical integer-count labels, so the inherited default
        # offset (1.35) sometimes isn't enough clearance -- push it out a bit
        # further so the title never overlaps the tick labels.
        hc.GetYaxis().SetTitleOffset(1.6)
        draw_opt = "PE" if is_data else "HIST"
        if i > 0:
            draw_opt += " SAME"
        hc.Draw(draw_opt)

    # ---------------------------------------------------------------------
    # Legend: one entry per process, from its SOLID (mH window/SR) curve
    # only -- dashed (sideband) curves are drawn but never get their own
    # legend entry, keeping one line per process.
    # ---------------------------------------------------------------------
    for hc, label, is_data, add_legend in drawn:
        if add_legend:
            legend.AddEntry(hc, label, "pe" if is_data else "l")
    legend.Draw()

    # ---------------------------------------------------------------------
    # Style key: a SEPARATE, general-purpose legend explaining the
    # solid/dashed line style itself (not tied to any process), so it's
    # visually distinct from the process legend above. The key TLine
    # objects are never drawn on the canvas -- only registered with this
    # legend for their line style.
    # ---------------------------------------------------------------------
    style_legend = ROOT.TLegend(*STYLE_KEY_LEGEND_POS)
    style_legend.SetFillStyle(0)
    style_legend.SetBorderSize(0)
    style_legend.SetTextSize(plot_style.LEGEND_TEXT_SIZE)
    style_legend.SetTextAlign(12)

    legend_keys = []  # keep these TLines alive (PyROOT) until c.SaveAs() below
    key_solid = ROOT.TLine()
    key_solid.SetLineColor(ROOT.kGray + 2)
    key_solid.SetLineWidth(3)
    key_solid.SetLineStyle(1)
    legend_keys.append(key_solid)
    style_legend.AddEntry(key_solid, f"m_{{H}} [{mass_lo:g}, {mass_hi:g}] GeV", "l")

    key_dashed = ROOT.TLine()
    key_dashed.SetLineColor(ROOT.kGray + 2)
    key_dashed.SetLineWidth(3)
    key_dashed.SetLineStyle(2)
    legend_keys.append(key_dashed)
    style_legend.AddEntry(key_dashed, "sideband", "l")

    style_legend.Draw()

    plot_style.cms_label(c, year)

    outpath = f"{myOutDir}{varname}_{category}{year}_Shape.png"
    c.SaveAs(outpath)
    print(f"   -> {outpath}")


def main():
    active_vars = get_active_vars(category)
    print(f"[shapes] category={category} year={year.lstrip('_')} "
          f"n_bkg={args.n_bkg} n_sig={args.n_sig} "
          f"scale={'log' if args.log else 'linear'} "
          f"data={'off' if args.no_data else 'on'} "
          f"unblind={args.unblind}")
    print(f"[shapes] output dir: {myOutDir}")
    print(f"[shapes] plotting {len(active_vars)} variables: {active_vars}")
    for v in active_vars:
        plot_shape(v)


if __name__ == "__main__":
    main()
