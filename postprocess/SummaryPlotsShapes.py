"""
SummaryPlotsShapes.py -- shape-comparison plots, restricted to the SR sideband
(m_mumu in [110,120] or [130,150] GeV).

Companion to SummaryPlots.py (which shows absolute stacked yields). This script
answers a different question -- "what does the SHAPE of each dominant process
look like, overlaid" -- by normalizing each selected histogram to unit integral
and drawing them as unfilled lines (no stack, no ratio panel).

Reads the histograms produced by makeHistos.py; does NO event processing. Run
makeHistos.py first:

    python makeHistos.py ggHcat 2024
    python SummaryPlotsShapes.py ggHcat 2024
    python SummaryPlotsShapes.py VBFcat 2024 --n-bkg 3 --n-sig 1

Backgrounds and data come from the SR_sideband region; SIGNAL comes from the
Unrestricted region, because a narrow resonance concentrates almost entirely in
[120,130] -- the exact core SR_sideband excludes -- so a region-restricted
signal shape would be a misleading sliver. Data is always shown unblinded here:
the SR-sideband restriction is itself the safety mechanism.

There are 6 backgrounds and 5 signals plus data; overlaying all 12 is
unreadable, so only the top --n-bkg backgrounds and top --n-sig signals BY
INTEGRAL for that specific variable/category are drawn (default 2 + 2), rather
than assuming DY always dominates -- dominance shifts between categories and
variables.
"""

import ROOT
import os
import sys
import argparse
import getpass

import plot_style
from plot_style import lumis
import plot_vars
import histo_config as cfg
from SummaryPlots import HistoFile

plot_style.setup_style()

DEFAULT_INDIR = f"/work/submit/{getpass.getuser()}/HmumuRun3/HISTOS/"
DEFAULT_OUTDIR = f"/home/submit/{getpass.getuser()}/public_html/HmumuRun3/"

# Region pairing this script works in: shapes from the sideband, signal from
# the unrestricted prediction.
SHAPE_REGION = "SR_sideband"
SIGNAL_REGION = cfg.SIGNAL_REGION_FOR[SHAPE_REGION]

# Colour palette local to this script -- SummaryPlots.py's stacked plots keep
# their own colours (plot_style.PROCESS_COLORS), unchanged.
#
# Colour is assigned by RANK (0 = most dominant by integral, for THIS
# variable/category), not by fixed process identity. With --n-bkg/--n-sig > 2,
# the extra entries reuse the last colour.
BKG_COLORS = [ROOT.kAzure + 7, ROOT.kAzure + 2]
SIG_COLORS = [ROOT.kRed + 2, ROOT.kRed - 4]

# Axis title offsets for THIS script's standalone canvas. Deliberately separate
# from plot_style.RATIO_X_TITLE_OFFSET (=1.0), which is tuned for
# SummaryPlots.py's ratio pad.
X_TITLE_OFFSET_SHAPES = 1.4
Y_TITLE_OFFSET_SHAPES = 2.0

# Legend position local to this script (not touching plot_style.LEGEND_POS,
# shared with the stacked plots). One row per process, so it can be compact.
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
    p.add_argument("category", choices=cfg.CATEGORIES)
    p.add_argument("year", help="data-taking year, e.g. 2024, 2025, 12022 ...")
    p.add_argument("-i", "--indir", default=DEFAULT_INDIR,
                   help="directory holding histos_<cat><year>.root "
                        "(default: %(default)s)")
    p.add_argument("-o", "--outdir", default=DEFAULT_OUTDIR,
                   help="output directory for PNGs (default: %(default)s)")
    p.add_argument("--n-bkg", type=int, default=2,
                   help="how many dominant backgrounds to show, ranked by "
                        "integral (default: 2)")
    p.add_argument("--n-sig", type=int, default=2,
                   help="how many dominant signals to show, ranked by "
                        "integral (default: 2)")
    p.add_argument("--log", action="store_true",
                   help="log y-axis for all plots (default: linear)")
    p.add_argument("--no-data", action="store_true",
                   help="don't overlay the data shape")
    p.add_argument("--vars", nargs="+", default=None,
                   help="only plot these variables (default: the category's "
                        "active set)")
    return p.parse_args()


def _pick_dominant(hists, proc_list, n):
    """Up to n (name, hist, integral) tuples from proc_list, ranked by integral
    descending, skipping missing/empty ones."""
    scored = []
    for name in proc_list:
        h = hists.get(name)
        if h is None:
            continue
        integral = h.Integral()
        if integral > 0:
            scored.append((name, h, integral))
    scored.sort(key=lambda t: t[2], reverse=True)
    return scored[:n]


def plot_shape(hfile, category, year, varname, outdir, args):

    if not hfile.has_variable(SHAPE_REGION, varname):
        print(f"   -> skipped '{varname}': not in histogram file")
        return

    titleX = plot_vars.get_xlabel(varname)

    bkg_hists = {p: hfile.get(SHAPE_REGION, varname, p) for p in cfg.BKG_PROCS}
    bkg_hists[cfg.DATA_PROCESS] = hfile.get(SHAPE_REGION, varname,
                                            cfg.DATA_PROCESS)
    sig_hists = {p: hfile.get(SIGNAL_REGION, varname, p) for p in cfg.SIG_PROCS}

    chosen_bkg = _pick_dominant(bkg_hists, cfg.BKG_PROCS, args.n_bkg)
    chosen_sig = _pick_dominant(sig_hists, cfg.SIG_PROCS, args.n_sig)
    print(f"[{varname}] dominant backgrounds: {[n for n, _, _ in chosen_bkg]}")
    print(f"[{varname}] dominant signals:     {[n for n, _, _ in chosen_sig]}")

    hData = bkg_hists.get(cfg.DATA_PROCESS)
    show_data = (not args.no_data) and hData is not None and hData.Integral() > 0

    # --- canvas ---
    c = ROOT.TCanvas("c", "", plot_style.CANVAS_W, plot_style.CANVAS_H)
    # cms_label() sizes its text as a fraction of THIS pad's height.
    # SummaryPlots.py draws it on pad1, only (1-RATIO_SPLIT) of the canvas
    # height; here we draw on the full canvas, so scale the top margin to keep
    # the label the same ABSOLUTE size in both scripts.
    c.SetTopMargin(plot_style.PAD_TOP_MARGIN * (1 - plot_style.RATIO_SPLIT))
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
        hc.SetDirectory(0)
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

    for rank, (name, h, _) in enumerate(chosen_bkg):
        _normalize(name, h, cfg.PROCESS_LABELS.get(name, name), _bkg_color(rank))
    for rank, (name, h, _) in enumerate(chosen_sig):
        _normalize(name, h, cfg.PROCESS_LABELS.get(name, name), _sig_color(rank))
    if show_data:
        _normalize(cfg.DATA_PROCESS, hData, cfg.PROCESS_LABELS["hData"],
                   ROOT.kBlack, is_data=True)

    if not drawn:
        print(f"[{varname}] nothing to draw (no non-empty histograms) -- skipping")
        return

    labels = plot_vars.get_bin_labels(varname)

    for i, (hc, label, is_data) in enumerate(drawn):
        if args.log:
            hc.SetMaximum(ymax * 10)
            hc.SetMinimum(max(ymax * 1e-6, 1e-8))
        else:
            hc.SetMaximum(ymax * 1.3)
            hc.SetMinimum(0.0)

        if labels:
            for b, lab in enumerate(labels, start=1):
                hc.GetXaxis().SetBinLabel(b, lab)

        hc.GetXaxis().SetTitle(titleX)
        # The inherited relative-mode axis title sizes (a fraction of pad
        # HEIGHT) were tuned for pad1 in the stacked plot, whose own x-axis is
        # suppressed. On this full-height standalone canvas they render
        # oversized and can push the title past the canvas edge. Reuse the
        # proven pixel-mode fix from the ratio pad instead.
        hc.GetXaxis().SetTitleFont(plot_style.RATIO_FONT_ABS)
        hc.GetXaxis().SetTitleSize(plot_style.RATIO_TITLE_SIZE_PX)
        hc.GetXaxis().SetTitleOffset(X_TITLE_OFFSET_SHAPES)
        hc.GetYaxis().SetTitle("Normalized Entries")
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

    # Mass-range label, matching SummaryPlots.py's SR_sideband label.
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.035)
    lines = [
        "H #rightarrow #mu#mu",
        f"{cfg.CATEGORY_LABELS.get(category, category)} ({year.lstrip('_')})",
        cfg.REGIONS[SHAPE_REGION]["label"],
    ]
    y0, dy = 0.85, 0.045
    for i, line in enumerate(lines):
        latex.DrawLatex(plot_style.PAD_LEFT_MARGIN + 0.05, y0 - i * dy, line)

    os.makedirs(outdir, exist_ok=True)
    outpath = f"{outdir}{varname}_{category}{year}_Shape.png"
    c.SaveAs(outpath)
    print(f"   -> {outpath}")


def main():
    args = parse_args()

    category = args.category
    year = "_" + args.year
    indir = args.indir if args.indir.endswith("/") else args.indir + "/"
    baseOutDir = args.outdir if args.outdir.endswith("/") else args.outdir + "/"

    if year not in lumis:
        sys.exit(f"ERROR: year '{args.year}' not found in lumis table "
                 f"(known: {[k.lstrip('_') for k in lumis if k.startswith('_')]})")

    hfile = HistoFile(indir + cfg.histo_filename(category, year))

    # Nested as a sibling of SummaryPlots.py's per-group trees under the same
    # category+year folder.
    outdir = os.path.join(baseOutDir, f"{category}{year}",
                          "shapes", SHAPE_REGION) + "/"

    active_vars = args.vars or cfg.get_active_vars(category)

    print(f"[shapes] category={category} year={args.year} "
          f"n_bkg={args.n_bkg} n_sig={args.n_sig} "
          f"scale={'log' if args.log else 'linear'} "
          f"data={'off' if args.no_data else 'on'}")
    print(f"[shapes] shapes from '{SHAPE_REGION}', signal from '{SIGNAL_REGION}'")
    print(f"[shapes] output dir: {outdir}")
    print(f"[shapes] plotting {len(active_vars)} variables")

    for varname in active_vars:
        plot_shape(hfile, category, year, varname, outdir, args)


if __name__ == "__main__":
    main()
