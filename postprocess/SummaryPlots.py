"""
SummaryPlots.py -- stacked data/MC comparison plots.

Reads the histograms produced by makeHistos.py, and make plots.
Run makeHistos.py first:

    python makeHistos.py ggHcat 2024
    python SummaryPlots.py ggHcat 2024

Output layout:
    <outdir>/<category>_<year>/<group>/<region>/<var>_<cat>_<year>_Stack.png
except the "mass" group, which skips the region split (the SR-sideband window
is just a gapped view of the same spectrum).
"""

import ROOT
import os
import sys
import argparse
import getpass

from utils import plot_style
from utils.plot_style import lumis
from utils import plot_vars
from utils import histo_config as cfg

plot_style.setup_style()

DEFAULT_INDIR = f"/work/submit/{getpass.getuser()}/HmumuRun3/HISTOS/"
DEFAULT_OUTDIR = f"/home/submit/{getpass.getuser()}/public_html/HmumuRun3/"


def parse_args():
    p = argparse.ArgumentParser(
        description="Stacked data/MC plots from a makeHistos.py histogram file."
    )
    p.add_argument("category", choices=cfg.CATEGORIES)
    p.add_argument("year", help="data-taking year, e.g. 2024, 2025, 12022 ...")
    p.add_argument("-i", "--indir", default=DEFAULT_INDIR,
                   help="directory holding histos_<cat><year>.root "
                        "(default: %(default)s)")
    p.add_argument("-o", "--outdir", default=DEFAULT_OUTDIR,
                   help="output directory for PNGs (default: %(default)s)")
    p.add_argument("--groups", nargs="+", default=None,
                   help="which variable groups to draw (default: all present)")
    p.add_argument("--no-SR-sideband", dest="sr_sideband",
                   action="store_false", default=True,
                   help="only produce the Inclusive region")
    p.add_argument("--linear", action="store_true",
                   help="linear y-axis (default: log)")
    return p.parse_args()


# ---------------------------------------------------------------------------
# Histogram access
# ---------------------------------------------------------------------------
class HistoFile:
    """Thin reader over makeHistos.py's output.

    Histograms are Clone()d on read so the caller never depends on the TFile
    staying open -- otherwise ROOT deletes them out from under the drawing
    code, which is the classic gotcha when moving from RDataFrame results to
    file-resident histograms.

    Every get() returns a FRESH clone, deliberately uncached. Callers mutate
    what they receive (SetLineColor, SetFillStyle, SetBinLabel, Scale), so a
    cached object handed to a second caller would arrive already styled as
    something else -- e.g. a signal histogram read once for the Inclusive
    region and again as the Unrestricted source for SR_sideband. Reads are
    cheap; the aliasing bug would not be.
    """

    def __init__(self, path):
        self.path = path
        self.f = ROOT.TFile.Open(path)
        if not self.f or self.f.IsZombie():
            sys.exit(f"ERROR: cannot open histogram file {path}\n"
                     f"       run makeHistos.py first")
        prov = self.f.Get("provenance")
        if prov:
            print(f"[plots] input: {prov.GetTitle()}")
        self._n = 0

    def get(self, region, varname, process):
        """Return a fresh detached TH1 clone, or None if absent."""
        h = self.f.Get(cfg.histo_path(region, varname, process))
        if not h:
            return None
        # Unique name per clone: ROOT keys objects by name, and two live
        # histograms sharing one would collide.
        self._n += 1
        clone = h.Clone(f"{region}_{varname}_{process}_c{self._n}")
        clone.SetDirectory(0)
        return clone

    def has_variable(self, region, varname):
        return any(self.get(region, varname, p) is not None
                   for p in cfg.ALL_PROCESSES)


def load_processes(hfile, region, varname):
    """Assemble {process: TH1} for one region+variable.

    Signals come from SIGNAL_REGION_FOR[region] -- the full, unrestricted
    prediction -- because a narrow resonance concentrates in [120,130], the
    exact core the SR-sideband window excludes. Backgrounds and data come from
    the region itself.
    """
    hists = {}
    for proc in cfg.BKG_PROCS + [cfg.DATA_PROCESS]:
        h = hfile.get(region, varname, proc)
        if h is not None:
            hists[proc] = h

    sig_region = cfg.SIGNAL_REGION_FOR[region]
    for proc in cfg.SIG_PROCS:
        h = hfile.get(sig_region, varname, proc)
        if h is not None:
            hists[proc] = h

    return hists


# ---------------------------------------------------------------------------
# Styling helpers
# ---------------------------------------------------------------------------
def _style_ratio_axis(axis, offset):
    """Shared absolute-pixel ratio-pad text style for one axis."""
    axis.SetTitleFont(plot_style.RATIO_FONT_ABS)
    axis.SetLabelFont(plot_style.RATIO_FONT_ABS)
    axis.SetTitleSize(plot_style.RATIO_TITLE_SIZE_PX)
    axis.SetLabelSize(plot_style.RATIO_LABEL_SIZE_PX)
    axis.SetTitleOffset(offset)


def _style_mc(hist, process):
    hist.SetLineWidth(3)
    col = plot_style.process_color(process)
    hist.SetLineColor(col)
    hist.SetFillColor(col)


def _style_data(hist):
    hist.SetMarkerStyle(plot_style.DATA_MARKER_STYLE)
    hist.SetMarkerSize(plot_style.DATA_MARKER_SIZE)
    hist.SetLineWidth(plot_style.DATA_LINE_WIDTH)
    hist.SetLineColor(ROOT.kBlack)


def _apply_bin_labels(hist, varname):
    labels = plot_vars.get_bin_labels(varname)
    if labels:
        for i, lab in enumerate(labels, start=1):
            hist.GetXaxis().SetBinLabel(i, lab)


def _draw_mass_range_label(category, year, region_label):
    ROOT.gPad.cd()
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.035)
    lines = [
        "H #rightarrow #mu#mu",
        f"{cfg.CATEGORY_LABELS.get(category, category)} ({year.lstrip('_')})",
        region_label,
    ]
    y0, dy = 0.85, 0.045
    for i, line in enumerate(lines):
        latex.DrawLatex(plot_style.PAD_LEFT_MARGIN + 0.05, y0 - i * dy, line)
    return latex        # returned so the caller keeps it alive


# ---------------------------------------------------------------------------
# The plot
# ---------------------------------------------------------------------------
def plot(hfile, category, year, varname, region_name, outdir, doLog=True):

    region = cfg.REGIONS[region_name]

    if not hfile.has_variable(region_name, varname):
        print(f"   -> skipped '{varname}' [{region_name}]: not in histogram file")
        return

    hists = load_processes(hfile, region_name, varname)
    titleX = plot_vars.get_xlabel(varname)

    hData = hists.get(cfg.DATA_PROCESS)
    hDY = hists.get("hDY")

    c, pad1, pad2 = plot_style.make_canvas_pads(doLog)

    # Backgrounds stacked + filled; signal drawn unstacked as solid outlines
    # at true scale, not summed into the background.
    BKGstack = ROOT.THStack()
    SIGstack = ROOT.THStack()

    for proc in cfg.BKG_PROCS:
        h = hists.get(proc)
        if h is None:
            continue
        print(f"   Integral {proc} = {h.Integral()}")
        _style_mc(h, proc)
        _apply_bin_labels(h, varname)
        BKGstack.Add(h)

    for proc in cfg.SIG_PROCS:
        h = hists.get(proc)
        if h is None:
            continue
        print(f"   Integral {proc} = {h.Integral()}")
        _style_mc(h, proc)
        _apply_bin_labels(h, varname)
        h.SetFillStyle(0)               # outline only, not filled/stacked
        SIGstack.Add(h)

    if BKGstack.GetNhists() == 0:
        print(f"   -> skipped '{varname}' [{region_name}]: no background histograms")
        return

    # --- y range ---
    rangeYax = 10 if doLog else 2
    if hData is not None and hDY is not None:
        BKGstack.SetMaximum(rangeYax * max(hData.GetMaximum(), hDY.GetMaximum()))
    if hDY is not None:
        floor = 1e9 if varname == "mva" else 1e6
        BKGstack.SetMinimum(hDY.GetMaximum() / floor)

    # --- upper pad ---
    # Draw order matters and matches the original exactly: data first, then
    # the stack WITHOUT "SAME". The stack then draws its own axes and owns the
    # visible frame, which is what makes the GetYaxis()/SetMaximum/SetMinimum
    # calls below take effect. Drawing the stack with "SAME" leaves the frame
    # belonging to hData, and the stack's y-axis title and range are silently
    # ignored.
    pad1.cd()
    if hData is not None:
        _style_data(hData)
        _apply_bin_labels(hData, varname)
        print(f"   Integral hData = {hData.Integral()}")
        hData.Draw("E")

    BKGstack.Draw("HIST")

    BKGstack.GetXaxis().SetLabelSize(0)   # glued layout: x-axis on ratio pad only
    BKGstack.GetXaxis().SetTitleSize(0)
    BKGstack.GetYaxis().SetTitle(
        "Events / GeV" if varname == "dimu_mass" else "Events")
    BKGstack.GetYaxis().SetTitleOffset(1.1)
    BKGstack.GetYaxis().SetLabelSize(0.04)
    BKGstack.GetYaxis().SetTitleSize(0.045)
    BKGstack.GetYaxis().ChangeLabel(1, -1, 0)

    if SIGstack.GetNhists() > 0:
        SIGstack.Draw("HIST NOSTACK SAME")
    if hData is not None:
        hData.Draw("E SAME")

    # --- ratio pad ---
    mcTOT = BKGstack.GetStack().Last()
    print(f"   ALL mcTOT integral = {mcTOT.Integral()}")

    if hData is not None:
        pad2.cd()
        print(f"   ALL data integral = {hData.Integral()}")
        ratio = hData.Clone("dataratio")
        ratio.SetDirectory(0)
        ratio.Divide(mcTOT)
        ratio.GetYaxis().SetTitle("data/MC")
        ratio.GetYaxis().SetRangeUser(*plot_vars.get_ratio_range(varname))
        ratio.GetXaxis().SetTitle(titleX)
        _style_ratio_axis(ratio.GetXaxis(), plot_style.RATIO_X_TITLE_OFFSET)
        _style_ratio_axis(ratio.GetYaxis(), plot_style.RATIO_Y_TITLE_OFFSET)
        ratio.Draw("pe")
        lineOne = ROOT.TLine(mcTOT.GetXaxis().GetXmin(), 1.,
                             mcTOT.GetXaxis().GetXmax(), 1.)
        lineOne.SetLineColor(11)
        lineOne.Draw("same")

    # --- legend ---
    pad1.cd()
    legend = ROOT.TLegend(*plot_style.LEGEND_POS)
    legend.SetNColumns(plot_style.LEGEND_NCOLUMNS)
    legend.SetColumnSeparation(plot_style.LEGEND_COLUMN_SEP)
    legend.SetFillStyle(0)
    legend.SetBorderSize(0)
    legend.SetTextSize(plot_style.LEGEND_TEXT_SIZE)
    legend.SetTextAlign(12)

    sig_order = (cfg.LEGEND_SIG_LEP_ORDER if category in cfg.LEP_CATEGORIES
                 else cfg.LEGEND_SIG_HAD_ORDER)
    for proc in cfg.LEGEND_BKG_ORDER + sig_order:
        h = hists.get(proc)
        if h is not None and h.Integral() > 0:
            style = "lep" if proc == cfg.DATA_PROCESS else (
                "f" if proc in cfg.BKG_PROCS else "l")
            legend.AddEntry(h, cfg.PROCESS_LABELS.get(proc, proc), style)
    legend.Draw()

    # --- CMS label + region / blinding annotation ---
    _cmslabel = plot_style.cms_label(pad1, year)

    # Matches the original branching:
    #   region-restricted  -> mass-range label, unconditionally for EVERY
    #                         variable (mva included), since the restriction
    #                         itself is what makes unblinded data safe
    #   Inclusive + blind window -> vertical lines marking it
    #   Inclusive, no blind window, not mva -> mass-range label
    #     (mva is excluded because it uses its own score-based cut instead)
    _keepalive = []
    blind_range = plot_vars.get_blind_range(varname)
    if region["filter"] is not None:
        _keepalive.append(
            _draw_mass_range_label(category, year, region["label"]))
    elif blind_range is not None:
        lo, hi = blind_range
        for x in (lo, hi):
            ln = ROOT.TLine(x, 0, x, 500000.)
            ln.SetLineColor(11)
            ln.Draw()
            _keepalive.append(ln)
    elif varname != "mva":
        _keepalive.append(
            _draw_mass_range_label(category, year, region["label"]))

    os.makedirs(outdir, exist_ok=True)
    outpath = f"{outdir}{varname}_{category}{year}_Stack.png"
    c.SaveAs(outpath)
    print(f"   -> {outpath}")


# ---------------------------------------------------------------------------
# Dispatch
# ---------------------------------------------------------------------------
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

    groups = args.groups or cfg.groups_for(category)
    doLog = not args.linear

    print(f"[plots] category={category} year={args.year} groups={groups}")
    print(f"[plots] SR-sideband: {'ON' if args.sr_sideband else 'OFF'}")
    print(f"[plots] base output dir: {baseOutDir}")

    for group in groups:
        variables = cfg.vars_in_group(category, group)
        if not variables:
            print(f"[plots] group '{group}': nothing to draw for {category}")
            continue

        if group in cfg.GROUPS_WITHOUT_REGION_SPLIT:
            regions = ["Inclusive"]
        else:
            regions = ["Inclusive"] + (["SR_sideband"] if args.sr_sideband else [])

        for region_name in regions:
            if group in cfg.GROUPS_WITHOUT_REGION_SPLIT:
                outdir = os.path.join(baseOutDir, f"{category}{year}", group) + "/"
            else:
                outdir = os.path.join(baseOutDir, f"{category}{year}",
                                      group, region_name) + "/"
            print(f"[plots] group '{group}' region '{region_name}' -> {outdir}")
            for varname in variables:
                plot(hfile, category, year, varname, region_name, outdir,
                     doLog=doLog)


if __name__ == "__main__":
    main()
