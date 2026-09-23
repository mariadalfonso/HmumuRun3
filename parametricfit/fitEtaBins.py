"""
fitEtaBins.py -- double-sided Crystal Ball fits of HiggsCandCorrMass split by
dimuon |eta| topology and by dimuon pT, for one signal production mode.

Depends only on sigFit.py and prepareFits.py:
  * sigFit.py supplies the shape, fit range, binning, the two-pass Minuit2
    fit (`fit_one`, which writes the ratio/pull canvases) and the plot style,
    so mu and sigma here are directly comparable with the signal workspaces.
  * prepareFits.py supplies the file lookup and the year handling, so the
    year may be a concrete era (12022, 22022, 12023, 22023, 2024) or a group
    (2022 = 12022 + 22022, 2023 = 12023 + 22023, Run3). A group is fitted as
    one combined sample, and every required snapshot must exist.

What it fits
------------
  inclusive   one reference fit over the whole sample

  eta         dimuon |eta| topology: BB, BE, EE from ETA_EDGES =
              [0, 1.4, 2.4]. Every unordered pair of regions is one
              exclusive category, so the three partition the sample.

  mu2d_BB     within each eta topology, a 2D grid in the two muon momenta:
  mu2d_BE     leading vs subleading pT, MU2D_WIDTH = 20 GeV cells from
  mu2d_EE     MU2D_MIN = 20 to MU2D_MAX = 120 GeV plus an open cell above.
              Each cell is fitted and the fitted sigma is drawn as the z
              axis of a TH2, giving a map of the mass resolution over the
              muon kinematics at fixed eta topology. mu and N_eff maps are
              produced alongside, on z scales shared across the three
              topologies so they can be compared directly.

The axes are LEADING and SUBLEADING pT rather than Muon1/Muon2, because the
pair is unordered: a Muon1-vs-Muon2 map would be symmetric about the
diagonal and half the fits redundant. Only the lower triangle is fitted.

Map cells are fitted WITHOUT writing per-cell canvases -- there are ~21 per
topology -- unless --cell-plots is given. Cells below --min-entries are
left unfilled and drawn blank, so an unfitted cell cannot be mistaken for
one with a small resolution.

The binning is set by the ETA_EDGES and MU2D_* constants at the top of this
file rather than by command-line options.

Usage
-----
    python fitEtaBins.py 2024
    python fitEtaBins.py 2024 --cell-plots
    python fitEtaBins.py Run3 --sig VH --category VHcat
    python fitEtaBins.py 22022 --file <snapshot> [<snapshot> ...]
"""

import argparse
from array import array
import getpass
import json
import os
from datetime import datetime, timezone

import ROOT

import sigFit as SF
from sigFit import fit_one, git_commit, XLOW, XHIGH, BINS_PER_GEV
from prepareFits import (ROOTFILES_DIR, SELMVA, SIGNAL_TAGS,
                         SIGNAL_YEAR_GROUPS, expand_year, get_signal_files,
                         normalize_year)

ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True)

TREE = "events"
MASSBINS = int(round((XHIGH - XLOW) * BINS_PER_GEV))
MASSCOL = "HiggsCandCorrMass"

# overlay colours; the first is sigFit's fit colour
PALETTE_RGB = [
    SF.FIT_COLOR_RGB,   # red
    (31, 78, 158),      # blue
    (34, 139, 34),      # green
    (230, 159, 0),      # orange
    (128, 64, 160),     # violet
    (0, 158, 176),      # teal
]

# gap from the left edge of "CMS" to "Internal", in units of the CMS text
# size ("CMS" in font 61 is ~2.2 wide). 600x420 pad -> 0.10 NDC, as in
# sigFit; 600x600 canvas -> 0.144.
CMS_EXTRA_OFFSET = 2.4

# z-axis colour palette for the 2D resolution maps. Any ROOT EColorPalette
# works: kPastel, kBird, kViridis, kCividis, kTemperatureMap, kRainBow, ...
# A SEQUENTIAL palette (kViridis, kCividis) is the safer choice for a
# resolution map, since the eye then reads "darker = better"; kPastel and
# kRainBow are not monotonic in lightness, so nearby cells can look more
# different than they are.
MAP_PALETTE = ROOT.kCandy

# overlay legend
OVERLAY_TEXT_SIZE = 0.026   # NDC; each entry is two lines of this
OVERLAY_GAP = 0.03          # NDC between the tallest curve and the legend

# ---------------------------------------------------------------------------
# the splits. Edit here; they are deliberately not command-line options.
# ---------------------------------------------------------------------------
# |eta| region edges -> regions -> every unordered pair is one category
ETA_EDGES = [0.0, 1.4, 2.4]

# 2D map of the fitted resolution over the two muon momenta, one map per
# |eta| topology. Cells are `MU2D_WIDTH` GeV wide in each muon, from
# MU2D_MIN up to MU2D_MAX, plus one open bin above MU2D_MAX.
#
# The axes are LEADING and SUBLEADING pT, not Muon1/Muon2: the pair is
# unordered, so a Muon1-vs-Muon2 map would be symmetric about the diagonal
# and half the fits would be redundant. Ordering by pT folds it into the
# lower triangle, halving the work and making every cell distinct.
MU2D_MIN = 20.0
MU2D_MAX = 120.0
MU2D_WIDTH = 20.0
MU2D_LEAD = "muLeadPt"
MU2D_SUB = "muSubPt"
MU2D_LEAD_AXIS = "leading muon p_{T} (GeV)"
MU2D_SUB_AXIS = "subleading muon p_{T} (GeV)"

ETA_AXIS = "|#eta| topology"

# short region names by number of |eta| regions
REGION_NAMES = {2: ["B", "E"], 3: ["B", "O", "E"]}

# columns defined on the fly when a snapshot does not have them
DERIVED = {
    "absEta1": "(float)std::abs(Muon1_eta)",
    "absEta2": "(float)std::abs(Muon2_eta)",
    MU2D_LEAD: "(float)std::max(Muon1_pt, Muon2_pt)",
    MU2D_SUB: "(float)std::min(Muon1_pt, Muon2_pt)",
}


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("year",
                   help="era (12022, 22022, 12023, 22023, 2024) or group "
                        "(2022, 2023, Run3)")
    p.add_argument("--sig", default="ggH", choices=sorted(SIGNAL_TAGS),
                   help="production mode (default: %(default)s)")
    p.add_argument("--category", default="ggHcat", choices=sorted(SELMVA),
                   help="snapshot category directory (default: %(default)s)")
    p.add_argument("--rootfiles-dir", default=ROOTFILES_DIR,
                   help="snapshot base dir (default: %(default)s)")
    p.add_argument("--file", nargs="+", default=None,
                   help="fit these snapshots instead of the looked-up ones. "
                        "The year is then only used for the label and "
                        "output names.")
    p.add_argument("--min-entries", type=int, default=500,
                   help="skip a map cell with fewer raw entries; a DCB has "
                        "six shape parameters (default: %(default)s)")
    p.add_argument("--cell-plots", action="store_true",
                   help="also write the ratio/pull canvases for every map "
                        "cell (dozens of files per topology)")
    p.add_argument("--mass", default=MASSCOL,
                   help="mass column to fit (default: %(default)s)")
    p.add_argument("--weight", default="w_allSF",
                   help="weight column, or 'none' for unweighted")
    p.add_argument("-o", "--plotdir",
                   default=f"/home/submit/{getpass.getuser()}/public_html/HmumuFits/EtaBinFits",
                   help="base output dir. Everything is written to a "
                        "<plotdir>/<year>/ subdirectory, so runs for different "
                        "years never overwrite each other.")
    p.add_argument("--no-inclusive", action="store_true",
                   help="skip the inclusive fit used as the delta reference")
    p.add_argument("--palette", default=None,
                   help="ROOT palette name for the 2D maps, e.g. kPastel, "
                        "kBird, kViridis, kCividis, kTemperatureMap "
                        f"(default: whatever MAP_PALETTE is set to)")
    p.add_argument("--pdf", action="store_true", help="also write .pdf")
    return p.parse_args()


# ---------------------------------------------------------------------------
# bins
# ---------------------------------------------------------------------------

def _region_range(lo, hi):
    return f"|#eta| < {hi:g}" if lo == 0 else f"{lo:g} #leq |#eta| < {hi:g}"


def _inside(col, edges, i):
    return f"{col} >= {edges[i]} && {col} < {edges[i + 1]}"


def eta_bins(edges=None):
    """Exclusive dimuon |eta| categories: BB, BE, EE.

    Region i is edges[i] <= |eta| < edges[i+1]. Each unordered pair (i, j)
    is one category; a mixed pair accepts either muon ordering, so every
    event lands in exactly one category and the categories partition the
    inclusive sample.
    """
    edges = list(ETA_EDGES if edges is None else edges)
    n = len(edges) - 1
    names = REGION_NAMES.get(n, [f"R{i}" for i in range(n)])

    bins = []
    for i in range(n):
        for j in range(i, n):
            short = names[i] + names[j]
            if i == j:
                filt = (f"({_inside('absEta1', edges, i)}) && "
                        f"({_inside('absEta2', edges, i)})")
                label = f"{short}: both {_region_range(edges[i], edges[i + 1])}"
            else:
                filt = (f"(({_inside('absEta1', edges, i)}) && "
                        f"({_inside('absEta2', edges, j)})) || "
                        f"(({_inside('absEta1', edges, j)}) && "
                        f"({_inside('absEta2', edges, i)}))")
                label = f"{short}: one {names[i]}, one {names[j]}"
            bins.append({"tag": f"eta_{short}", "short": short,
                         "label": label, "lo": None, "hi": None,
                         "filt": filt, "group": "eta"})
    return bins


def mu2d_edges():
    """Cell edges in one muon pT, plus an open bin above the last one."""
    n = int(round((MU2D_MAX - MU2D_MIN) / MU2D_WIDTH))
    return [MU2D_MIN + i * MU2D_WIDTH for i in range(n + 1)]


def _cell_range(edges, i):
    """(lo, hi, label, filter-upper) for cell i. i == len(edges)-1 is open."""
    lo = edges[i]
    if i == len(edges) - 1:
        return lo, lo + MU2D_WIDTH, f"#geq {lo:g}", None
    return lo, edges[i + 1], f"{lo:g}-{edges[i + 1]:g}", edges[i + 1]


def _cell_cut(col, edges, i):
    lo, _, _, hi = _cell_range(edges, i)
    return f"{col} >= {lo}" if hi is None else f"{col} >= {lo} && {col} < {hi}"


def mu2d_cells(eta_bin):
    """Fit cells for the (leading, subleading) pT map inside one eta bin.

    Only the lower triangle is generated: the subleading muon cannot be in
    a higher pT cell than the leading one, so cells above the diagonal are
    empty by construction and are not fitted.
    """
    edges = mu2d_edges()
    n = len(edges)
    out = []
    for il in range(n):                       # leading
        for isub in range(il + 1):            # subleading <= leading
            lo_l, hi_l, lab_l, _ = _cell_range(edges, il)
            lo_s, hi_s, lab_s, _ = _cell_range(edges, isub)
            filt = (f"({eta_bin['filt']}) && "
                    f"({_cell_cut(MU2D_LEAD, edges, il)}) && "
                    f"({_cell_cut(MU2D_SUB, edges, isub)})")
            out.append({
                "tag": (f"{eta_bin['tag']}_lead{lo_l:g}_sub{lo_s:g}"
                        .replace(".", "p")),
                "short": f"{lab_l} / {lab_s}",
                "label": (f"{eta_bin['short']}, lead {lab_l}, "
                          f"sub {lab_s} GeV"),
                "filt": filt,
                "group": f"mu2d_{eta_bin['short']}",
                "eta_short": eta_bin["short"],
                "ix": il + 1, "iy": isub + 1,   # TH2 bin indices
                "open_x": il == n - 1, "open_y": isub == n - 1,
            })
    return out


def coverage_filters():
    """Events some bin can accept, for the overflow report.

    The muon pT map has an open top cell, so the only way to fall outside
    is a muon below MU2D_MIN -- reported separately since the map, unlike
    the eta categories, then does not add up to the inclusive fit.
    """
    e0, e1 = ETA_EDGES[0], ETA_EDGES[-1]
    return {
        "eta": (f"absEta1 >= {e0} && absEta1 < {e1} && "
                f"absEta2 >= {e0} && absEta2 < {e1}"),
        "muon pT map": f"{MU2D_SUB} >= {MU2D_MIN}",
    }


# ---------------------------------------------------------------------------
# plotting (style from sigFit.py)
# ---------------------------------------------------------------------------

def color(i):
    return ROOT.TColor.GetColor(*PALETTE_RGB[i % len(PALETTE_RGB)])


def make_canvas(name):
    """Single-pad canvas with sigFit's geometry."""
    c = ROOT.TCanvas(name, "", SF.CANVAS_W, SF.CANVAS_H)
    c.SetLeftMargin(SF.PAD_LEFT_MARGIN)
    c.SetRightMargin(SF.PAD_RIGHT_MARGIN)
    c.SetTopMargin(SF.PAD_TOP_MARGIN)
    c.SetBottomMargin(SF.PAD_BOTTOM_MARGIN)
    return c


def cms_label(pad, year):
    """sigFit.cms_label, with the "Internal" offset scaled to this pad.

    sigFit puts the extra text a fixed 0.10 NDC after "CMS", which suits its
    600x420 upper pad. ROOT scales text by a pad's SHORTER side while NDC x
    is a fraction of its width, so on the square canvases used here "CMS" is
    wider in NDC and the two overlap. Text, sizes and luminosity still come
    from sigFit.
    """
    pad.cd()
    t = pad.GetTopMargin()
    l = pad.GetLeftMargin()
    r = pad.GetRightMargin()
    y = 1 - t + 0.2 * t
    cms_size = SF.CMS_TEXT_SIZE_FRAC * t

    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextAngle(0)
    latex.SetTextColor(ROOT.kBlack)

    latex.SetTextFont(61)
    latex.SetTextAlign(11)
    latex.SetTextSize(cms_size)
    latex.DrawLatex(l, y, SF.CMS_TEXT)

    if SF.CMS_EXTRA_TEXT:
        w_px = pad.GetWw() * pad.GetAbsWNDC()
        h_px = pad.GetWh() * pad.GetAbsHNDC()
        dx = (CMS_EXTRA_OFFSET * cms_size * min(w_px, h_px) / w_px
              if w_px > 0 else 0.10)
        latex.SetTextFont(52)
        latex.SetTextSize(SF.EXTRA_OVER_CMS_TEXT_SIZE * cms_size)
        latex.DrawLatex(l + dx, y, SF.CMS_EXTRA_TEXT)

    lumi = SF.get_lumi(year) if year else 0.0
    if lumi <= 0:
        print(f"warning: no luminosity found for year tag '{year}'")
    latex.SetTextFont(42)
    latex.SetTextAlign(31)
    latex.SetTextSize(SF.LUMI_TEXT_SIZE_FRAC * t)
    if lumi > 0:
        txt = "%.1f fb^{-1} (%.1f TeV)" % (lumi, SF.SQRT_S_TEV)
    else:
        txt = "(%.1f TeV)" % SF.SQRT_S_TEV
    latex.DrawLatex(1 - r, y, txt)
    return latex


def save(c, out, save_pdf):
    c.SaveAs(out + ".png")
    if save_pdf:
        c.SaveAs(out + ".pdf")
    c.Close()


def overlay_pdfs(x, fits, plotdir, year, group, axis, save_pdf=False):
    """Every bin's fitted PDF on one frame, each normalised to unity.

    All PDFs are built on the same observable x; a PDF plotted on a frame of
    a different variable would be drawn as a constant.

    Layout: one legend band across the top of the frame, each entry giving
    the bin on the first line and its mu, sigma on the second (two columns
    above three bins). The y axis is raised so the curves end below it.
    """
    legs = [f for f in fits if f.get("group") == group]
    if len(legs) < 2:
        return
    frame = x.frame(ROOT.RooFit.Title(""))
    ymax = 0.0
    for i, f in enumerate(legs):
        f["pdf"].plotOn(frame, ROOT.RooFit.Name(f["tag"]),
                        ROOT.RooFit.LineColor(color(i)),
                        ROOT.RooFit.LineWidth(SF.DATA_LINE_WIDTH))
        curve = frame.findObject(f["tag"])
        if curve:
            ys = curve.GetY()
            ymax = max([ymax] + [ys[k] for k in range(curve.GetN())])

    # legend band
    ncol = 1 if len(legs) <= 3 else 2
    ts = OVERLAY_TEXT_SIZE * (1.0 if ncol == 1 else 0.85)   # columns are narrow
    nrow = -(-len(legs) // ncol)
    row_h = 2 * SF.LINE_SPACING * ts            # two text lines per entry
    left = SF.PAD_LEFT_MARGIN + 0.02
    right = 1 - SF.PAD_RIGHT_MARGIN - 0.02
    top = 1 - SF.PAD_TOP_MARGIN - 0.02
    bottom = top - nrow * row_h

    leg = ROOT.TLegend(left, bottom, right, top)
    leg.SetNColumns(ncol)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(ts)
    leg.SetMargin(0.08 * ncol)                  # line sample ~same NDC width
    for f in legs:
        curve = frame.findObject(f["tag"])
        if curve:
            leg.AddEntry(curve,
                         f"#splitline{{{f['label']}}}"
                         f"{{#mu = {f['mu']:.3f}, #sigma = {f['sigma']:.3f} GeV}}",
                         "l")

    # headroom: the tallest curve must end OVERLAY_GAP below the legend
    plot_h = 1 - SF.PAD_TOP_MARGIN - SF.PAD_BOTTOM_MARGIN
    usable = (bottom - OVERLAY_GAP - SF.PAD_BOTTOM_MARGIN) / plot_h
    if ymax > 0:
        frame.SetMinimum(0.0)
        frame.SetMaximum(ymax / max(usable, 0.3))
    frame.GetYaxis().SetTitle("a.u.")
    frame.GetYaxis().SetTitleOffset(1.4)

    c = make_canvas(f"c_overlay_{group}")
    frame.Draw()
    leg.Draw()
    cms_label(c, year)
    save(c, f"{plotdir}/overlay_{group}_{year}", save_pdf)


def summary_graph(fits, plotdir, year, group, axis, save_pdf=False):
    """mu and sigma vs bin, the plot that actually answers the question."""
    legs = [f for f in fits if f.get("group") == group]
    if len(legs) < 2:
        return
    n = len(legs)
    # pT bins carry numeric edges; the eta categories do not
    categorical = legs[0]["lo"] is None
    ax = axis

    for ytit, key in [("fitted #mu (GeV)", "mu"),
                      ("fitted #sigma (GeV)", "sigma")]:
        g = ROOT.TGraphErrors(n)
        for i, f in enumerate(legs):
            if categorical:
                xc, xe = i + 0.5, 0.0
            else:
                # an open bin has no real upper edge; hi is a drawing value
                xc = 0.5 * (f["lo"] + f["hi"])
                xe = 0.5 * (f["hi"] - f["lo"])
            g.SetPoint(i, xc, f[key])
            g.SetPointError(i, xe, f[key + "Err"])
        g.SetTitle(f";{ax};{ytit}")
        g.SetMarkerStyle(SF.DATA_MARKER_STYLE)
        g.SetMarkerSize(SF.DATA_MARKER_SIZE)
        g.SetLineWidth(SF.DATA_LINE_WIDTH)

        c = make_canvas(f"c_{key}_{group}")
        if categorical:
            # labelled frame: one bin per category
            lo_y = min(f[key] - f[key + "Err"] for f in legs)
            hi_y = max(f[key] + f[key + "Err"] for f in legs)
            pad_y = 0.15 * (hi_y - lo_y) or 0.05 * abs(hi_y) or 1.0
            frame = ROOT.TH1F(f"frame_{key}_{group}", f";{ax};{ytit}",
                              n, 0, n)
            frame.SetDirectory(0)
            for i, f in enumerate(legs):
                frame.GetXaxis().SetBinLabel(i + 1, f["short"])
            frame.GetXaxis().SetLabelSize(
                0.06 if n <= 6 else 0.035)
            frame.SetMinimum(lo_y - pad_y)
            frame.SetMaximum(hi_y + pad_y)
            frame.GetYaxis().SetTitleOffset(1.5)
            frame.Draw()
            g.Draw("P SAME")
        else:
            g.GetYaxis().SetTitleOffset(1.5)
            g.Draw("AP")
        cms_label(c, year)
        save(c, f"{plotdir}/{key}_vs_{group}_{year}", save_pdf)


# ---------------------------------------------------------------------------

def resolution_map(cells, plotdir, year, eta_short, label, key, ytit,
                   zrange, save_pdf=False, fmt="4.2f"):
    """One 2D map: leading vs subleading muon pT, z = the fitted quantity.

    Empty and skipped cells are left unfilled and drawn blank (the "0"
    suffix on the COLZ option), so a cell with no fit cannot be mistaken
    for one with a small value.

    `zrange` is passed in and shared across the eta regions, otherwise each
    map auto-scales to its own range and the three cannot be compared.
    """
    edges = mu2d_edges()
    n = len(edges)
    # the open top cell is drawn one width wide; its upper edge is nominal
    axis_edges = array("d", edges + [edges[-1] + MU2D_WIDTH])

    h = ROOT.TH2D(f"map_{key}_{eta_short}_{year}", "",
                  n, axis_edges, n, axis_edges)
    h.SetDirectory(0)
    h.SetStats(0)

    filled = 0
    for c in cells:
        if c.get(key) is None:
            continue
        h.SetBinContent(c["ix"], c["iy"], c[key])
        filled += 1
    if not filled:
        return

    h.SetTitle(f";{MU2D_LEAD_AXIS};{MU2D_SUB_AXIS};{ytit}")
    h.GetZaxis().SetTitleOffset(1.25)
    for ax in (h.GetXaxis(), h.GetYaxis(), h.GetZaxis()):
        ax.SetTitleFont(43)
        ax.SetLabelFont(43)
        ax.SetTitleSize(24)
        ax.SetLabelSize(20)
    h.GetXaxis().SetTitleOffset(1.15)
    h.GetYaxis().SetTitleOffset(1.45)
    h.SetMinimum(zrange[0])
    h.SetMaximum(zrange[1])
    h.SetMarkerSize(1.3)

    c1 = make_canvas(f"c_map_{key}_{eta_short}")
    c1.SetRightMargin(0.17)
    ROOT.gStyle.SetPalette(MAP_PALETTE)
    ROOT.gStyle.SetPaintTextFormat(fmt)
    h.Draw("COLZ0 TEXT")

    # Inside the frame, top-left: the CMS header occupies
    # y = 1 - topMargin + 0.2*topMargin (0.936 by default), so anything
    # above ~0.90 collides with it. The upper-left of the map is empty by
    # construction (subleading <= leading), so the text has that space.
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(0.032)
    x0 = SF.PAD_LEFT_MARGIN + 0.04
    latex.DrawLatex(x0, 0.86, f"{label}")
    latex.DrawLatex(x0, 0.81, f"|#eta| topology: {eta_short}")
    if any(c.get("open_x") or c.get("open_y") for c in cells):
        latex.SetTextSize(0.024)
        latex.DrawLatex(x0, 0.76,
                        f"overflow bin: #geq {MU2D_MAX:g} GeV")

    cms_label(c1, year)
    save(c1, f"{plotdir}/{key}_map_{eta_short}_{year}", save_pdf)


def resolution_maps(fits, plotdir, year, label, save_pdf=False):
    """The three maps (sigma, mu, N_eff) for every eta topology, on a
    common z scale per quantity."""
    shorts = [b["short"] for b in eta_bins()]
    by_eta = {sh: [f for f in fits if f.get("group") == f"mu2d_{sh}"]
              for sh in shorts}
    by_eta = {sh: v for sh, v in by_eta.items() if v}
    if not by_eta:
        return

    panels = [
        ("sigma", "fitted #sigma (GeV)", "4.2f"),
        ("mu", "fitted #mu (GeV)", "6.1f"),
        ("n_eff", "N_{eff}", "4.0f"),
    ]
    for key, ytit, fmt in panels:
        vals = [f[key] for v in by_eta.values() for f in v
                if f.get(key) is not None]
        if not vals:
            continue
        lo, hi = min(vals), max(vals)
        pad = 0.05 * (hi - lo) if hi > lo else (0.05 * abs(hi) or 1.0)
        zrange = (lo - pad, hi + pad)
        for sh, v in by_eta.items():
            resolution_map(v, plotdir, year, sh, label, key, ytit,
                           zrange, save_pdf, fmt)
        print(f"  {key} map: z range {zrange[0]:.3g} .. {zrange[1]:.3g} "
              f"(shared across {', '.join(by_eta)})")


def fit_cell(x, hist, tag, min_entries):
    """Fit one map cell WITHOUT writing canvases.

    sigFit.fit_one always draws the ratio and pull panels, which is right
    for a handful of bins and wrong for ~60 map cells per eta region. This
    does the same two-pass Minuit2 fit with the same pinned nuisances and
    returns just the numbers.
    """
    norm = SF.range_integral(hist)
    n_eff = SF.effective_entries(hist)
    if hist.GetEntries() < min_entries or norm <= 0:
        return None

    data = ROOT.RooDataHist(f"dh_{tag}", "data", ROOT.RooArgList(x), hist)
    pdf, bundle = SF.make_signal_pdf(x, tag)
    nom, nuis = bundle["nom"], bundle["nuis"]
    for v in nuis.values():
        v.setVal(0.0)
        v.setConstant(True)

    opts = [ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.Strategy(2),
            ROOT.RooFit.Save(True), ROOT.RooFit.Range("full"),
            ROOT.RooFit.SumW2Error(True), ROOT.RooFit.PrintLevel(-1)]
    pdf.fitTo(data, *opts)
    res = pdf.fitTo(data, *opts)

    parked = [k for k, v in nom.items() if SF.near_bound(v)]
    ok = res.status() == 0 and res.covQual() == 3 and not parked
    return {"mu": nom["mu"].getVal(), "muErr": nom["mu"].getError(),
            "sigma": nom["sigma"].getVal(),
            "sigmaErr": nom["sigma"].getError(),
            "yield": norm, "n_eff": n_eff,
            "status": res.status(), "cov_qual": res.covQual(),
            "params_at_bound": parked, "ok": ok,
            "_keep": (pdf, bundle, data)}


def input_files(args, year):
    if args.file:
        missing = [f for f in args.file if not os.path.exists(f)]
        if missing:
            raise SystemExit(f"missing snapshot(s): {', '.join(missing)}")
        return list(args.file)
    try:
        return get_signal_files(args.sig, args.category, year,
                                args.rootfiles_dir)
    except (FileNotFoundError, ValueError) as e:
        raise SystemExit(str(e))


def main():
    global MAP_PALETTE
    args = parse_args()
    if args.palette:
        if not hasattr(ROOT, args.palette):
            raise SystemExit(f"ROOT has no palette named {args.palette}")
        MAP_PALETTE = getattr(ROOT, args.palette)
    SF.setup_style()
    year = normalize_year(args.year)

    files = input_files(args, year)
    eras = expand_year(year, SIGNAL_YEAR_GROUPS)

    df = ROOT.RDataFrame(TREE, files)
    cols = set(str(c) for c in df.GetColumnNames())

    need = ["Muon1_eta", "Muon2_eta", "Muon1_pt", "Muon2_pt", args.mass]
    missing = [c for c in need if c not in cols]
    if missing:
        raise SystemExit(f"input has no {', '.join(missing)}")

    for name, expr in DERIVED.items():
        if name not in cols:
            df = df.Define(name, expr)

    weight = None if args.weight == "none" else args.weight
    if weight and weight not in cols:
        raise SystemExit(f"input has no weight column {weight}")

    # same sanity cut as prepareFits.getHisto
    df = df.Filter(f"!std::isnan({args.mass})", "valid mass")

    # every output goes under <plotdir>/<year>/, and every fit tag carries the
    # year too, so runs for different years cannot overwrite one another
    outdir = os.path.join(args.plotdir, args.category, year)
    os.makedirs(outdir, exist_ok=True)

    # eta categories (inclusive in pT), pT bins (inclusive in eta), then a
    # pT scan inside each eta category
    eta_shorts = [b["short"] for b in eta_bins()]
    groups = [("eta", ETA_AXIS)]

    print(f"sample  : {args.sig} ({args.category})"
          f"\nyear    : {year}"
          + ("" if args.file else f"  (eras: {', '.join(eras)})")
          + f"\nfiles   : {len(files)}"
          + "".join(f"\n          {f}" for f in files)
          + f"\nmass    : {args.mass}"
          f"\neta     : {ETA_EDGES}  -> "
          f"{', '.join(b['short'] for b in eta_bins())}"
          f"\nmuon map: {MU2D_MIN:g}-{MU2D_MAX:g} GeV in "
          f"{MU2D_WIDTH:g} GeV cells + open, lead vs sub "
          f"({len(mu2d_cells(eta_bins()[0]))} cells per topology)"
          + f"\ngroups  : {', '.join(g for g, _ in groups)}"
          f"\nweight  : {weight or 'unweighted'}"
          f"\nwindow  : {XLOW}-{XHIGH} GeV, {MASSBINS} bins"
          f"\nplots   : {outdir}")
    if args.file:
        print("note: --file given; the luminosity label is still taken "
              f"from year '{year}'")

    # --- book everything, so the files are read once ----------------------
    def histo(tag, filt):
        d = df.Filter(filt) if filt else df
        spec = (f"h_{tag}", f";m_{{#mu#mu}} (GeV);Events / "
                            f"{1. / BINS_PER_GEV:.2f} GeV",
                MASSBINS, XLOW, XHIGH)
        return d.Histo1D(spec, args.mass, weight) if weight \
            else d.Histo1D(spec, args.mass)

    booked = []
    if not args.no_inclusive:
        booked.append({"tag": f"inclusive_{year}", "short": "incl",
                       "label": "inclusive", "lo": None, "hi": None,
                       "filt": None, "group": None,
                       "h": histo("inclusive", None)})

    all_bins = eta_bins()
    for b in all_bins:
        tag = f"{b['tag']}_{year}"
        booked.append({**b, "tag": tag, "h": histo(tag, b["filt"])})

    # report what falls outside the edges (pT is unbounded above) instead of
    # letting the bins silently not add up to the inclusive fit
    n_all = df.Count()
    n_in = {k: df.Filter(v).Count() for k, v in coverage_filters().items()}

    # --- fit --------------------------------------------------------------
    # one observable for every fit, so the PDFs can share the overlay frame
    x = ROOT.RooRealVar(f"mh_etapt_{year}", "m_{#mu#mu}", XLOW, XHIGH)
    x.setRange("full", XLOW, XHIGH)

    fits, dump = [], {}
    for b in booked:
        print(f"\n--- fit {b['tag']} ---")
        h = b["h"].GetValue()
        h.SetDirectory(0)
        if h.GetEntries() < 100:
            print(f"  !! only {h.GetEntries():.0f} entries, skipping "
                  "(a DCB has six shape parameters)")
            continue
        res = fit_one(x, h, b["tag"], b["label"], outdir,
                      freeze=False, year=year, save_pdf=args.pdf)
        if res is None:
            continue
        pdf, bundle, norm_var, record = res
        par = record["parameters"]
        fits.append({**b, "pdf": pdf, "bundle": bundle, "norm_var": norm_var,
                     "mu": par["mu"]["value"], "muErr": par["mu"]["error"],
                     "sigma": par["sigma"]["value"],
                     "sigmaErr": par["sigma"]["error"],
                     "chi2": record["fit_quality"]["chi2_ndf"]})
        dump[b["tag"]] = {"label": b["label"], "short": b["short"],
                          "group": b["group"], "filter": b["filt"],
                          "lo": b["lo"], "hi": b["hi"], **record}

    if not fits:
        raise SystemExit("no successful fits")

    # --- report -----------------------------------------------------------
    ref = next((f for f in fits if f["label"] == "inclusive"), fits[0])
    wlab = max([34] + [len(f["label"]) + 2 for f in fits])
    width = wlab + 70
    for group, axis in [(None, "reference")] + groups:
        rows = [f for f in fits if f.get("group") == group]
        if not rows:
            continue
        print("\n" + "=" * width)
        print(f"{args.mass} split by {group or 'nothing'}, {args.sig} {year}"
              f"   (reference: {ref['label']})")
        print("-" * width)
        print(f"{'bin':<{wlab}}{'mu [GeV]':>18}{'sigma [GeV]':>18}"
              f"{'chi2/ndf':>10}{'d(mu)/mu':>11}{'d(sig)/sig':>12}")
        for f in rows:
            dmu = 0.0 if f is ref else (f["mu"] - ref["mu"]) / ref["mu"]
            dsg = 0.0 if f is ref else (f["sigma"] - ref["sigma"]) / ref["sigma"]
            print(f"{f['label']:<{wlab}}"
                  f"{f['mu']:>10.3f} +/-{f['muErr']:<6.3f}"
                  f"{f['sigma']:>10.3f} +/-{f['sigmaErr']:<6.3f}"
                  f"{f['chi2']:>10.3f}{dmu:>+11.5f}{dsg:>+12.4f}")
        print("=" * width)

    tot = n_all.GetValue()
    for k, cnt in n_in.items():
        inb = cnt.GetValue()
        if tot > inb:
            print(f"note: {tot - inb} of {tot} events "
                  f"({100. * (tot - inb) / tot:.2f} %) fall outside the "
                  f"{k} edges and are in no {k} bin")

    for group, axis in groups:
        overlay_pdfs(x, fits, outdir, year, group, axis, args.pdf)
        summary_graph(fits, outdir, year, group, axis, args.pdf)
    # --- the 2D resolution maps -------------------------------------
    print(f"\nfitting the leading-vs-subleading pT map "
          f"({args.min_entries}+ entries per cell)")
    map_fits, keep_cells = [], []
    for eb in eta_bins():
        cells = mu2d_cells(eb)
        booked_cells = [(c, histo(f"{c['tag']}_{year}", c["filt"]))
                        for c in cells]
        nok = 0
        for c, hh in booked_cells:
            h = hh.GetValue()
            h.SetDirectory(0)
            if args.cell_plots:
                r = fit_one(x, h, f"{c['tag']}_{year}", c["label"], outdir,
                            freeze=False, year=year, save_pdf=args.pdf)
                if r is None:
                    continue
                _, _, _, rec = r
                out = {"mu": rec["parameters"]["mu"]["value"],
                       "muErr": rec["parameters"]["mu"]["error"],
                       "sigma": rec["parameters"]["sigma"]["value"],
                       "sigmaErr": rec["parameters"]["sigma"]["error"],
                       "yield": rec["normalization"]["integral"],
                       "n_eff": rec["normalization"]["effective_entries"],
                       "status": rec["fit_quality"]["status"],
                       "cov_qual": rec["fit_quality"]["cov_qual"],
                       "params_at_bound": rec["fit_quality"]["params_at_bound"],
                       "ok": rec["fit_quality"]["ok"]}
            else:
                out = fit_cell(x, h, f"{c['tag']}_{year}", args.min_entries)
            if out is None:
                map_fits.append({**c, "sigma": None, "mu": None,
                                 "n_eff": None, "skipped": True})
                dump[c["tag"]] = {
                    "label": c["label"], "short": c["short"],
                    "group": c["group"], "filter": c["filt"],
                    "skipped": f"fewer than {args.min_entries} entries",
                    "fit_quality": {"ok": True},   # skipped, not failed
                }
                continue
            keep_cells.append(out.pop("_keep", None))
            nok += 1
            map_fits.append({**c, **out})
            # Same nested schema as the fit_one records, so the JSON is
            # uniform and the final "which fits need checking" scan does
            # not have to know where an entry came from.
            dump[c["tag"]] = {
                "label": c["label"], "short": c["short"],
                "group": c["group"], "filter": c["filt"],
                "parameters": {
                    "mu": {"value": out["mu"], "error": out["muErr"]},
                    "sigma": {"value": out["sigma"],
                              "error": out["sigmaErr"]},
                },
                "normalization": {"integral": out["yield"],
                                  "effective_entries": out["n_eff"]},
                "fit_quality": {"status": out["status"],
                                "cov_qual": out["cov_qual"],
                                "params_at_bound": out["params_at_bound"],
                                "ok": out["ok"]},
            }
        print(f"  {eb['short']}: {nok} of {len(cells)} cells fitted")

    resolution_maps(map_fits, outdir, year,
                    f"{args.sig} {year}", args.pdf)

    out = {"provenance": {
        "written_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_commit": git_commit(), "files": files, "year": year,
        "eras": None if args.file else eras,
        "sig": args.sig, "category": args.category,
        "mass_column": args.mass,
        "eta_edges": ETA_EDGES,
        "mu2d": {"min": MU2D_MIN, "max": MU2D_MAX, "width": MU2D_WIDTH,
                 "lead": MU2D_LEAD, "sub": MU2D_SUB,
                 "min_entries": args.min_entries},
        "groups": [g for g, _ in groups], "weight": weight,
        "fit_range": [XLOW, XHIGH], "n_bins": MASSBINS}, "fits": dump}
    jf = os.path.join(outdir, f"fitEtaBins_{args.sig}_{year}.json")
    with open(jf, "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
    print(f"\nwrote {jf}")

    bad = [t for t, r in dump.items() if not r["fit_quality"]["ok"]]
    if bad:
        print(f"WARNING: check these fits: {', '.join(bad)}")


if __name__ == "__main__":
    main()
