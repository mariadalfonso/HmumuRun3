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
Three sets, all against the same inclusive reference fit:

  eta    dimuon |eta| topology: BB, BE, EE from the region edges
         ETA_EDGES = [0, 1.4, 2.4]. Every unordered pair of regions is one
         exclusive category, so the three partition the inclusive sample.
         INCLUSIVE IN pT -- this is the original behaviour.

  pt     dimuon pT bins from PT_EDGES = [0, 50, 100, 200, 400] GeV, using
         HiggsCandCorrPt. Inclusive in eta. Resolution improves with pT for
         a fixed eta (the muons are straighter but better measured in the
         tracker), so this separates the two effects.

  pt_BB  a pT scan WITHIN each eta category, one group per region
  pt_BE  (pt_BB, pt_BE, pt_EE). This is the point of the script: it answers
  pt_EE  whether the pT dependence of the resolution is the same in BB as in
         EE, which the two inclusive scans above cannot.

A combined plot puts all three regions' pT scans on one set of axes, with
the eta-inclusive scan as a black reference. The open pT bin is omitted
there -- it has no upper edge, so no honest x position on a numeric axis --
but it is fitted and reported like every other bin.

Both splits fold the two hemispheres together, so neither can see a
forward/backward asymmetry -- they are for resolution, not for scale.

The binning is set by the ETA_EDGES / PT_EDGES constants at the top of this
file rather than by command-line options.

Usage
-----
    python fitEtaBins.py 2024
    python fitEtaBins.py 2024 --cross
    python fitEtaBins.py Run3 --sig VH --category VHcat
    python fitEtaBins.py 22022 --file <snapshot> [<snapshot> ...]
"""

import argparse
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

# overlay legend
OVERLAY_TEXT_SIZE = 0.026   # NDC; each entry is two lines of this
OVERLAY_GAP = 0.03          # NDC between the tallest curve and the legend

# ---------------------------------------------------------------------------
# the two splits. Edit here; they are deliberately not command-line options.
# ---------------------------------------------------------------------------
# |eta| region edges -> regions -> every unordered pair is one category
ETA_EDGES = [0.0, 1.4, 2.4]

# dimuon pT bin edges in GeV. With PT_OVERFLOW the last edge opens into a
# further bin holding everything above it, so the pT bins then partition the
# whole sample the way the eta categories already do.
PT_EDGES = [0.0, 20.0, 50.0, 100.0, 200.0, 400.0]
PT_OVERFLOW = True
# x position and half-width used to DRAW the open bin on the per-group pT
# summary graphs. Nominal only -- the bin has no upper edge, and the
# by-eta summary leaves it out entirely rather than drawing a made-up
# width there.
PT_OVERFLOW_PLOT_HI = 600.0
PTCOL = "HiggsCandCorrPt"

ETA_AXIS = "|#eta| topology"
PT_AXIS = "p_{T}^{#mu#mu} (GeV)"

# short region names by number of |eta| regions
REGION_NAMES = {2: ["B", "E"], 3: ["B", "O", "E"]}

# columns defined on the fly when a snapshot does not have them
DERIVED = {
    "absEta1": "(float)std::abs(Muon1_eta)",
    "absEta2": "(float)std::abs(Muon2_eta)",
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
    p.add_argument("--mass", default=MASSCOL,
                   help="mass column to fit (default: %(default)s)")
    p.add_argument("--weight", default="w_allSF",
                   help="weight column, or 'none' for unweighted")
    p.add_argument("-o", "--plotdir",
                   default=f"/home/submit/{getpass.getuser()}/public_html/HmumuRun3/EtaBinFits",
                   help="base output dir. Everything is written to a "
                        "<plotdir>/<year>/ subdirectory, so runs for different "
                        "years never overwrite each other.")
    p.add_argument("--no-inclusive", action="store_true",
                   help="skip the inclusive fit used as the delta reference")
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


def _pt_tag(lo, hi):
    return f"{lo:g}_{hi:g}".replace(".", "p")


def pt_bins(edges=None, overflow=None):
    """Dimuon pT bins, inclusive in eta.

    With `overflow` the final bin is open above the last edge. Without it,
    high-pT events land in no bin and are reported as uncovered.
    """
    edges = list(PT_EDGES if edges is None else edges)
    overflow = PT_OVERFLOW if overflow is None else overflow
    bins = []
    for lo, hi in zip(edges, edges[1:]):
        bins.append({"tag": f"pt_{_pt_tag(lo, hi)}",
                     "short": f"{lo:g}-{hi:g}",
                     "label": f"{lo:g} #leq {PT_AXIS} < {hi:g}",
                     "lo": lo, "hi": hi,
                     "filt": f"{PTCOL} >= {lo} && {PTCOL} < {hi}",
                     "group": "pt"})
    if overflow:
        lo = edges[-1]
        bins.append({"tag": f"pt_{_pt_tag(lo, 0)}".replace("_0", "_inf"),
                     "short": f">{lo:g}",
                     "label": f"{PT_AXIS} #geq {lo:g}",
                     # no upper edge: these are for drawing only
                     "lo": lo, "hi": PT_OVERFLOW_PLOT_HI, "open": True,
                     "filt": f"{PTCOL} >= {lo}",
                     "group": "pt"})
    return bins


def cross_bins(eta_edges=None, pt_edges=None):
    """A pT scan WITHIN each |eta| topology.

    Each eta region gets its own group ("pt_BB", "pt_BE", "pt_EE"), so every
    region produces its own overlay and its own mu/sigma-versus-pT graphs
    rather than all 15 bins landing on one categorical axis.
    """
    out = []
    for e in eta_bins(eta_edges):
        for p in pt_bins(pt_edges):
            out.append({
                "tag": f"{e['tag']}_{p['tag']}",
                "short": p["short"],          # the pT bin, within this region
                "label": f"{e['short']}, {p['label']}",
                "lo": p["lo"], "hi": p["hi"], "open": p.get("open", False),
                "filt": f"({e['filt']}) && ({p['filt']})",
                "group": f"pt_{e['short']}", "eta_short": e["short"]})
    return out


def coverage_filters():
    """Events some bin of each split can accept, for the overflow report.

    The pT entry is dropped when PT_OVERFLOW is on, because the open bin
    then leaves nothing uncovered above the last edge.
    """
    e0, e1 = ETA_EDGES[0], ETA_EDGES[-1]
    out = {"eta": (f"absEta1 >= {e0} && absEta1 < {e1} && "
                   f"absEta2 >= {e0} && absEta2 < {e1}")}
    if not PT_OVERFLOW:
        out["pt"] = (f"{PTCOL} >= {PT_EDGES[0]} && "
                     f"{PTCOL} < {PT_EDGES[-1]}")
    return out


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

def summary_by_eta(fits, plotdir, year, save_pdf=False):
    """mu and sigma versus pT, one series per |eta| category.

    The single most useful output: three curves on one set of axes show
    directly whether the pT dependence of the resolution is the same in the
    barrel as in the endcap. The eta-inclusive scan is drawn in black as a
    reference.
    """
    # The open bin has no upper edge, so it has no honest x position or
    # width on a numeric axis: left out here. It is still fitted, and still
    # appears in the tables, the JSON and its own ratio/pull canvases.
    def closed(rows):
        return [f for f in rows if not f.get("open")]

    shorts = [b["short"] for b in eta_bins()]
    series = [(sh, closed([f for f in fits
                           if f.get("group") == f"pt_{sh}"]))
              for sh in shorts]
    series = [(sh, rows) for sh, rows in series if len(rows) >= 2]
    if not series:
        return
    ref_rows = closed([f for f in fits if f.get("group") == "pt"])

    for ytit, key in [("fitted #mu (GeV)", "mu"),
                      ("fitted #sigma (GeV)", "sigma")]:
        graphs, keep = [], []

        def build(rows, col, style):
            g = ROOT.TGraphErrors(len(rows))
            for i, f in enumerate(rows):
                g.SetPoint(i, 0.5 * (f["lo"] + f["hi"]), f[key])
                g.SetPointError(i, 0.5 * (f["hi"] - f["lo"]),
                                f[key + "Err"])
            g.SetMarkerStyle(style)
            g.SetMarkerSize(SF.DATA_MARKER_SIZE)
            g.SetMarkerColor(col)
            g.SetLineColor(col)
            g.SetLineWidth(SF.DATA_LINE_WIDTH)
            return g

        if ref_rows:
            graphs.append(("inclusive in #eta",
                           build(ref_rows, ROOT.kBlack, 24)))
        for i, (sh, rows) in enumerate(series):
            graphs.append((sh, build(rows, color(i),
                                     SF.DATA_MARKER_STYLE)))

        # frame spanning every point, with headroom for the legend
        vals = [f[key] for _, rows in series for f in rows] + \
               [f[key] for f in ref_rows]
        errs = [f[key + "Err"] for _, rows in series for f in rows] + \
               [f[key + "Err"] for f in ref_rows]
        lo_y = min(v - e for v, e in zip(vals, errs))
        hi_y = max(v + e for v, e in zip(vals, errs))
        span = (hi_y - lo_y) or 0.05 * abs(hi_y) or 1.0

        xmax = max(f["hi"] for _, rows in series for f in rows)
        frame = ROOT.TH1F(f"frame_{key}_byeta", f";{PT_AXIS};{ytit}",
                          1, PT_EDGES[0], xmax)
        frame.SetDirectory(0)
        frame.SetMinimum(lo_y - 0.10 * span)
        frame.SetMaximum(hi_y + 0.45 * span)     # room for the legend
        frame.GetYaxis().SetTitleOffset(1.5)

        c = make_canvas(f"c_{key}_byeta")
        frame.Draw()
        for _, g in graphs:
            g.Draw("P SAME")
        keep += [frame] + [g for _, g in graphs]

        leg = ROOT.TLegend(SF.PAD_LEFT_MARGIN + 0.04,
                           1 - SF.PAD_TOP_MARGIN - 0.04 - 0.055 * len(graphs),
                           SF.PAD_LEFT_MARGIN + 0.34,
                           1 - SF.PAD_TOP_MARGIN - 0.04)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(42)
        leg.SetTextSize(0.030)
        leg.SetHeader("|#eta| topology")
        for lab, g in graphs:
            leg.AddEntry(g, lab, "lp")
        leg.Draw()
        keep.append(leg)

        cms_label(c, year)
        save(c, f"{plotdir}/{key}_vs_pt_byeta_{year}", save_pdf)


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
    args = parse_args()
    SF.setup_style()
    year = normalize_year(args.year)

    files = input_files(args, year)
    eras = expand_year(year, SIGNAL_YEAR_GROUPS)

    df = ROOT.RDataFrame(TREE, files)
    cols = set(str(c) for c in df.GetColumnNames())

    need = ["Muon1_eta", "Muon2_eta", PTCOL, args.mass]
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
    outdir = os.path.join(args.plotdir, year)
    os.makedirs(outdir, exist_ok=True)

    # eta categories (inclusive in pT), pT bins (inclusive in eta), then a
    # pT scan inside each eta category
    eta_shorts = [b["short"] for b in eta_bins()]
    groups = [("eta", ETA_AXIS), ("pt", PT_AXIS)]
    groups += [(f"pt_{sh}", f"{PT_AXIS}, {sh}") for sh in eta_shorts]

    print(f"sample  : {args.sig} ({args.category})"
          f"\nyear    : {year}"
          + ("" if args.file else f"  (eras: {', '.join(eras)})")
          + f"\nfiles   : {len(files)}"
          + "".join(f"\n          {f}" for f in files)
          + f"\nmass    : {args.mass}"
          f"\neta     : {ETA_EDGES}  -> "
          f"{', '.join(b['short'] for b in eta_bins())}"
          f"\npT      : {PT_EDGES} GeV ({PTCOL})"
          + (f" + open bin >{PT_EDGES[-1]:g}" if PT_OVERFLOW else "")
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

    all_bins = eta_bins() + pt_bins() + cross_bins()
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
    summary_by_eta(fits, outdir, year, args.pdf)

    out = {"provenance": {
        "written_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_commit": git_commit(), "files": files, "year": year,
        "eras": None if args.file else eras,
        "sig": args.sig, "category": args.category,
        "mass_column": args.mass, "pt_column": PTCOL,
        "eta_edges": ETA_EDGES, "pt_edges": PT_EDGES,
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
