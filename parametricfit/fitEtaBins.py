"""
fitEtaBins.py -- double-sided Crystal Ball fits of HiggsCandCorrMass in bins
of muon eta or dimuon rapidity, for one signal production mode.

Depends only on sigFit.py and prepareFits.py:
  * sigFit.py supplies the shape, fit range, binning, the two-pass Minuit2
    fit (`fit_one`, which writes the ratio/pull canvases) and the plot style,
    so mu and sigma here are directly comparable with the signal workspaces.
  * prepareFits.py supplies the file lookup and the year handling, so the
    year may be a concrete era (12022, 22022, 12023, 22023, 2024) or a group
    (2022 = 12022 + 22022, 2023 = 12023 + 22023, Run3). A group is fitted as
    one combined sample, and every required snapshot must exist.

Binning variable (--binvar)
---------------------------
The dimuon mass depends on both muons, so "eta" has to be made specific:

  mu1eta      leading muon eta. Direct, but each bin mixes topologies because
              the subleading muon is unconstrained. (default)
  mu2eta      subleading muon eta.
  rapidity    HiggsCandCorrRapidity, a genuine per-event quantity and the
              cleanest choice for a forward/backward comparison. Note y_mumu
              is not bounded by 2.4, so events outside the edges are counted
              and reported as overflow rather than silently dropped.
  bothsame    both muons required in the same bin. Unambiguous, but discards
              the mixed-topology events, so the statistics drop.
  absetamax   max(|eta_1|,|eta_2|), the barrel/overlap/endcap split. Folds the
              two hemispheres together, so it cannot see a forward/backward
              asymmetry -- use it for resolution, not for scale.
  topology    dimuon |eta| topology. --edges are |eta| edges defining regions
              (default 0 1.4 2.4 -> B, E), and every unordered pair of regions
              is one exclusive category: BB, BE, EE. Three regions give
              B, O, E and six categories (BB BO BE OO OE EE). Each event lands
              in exactly one category, so the categories partition the
              inclusive sample. Like absetamax it folds the hemispheres.

Default edges for the signed variables are -2.4 -1.2 1.2 2.4, which tests
whether the calibration behaves the same in the +z and -z halves of the
detector.

Usage
-----
    python fitEtaBins.py 22022
    python fitEtaBins.py 2022 --binvar topology                  # BB, BE, EE
    python fitEtaBins.py 2023 --binvar topology --edges 0 0.9 1.8 2.4
    python fitEtaBins.py 2024 --sig VH --category VHcat --binvar rapidity
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

# binvar -> (expression or column, the branches it needs, axis label)
BINVAR = {
    "mu1eta":    ("Muon1_eta", ["Muon1_eta"], "#eta^{#mu1}"),
    "mu2eta":    ("Muon2_eta", ["Muon2_eta"], "#eta^{#mu2}"),
    "rapidity":  ("HiggsCandCorrRapidity", ["HiggsCandCorrRapidity"], "y_{#mu#mu}"),
    "bothsame":  (None, ["Muon1_eta", "Muon2_eta"], "#eta^{#mu1,#mu2}"),
    "absetamax": ("absEtaMax", ["Muon1_eta", "Muon2_eta"], "max|#eta^{#mu}|"),
    "topology":  (None, ["Muon1_eta", "Muon2_eta"], "|#eta| topology"),
}

DEFAULT_EDGES = {"topology": [0.0, 1.4, 2.4]}
SIGNED_DEFAULT_EDGES = [-2.4, -1.2, 1.2, 2.4]

# short region names for topology mode, by number of |eta| regions
REGION_NAMES = {2: ["B", "E"], 3: ["B", "O", "E"]}

# columns defined on the fly when a snapshot does not have them
DERIVED = {
    "absEtaMax": "(float)std::max(std::abs(Muon1_eta),std::abs(Muon2_eta))",
    "absEta1":   "(float)std::abs(Muon1_eta)",
    "absEta2":   "(float)std::abs(Muon2_eta)",
}
NEEDED_DERIVED = {"absetamax": ["absEtaMax"],
                  "topology": ["absEta1", "absEta2"]}


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
    p.add_argument("--binvar", default="mu1eta", choices=sorted(BINVAR))
    p.add_argument("--edges", type=float, nargs="+", default=None,
                   help="bin edges, ascending. For topology these are |eta| "
                        "region edges (default 0 1.4 2.4); otherwise "
                        "default -2.4 -1.2 1.2 2.4")
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

def bin_filter(binvar, lo, hi):
    if binvar == "bothsame":
        return (f"Muon1_eta >= {lo} && Muon1_eta < {hi} && "
                f"Muon2_eta >= {lo} && Muon2_eta < {hi}")
    col = BINVAR[binvar][0]
    return f"{col} >= {lo} && {col} < {hi}"


def bin_label(binvar, lo, hi):
    ax = BINVAR[binvar][2]
    if binvar == "bothsame":
        return f"{lo:g} #leq {ax} < {hi:g} (both)"
    return f"{lo:g} #leq {ax} < {hi:g}"


def bin_tag(lo, hi):
    def f(v):
        return f"{'m' if v < 0 else 'p'}{abs(v):g}".replace(".", "p")
    return f"{f(lo)}_{f(hi)}"


def _region_range(lo, hi):
    return f"|#eta| < {hi:g}" if lo == 0 else f"{lo:g} #leq |#eta| < {hi:g}"


def topology_bins(edges):
    """Exclusive dimuon categories from |eta| region edges.

    Region i is edges[i] <= |eta| < edges[i+1]. Each unordered pair (i, j)
    is one category; a mixed pair accepts either muon ordering.
    """
    n = len(edges) - 1
    names = REGION_NAMES.get(n, [f"R{i}" for i in range(n)])

    def inside(col, i):
        return f"{col} >= {edges[i]} && {col} < {edges[i + 1]}"

    bins = []
    for i in range(n):
        for j in range(i, n):
            short = names[i] + names[j]
            if i == j:
                filt = f"({inside('absEta1', i)}) && ({inside('absEta2', i)})"
                label = f"{short}: both {_region_range(edges[i], edges[i + 1])}"
            else:
                filt = (f"(({inside('absEta1', i)}) && ({inside('absEta2', j)}))"
                        f" || (({inside('absEta1', j)}) && ({inside('absEta2', i)}))")
                label = f"{short}: one {names[i]}, one {names[j]}"
            bins.append({"tag": f"topology_{short}", "short": short,
                         "label": label, "lo": None, "hi": None, "filt": filt})
    return bins


def make_bins(binvar, edges):
    """List of {tag, short, label, lo, hi, filt} (tags without the year)."""
    if binvar == "topology":
        return topology_bins(edges)
    return [{"tag": f"{binvar}_{bin_tag(lo, hi)}", "short": None,
             "label": bin_label(binvar, lo, hi), "lo": lo, "hi": hi,
             "filt": bin_filter(binvar, lo, hi)}
            for lo, hi in zip(edges, edges[1:])]


def coverage_filter(binvar, edges):
    """Events that some bin can accept, or None if not meaningful."""
    lo, hi = edges[0], edges[-1]
    if binvar == "topology":
        return (f"absEta1 >= {lo} && absEta1 < {hi} && "
                f"absEta2 >= {lo} && absEta2 < {hi}")
    col = BINVAR[binvar][0]
    return f"{col} >= {lo} && {col} < {hi}" if col else None


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


def overlay_pdfs(x, fits, plotdir, year, binvar, save_pdf=False):
    """Every bin's fitted PDF on one frame, each normalised to unity.

    All PDFs are built on the same observable x; a PDF plotted on a frame of
    a different variable would be drawn as a constant.

    Layout: one legend band across the top of the frame, each entry giving
    the bin on the first line and its mu, sigma on the second (two columns
    above three bins). The y axis is raised so the curves end below it.
    """
    legs = [f for f in fits if f["label"] != "inclusive"]
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

    c = make_canvas(f"c_overlay_{binvar}")
    frame.Draw()
    leg.Draw()
    cms_label(c, year)
    save(c, f"{plotdir}/overlay_{binvar}_{year}", save_pdf)


def summary_graph(fits, plotdir, year, binvar, save_pdf=False):
    """mu and sigma vs bin, the plot that actually answers the question."""
    legs = [f for f in fits if f["label"] != "inclusive"]
    if len(legs) < 2:
        return
    n = len(legs)
    categorical = legs[0]["lo"] is None
    ax = BINVAR[binvar][2]

    for ytit, key in [("fitted #mu (GeV)", "mu"),
                      ("fitted #sigma (GeV)", "sigma")]:
        g = ROOT.TGraphErrors(n)
        for i, f in enumerate(legs):
            if categorical:
                xc, xe = i + 0.5, 0.0
            else:
                xc = 0.5 * (f["lo"] + f["hi"])
                xe = 0.5 * (f["hi"] - f["lo"])
            g.SetPoint(i, xc, f[key])
            g.SetPointError(i, xe, f[key + "Err"])
        g.SetTitle(f";{ax};{ytit}")
        g.SetMarkerStyle(SF.DATA_MARKER_STYLE)
        g.SetMarkerSize(SF.DATA_MARKER_SIZE)
        g.SetLineWidth(SF.DATA_LINE_WIDTH)

        c = make_canvas(f"c_{key}_{binvar}")
        if categorical:
            # labelled frame: one bin per category
            lo_y = min(f[key] - f[key + "Err"] for f in legs)
            hi_y = max(f[key] + f[key + "Err"] for f in legs)
            pad_y = 0.15 * (hi_y - lo_y) or 0.05 * abs(hi_y) or 1.0
            frame = ROOT.TH1F(f"frame_{key}_{binvar}", f";{ax};{ytit}",
                              n, 0, n)
            frame.SetDirectory(0)
            for i, f in enumerate(legs):
                frame.GetXaxis().SetBinLabel(i + 1, f["short"])
            frame.GetXaxis().SetLabelSize(0.06)
            frame.SetMinimum(lo_y - pad_y)
            frame.SetMaximum(hi_y + pad_y)
            frame.GetYaxis().SetTitleOffset(1.5)
            frame.Draw()
            g.Draw("P SAME")
        else:
            g.GetYaxis().SetTitleOffset(1.5)
            g.Draw("AP")
        cms_label(c, year)
        save(c, f"{plotdir}/{key}_vs_{binvar}_{year}", save_pdf)


# ---------------------------------------------------------------------------

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

    edges = list(args.edges if args.edges is not None
                 else DEFAULT_EDGES.get(args.binvar, SIGNED_DEFAULT_EDGES))
    if len(edges) < 2 or any(b <= a for a, b in zip(edges, edges[1:])):
        raise SystemExit("--edges must be ascending and give at least one bin")
    if args.binvar in ("topology", "absetamax") and edges[0] < 0:
        raise SystemExit(f"--binvar {args.binvar} uses |eta|: edges must be >= 0")

    files = input_files(args, year)
    eras = expand_year(year, SIGNAL_YEAR_GROUPS)

    df = ROOT.RDataFrame(TREE, files)
    cols = set(str(c) for c in df.GetColumnNames())

    need = BINVAR[args.binvar][1] + [args.mass]
    missing = [c for c in need if c not in cols]
    if missing:
        raise SystemExit(f"input has no {', '.join(missing)}")

    for name in NEEDED_DERIVED.get(args.binvar, []):
        if name not in cols:
            df = df.Define(name, DERIVED[name])

    weight = None if args.weight == "none" else args.weight
    if weight and weight not in cols:
        raise SystemExit(f"input has no weight column {weight}")

    # same sanity cut as prepareFits.getHisto
    df = df.Filter(f"!std::isnan({args.mass})", "valid mass")

    # every output goes under <plotdir>/<year>/, and every fit tag carries the
    # year too, so runs for different years cannot overwrite one another
    outdir = os.path.join(args.plotdir, year)
    os.makedirs(outdir, exist_ok=True)

    print(f"sample  : {args.sig} ({args.category})"
          f"\nyear    : {year}"
          + ("" if args.file else f"  (eras: {', '.join(eras)})")
          + f"\nfiles   : {len(files)}"
          + "".join(f"\n          {f}" for f in files)
          + f"\nmass    : {args.mass}\nbinvar  : {args.binvar}"
          f"\nedges   : {edges}\nweight  : {weight or 'unweighted'}"
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
        booked.append({"tag": f"inclusive_{year}", "short": None,
                       "label": "inclusive", "lo": None, "hi": None,
                       "filt": None, "h": histo("inclusive", None)})
    for b in make_bins(args.binvar, edges):
        tag = f"{b['tag']}_{year}"
        booked.append({**b, "tag": tag, "h": histo(tag, b["filt"])})

    # report what falls outside the edges (y_mumu is unbounded, and topology
    # needs both muons inside) instead of letting the bins silently not add
    # up to the inclusive fit
    n_all = df.Count()
    cover = coverage_filter(args.binvar, edges)
    n_in = df.Filter(cover).Count() if cover else None

    # --- fit --------------------------------------------------------------
    # one observable for every fit, so the PDFs can share the overlay frame
    x = ROOT.RooRealVar(f"mh_{args.binvar}_{year}", "m_{#mu#mu}", XLOW, XHIGH)
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
                          "filter": b["filt"], "lo": b["lo"], "hi": b["hi"],
                          **record}

    if not fits:
        raise SystemExit("no successful fits")

    # --- report -----------------------------------------------------------
    ref = next((f for f in fits if f["label"] == "inclusive"), fits[0])
    wlab = max([34] + [len(f["label"]) + 2 for f in fits])
    width = wlab + 70
    print("\n" + "=" * width)
    print(f"{args.mass} in bins of {args.binvar}, {args.sig} {year}   "
          f"(reference: {ref['label']})")
    print("-" * width)
    print(f"{'bin':<{wlab}}{'mu [GeV]':>18}{'sigma [GeV]':>18}"
          f"{'chi2/ndf':>10}{'d(mu)/mu':>11}{'d(sig)/sig':>12}")
    for f in fits:
        dmu = 0.0 if f is ref else (f["mu"] - ref["mu"]) / ref["mu"]
        dsg = 0.0 if f is ref else (f["sigma"] - ref["sigma"]) / ref["sigma"]
        print(f"{f['label']:<{wlab}}"
              f"{f['mu']:>10.3f} +/-{f['muErr']:<6.3f}"
              f"{f['sigma']:>10.3f} +/-{f['sigmaErr']:<6.3f}"
              f"{f['chi2']:>10.3f}{dmu:>+11.5f}{dsg:>+12.4f}")
    print("=" * width)

    # forward/backward: the point of a signed binning
    outer = [f for f in fits if f["lo"] is not None and
             (f["lo"] <= edges[0] or f["hi"] >= edges[-1])]
    if len(outer) == 2 and outer[0]["lo"] < 0 < outer[1]["hi"]:
        a, b_ = outer[0], outer[1]
        d = b_["mu"] - a["mu"]
        e = (a["muErr"] ** 2 + b_["muErr"] ** 2) ** 0.5
        print(f"forward/backward: mu({b_['label']}) - mu({a['label']}) = "
              f"{d:+.4f} +/- {e:.4f} GeV  ({abs(d)/e if e else 0:.1f} sigma)")

    if n_in is not None:
        tot, inb = n_all.GetValue(), n_in.GetValue()
        if tot > inb:
            print(f"note: {tot - inb} of {tot} events "
                  f"({100. * (tot - inb) / tot:.2f} %) fall outside "
                  f"[{edges[0]}, {edges[-1]}] and are in no bin")

    overlay_pdfs(x, fits, outdir, year, args.binvar, args.pdf)
    summary_graph(fits, outdir, year, args.binvar, args.pdf)

    out = {"provenance": {
        "written_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_commit": git_commit(), "files": files, "year": year,
        "eras": None if args.file else eras,
        "sig": args.sig, "category": args.category,
        "mass_column": args.mass, "binvar": args.binvar,
        "edges": edges, "weight": weight,
        "fit_range": [XLOW, XHIGH], "n_bins": MASSBINS}, "fits": dump}
    jf = os.path.join(outdir, f"fitEtaBins_{args.binvar}_{args.sig}_{year}.json")
    with open(jf, "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
    print(f"\nwrote {jf}")

    bad = [t for t, r in dump.items() if not r["fit_quality"]["ok"]]
    if bad:
        print(f"WARNING: check these fits: {', '.join(bad)}")


if __name__ == "__main__":
    main()
