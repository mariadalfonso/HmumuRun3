"""
massResVsMass.py -- relative mass resolution versus invariant mass.

Plots the per-event RELATIVE mass error, sigma_m/m, against the invariant
mass, for one or more signal samples and for both mass definitions in the
snapshots:

    preFSR   HiggsCandMass      / HiggsCandMassErr
    postFSR  HiggsCandCorrMass  / HiggsCandCorrMassErr

A variant is a PAIR: mixing the post-FSR mass with the pre-FSR error would
be wrong, so they are never selected independently.

    sigma_m / m = 0.5 * sqrt((sigma_1/pT1)^2 + (sigma_2/pT2)^2)

The 0.5 is there because each momentum enters the mass under a square
root. Multiplying by m gives the absolute error in GeV, which is what is
directly comparable with a fitted cb_sigma, so both are available here:
--absolute switches the y axis over.

What it draws
-------------
  <tag>_2d        2D map of sigma_m/m versus m, one per sample, with the
                  median and the 16-84% band drawn over it. The map shows
                  the correlation; the band is what you would quote.

  profile         median sigma_m/m versus m, one curve per sample, on one
                  set of axes. Written only when there is MORE THAN ONE
                  sample (or with --force-profile): for a single sample the
                  same median and band are already drawn on its 2D map, and
                  the pre/post-FSR comparison is covered by fsr_compare.

  <tag>_1d        distribution of sigma_m/m integrated over the window,
                  one curve per sample, with the median marked.

  fsr_compare     pre- versus post-FSR median for each sample (dashed vs
                  solid). With --variant both, this is the plot the whole
                  exercise is for: FSR recovery puts the radiated photon
                  back into the mass, so the solid curve should sit BELOW
                  the dashed one. The console also prints the change in
                  the integrated median per sample.

Why the mass dependence matters
-------------------------------
A signal shape fitted over 110-150 GeV assumes the resolution is constant
across the window. sigma_m/m is roughly flat in m by construction (it is a
ratio), but the ABSOLUTE error is not -- it grows with m. If the relative
error also slopes, the constant-width assumption is worse than it looks,
and the residual oscillation in the single-DSCB pull may be partly this.

Usage
-----
    python massResVsMass.py 2024                      # both variants
    python massResVsMass.py 2024 --sig ggH qqH VH ttH
    python massResVsMass.py 2024 --variant postFSR
    python massResVsMass.py Run3 --sig ggH --category VBFcat
    python massResVsMass.py 2024 --absolute
    python massResVsMass.py 2024 --file <snapshot> [<snapshot> ...]
"""

import argparse
from array import array
import json
import os
from datetime import datetime, timezone

import ROOT

import sigFit as SF
from sigFit import git_commit, XLOW, XHIGH
from prepareFits import (ROOTFILES_DIR, SELMVA, SIGNAL_TAGS,
                         SIGNAL_YEAR_GROUPS, expand_year, get_signal_files,
                         normalize_year)

ROOT.gROOT.SetBatch(True)
ROOT.TH1.SetDefaultSumw2(True)

TREE = "events"

# The two mass definitions in the snapshots. Each has its own relative-error
# branch, so a variant is a (mass, relative error) PAIR -- mixing the
# post-FSR mass with the pre-FSR error, or vice versa, would be wrong.
VARIANTS = {
    "preFSR":  ("HiggsCandMass", "HiggsCandMassErr",
                "pre FSR corr."),
    "postFSR": ("HiggsCandCorrMass", "HiggsCandCorrMassErr",
                "post FSR corr."),
}
VARIANT_ORDER = ["preFSR", "postFSR"]

# binning of the 2D map
MASS_NBINS = 80                      # 0.5 GeV over 110-150
REL_NBINS = 100
REL_MAX = 0.05                       # sigma_m/m; ~0.012 is typical
ABS_MAX = 6.0                        # GeV, when --absolute

# quantiles drawn as the band, and the one drawn as the central curve
Q_LO, Q_MID, Q_HI = 0.16, 0.50, 0.84

PALETTE_RGB = [
    SF.FIT_COLOR_RGB,
    (31, 78, 158),
    (34, 139, 34),
    (230, 159, 0),
    (128, 64, 160),
    (0, 158, 176),
]

CMS_EXTRA_OFFSET = 2.4

# cms_label also reports the NDC x just past "Internal", this many
# CMS-text-widths from the left margin. Not currently used -- the 2D
# annotation sits inside the frame, because the header line has room only
# for "CMS Internal" and the luminosity.
ANN_OFFSET = 6.2

# The header is sized from the pad top margin, and ROOT scales NDC text by
# the pad's SHORTER side. Every canvas here is square, which makes the
# header noticeably larger than on sigFit's 600x420 fit pads, so scale it
# down. Applied to the 1D plots, which keep sigFit's one-line header.
CMS_SCALE = 0.72

# The 2D map puts EVERYTHING above the frame instead, in two lines with a
# left- and a right-aligned item each (the ATLAS convention). Nothing then
# sits on top of the colour map, which is the only reliable way to keep it
# legible whatever the palette and wherever the events fall. The price is a
# taller top margin.
# Mass range DRAWN on the 2D map. The histograms are still booked over the
# full [XLOW, XHIGH] window, so the quantiles, the profiles and the JSON are
# unaffected -- only the view zooms. Outside ~115-135 the map is a sparse
# scatter that compresses the interesting region.
MAP_MASS_RANGE = (115.0, 135.0)

# Log z on the 2D map. The occupancy spans orders of magnitude -- the peak
# cell holds a few hundred times what a tail cell does -- so on a linear
# scale everything outside the core reads as the lowest colour. Log z shows
# the tail structure instead, at the cost of making the core look less
# dominant than it is.
MAP_LOGZ = False

HDR_TOP_MARGIN = 0.135
HDR_CMS_SIZE = 0.046
HDR_SUB_SIZE = 0.032

# z-axis palette for the 2D map. Any ROOT EColorPalette: kLightTemperature,
# kBird, kViridis, kCividis, ...
MAP_PALETTE = ROOT.kLightTemperature


def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("year",
                   help="era (12022, 22022, 12023, 22023, 2024) or group "
                        "(2022, 2023, Run3)")
    p.add_argument("--sig", nargs="+", default=["ggH"],
                   choices=sorted(SIGNAL_TAGS),
                   help="one or more production modes; several are overlaid "
                        "(default: %(default)s)")
    p.add_argument("--category", default="ggHcat", choices=sorted(SELMVA),
                   help="snapshot category directory (default: %(default)s)")
    p.add_argument("--rootfiles-dir", default=ROOTFILES_DIR)
    p.add_argument("--file", nargs="+", default=None,
                   help="fit these files instead of the lookup; only valid "
                        "with a single --sig, since the files cannot be "
                        "attributed to a mode")
    p.add_argument("--variant", default="both",
                   choices=VARIANT_ORDER + ["both"],
                   help="which mass definition: preFSR "
                        "(HiggsCandMass / HiggsCandMassErr), postFSR "
                        "(HiggsCandCorrMass / HiggsCandCorrMassErr), or "
                        "both plus the comparison (default: %(default)s)")
    p.add_argument("--mass", default=None,
                   help="override the mass branch; requires --relerr too, "
                        "and replaces --variant with a single custom one")
    p.add_argument("--relerr", default=None,
                   help="override the relative-error branch")
    p.add_argument("--absolute", action="store_true",
                   help="plot m * sigma_m/m in GeV instead of the ratio")
    p.add_argument("--weight", default="w_allSF",
                   help="weight column, or 'none'")
    p.add_argument("--cut", default=None,
                   help="extra selection, e.g. 'discrMVA>0.8'")
    p.add_argument("-o", "--plotdir",
                   default=os.path.expanduser(
                       "~/public_html/HmumuFits/massResVsMass"))
    p.add_argument("--logz", action="store_true",
                   help="log z axis on the 2D map, to see the tail "
                        "occupancy rather than only the core")
    p.add_argument("--map-mass-range", nargs=2, type=float, default=None,
                   metavar=("LO", "HI"),
                   help=f"mass range drawn on the 2D map (default "
                        f"{MAP_MASS_RANGE[0]:g} {MAP_MASS_RANGE[1]:g}); the "
                        f"histograms are still filled over the full window, "
                        f"so only the view changes. Pass 0 0 for no zoom.")
    p.add_argument("--force-profile", action="store_true",
                   help="write the standalone median-vs-mass plot even for a "
                        "single sample, where it duplicates the curves "
                        "already drawn on the 2D map")
    p.add_argument("--palette", default=None,
                   help="ROOT palette for the 2D map, e.g. kLightTemperature, "
                        "kBird, kViridis (default: MAP_PALETTE)")
    p.add_argument("--pdf", action="store_true", help="also write .pdf")
    return p.parse_args()


# ---------------------------------------------------------------------------
# style (from sigFit, as in fitEtaBins.py)
# ---------------------------------------------------------------------------

def color(i):
    return ROOT.TColor.GetColor(*PALETTE_RGB[i % len(PALETTE_RGB)])


def make_canvas(name, width=None, right=None):
    c = ROOT.TCanvas(name, "", width or SF.CANVAS_W, SF.CANVAS_H)
    c.SetTopMargin(SF.PAD_TOP_MARGIN)
    c.SetBottomMargin(SF.PAD_BOTTOM_MARGIN)
    c.SetLeftMargin(SF.PAD_LEFT_MARGIN)
    c.SetRightMargin(SF.PAD_RIGHT_MARGIN if right is None else right)
    c.SetTickx()
    c.SetTicky()
    return c


def cms_label(pad, year, scale=1.0):
    """sigFit.cms_label with the "Internal" offset scaled to this pad.

    `scale` shrinks the whole header, to free room on the same line for
    another label. Returns (latex, x_end, y) so the caller can place text
    just after "Internal" without guessing.
    """
    pad.cd()
    t, l, r = pad.GetTopMargin(), pad.GetLeftMargin(), pad.GetRightMargin()
    y = 1 - t + 0.2 * t
    size = SF.CMS_TEXT_SIZE_FRAC * t * scale

    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextColor(ROOT.kBlack)
    tex.SetTextFont(61)
    tex.SetTextAlign(11)
    tex.SetTextSize(size)
    tex.DrawLatex(l, y, SF.CMS_TEXT)

    if SF.CMS_EXTRA_TEXT:
        w = pad.GetWw() * pad.GetAbsWNDC()
        h = pad.GetWh() * pad.GetAbsHNDC()
        dx = CMS_EXTRA_OFFSET * size * min(w, h) / w if w > 0 else 0.10
        tex.SetTextFont(52)
        tex.SetTextSize(SF.EXTRA_OVER_CMS_TEXT_SIZE * size)
        tex.DrawLatex(l + dx, y, SF.CMS_EXTRA_TEXT)

    lumi = SF.get_lumi(year) if year else 0.0
    tex.SetTextFont(42)
    tex.SetTextAlign(31)
    tex.SetTextSize(SF.LUMI_TEXT_SIZE_FRAC * t * scale)
    tex.DrawLatex(1 - r, y,
                  ("%.1f fb^{-1} (%.1f TeV)" % (lumi, SF.SQRT_S_TEV))
                  if lumi > 0 else "(%.1f TeV)" % SF.SQRT_S_TEV)

    w = pad.GetWw() * pad.GetAbsWNDC()
    h = pad.GetWh() * pad.GetAbsHNDC()
    scl = min(w, h) / w if w > 0 else 1.0
    return tex, l + ANN_OFFSET * size * scl, y


def header_2d(pad, year, left2="", right2=""):
    """Two header lines ABOVE the frame.

        CMS Internal                        <lumi> (13.6 TeV)
        <left2>                                        <right2>

    Returns the TLatex so the caller keeps it alive.
    """
    pad.cd()
    t, l, r = pad.GetTopMargin(), pad.GetLeftMargin(), pad.GetRightMargin()
    top = 1.0 - t                       # NDC y of the frame's top edge
    # The two lines sit close together: the second is a caption for the
    # first, not a separate block. With TextAlign 11 the y is the baseline,
    # so line 1 rises HDR_CMS_SIZE above y1 -- the gap must stay above that
    # or the two will touch.
    y1 = top + 0.42 * t
    y2 = top + 0.10 * t

    w = pad.GetWw() * pad.GetAbsWNDC()
    h = pad.GetWh() * pad.GetAbsHNDC()
    scl = min(w, h) / w if w > 0 else 1.0

    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextColor(ROOT.kBlack)

    tex.SetTextFont(61)
    tex.SetTextAlign(11)
    tex.SetTextSize(HDR_CMS_SIZE)
    tex.DrawLatex(l, y1, SF.CMS_TEXT)
    if SF.CMS_EXTRA_TEXT:
        tex.SetTextFont(52)
        tex.SetTextSize(SF.EXTRA_OVER_CMS_TEXT_SIZE * HDR_CMS_SIZE)
        tex.DrawLatex(l + CMS_EXTRA_OFFSET * HDR_CMS_SIZE * scl, y1,
                      SF.CMS_EXTRA_TEXT)

    lumi = SF.get_lumi(year) if year else 0.0
    tex.SetTextFont(42)
    tex.SetTextAlign(31)
    tex.SetTextSize(HDR_SUB_SIZE)
    tex.DrawLatex(1 - r, y1,
                  ("%.1f fb^{-1} (%.1f TeV)" % (lumi, SF.SQRT_S_TEV))
                  if lumi > 0 else "(%.1f TeV)" % SF.SQRT_S_TEV)

    if left2:
        tex.SetTextAlign(11)
        tex.DrawLatex(l, y2, left2)
    if right2:
        tex.SetTextAlign(31)
        tex.DrawLatex(1 - r, y2, right2)
    return tex


def style_axes(obj, xoff=1.15, yoff=1.55):
    for ax, off in ((obj.GetXaxis(), xoff), (obj.GetYaxis(), yoff)):
        ax.SetTitleFont(43)
        ax.SetLabelFont(43)
        ax.SetTitleSize(24)
        ax.SetLabelSize(20)
        ax.SetTitleOffset(off)


def save(c, out, save_pdf):
    os.makedirs(os.path.dirname(out), exist_ok=True)
    c.SaveAs(out + ".png")
    if save_pdf:
        c.SaveAs(out + ".pdf")


# ---------------------------------------------------------------------------
# quantile profile
# ---------------------------------------------------------------------------

def quantile_profile(h2):
    """Per-mass-column quantiles of the y distribution.

    Returns (centres, lo, mid, hi, n) with one entry per x bin that has
    entries. Quantiles rather than a TProfile mean: the relative-error
    distribution has a long high tail, so the mean sits above the bulk and
    the median is the honest "typical" value.
    """
    xs, los, mids, his, ns = [], [], [], [], []
    probs = array("d", [Q_LO, Q_MID, Q_HI])
    qs = array("d", [0.0, 0.0, 0.0])
    for ix in range(1, h2.GetNbinsX() + 1):
        proj = h2.ProjectionY(f"_py{ix}", ix, ix)
        n = proj.GetEntries()
        if n < 20 or proj.Integral() <= 0:
            proj.Delete()
            continue
        proj.GetQuantiles(3, qs, probs)
        xs.append(h2.GetXaxis().GetBinCenter(ix))
        los.append(qs[0])
        mids.append(qs[1])
        his.append(qs[2])
        ns.append(n)
        proj.Delete()
    return xs, los, mids, his, ns


def band_graph(xs, los, his, col):
    """Shaded 16-84% band as a closed polygon."""
    n = len(xs)
    g = ROOT.TGraph(2 * n + 1)
    for i in range(n):
        g.SetPoint(i, xs[i], his[i])
    for i in range(n):
        g.SetPoint(n + i, xs[n - 1 - i], los[n - 1 - i])
    g.SetPoint(2 * n, xs[0], his[0])
    g.SetFillColorAlpha(col, 0.25)
    g.SetLineWidth(0)
    return g


def median_graph(xs, mids, col, style=20):
    g = ROOT.TGraph(len(xs))
    for i, (xx, mm) in enumerate(zip(xs, mids)):
        g.SetPoint(i, xx, mm)
    g.SetLineColor(col)
    g.SetLineWidth(3)
    g.SetMarkerColor(col)
    g.SetMarkerStyle(style)
    g.SetMarkerSize(0.7)
    return g


# ---------------------------------------------------------------------------
# plots
# ---------------------------------------------------------------------------

def plot_2d(entry, plotdir, year, label, ytit, save_pdf):
    """The 2D map for one sample, with its median and band drawn over."""
    h2 = entry["h2"]
    h2.SetStats(0)
    h2.SetTitle(f";m_{{#mu#mu}} (GeV);{ytit};events")
    style_axes(h2)
    h2.GetZaxis().SetLabelFont(43)
    h2.GetZaxis().SetLabelSize(18)

    if MAP_MASS_RANGE:
        h2.GetXaxis().SetRangeUser(*MAP_MASS_RANGE)

    c = make_canvas(f"c2d_{entry['sig']}_{entry['variant']}", right=0.16)
    c.SetTopMargin(HDR_TOP_MARGIN)
    if MAP_LOGZ:
        # An empty bin is 0, which has no place on a log axis: ROOT then
        # drops it to the underflow colour, indistinguishable from a bin
        # holding one event. Setting the floor just below 1 keeps the two
        # apart -- empty stays blank, one event gets the lowest colour.
        c.SetLogz()
        h2.SetMinimum(0.5)
    ROOT.gStyle.SetPalette(MAP_PALETTE)
    h2.Draw("COLZ")

    keep = []
    xs, los, mids, his, _ = entry["prof"]
    if xs:
        band = band_graph(xs, los, his, ROOT.kBlack)
        band.SetFillColorAlpha(ROOT.kBlack, 0.0)   # outline only over COLZ
        med = median_graph(xs, mids, ROOT.kBlack)
        med.SetLineStyle(ROOT.kSolid)
        med.Draw("L SAME")
        for vals, sty in ((los, ROOT.kDashed), (his, ROOT.kDashed)):
            g = median_graph(xs, vals, ROOT.kBlack)
            g.SetLineStyle(sty)
            g.SetLineWidth(2)
            g.SetMarkerSize(0)
            g.Draw("L SAME")
            keep.append(g)
        keep += [band, med]

    # The header line only fits "CMS Internal" and the luminosity -- the
    # gap between them is about 0.19 NDC, not enough for the sample too.
    # Everything else goes INSIDE the frame, upper left, where the
    # sigma_m/m distribution is sparse.
    # Sample only: the category is already in the output path, and a
    # shorter caption keeps the whole second line on one row.
    keep.append(header_2d(
        c, year,
        left2=f"{entry['sig']}, {entry['desc']}",
        # "+/-1 sigma" is the usual shorthand for the 16-84% range. It is
        # exact only for a Gaussian, and this distribution has a long high
        # tail, so the band is not symmetric about the median.
        right2="solid: median,  dashed: #pm1#sigma"))

    sfx = "_logz" if MAP_LOGZ else ""
    save(c, f"{plotdir}/{entry['sig']}_2d_{entry['variant']}{sfx}_{year}",
         save_pdf)


def plot_profile(entries, plotdir, year, label, ytit, save_pdf,
                 suffix=""):
    """Median versus mass, all samples on one set of axes."""
    live = [e for e in entries if e["prof"][0]]
    if not live:
        return

    allv = [v for e in live for v in e["prof"][1] + e["prof"][3]]
    lo, hi = min(allv), max(allv)
    span = (hi - lo) or 0.1 * abs(hi) or 1.0

    frame = ROOT.TH1F(f"prof_frame{suffix}", f";m_{{#mu#mu}} (GeV);{ytit}",
                      1, XLOW, XHIGH)
    frame.SetDirectory(0)
    frame.SetStats(0)
    frame.SetMinimum(lo - 0.15 * span)
    frame.SetMaximum(hi + 0.35 * span)
    style_axes(frame)

    c = make_canvas(f"cprof{suffix}")
    frame.Draw()

    keep = [frame]
    leg = ROOT.TLegend(SF.PAD_LEFT_MARGIN + 0.04,
                       1 - SF.PAD_TOP_MARGIN - 0.04 - 0.05 * len(live),
                       SF.PAD_LEFT_MARGIN + 0.40,
                       1 - SF.PAD_TOP_MARGIN - 0.04)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.028)

    for i, e in enumerate(live):
        xs, los, mids, his, _ = e["prof"]
        col = color(i)
        band = band_graph(xs, los, his, col)
        band.Draw("F SAME")
        med = median_graph(xs, mids, col)
        med.Draw("L SAME")
        keep += [band, med]
        # slope over the window, from the end points of the median
        slope = (mids[-1] - mids[0]) / (xs[-1] - xs[0]) if len(xs) > 1 else 0.0
        leg.AddEntry(med, f"{e['sig']}  (median "
                          f"{mids[len(mids) // 2]:.4g})", "l")
        e["slope_per_gev"] = slope

    leg.Draw()
    keep.append(leg)

    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextFont(42)
    tex.SetTextSize(0.026)
    tex.DrawLatex(1 - SF.PAD_RIGHT_MARGIN - 0.40,
                  SF.PAD_BOTTOM_MARGIN + 0.05,
                  "line: median   band: 16-84%")
    tex.DrawLatex(1 - SF.PAD_RIGHT_MARGIN - 0.40,
                  SF.PAD_BOTTOM_MARGIN + 0.02, label)
    keep.append(tex)

    keep.append(cms_label(c, year, scale=CMS_SCALE)[0])
    save(c, f"{plotdir}/profile{suffix}_{year}", save_pdf)


def plot_1d(entries, plotdir, year, label, xtit, save_pdf,
            suffix=""):
    """Distribution integrated over the mass window, all samples overlaid."""
    live = [e for e in entries if e["h1"].Integral() > 0]
    if not live:
        return

    c = make_canvas(f"c1d{suffix}")
    keep = []
    leg = ROOT.TLegend(1 - SF.PAD_RIGHT_MARGIN - 0.40,
                       1 - SF.PAD_TOP_MARGIN - 0.06 - 0.05 * len(live),
                       1 - SF.PAD_RIGHT_MARGIN - 0.02,
                       1 - SF.PAD_TOP_MARGIN - 0.06)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.028)

    ymax = 0.0
    for e in live:
        h = e["h1"]
        if h.Integral() > 0:
            h.Scale(1.0 / h.Integral())
        ymax = max(ymax, h.GetMaximum())

    for i, e in enumerate(live):
        h = e["h1"]
        h.SetStats(0)
        h.SetLineColor(color(i))
        h.SetLineWidth(3)
        h.SetTitle(f";{xtit};fraction of events")
        style_axes(h)
        h.SetMaximum(1.35 * ymax)
        h.Draw("HIST" if i == 0 else "HIST SAME")
        leg.AddEntry(h, f"{e['sig']}  (median {e['median']:.4g})", "l")
        keep.append(h)

        ln = ROOT.TLine(e["median"], 0.0, e["median"], 1.05 * ymax)
        ln.SetLineColor(color(i))
        ln.SetLineStyle(ROOT.kDashed)
        ln.SetLineWidth(2)
        ln.Draw("SAME")
        keep.append(ln)

    leg.Draw()
    keep.append(leg)

    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextFont(42)
    tex.SetTextSize(0.028)
    tex.DrawLatex(SF.PAD_LEFT_MARGIN + 0.04, 0.86, label)
    tex.SetTextSize(0.024)
    tex.DrawLatex(SF.PAD_LEFT_MARGIN + 0.04, 0.82,
                  "dashed: median")
    keep.append(tex)

    keep.append(cms_label(c, year, scale=CMS_SCALE)[0])
    save(c, f"{plotdir}/dist{suffix}_{year}", save_pdf)


def plot_fsr_compare(entries, plotdir, year, label, ytit, save_pdf):
    """Pre- versus post-FSR median resolution, one colour per sample.

    Dashed is before FSR recovery, solid after. This is the plot the two
    variants exist for: FSR recovery adds the radiated photon back into the
    mass, so it should pull the low-mass tail back under the peak and
    IMPROVE the resolution. If the solid curve does not sit below the
    dashed one, the recovery is not doing what it should.

    Only the medians are drawn -- overlaying four bands would be unreadable.
    """
    samples = []
    for e in entries:
        if e["sig"] not in samples:
            samples.append(e["sig"])
    pairs = []
    for sig in samples:
        pre = next((e for e in entries
                    if e["sig"] == sig and e["variant"] == "preFSR"
                    and e["prof"][0]), None)
        post = next((e for e in entries
                     if e["sig"] == sig and e["variant"] == "postFSR"
                     and e["prof"][0]), None)
        if pre and post:
            pairs.append((sig, pre, post))
    if not pairs:
        return

    allv = [v for _, a, b in pairs for e in (a, b) for v in e["prof"][2]]
    lo, hi = min(allv), max(allv)
    span = (hi - lo) or 0.1 * abs(hi) or 1.0

    frame = ROOT.TH1F("fsr_frame", f";m_{{#mu#mu}} (GeV);{ytit}",
                      1, XLOW, XHIGH)
    frame.SetDirectory(0)
    frame.SetStats(0)
    frame.SetMinimum(lo - 0.15 * span)
    frame.SetMaximum(hi + 0.40 * span)
    style_axes(frame)

    c = make_canvas("cfsr")
    frame.Draw()
    keep = [frame]

    leg = ROOT.TLegend(SF.PAD_LEFT_MARGIN + 0.04,
                       1 - SF.PAD_TOP_MARGIN - 0.04 - 0.048 * (2 * len(pairs)),
                       SF.PAD_LEFT_MARGIN + 0.46,
                       1 - SF.PAD_TOP_MARGIN - 0.04)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.026)

    for i, (sig, pre, post) in enumerate(pairs):
        col = color(i)
        for e, style in ((pre, ROOT.kDashed), (post, ROOT.kSolid)):
            xs, _, mids, _, _ = e["prof"]
            g = median_graph(xs, mids, col)
            g.SetLineStyle(style)
            g.SetMarkerSize(0)
            g.Draw("L SAME")
            keep.append(g)
            med = mids[len(mids) // 2]
            leg.AddEntry(g, f"{sig}, {e['variant']}  ({med:.4g})", "l")

        # improvement at the peak, from the middle of each profile
        m_pre = pre["prof"][2][len(pre["prof"][2]) // 2]
        m_post = post["prof"][2][len(post["prof"][2]) // 2]
        pre["fsr_gain"] = None
        post["fsr_gain"] = (m_post - m_pre) / m_pre if m_pre else None

    leg.Draw()
    keep.append(leg)

    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextFont(42)
    tex.SetTextSize(0.026)
    tex.DrawLatex(1 - SF.PAD_RIGHT_MARGIN - 0.44,
                  SF.PAD_BOTTOM_MARGIN + 0.05,
                  "dashed: before FSR   solid: after FSR")
    tex.DrawLatex(1 - SF.PAD_RIGHT_MARGIN - 0.44,
                  SF.PAD_BOTTOM_MARGIN + 0.02, label)
    keep.append(tex)

    keep.append(cms_label(c, year, scale=CMS_SCALE)[0])
    save(c, f"{plotdir}/fsr_compare_{year}", save_pdf)


# ---------------------------------------------------------------------------

def main():
    global MAP_PALETTE, MAP_MASS_RANGE, MAP_LOGZ
    args = parse_args()
    MAP_LOGZ = args.logz or MAP_LOGZ
    if args.map_mass_range:
        lo, hi = args.map_mass_range
        MAP_MASS_RANGE = None if hi <= lo else (lo, hi)
    if args.palette:
        if not hasattr(ROOT, args.palette):
            raise SystemExit(f"ROOT has no palette named {args.palette}")
        MAP_PALETTE = getattr(ROOT, args.palette)
    SF.setup_style()
    year = normalize_year(args.year)

    if args.file and len(args.sig) > 1:
        raise SystemExit("--file cannot be combined with several --sig: the "
                         "files cannot be attributed to a production mode")

    # which (mass, relerr) pairs to run
    if args.mass or args.relerr:
        if not (args.mass and args.relerr):
            raise SystemExit("--mass and --relerr must be given together: a "
                             "mass definition and its own error branch")
        variants = [("custom", args.mass, args.relerr, "custom")]
    elif args.variant == "both":
        variants = [(k, *VARIANTS[k]) for k in VARIANT_ORDER]
    else:
        variants = [(args.variant, *VARIANTS[args.variant])]

    absolute = args.absolute
    ytit = "#sigma_{m} (GeV)" if absolute else "#sigma_{m} / m"
    ymax = ABS_MAX if absolute else REL_MAX

    outdir = os.path.join(args.plotdir, args.category, year)
    os.makedirs(outdir, exist_ok=True)
    label = f"{args.category}, {year}"

    print(f"category : {args.category}"
          f"\nyear     : {year}"
          f"\nsamples  : {', '.join(args.sig)}"
          f"\nvariants : "
          + "; ".join(f"{k} ({m} / {e})" for k, m, e, _ in variants)
          + f"\ny axis   : {'absolute (GeV)' if absolute else 'relative'}"
            f"\ncut      : {args.cut or 'none'}"
            f"\nplots    : {outdir}")

    # ---- book every (sample, variant) before triggering any loop -------
    entries = []
    for sig in args.sig:
        files = (list(args.file) if args.file else
                 get_signal_files(sig, args.category, year,
                                  args.rootfiles_dir))
        if not files:
            print(f"  !! no files for {sig}, skipping")
            continue
        print(f"\n--- {sig}: {len(files)} file(s) ---")
        for f in files:
            print(f"      {f}")

        df0 = ROOT.RDataFrame(TREE, files)
        cols = set(str(c) for c in df0.GetColumnNames())
        weight = None if args.weight == "none" else args.weight
        if weight and weight not in cols:
            raise SystemExit(f"{sig}: input has no weight column {weight}")

        for vkey, mcol, ecol, vdesc in variants:
            missing = [c for c in (mcol, ecol) if c not in cols]
            if missing:
                print(f"  !! {sig}: no {', '.join(missing)}, "
                      f"skipping variant {vkey}")
                continue

            df = (df0
                  .Filter(f"!std::isnan({mcol})", "valid mass")
                  .Filter(f"{ecol} > 0", "valid resolution")
                  .Filter(f"{mcol} >= {XLOW} && {mcol} < {XHIGH}",
                          "fit window"))
            if args.cut:
                df = df.Filter(args.cut, "user cut")

            yexpr = (f"(float)({mcol} * {ecol})" if absolute
                     else f"(float){ecol}")
            df = df.Define("_yval", yexpr)

            sp2 = (f"h2_{sig}_{vkey}", "", MASS_NBINS, XLOW, XHIGH,
                   REL_NBINS, 0.0, ymax)
            sp1 = (f"h1_{sig}_{vkey}", "", REL_NBINS, 0.0, ymax)
            h2 = (df.Histo2D(sp2, mcol, "_yval", weight) if weight
                  else df.Histo2D(sp2, mcol, "_yval"))
            h1 = (df.Histo1D(sp1, "_yval", weight) if weight
                  else df.Histo1D(sp1, "_yval"))
            entries.append({"sig": sig, "variant": vkey, "desc": vdesc,
                            "mass_col": mcol, "relerr_col": ecol,
                            "files": files,
                            "_h2": h2, "_h1": h1, "_n": df.Count()})

    if not entries:
        raise SystemExit("nothing to plot")

    # ---- trigger, then profile ----------------------------------------
    probs = array("d", [Q_LO, Q_MID, Q_HI])
    for e in entries:
        e["h2"] = e.pop("_h2").GetValue()
        e["h2"].SetDirectory(0)
        e["h1"] = e.pop("_h1").GetValue()
        e["h1"].SetDirectory(0)
        e["n"] = e.pop("_n").GetValue()
        e["prof"] = quantile_profile(e["h2"])
        qs = array("d", [0.0, 0.0, 0.0])
        e["h1"].GetQuantiles(3, qs, probs)
        e["q16"], e["median"], e["q84"] = qs[0], qs[1], qs[2]
        xs, _, mids, _, _ = e["prof"]
        e["slope_per_gev"] = ((mids[-1] - mids[0]) / (xs[-1] - xs[0])
                              if len(xs) > 1 else float("nan"))

    # ---- report --------------------------------------------------------
    print(f"\n{'sample':<7}{'variant':<10}{'events':>11}{'q16':>11}"
          f"{'median':>11}{'q84':>11}{'slope/GeV':>12}")
    for e in entries:
        print(f"{e['sig']:<7}{e['variant']:<10}{e['n']:>11d}"
              f"{e['q16']:>11.5g}{e['median']:>11.5g}{e['q84']:>11.5g}"
              f"{e['slope_per_gev']:>12.3g}")

    # FSR gain per sample, from the integrated medians
    if len(variants) > 1:
        print(f"\n{'sample':<7}{'pre':>11}{'post':>11}{'change':>10}")
        for sig in args.sig:
            pre = next((e for e in entries if e["sig"] == sig
                        and e["variant"] == "preFSR"), None)
            post = next((e for e in entries if e["sig"] == sig
                         and e["variant"] == "postFSR"), None)
            if not (pre and post and pre["median"]):
                continue
            rel = (post["median"] - pre["median"]) / pre["median"]
            print(f"{sig:<7}{pre['median']:>11.5g}{post['median']:>11.5g}"
                  f"{rel:>+9.1%}")
        print("a NEGATIVE change means FSR recovery improved the resolution")

    # ---- plots ---------------------------------------------------------
    for e in entries:
        plot_2d(e, outdir, year, label, ytit, args.pdf)
    n_samples = len({e["sig"] for e in entries})
    for vkey, _, _, _ in variants:
        sub = [e for e in entries if e["variant"] == vkey]
        sfx = f"_{vkey}" if len(variants) > 1 else ""
        # The standalone profile only adds something when several samples
        # share the axes: with one sample its median and band are already
        # drawn on that sample's 2D map, and the pre/post comparison is
        # covered by fsr_compare.
        if n_samples > 1 or args.force_profile:
            plot_profile(sub, outdir, year, f"{label}, {vkey}", ytit,
                         args.pdf, sfx)
        plot_1d(sub, outdir, year, f"{label}, {vkey}", ytit, args.pdf, sfx)
    if len(variants) > 1:
        plot_fsr_compare(entries, outdir, year, label, ytit, args.pdf)

    # ---- json ----------------------------------------------------------
    dump = {}
    for e in entries:
        xs, los, mids, his, ns = e["prof"]
        dump.setdefault(e["sig"], {})[e["variant"]] = {
            "mass_column": e["mass_col"], "relerr_column": e["relerr_col"],
            "n_events": e["n"], "q16": e["q16"], "median": e["median"],
            "q84": e["q84"], "slope_per_gev": e["slope_per_gev"],
            "profile": {"mass": xs, "q16": los, "median": mids,
                        "q84": his, "n": ns}}

    out = {"provenance": {
        "written_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "git_commit": git_commit(), "year": year,
        "eras": None if args.file else expand_year(year, SIGNAL_YEAR_GROUPS),
        "category": args.category, "samples": args.sig,
        "variants": {k: {"mass": m, "relerr": er} for k, m, er, _ in variants},
        "absolute": absolute, "cut": args.cut, "weight": args.weight,
        "map_mass_range": MAP_MASS_RANGE, "map_logz": MAP_LOGZ,
        "mass_window": [XLOW, XHIGH],
        "quantiles": [Q_LO, Q_MID, Q_HI]}, "samples": dump}
    jf = os.path.join(outdir, f"massResVsMass_{year}.json")
    with open(jf, "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
    print(f"\nwrote {jf}")
    print(f"plots {outdir}")


if __name__ == "__main__":
    main()
