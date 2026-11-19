"""
bkgTable.py -- the Run 2 note's Table 6.13 from bkgFit.py's output; bkgFit.py
redraws it after every bin it fits.

One row per BDT bin, one column per family: the shape degrees of freedom
of the family's F-test pick (its free parameters, the normalisation aside)
and the chi2/ndf of its background-only fit to the data sideband, as
bkgFit.py computed it (110-120 and 130-150 GeV, low-count bins merged).
A family with no order that passes the goodness of fit is "fails GOF",
with the best p-value its orders reached; one not fitted in a bin is "-".
The default member (pdfindex 0, lowest penalised NLL) has a blue cell.

The rows are the bins of the BDT split in force (bdtBins.SELMVA) with a
Bkg_<tag>_params.json in the wsdir made with today's cut, so the table is
complete once the last bin is fitted. Written to
<plotdir>/<cat>/bkg_table_<cat>_<year>.png.
"""

import json
import os

import ROOT

from . import bdtBins, style
from . import pdfDefinitions as pdefs

ROOT.gROOT.SetBatch(True)

# pixel sizes of the drawn table
ROW_PX, HEAD_PX, LABEL_PX, CELL_PX = 30, 34, 120, 196
TOP_PX, BOTTOM_PX, SIDE_PX = 40, 58, 16
TEXT_PX, SMALL_PX = 15, 13
GRID_RGB = (190, 190, 190)
BEST_RGB = (232, 240, 252)          # light blue behind the default member
FAIL_RGB = (120, 120, 120)          # grey text: no order passes the GOF


def bin_number(b):
    digits = b[len(b.rstrip("0123456789")):]
    return int(digits) if digits else -1


def read_bins(wsdir, category, year):
    """{bin: bkgFit.py record} for the BDT bins of the split in force
    (bdtBins.SELMVA) that have a fit in wsdir made with today's cut: fits
    of another split left in the same wsdir are not shown."""
    recs = {}
    cuts = bdtBins.SELMVA.get(category, {})
    for b in sorted((b for b in cuts if b not in bdtBins.INCLUSIVE),
                    key=bin_number):
        path = os.path.join(wsdir, f"Bkg_{category}_{b}_{year}_params.json")
        if not os.path.isfile(path):
            continue
        with open(path) as fh:
            rec = json.load(fh)
        sel = rec["provenance"].get("selection")
        if sel is not None and sel != cuts[b]:
            continue
        recs[b] = rec
    return recs


def families(recs):
    """The families any bin fitted, DEFAULT_FAMILIES first"""
    fam_of = {m: f for f, ms in pdefs.BKG_FAMILIES.items() for m in ms}
    seen = {fam_of[n] for r in recs.values() for n in
            list(r["candidates"]) + list(r.get("failed_fits", {}))
            if n in fam_of}
    order = pdefs.DEFAULT_FAMILIES + [f for f in pdefs.BKG_FAMILIES
                                      if f not in pdefs.DEFAULT_FAMILIES]
    return [f for f in order if f in seen]


def cell(rec, fam):
    """(kind, text, plain): kind ok / best / fail / none"""
    pick = rec["family_choice"].get(fam)
    if pick:
        c = rec["candidates"][pick]
        kind = "best" if rec.get("best") == pick else "ok"
        return (kind, "%d-dof, #chi^{2}/ndf = %.2f" % (c["npar"], c["chi2_ndf"]),
                "%d-dof, chi2/ndf=%.2f" % (c["npar"], c["chi2_ndf"]))
    tried = {n: c for n, c in rec["candidates"].items()
             if n in pdefs.BKG_FAMILIES[fam]}
    if tried:
        p = max(c["pvalue"] for c in tried.values())
        return ("fail", "fails GOF (best p = %.0e)" % p if p > 0 else
                "fails GOF (p = 0)", "fails GOF (best p=%.1e)" % p)
    return ("none", "-", "-")


def header(fam):
    return pdefs.FAMILY_LABEL.get(fam, fam).format(n="")


def draw(recs, fams, category, year, outbase, pdf):
    bins = list(recs)
    w = 2 * SIDE_PX + LABEL_PX + CELL_PX * len(fams)
    h = TOP_PX + HEAD_PX + ROW_PX * len(bins) + BOTTOM_PX
    style.setup_style()
    c = ROOT.TCanvas("c_table", "", w, h)
    c.SetMargin(0, 0, 0, 0)
    keep = []

    def X(px):
        return px / float(w)

    def Y(px):                       # pixels from the top
        return 1 - px / float(h)

    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(43)
    keep.append(lat)
    # title band: the category and data, what the table is
    lat.SetTextSize(14)
    lat.SetTextAlign(33)
    lat.DrawLatex(X(w - SIDE_PX), Y(12), "%s %s (%.1f TeV)" % (
        category, year, style.SQRT_S_TEV))
    lat.SetTextAlign(13)
    lat.SetTextSize(TEXT_PX)
    lat.DrawLatex(X(SIDE_PX), Y(12), "Shape degrees of freedom selected by "
                  "the F-test, and #chi^{2}/ndf of the background-only fit "
                  "to the data sideband")

    top = TOP_PX
    x0 = SIDE_PX
    # header row
    grid = ROOT.TColor.GetColor(*GRID_RGB)
    lat.SetTextAlign(22)
    lat.SetTextSize(TEXT_PX)
    lat.DrawLatex(X(x0 + LABEL_PX / 2.0), Y(top + HEAD_PX / 2.0),
                  "#bf{Category}")
    for k, fam in enumerate(fams):
        lat.DrawLatex(X(x0 + LABEL_PX + CELL_PX * (k + 0.5)),
                      Y(top + HEAD_PX / 2.0), f"#bf{{{header(fam)}}}")
    # rows
    best_col = ROOT.TColor.GetColor(*BEST_RGB)
    fail_col = ROOT.TColor.GetColor(*FAIL_RGB)
    for i, b in enumerate(bins):
        ytop = top + HEAD_PX + ROW_PX * i
        yc = Y(ytop + ROW_PX / 2.0)
        lat.SetTextColor(ROOT.kBlack)
        lat.SetTextSize(TEXT_PX)
        lat.DrawLatex(X(x0 + LABEL_PX / 2.0), yc, style.BIN_LABEL.get(b, b))
        for k, fam in enumerate(fams):
            kind, text, _ = cell(recs[b], fam)
            xl = x0 + LABEL_PX + CELL_PX * k
            if kind == "best":
                box = ROOT.TBox(X(xl + 1), Y(ytop + ROW_PX - 1),
                                X(xl + CELL_PX - 1), Y(ytop + 1))
                box.SetFillColor(best_col)
                box.Draw()
                keep.append(box)
            lat.SetTextColor(fail_col if kind in ("fail", "none")
                             else ROOT.kBlack)
            lat.SetTextSize(SMALL_PX if kind == "fail" else TEXT_PX)
            lat.DrawLatex(X(xl + CELL_PX / 2.0), yc, text)
    lat.SetTextColor(ROOT.kBlack)
    # lines: under the header, between the label column and the cells,
    # and the table's top and bottom
    ybot = top + HEAD_PX + ROW_PX * len(bins)
    xr = x0 + LABEL_PX + CELL_PX * len(fams)
    for y_px, width in ((top, 2), (top + HEAD_PX, 1), (ybot, 2)):
        ln = ROOT.TLine(X(x0), Y(y_px), X(xr), Y(y_px))
        ln.SetNDC()
        ln.SetLineWidth(width)
        ln.Draw()
        keep.append(ln)
    ln = ROOT.TLine(X(x0 + LABEL_PX), Y(top), X(x0 + LABEL_PX), Y(ybot))
    ln.SetNDC()
    ln.SetLineColor(grid)
    ln.Draw()
    keep.append(ln)
    for i in range(1, len(bins)):
        y_px = top + HEAD_PX + ROW_PX * i
        ln = ROOT.TLine(X(x0), Y(y_px), X(xr), Y(y_px))
        ln.SetNDC()
        ln.SetLineColor(grid)
        ln.SetLineStyle(3)
        ln.Draw()
        keep.append(ln)
    # footnote
    rec = next(iter(recs.values()))
    prov = rec["provenance"]
    lat.SetTextAlign(13)
    lat.SetTextSize(SMALL_PX)
    lat.SetTextColor(ROOT.kGray + 2)
    lat.DrawLatex(X(SIDE_PX), Y(ybot + 10),
                  "dof: free shape parameters of the pick (normalisation "
                  "aside); #chi^{2} over the sideband %g-%g and %g-%g GeV, "
                  "low-count bins merged." % (
                      prov["window"][0], prov["blind_region"][0],
                      prov["blind_region"][1], prov["window"][1]))
    lat.DrawLatex(X(SIDE_PX), Y(ybot + 30),
                  "F-test at p < %g over the orders with GOF p > %g; "
                  "blue: the default member (lowest penalised NLL)." % (
                      prov["f_threshold"], prov["gof_threshold"]))
    keep.append(c)
    c.SaveAs(f"{outbase}.png")
    if pdf:
        c.SaveAs(f"{outbase}.pdf")


def write_table(wsdir, category, year, plotdir, pdf=False):
    """Print and draw the table of every fitted bin; returns the PNG path,
    or None if no BDT bin has a fit yet."""
    recs = read_bins(wsdir, category, year)
    if not recs:
        return None
    fams = families(recs)
    wid = 29
    print("\n" + "category".ljust(12) + "".join(
        header(f).replace("#times", "x").ljust(wid) for f in fams))
    for b, rec in recs.items():
        print(b.ljust(12) + "".join(
            (cell(rec, f)[2] + (" *" if cell(rec, f)[0] == "best" else ""))
            .ljust(wid) for f in fams))
    print("* the default member (pdfindex 0)")
    os.makedirs(os.path.join(plotdir, category), exist_ok=True)
    outbase = os.path.join(plotdir, category, f"bkg_table_{category}_{year}")
    draw(recs, fams, category, year, outbase, pdf)
    return f"{outbase}.png"
