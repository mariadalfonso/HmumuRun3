#!/usr/bin/env python3
"""
bkgMassInBDT.py -- background dimuon mass in slices of the BDT discriminant,
with the control that tests whether the muon pairing sculpts the mass.

    python bkgMassInBDT.py 2024 -d /work/submit/kbai/HmumuRun3/bkg
    python bkgMassInBDT.py 2024 -d ... --mass-range 70 200 --cat VLcat
    python bkgMassInBDT.py 2024 -d ... -s 204,205 -o ~/public_html/HmumuRun3/bkgshape

Samples are summed with w_allSF, so the yields are expected events.

Why the n==2 / n>=3 split
-------------------------
Slicing the background mass in BDT bins shows whether the BDT shapes the mass,
which it does by construction -- the two are correlated. It cannot by itself
show whether the PAIRING shapes the mass, because both effects act at once.

The control is the good-muon multiplicity. With exactly two good muons the
category vetoes leave no pairing freedom at all: the pair is both muons, and
no criterion can change that. With three or more the algorithm chooses. Both
populations sit in the same category and the same BDT bin and differ only in
whether a choice was made, so any excess near m_H in the n>=3 shape relative
to n==2 is the sculpting, directly.

The script reports, per BDT bin and per multiplicity:

    yield, entries, and the fraction of the mass range inside 120-130 GeV

If the 120-130 fraction agrees between n==2 and n>=3 within its uncertainty,
the pairing is not building structure under the peak. If n>=3 is high, it is.

Statistics are the limiting factor: category 3 needs three tight muons above
20 GeV, so a VV sample may leave only a handful of events in the fit range.
The yields are printed first for that reason, and --no-cat drops the category
requirement (keeping nGoodMuons) to recover events when the test would
otherwise be unanswerable.
"""

import argparse
import getpass
import os
import sys

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

TREE = "events"
DEFAULT_SAMPLES = [201, 202, 203, 204, 205, 206, 213, 214, 215, 216]
SAMPLE_LABEL = {201: "WZto2L2Q", 202: "ZZto2L2Nu", 203: "ZZto2L2Q",
                204: "WZto3LNu", 205: "ZZto4L", 206: "WWto2L2Nu",
                213: "WZZ", 214: "ZZZ", 215: "WWW", 216: "WWZ"}
PEAK = (120.0, 130.0)

LAB, TIT, HDR = 20, 23, 26


def find_mva(cols):
    """The discriminant branch, whose name depends on the snapshot vintage."""
    for c in ("discrMVA", "discrMVA0"):
        if c in cols:
            return c
    return None


def collect(paths, args):
    """Sum the histograms over samples, per BDT bin and multiplicity class.

    Returns hists[(bin index, class)] -> TH1D and a per-cell yield table.
    Everything is booked before the first GetValue, so each file is read once.
    """
    lo, hi = args.mass_range
    edges = args.bdt_bins
    classes = [("n=2", "nGoodMuons == 2"), ("n>=3", "nGoodMuons >= 3")]

    hists, booked, used, persample = {}, [], [], []
    for k in range(len(edges) - 1):
        for cname, _ in classes:
            h = ROOT.TH1D(f"h_{k}_{cname}", "", args.nbins, lo, hi)
            h.Sumw2()
            hists[(k, cname)] = h
    # BDT-inclusive: the pairing acting on its own, with no BDT selection to
    # confound it. Summed over samples, and kept per sample as well.
    hincl, hsamp, booked_x = {}, {}, []
    # exact multiplicities per sample, for the dominant-sample plot: the two
    # classes above lump n=3 with n=4, and in ZZ->4l those are different
    # topologies rather than different amounts of the same one
    MULT = [("n=2", "nGoodMuons == 2"), ("n=3", "nGoodMuons == 3"),
            ("n=4", "nGoodMuons == 4"), ("n>=5", "nGoodMuons >= 5")]
    hmult, booked_m = {}, []
    for cname, _ in classes:
        h = ROOT.TH1D(f"hincl_{cname}", "", args.nbins, lo, hi)
        h.Sumw2()
        hincl[cname] = h

    for sid, path in paths:
        df = ROOT.RDataFrame(TREE, path)
        cols = {str(c) for c in df.GetColumnNames()}
        mva = args.mva or find_mva(cols)
        if mva is None:
            print(f"  {SAMPLE_LABEL.get(sid, sid)}: no discrMVA branch, skipped")
            continue
        for need in ("HiggsCandCorrMass", "w_allSF", "nGoodMuons"):
            if need not in cols:
                print(f"  {SAMPLE_LABEL.get(sid, sid)}: no {need}, skipped")
                break
        else:
            sel = df
            if args.cat_filter and "category" in cols:
                sel = sel.Filter(args.cat_filter, "category")
            # per-sample totals inside the mass range, so the contribution of
            # each sample to the sum is visible rather than assumed
            inr = sel.Filter(f"HiggsCandCorrMass > {lo} && "
                             f"HiggsCandCorrMass < {hi}", "mass range")
            persample.append((
                sid, df.Count(), inr.Count(), inr.Sum("w_allSF"),
                inr.Filter("nGoodMuons == 2").Count(),
                inr.Filter("nGoodMuons >= 3").Count(),
                inr.Filter("nGoodMuons >= 3").Sum("w_allSF")))
            for cname, cut in classes:
                hh = inr.Filter(cut).Histo1D(
                    (f"x_{sid}_{cname}", "", args.nbins, lo, hi),
                    "HiggsCandCorrMass", "w_allSF")
                booked_x.append((sid, cname, hh))
            for mname, cut in MULT:
                hh = inr.Filter(cut).Histo1D(
                    (f"m_{sid}_{mname}", "", args.nbins, lo, hi),
                    "HiggsCandCorrMass", "w_allSF")
                booked_m.append((sid, mname, hh))
            for k in range(len(edges) - 1):
                a, b = edges[k], edges[k + 1]
                inbin = sel.Filter(f"{mva} >= {a} && {mva} < {b}", f"bdt{k}")
                for cname, cut in classes:
                    hh = inbin.Filter(cut).Histo1D(
                        (f"t_{sid}_{k}_{cname}", "", args.nbins, lo, hi),
                        "HiggsCandCorrMass", "w_allSF")
                    booked.append((k, cname, hh))
            used.append((sid, mva))

    if not booked:
        return None, None, None, None, None, None, None

    for k, cname, hh in booked:            # triggers one loop per file
        hists[(k, cname)].Add(hh.GetPtr())
    for sid, cname, hh in booked_x:
        h = hh.GetPtr().Clone(f"s_{sid}_{cname}")
        hsamp[(sid, cname)] = h
        hincl[cname].Add(h)
    for sid, mname, hh in booked_m:
        hmult[(sid, mname)] = hh.GetPtr().Clone(f"mm_{sid}_{mname}")

    rows = [dict(sid=sid, n_tot=a.GetValue(), n=b.GetValue(),
                 y=c.GetValue(), n2=d.GetValue(), n3=e.GetValue(),
                 y3=f.GetValue())
            for sid, a, b, c, d, e, f in persample]

    table = {}
    for (k, cname), h in hists.items():
        n = h.GetEntries()
        y = h.Integral()
        b1 = h.FindFixBin(PEAK[0] + 1e-6)
        b2 = h.FindFixBin(PEAK[1] - 1e-6)
        ypk = h.Integral(b1, b2)
        err = ROOT.Double_t(0) if hasattr(ROOT, "Double_t") else 0.0
        table[(k, cname)] = dict(entries=n, yield_=y, peak=ypk,
                                 frac=(ypk / y if y else 0.0))
    return hists, table, used, rows, hincl, hsamp, hmult


def report_samples(rows, args):
    lo, hi = args.mass_range
    print(f"\nper-sample contribution to the sum   "
          f"(mass {lo:.0f}-{hi:.0f} GeV)")
    print(f"  {'sample':<12}{'snapshot':>12}{'in range':>11}{'yield':>12}"
          f"{'n=2':>9}{'n>=3':>9}{'yield n>=3':>12}")
    print("  " + "-" * 75)
    tot = dict(n_tot=0, n=0, y=0.0, n2=0, n3=0, y3=0.0)
    for r in sorted(rows, key=lambda r: -r["n3"]):
        print(f"  {SAMPLE_LABEL.get(r['sid'], r['sid']):<12}"
              f"{r['n_tot']:>12,}{r['n']:>11,}{r['y']:>12.4g}"
              f"{r['n2']:>9,}{r['n3']:>9,}{r['y3']:>12.4g}")
        for k in tot:
            tot[k] += r[k]
    print("  " + "-" * 75)
    print(f"  {'SUM':<12}{tot['n_tot']:>12,}{tot['n']:>11,}{tot['y']:>12.4g}"
          f"{tot['n2']:>9,}{tot['n3']:>9,}{tot['y3']:>12.4g}")
    print("  all samples above are added into one histogram, weighted by "
          "w_allSF")
    if tot["n3"] == 0:
        print("  WARNING: no event has three or more good muons, so the "
              "pairing control\n           cannot be formed at all")
    return tot


def report(table, args):
    edges = args.bdt_bins
    print(f"\n{'BDT bin':<14}{'class':>7}{'entries':>10}{'yield':>12}"
          f"{'in 120-130':>12}{'fraction':>11}")
    print("-" * 66)
    for k in range(len(edges) - 1):
        for cname in ("n=2", "n>=3"):
            t = table[(k, cname)]
            if t["entries"] == 0:
                continue
            # binomial error on the fraction, from the raw entry count
            f, n = t["frac"], max(t["entries"], 1)
            ef = np.sqrt(max(f * (1 - f), 0) / n)
            print(f"{edges[k]:.1f}-{edges[k+1]:<9.1f}{cname:>7}"
                  f"{t['entries']:>10.0f}{t['yield_']:>12.4g}"
                  f"{t['peak']:>12.4g}{100*f:>9.1f}%"
                  f" +-{100*ef:.1f}")
    print(f"\nfraction = yield in {PEAK[0]:.0f}-{PEAK[1]:.0f} over the whole "
          f"{args.mass_range[0]:.0f}-{args.mass_range[1]:.0f} range")
    print("the n=2 rows are the control: no pairing freedom, so their fraction")
    print("is what the mass shape looks like with the pairing switched off")

    print(f"\n{'BDT bin':<14}{'n=2 frac':>12}{'n>=3 frac':>12}{'difference':>14}")
    print("-" * 52)
    for k in range(len(edges) - 1):
        a, b = table[(k, "n=2")], table[(k, "n>=3")]
        if a["entries"] == 0 or b["entries"] == 0:
            print(f"{edges[k]:.1f}-{edges[k+1]:<9.1f}"
                  f"{'(too few events)':>38}")
            continue
        fa, fb = a["frac"], b["frac"]
        ea = np.sqrt(max(fa*(1-fa), 0)/max(a["entries"], 1))
        eb = np.sqrt(max(fb*(1-fb), 0)/max(b["entries"], 1))
        d, ed = fb - fa, np.sqrt(ea*ea + eb*eb)
        sig = abs(d) / ed if ed else 0.0
        print(f"{edges[k]:.1f}-{edges[k+1]:<9.1f}{100*fa:>11.1f}%"
              f"{100*fb:>11.1f}%{100*d:>+9.1f} +-{100*ed:.1f}"
              f"  ({sig:.1f} sigma)")
    print("\na difference compatible with zero means the pairing is not "
          "building\nstructure under the peak in this sample")


def peak_frac(h):
    """(fraction of the range inside the peak window, binomial error)."""
    y = h.Integral()
    if y <= 0:
        return 0.0, 0.0
    b1 = h.FindFixBin(PEAK[0] + 1e-6)
    b2 = h.FindFixBin(PEAK[1] - 1e-6)
    f = h.Integral(b1, b2) / y
    n = max(h.GetEntries(), 1)
    return f, np.sqrt(max(f * (1 - f), 0) / n)


def overlay(entries, xtitle, header, stem, args, legend_right=True,
            normalize=True, logy=False):
    """entries: list of (TH1D, legend label).

    normalize=True compares shapes, each scaled to unit area. False keeps the
    w_allSF yields, so samples of very different size can be read against each
    other on a log scale -- which is the only way the per-sample plots are
    legible, since their statistics differ by orders of magnitude.
    """
    lo, hi = args.mass_range
    live = [(h, lab) for h, lab in entries if h.Integral() > 0]
    if not live:
        print(f"    {os.path.basename(stem)}: nothing to draw")
        return None
    c = ROOT.TCanvas("c_" + os.path.basename(stem), "", 880, 650)
    c.SetLeftMargin(0.14)
    c.SetRightMargin(0.04)
    c.SetBottomMargin(0.14)
    c.SetTopMargin(0.09)
    cols = [ROOT.kAzure + 2, ROOT.kRed + 1, ROOT.kGreen + 3, ROOT.kOrange + 7,
            ROOT.kViolet + 6, ROOT.kTeal + 3, ROOT.kMagenta + 2,
            ROOT.kSpring - 1, ROOT.kGray + 2, ROOT.kPink + 7]
    x0 = 0.60 if legend_right else 0.17
    leg = ROOT.TLegend(x0, 0.88 - 0.055 * len(live), x0 + 0.34, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(43)
    leg.SetTextSize(LAB)
    hs = []
    for k, (h0, lab) in enumerate(live):
        h = h0.Clone(h0.GetName() + ("_nrm" if normalize else "_abs"))
        if normalize:
            h.Scale(1.0 / h.Integral())
        h.SetLineColor(cols[k % len(cols)])
        h.SetLineWidth(2)
        h.SetMarkerColor(cols[k % len(cols)])
        h.SetMarkerStyle(20)
        h.SetMarkerSize(0.7)
        hs.append(h)
        f, ef = peak_frac(h0)
        leg.AddEntry(h, f"{lab}  ({h0.GetEntries():.0f} ev, "
                        f"{100*f:.0f}#pm{100*ef:.0f}% in peak)", "lep")
    ymax = max(h.GetMaximum() + h.GetBinError(h.GetMaximumBin()) for h in hs)
    if logy:
        c.SetLogy(1)
        # a log axis needs a positive floor; take half the smallest positive
        # bin, so the thinnest sample is still on the canvas. Bins that are
        # zero or negative (amcatnlo weights can go negative) fall off the
        # bottom rather than breaking the axis.
        pos = [h.GetBinContent(b) for h in hs
               for b in range(1, h.GetNbinsX() + 1)
               if h.GetBinContent(b) > 0]
        ylo = 0.5 * min(pos) if pos else 1e-6
        ytop = 20.0 * ymax
    else:
        ylo, ytop = 0.0, 1.45 * ymax
    for k, h in enumerate(hs):
        h.SetMaximum(ytop)
        h.SetMinimum(ylo)
        for ax, t in ((h.GetXaxis(), xtitle),
                      (h.GetYaxis(), "Normalized Events" if normalize
                       else "Events (w_{allSF})")):
            ax.SetTitle(t)
            ax.SetTitleFont(43)
            ax.SetTitleSize(TIT)
            ax.SetLabelFont(43)
            ax.SetLabelSize(LAB)
        h.GetXaxis().SetTitleOffset(1.4)
        h.GetYaxis().SetTitleOffset(1.6)
        h.Draw("HIST E" if k == 0 else "HIST E SAME")
    leg.Draw()
    lines = []
    for x in PEAK:
        if lo < x < hi:
            ln = ROOT.TLine(x, ylo, x, ytop)
            ln.SetLineColor(ROOT.kGray + 2)
            ln.SetLineStyle(2)
            ln.Draw()
            lines.append(ln)
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(43)
    lat.SetTextSize(HDR)
    lat.SetTextAlign(12)
    lat.DrawLatex(0.14, 0.955, header)
    c.SaveAs(f"{stem}.png")
    print(f"    wrote {stem}.png")
    return c, hs, leg, lat, lines


def draw_dominant(hmult, args, outdir):
    """The two samples that carry the multi-muon events, by exact multiplicity.

    One pad per sample, so colour no longer has to encode the sample and can
    encode the multiplicity instead. The three classes are drawn as markers
    with error bars, each with its own colour and marker shape, and are nudged
    apart horizontally within each bin so their error bars do not sit on top
    of one another. No class is treated as the reference: the multiplicities
    are different topologies and their shapes differ for reasons other than
    the pairing.

    n=3 and n=4 are kept apart rather than lumped into n>=3: in ZZ->4l they
    are different topologies, the four-lepton final state against one that
    lost a lepton, and in WZ->3lnu n=4 should barely exist at all.
    """
    os.makedirs(outdir, exist_ok=True)
    lo, hi = args.mass_range
    CLASSES = ("n=2", "n=3", "n=4")
    # (colour, marker, marker size, horizontal nudge as a fraction of a bin)
    STYLE = {"n=2": (ROOT.kBlack, 24, 1.1, -0.22),
             "n=3": (ROOT.kRed + 1, 20, 1.0, 0.0),
             "n=4": (ROOT.kAzure + 2, 25, 1.2, +0.22)}
    # (sample, class) pairs left off: WZ->3lnu barely populates n=4
    DOM_SKIP = {(204, "n=4")}
    samples = [sid for sid in (204, 205)
               if any((sid, m) in hmult for m in CLASSES)]
    if not samples:
        print("    neither WZto3LNu nor ZZto4L is present, "
              "so the dominant-sample plot is skipped")
        return None

    rebin = args.dom_rebin
    if rebin > 1 and args.nbins % rebin:
        print(f"    --dom-rebin {rebin} does not divide --nbins {args.nbins}, "
              "using the native binning")
        rebin = 1

    def to_graph(h, nudge):
        """TH1 -> TGraphErrors with points shifted by nudge * bin width."""
        g = ROOT.TGraphErrors()
        for b in range(1, h.GetNbinsX() + 1):
            if h.GetBinContent(b) <= 0:
                continue
            w = h.GetBinWidth(b)
            n = g.GetN()
            g.SetPoint(n, h.GetBinCenter(b) + nudge * w, h.GetBinContent(b))
            g.SetPointError(n, 0.0, h.GetBinError(b))
        return g

    W, H = 760, 640
    c = ROOT.TCanvas("c_dom", "", W * len(samples), H)
    keep, skipped, drawn = [c], [], 0

    for i, sid in enumerate(samples):
        label = SAMPLE_LABEL.get(sid, str(sid))
        hs = {}
        for m in CLASSES:
            if (sid, m) in DOM_SKIP:
                continue
            h0 = hmult.get((sid, m))
            if h0 is None or h0.Integral() <= 0:
                continue
            if h0.GetEntries() < args.min_entries:
                skipped.append(f"{label} {m} ({h0.GetEntries():.0f})")
                continue
            h = h0.Clone(f"dom_{sid}_{m}")
            if rebin > 1:
                h.Rebin(rebin)
            h.Scale(1.0 / h.Integral())
            hs[m] = (h, h0)
        if not hs:
            print(f"    {label}: no usable multiplicity class, pad skipped")
            continue

        c.cd()
        pad = ROOT.TPad(f"pad_{sid}", "", i / len(samples), 0.0,
                        (i + 1) / len(samples), 1.0)
        pad.SetLeftMargin(0.16)
        pad.SetRightMargin(0.04)
        pad.SetBottomMargin(0.13)
        pad.SetTopMargin(0.08)
        pad.Draw()
        pad.cd()
        keep.append(pad)

        ymax = max(h.GetMaximum() + h.GetBinError(h.GetMaximumBin())
                   for h, _ in hs.values())
        ytop = 1.6 * ymax
        frame = next(iter(hs.values()))[0].Clone(f"domframe_{sid}")
        frame.Reset()
        frame.SetMaximum(ytop)
        frame.SetMinimum(0.0)
        for ax, t in ((frame.GetXaxis(), "m_{#mu#mu} [GeV]"),
                      (frame.GetYaxis(), "Normalized Events")):
            ax.SetTitle(t)
            ax.SetTitleFont(43)
            ax.SetTitleSize(TIT)
            ax.SetLabelFont(43)
            ax.SetLabelSize(LAB)
        frame.GetXaxis().SetTitleOffset(1.4)
        frame.GetYaxis().SetTitleOffset(1.7)
        frame.GetYaxis().SetNdivisions(505)
        frame.Draw("AXIS")

        box = ROOT.TBox(max(PEAK[0], lo), 0.0, min(PEAK[1], hi), ytop)
        box.SetFillColorAlpha(ROOT.kGray, 0.6)
        box.SetFillStyle(1001)
        box.SetLineWidth(0)
        box.Draw()
        keep += [frame, box]

        leg = ROOT.TLegend(0.58, 0.66, 0.96, 0.90)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(43)
        leg.SetTextSize(LAB - 2)
        leg.AddEntry(box, f"Higgs SR, {PEAK[0]:.0f}-{PEAK[1]:.0f} GeV", "f")
        for m in CLASSES:
            if m not in hs:
                continue
            h, h0 = hs[m]
            col, mk, mks, nudge = STYLE[m]
            if args.dom_markers:
                obj = to_graph(h, nudge)
                obj.SetMarkerStyle(mk)
                obj.SetMarkerSize(mks)
                opt, lopt = "P Z SAME", "pe"
            else:
                # step lines: the eye follows the shape of each class, so a
                # bump or shoulder under the SR stands out against the trend
                obj = h
                obj.SetMarkerSize(0)
                opt, lopt = "HIST E SAME", "le"
            obj.SetLineColor(col)
            obj.SetLineWidth(2)
            obj.SetMarkerColor(col)
            obj.Draw(opt)
            # h0 is unnormalised and filled with w_allSF, so its integral is
            # the expected yield; GetEntries is the raw MC count, which is
            # what sets the statistical precision, so both are shown
            leg.AddEntry(obj, f"{m}  ({h0.Integral():.3g} exp, "
                              f"{h0.GetEntries():.0f} MC)", lopt)
            keep.append(obj)
        leg.Draw()
        pad.RedrawAxis()

        lat = ROOT.TLatex()
        lat.SetNDC()
        lat.SetTextFont(43)
        lat.SetTextSize(HDR)
        lat.SetTextAlign(12)
        lat.DrawLatex(0.16, 0.955, f"{label}, {args.cat}")
        keep += [leg, lat]
        drawn += 1

    if skipped:
        print("    too few entries to draw: " + ", ".join(skipped))
    if not drawn:
        print("    the dominant samples have no usable multiplicity classes")
        return None
    stem = os.path.join(outdir, "bkgMass_noBDT_dominant")
    c.SaveAs(f"{stem}.png")
    print(f"    wrote {stem}.png")
    return keep


def draw_nobdt(hincl, hsamp, args, outdir):
    """The pairing on its own: no BDT cut anywhere in these plots."""
    os.makedirs(outdir, exist_ok=True)
    keep = []

    # 1. the clean test: n=2 against n>=3, all samples, no BDT selection
    keep.append(overlay(
        [(hincl["n=2"], "n=2 good muons (no pairing freedom)"),
         (hincl["n>=3"], "n#geq3 good muons (pairing chooses)")],
        "m_{#mu#mu} [GeV]",
        f"{args.cat} background, NO BDT cut -- pairing only",
        os.path.join(outdir, "bkgMass_noBDT_pairing"), args))

    # 2. per sample, so one sample's topology cannot look like sculpting
    for cname, tag in (("n=2", "n2"), ("n>=3", "n3ormore")):
        ents = []
        for sid in sorted({s for s, _ in hsamp},
                          key=lambda s: -hsamp[(s, cname)].GetEntries()):
            h = hsamp[(sid, cname)]
            if h.GetEntries() < args.min_entries:
                continue
            ents.append((h, SAMPLE_LABEL.get(sid, str(sid))))
        keep.append(overlay(
            ents, "m_{#mu#mu} [GeV]",
            f"{args.cat}, NO BDT cut, {cname} good muons, by sample "
            f"(w_allSF yields, log scale)",
            os.path.join(outdir, f"bkgMass_noBDT_bySample_{tag}"), args,
            normalize=False, logy=True))
    return keep


def report_nobdt(hincl, hsamp, args):
    print("\nNO BDT cut -- the pairing acting on its own")
    print(f"  {'class':<10}{'entries':>10}{'yield':>12}"
          f"{'in 120-130':>12}{'fraction':>14}")
    print("  " + "-" * 58)
    fr = {}
    for cname in ("n=2", "n>=3"):
        h = hincl[cname]
        f, ef = peak_frac(h)
        fr[cname] = (f, ef)
        b1 = h.FindFixBin(PEAK[0] + 1e-6)
        b2 = h.FindFixBin(PEAK[1] - 1e-6)
        print(f"  {cname:<10}{h.GetEntries():>10.0f}{h.Integral():>12.4g}"
              f"{h.Integral(b1, b2):>12.4g}"
              f"{100*f:>10.1f}% +-{100*ef:.1f}")
    (fa, ea), (fb, eb) = fr["n=2"], fr["n>=3"]
    d, ed = fb - fa, np.sqrt(ea*ea + eb*eb)
    print(f"  difference (n>=3 minus n=2): {100*d:+.1f} +-{100*ed:.1f} points"
          f"  ({abs(d)/ed if ed else 0:.1f} sigma)")
    print("  This is the pairing test with the BDT removed entirely. Compare")
    print("  it with the per-BDT-bin table: if the difference is flat in BDT")
    print("  and matches this, the effect is the pairing; if it appears only")
    print("  at high BDT, it is the BDT reshaping an already-selected sample.")

    print("\n  per sample, no BDT cut")
    print(f"    {'sample':<12}{'n=2 frac':>11}{'n>=3 frac':>12}{'difference':>16}")
    print("    " + "-" * 51)
    for sid in sorted({s for s, _ in hsamp}):
        h2, h3 = hsamp[(sid, "n=2")], hsamp[(sid, "n>=3")]
        if h3.GetEntries() < args.min_entries:
            continue
        f2, e2 = peak_frac(h2)
        f3, e3 = peak_frac(h3)
        dd, edd = f3 - f2, np.sqrt(e2*e2 + e3*e3)
        print(f"    {SAMPLE_LABEL.get(sid, sid):<12}{100*f2:>10.1f}%"
              f"{100*f3:>11.1f}%{100*dd:>+11.1f} +-{100*edd:.1f}")
    print(f"    (samples with fewer than {args.min_entries} n>=3 entries omitted)")


def draw(hists, args, outdir, nmva):
    os.makedirs(outdir, exist_ok=True)
    edges = args.bdt_bins
    lo, hi = args.mass_range
    cols = [ROOT.kViolet + 6, ROOT.kBlue + 1, ROOT.kAzure + 1,
            ROOT.kTeal + 3, ROOT.kGreen + 2, ROOT.kSpring - 1,
            ROOT.kYellow + 2, ROOT.kOrange + 1, ROOT.kRed + 1,
            ROOT.kPink + 7]
    keep = []

    # --- all BDT slices overlaid, shapes -------------------------------
    for cname, tag in (("n=2", "n2"), ("n>=3", "n3ormore")):
        c = ROOT.TCanvas("c_" + tag, "", 860, 640)
        c.SetLeftMargin(0.13)
        c.SetRightMargin(0.04)
        c.SetBottomMargin(0.14)
        c.SetTopMargin(0.09)
        leg = ROOT.TLegend(0.68, 0.45, 0.96, 0.88)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(43)
        leg.SetTextSize(int(0.8 * LAB))
        drawn, hs = False, []
        for k in range(len(edges) - 1):
            h = hists[(k, cname)].Clone(f"s_{tag}_{k}")
            if h.Integral() <= 0:
                continue
            h.Scale(1.0 / h.Integral())
            h.SetLineColor(cols[k % len(cols)])
            h.SetLineWidth(2)
            hs.append(h)
            leg.AddEntry(h, f"{edges[k]:.1f}-{edges[k+1]:.1f} "
                            f"({hists[(k,cname)].GetEntries():.0f})", "l")
        if not hs:
            continue
        ymax = max(h.GetMaximum() for h in hs)
        for h in hs:
            h.SetMaximum(1.45 * ymax)
            h.SetMinimum(0.0)
            for ax, t in ((h.GetXaxis(), "m_{#mu#mu} [GeV]"),
                          (h.GetYaxis(), "Normalized Events")):
                ax.SetTitle(t)
                ax.SetTitleFont(43)
                ax.SetTitleSize(TIT)
                ax.SetLabelFont(43)
                ax.SetLabelSize(LAB)
            h.GetXaxis().SetTitleOffset(1.4)
            h.GetYaxis().SetTitleOffset(1.6)
            h.SetMarkerStyle(1)          # no marker, error bars only
            h.Draw("HIST E" if not drawn else "HIST E SAME")
            drawn = True
        leg.Draw()
        for x in PEAK:
            if lo < x < hi:
                ln = ROOT.TLine(x, 0, x, 1.45 * ymax)
                ln.SetLineColor(ROOT.kGray + 2)
                ln.SetLineStyle(2)
                ln.Draw()
                keep.append(ln)
        lat = ROOT.TLatex()
        lat.SetNDC()
        lat.SetTextFont(43)
        lat.SetTextSize(HDR)
        lat.SetTextAlign(12)
        lat.DrawLatex(0.13, 0.955,
                      f"{args.cat} background, {nmva} slices, {cname} good muons")
        stem = os.path.join(outdir, f"bkgMass_BDTslices_{tag}")
        c.SaveAs(f"{stem}.png")
        print(f"    wrote {stem}.png")
        keep.append((c, hs, leg, lat))

    # --- the control: n=2 against n>=3, per BDT bin --------------------
    live = [k for k in range(len(edges) - 1)
            if hists[(k, "n=2")].Integral() > 0
            and hists[(k, "n>=3")].Integral() > 0]
    if not live:
        print("    no BDT bin has both multiplicity classes populated, "
              "so the control plot is skipped")
        return keep

    nx = min(4, len(live))
    ny = int(np.ceil(len(live) / nx))

    # Divided pads are small, and TitleOffset with pixel fonts (font 43) puts
    # the axis titles a fixed number of pixels from the axis regardless of pad
    # size, which is what pushed them outside the margin. The titles are drawn
    # by hand in pad-local NDC instead, and the margins are set from the pixel
    # widths that have to fit inside them.
    PW, PH = 430, 370                  # pad size in pixels
    lab_px, tit_px, leg_px = 15, 17, 14
    mL, mR, mB, mT = 0.20, 0.05, 0.19, 0.13

    c = ROOT.TCanvas("c_ctrl", "", PW * nx, PH * ny)
    c.Divide(nx, ny, 0.004, 0.004)
    for p, k in enumerate(live, 1):
        c.cd(p)
        ROOT.gPad.SetLeftMargin(mL)
        ROOT.gPad.SetRightMargin(mR)
        ROOT.gPad.SetBottomMargin(mB)
        ROOT.gPad.SetTopMargin(mT)
        h2 = hists[(k, "n=2")].Clone(f"c2_{k}")
        h3 = hists[(k, "n>=3")].Clone(f"c3_{k}")
        for h, col in ((h2, ROOT.kAzure + 2), (h3, ROOT.kRed + 1)):
            if h.Integral() > 0:
                h.Scale(1.0 / h.Integral())
            h.SetLineColor(col)
            h.SetLineWidth(2)
            h.SetMarkerColor(col)
            h.SetMarkerStyle(20)
            h.SetMarkerSize(0.6)
        ymax = max(h2.GetMaximum() + h2.GetBinError(h2.GetMaximumBin()),
                   h3.GetMaximum() + h3.GetBinError(h3.GetMaximumBin()))
        h2.SetMaximum(1.45 * ymax)
        h2.SetMinimum(0.0)
        # titles are hand-drawn below, so the axes carry labels only
        h2.GetXaxis().SetTitle("")
        h2.GetYaxis().SetTitle("")
        for ax in (h2.GetXaxis(), h2.GetYaxis()):
            ax.SetLabelFont(43)
            ax.SetLabelSize(lab_px)
        h2.GetXaxis().SetNdivisions(505)   # fewer labels, they crowd at this size
        h2.GetYaxis().SetNdivisions(505)
        h2.Draw("HIST E")
        h3.Draw("HIST E SAME")

        span_x, span_y = 1.0 - mL - mR, 1.0 - mB - mT
        xt = ROOT.TLatex()
        xt.SetNDC()
        xt.SetTextFont(43)
        xt.SetTextSize(tit_px)
        xt.SetTextAlign(22)
        xt.DrawLatex(mL + 0.5 * span_x, 0.7 * tit_px / PH, "m_{#mu#mu} [GeV]")
        yt = ROOT.TLatex()
        yt.SetNDC()
        yt.SetTextFont(43)
        yt.SetTextSize(tit_px)
        yt.SetTextAlign(22)
        yt.SetTextAngle(90)
        yt.DrawLatex(0.7 * tit_px / PW, mB + 0.5 * span_y, "Normalized Events")

        lg = ROOT.TLegend(mL + 0.42 * span_x, 1.0 - mT - 0.20 * span_y,
                          1.0 - mR - 0.01, 1.0 - mT - 0.01)
        lg.SetBorderSize(0)
        lg.SetFillStyle(0)
        lg.SetTextFont(43)
        lg.SetTextSize(leg_px)
        lg.AddEntry(h2, f"n=2 ({hists[(k,'n=2')].GetEntries():.0f})", "lep")
        lg.AddEntry(h3, f"n#geq3 ({hists[(k,'n>=3')].GetEntries():.0f})", "lep")
        lg.Draw()

        tx = ROOT.TLatex()
        tx.SetNDC()
        tx.SetTextFont(43)
        tx.SetTextSize(tit_px)
        tx.SetTextAlign(12)
        tx.DrawLatex(mL, 1.0 - 0.5 * mT, f"BDT {edges[k]:.1f}-{edges[k+1]:.1f}")
        keep.append((h2, h3, lg, tx, xt, yt))
    stem = os.path.join(outdir, "bkgMass_pairingControl")
    c.SaveAs(f"{stem}.png")
    print(f"    wrote {stem}.png")
    keep.append(c)
    return keep


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("year")
    p.add_argument("-d", "--indir",
                   default=f"/work/submit/{getpass.getuser()}/HmumuRun3/bkg")
    p.add_argument("--cat", default="VLcat", help="category dir (default: %(default)s)")
    p.add_argument("-s", "--samples", default=None,
                   help="comma-separated sample ids (default: the vv and vvv groups)")
    p.add_argument("--bdt-bins", default="0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0",
                   help="slice edges (default: %(default)s)")
    p.add_argument("--mass-range", nargs=2, type=float, default=[110.0, 150.0],
                   metavar=("LO", "HI"),
                   help="HiggsCandCorrMass range (default: %(default)s). "
                        "Use 70 200 to see the Z as well.")
    p.add_argument("--nbins", type=int, default=40,
                   help="mass bins (default: %(default)s; widen when statistics are thin)")
    p.add_argument("--mva", default=None,
                   help="discriminant branch; autodetected as discrMVA or discrMVA0")
    p.add_argument("--cat-filter", default="category>0",
                   help="extra filter on the category column (default: %(default)s)")
    p.add_argument("--no-cat", action="store_true",
                   help="drop the category filter entirely, to recover statistics: "
                        "the pairing does not care about the lepton counting")
    p.add_argument("-o", "--plotdir",
                   default=f"/home/submit/{getpass.getuser()}/public_html/"
                           f"HmumuRun3/bkgshape")
    p.add_argument("--min-entries", type=int, default=20,
                   help="omit a sample from the per-sample plots and table "
                        "when it has fewer n>=3 entries than this "
                        "(default: %(default)s)")
    p.add_argument("--dom-rebin", type=int, default=2,
                   help="merge this many mass bins in the dominant-sample "
                        "plot only, where n=4 is sparse (default: %(default)s)")
    p.add_argument("--dom-markers", action="store_true",
                   help="draw the dominant-sample plot as offset markers "
                        "instead of step lines")
    p.add_argument("-j", "--ncores", type=int, default=0)
    args = p.parse_args()

    args.bdt_bins = [float(x) for x in args.bdt_bins.split(",")]
    if len(args.bdt_bins) < 2:
        sys.exit("need at least two BDT edges")
    if args.no_cat:
        args.cat_filter = None
    if args.ncores != 1:
        ROOT.ROOT.EnableImplicitMT(args.ncores)

    ids = ([int(x) for x in args.samples.split(",")] if args.samples
           else DEFAULT_SAMPLES)
    paths = []
    for sid in ids:
        pth = os.path.join(args.indir, args.cat,
                           f"snapshot_mc_{sid}_{args.year}_{args.cat}.root")
        if os.path.exists(pth):
            paths.append((sid, pth))
        else:
            print(f"  missing: {pth}")
    if not paths:
        sys.exit("no snapshots found")

    print(f"{args.cat}, {len(paths)} sample(s), "
          f"mass {args.mass_range[0]:.0f}-{args.mass_range[1]:.0f} GeV, "
          f"{len(args.bdt_bins)-1} BDT slices")
    print("samples: " + ", ".join(SAMPLE_LABEL.get(s, str(s)) for s, _ in paths))
    if args.cat_filter:
        print(f"category filter: {args.cat_filter}")
    else:
        print("category filter: none (--no-cat)")

    hists, table, used, rows, hincl, hsamp, hmult = collect(paths, args)
    if hists is None:
        sys.exit("nothing to plot")
    mvas = sorted({m for _, m in used})
    print(f"discriminant: {', '.join(mvas)}")

    report_samples(rows, args)
    report_nobdt(hincl, hsamp, args)
    report(table, args)
    print()
    _keep = draw_nobdt(hincl, hsamp, args, args.plotdir)
    _keepd = draw_dominant(hmult, args, args.plotdir)
    _keep2 = draw(hists, args, args.plotdir, len(args.bdt_bins) - 1)


if __name__ == "__main__":
    main()
