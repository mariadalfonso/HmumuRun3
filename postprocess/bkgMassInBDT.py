#!/usr/bin/env python3
"""
bkgMassInBDT.py -- background dimuon mass in slices of the BDT discriminant,
for any category, with the controls that say WHY the mass shape moves.

    python bkgMassInBDT.py 2024 --cat VLcat
    python bkgMassInBDT.py 2024 --cat ggHcat -d ... --mass-range 70 200
    python bkgMassInBDT.py 2024 --cat VLcat -d ... -s 204,205 -o ~/public_html/HmumuRun3/bkgshape

The samples are whatever snapshots are in <indir>/<cat>/: the signal_hmm ids
are signal (or from --sig-indir), every other MC snapshot is background, so
the plots use whatever has been run; -s and --sig give explicit lists instead.
Samples are summed with w_allSF (times boson_ptWeight for the DY samples, as
in MACRO/prepareHisto.py), so the yields are expected events. Plots go to
<plotdir>/<cat>/.

Slicing the background mass in BDT bins shows whether the BDT shapes the mass,
which it does by construction -- the two are correlated. On its own it cannot
say why. Two controls separate the causes.

Control 1, every category: the process composition
--------------------------------------------------
The summed shape can move with the BDT score because the BDT reshapes the mass
within one process, through inputs correlated with m_mumu, or because it
changes the mix of processes, whose mass shapes differ (DY falls steeply,
ttbar is broad and flat). Each process group is sliced on its own: if a
group's shape is the same in every BDT slice and only the composition moves,
the summed shape changes because the mix changes. If the group's shape itself
moves, the BDT is sculpting that background.

Control 2, VLcat and TTLcat only: the good-muon multiplicity
------------------------------------------------------------
Whether the PAIRING shapes the mass. With exactly two good muons the category
vetoes leave no pairing freedom at all: the pair is both muons, and no
criterion can change that. With three or more the algorithm chooses. Both
populations sit in the same category and the same BDT bin and differ only in
whether a choice was made, so any excess near m_H in the n>=3 shape relative
to n==2 is the sculpting, directly.

The other categories veto extra muons (exactly two loose and two good), so
every event has n==2 and this control does not exist there; it is switched
off by default and --mult forces it on.

The script reports, per BDT bin and per process group or multiplicity:

    yield, entries, and the fraction of the mass range inside 120-130 GeV

If the 120-130 fraction agrees between n==2 and n>=3 within its uncertainty,
the pairing is not building structure under the peak. If n>=3 is high, it is.

Statistics are the limiting factor for control 2: category 3 needs three tight
muons above 20 GeV, so a VV sample may leave only a handful of events in the
fit range. The yields are printed first for that reason, and --no-cat drops
the category requirement (keeping nGoodMuons) to recover events when the test
would otherwise be unanswerable.
"""

import argparse
import getpass
import glob
import os
import sys

import numpy as np
import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

TREE = "events"
PEAK = (120.0, 130.0)

SAMPLE_LABEL = {10: "VBFH", 11: "ggH", 12: "WminusH", 13: "WplusH",
                14: "ZH", 15: "ttH", 17: "bbH",
                20: "VBFH_Zg", 21: "ggH_Zg", 22: "WminusH_Zg", 23: "WplusH_Zg",
                24: "ZH_Zg", 25: "ttH_Zg", 26: "ttH_Zg",
                130: "WtoMuNu",
                100: "DYto2L", 103: "DYto2Mu", 104: "DYto2Tau",
                109: "DYto2Mu_105-160", 126: "DYto2Mu_50-130",
                127: "DYto2Mu_130-200", 128: "DYto2Mu_200-400",
                101: "EWK_Z2J", 98: "EWK_Z2J_105-160_hw", 99: "EWK_Z2J_105-160_py",
                140: "TTto2L2Nu",
                105: "TWminus", 106: "TbarWplus", 107: "TTtoLNu2Q",
                201: "WZto2L2Q", 202: "ZZto2L2Nu", 203: "ZZto2L2Q",
                204: "WZto3LNu", 205: "ZZto4L", 206: "WWto2L2Nu",
                213: "WZZ", 214: "ZZZ", 215: "WWW", 216: "WWZ",
                221: "TTLNu-1J", 222: "TTLNu-EWK", 223: "TTLL_4-50",
                224: "TTLL_50", 225: "TTZ-QQ", 226: "TZQB",
                227: "TWZ_T2Q", 228: "TWZ_W2Q", 229: "TWZ_3L", 230: "TTTT",
                231: "TTWZ", 232: "TTWW", 233: "TTZZ", 234: "TTW-QQ",
                235: "TTTWminus", 236: "TTTWplus", 237: "TTNuNu",
                249: "TTG-1J"}

SIGNAL = [10, 11, 12, 13, 14, 15, 17]          # signal_hmm in samples.yaml
VV = [201, 202, 203, 204, 205, 206]
VVV = [213, 214, 215, 216]
# Categories whose muon selection leaves room for a third good muon, so the
# multiplicity control (control 2) can be formed. The others apply mu_veto.
MULT_CATS = {"VLcat", "TTLcat"}
CATS = ["ggHcat", "VBFcat", "Zinvcat", "VLcat", "VHcat", "TTLcat", "TTHcat"]


# Alternatives to a nominal sample: other DY QCD generators and mass bins
# (overlapping 100/103), the EWK Z samples generated in 105-160 GeV only
# (overlapping 101), the amcatnlo VV set, and systematic variations.
# Left out of a discovered set so they are not added on top of the nominal;
# -s still picks them deliberately.
ALTERNATIVES = {98, 99, 109, 126, 127, 128, *range(141, 150), *range(207, 213),
                *range(242, 246)}


def discover(indir, cat, year, data=False):
    """(id, path) of every MC snapshot in <indir>/<cat>/ for this year, or
    of every data snapshot with data=True (data carries negative ids)."""
    pre, post = "snapshot_mc_", f"_{year}_{cat}.root"
    out = []
    for path in glob.glob(os.path.join(indir, cat, f"{pre}*{post}")):
        sid = os.path.basename(path)[len(pre):-len(post)]
        if sid.lstrip("-").isdigit() and sid.startswith("-") == data:
            out.append((int(sid), path))
    return sorted(out)


def explicit(ids, indir, cat, year, what):
    """(id, path) for a given id list, reporting the ones not on disk."""
    out = []
    for sid in ids:
        pth = os.path.join(indir, cat, f"snapshot_mc_{sid}_{year}_{cat}.root")
        if os.path.exists(pth):
            out.append((sid, pth))
        else:
            print(f"  missing {what}: {pth}")
    return out


# process groups, in drawing order; a sample in no group falls into "other"
GROUPS = [("DY", {100, 103, 104, 109, 126, 127, 128}),
          ("EWK Z", {98, 99, 101}),
          ("ttbar", {140}),
          ("single top", {105, 106, 107}),
          ("ttV", {*range(221, 240), *range(246, 250)}),
          ("VV/VVV", set(VV + VVV)),
          ("W+jets", {129, 130}),
          ("H to Zgamma", set(range(20, 27)))]
# DY pT reweighting, applied to the same ids as MACRO/prepareHisto.py
DY_PTWEIGHT = {100, 103, 104, 109}


def group_of(sid):
    for name, ids in GROUPS:
        if sid in ids:
            return name
    return "other"


def tag_of(name):
    return name.replace(" ", "").replace("/", "")


LAB, TIT, HDR = 20, 23, 26


def find_mva(cols):
    """The discriminant branch, whose name depends on the snapshot vintage."""
    for c in ("discrMVA", "discrMVA0"):
        if c in cols:
            return c
    return None


def open_sample(sid, path, args, need=("HiggsCandCorrMass", "w_allSF")):
    """RDataFrame for one snapshot, weighted and cut to the mass range.

    Returns (df, inr, mva): the raw frame, the frame after the category,
    extra and mass-range filters with the weight defined as wgt, and the
    discriminant branch. None when the snapshot lacks a needed column.
    """
    lo, hi = args.mass_range
    label = SAMPLE_LABEL.get(sid, sid)
    # a job that is still writing, or died while writing, leaves a file ROOT
    # cannot open; skip it rather than stop the whole run
    try:                                   # newer ROOT raises, older returns null
        f = ROOT.TFile.Open(path)
    except OSError:
        f = None
    ok = bool(f) and not f.IsZombie() and bool(f.Get(TREE))
    if f:
        f.Close()
    if not ok:
        print(f"  {label}: {os.path.basename(path)} is unreadable or has no "
              f"'{TREE}' tree (job still running, or crashed?), skipped")
        return None
    df = ROOT.RDataFrame(TREE, path)
    cols = {str(c) for c in df.GetColumnNames()}
    mva = args.mva or find_mva(cols)
    if mva is None:
        print(f"  {label}: no discrMVA branch, skipped")
        return None
    missing = [c for c in need if c not in cols]
    if missing:
        print(f"  {label}: no {missing[0]}, skipped")
        return None
    wexpr = "w_allSF" if "w_allSF" in cols else "1.0"
    if (sid in DY_PTWEIGHT and "boson_ptWeight" in cols
            and not args.no_ptweight):
        wexpr = "w_allSF * boson_ptWeight"
    sel = df.Define("wgt", wexpr)
    if args.cat_filter and "category" in cols:
        sel = sel.Filter(args.cat_filter, "category")
    if args.filter:
        sel = sel.Filter(args.filter, "extra filter")
    inr = sel.Filter(f"HiggsCandCorrMass > {lo} && "
                     f"HiggsCandCorrMass < {hi}", "mass range")
    return df, inr, mva


def collect(paths, args, sig_paths=(), data_paths=()):
    """Book and fill every histogram, so each file is read once.

    Always, for control 1:
        incl[sid], slice[(sid, k)]   mass per sample, no BDT cut and per BDT bin
    With args.mult, for control 2:
        hists[(k, class)], table      summed over samples, per BDT bin
        cls_incl[class], hsamp[(sid, class)], hmult[(sid, exact n)]
    For the score plot, inside the mass range:
        score_bkg[sid], score_sig[sid]   BDT score per sample
        sig_slice[k]                     signal yield per BDT bin
    and in the sidebands only, outside the blinded window:
        score_bkg_sb[sid], score_data[sid], bkg_sb_slice[k], data_slice[k]
    Data is never booked inside the blinded window.
    Returns None when no sample could be used.
    """
    lo, hi = args.mass_range
    edges = args.bdt_bins
    classes = [("n=2", "nGoodMuons == 2"), ("n>=3", "nGoodMuons >= 3")]
    # exact multiplicities per sample, for the dominant-sample plot: the two
    # classes above lump n=3 with n=4, and in ZZ->4l those are different
    # topologies rather than different amounts of the same one
    MULT = [("n=2", "nGoodMuons == 2"), ("n=3", "nGoodMuons == 3"),
            ("n=4", "nGoodMuons == 4"), ("n>=5", "nGoodMuons >= 5")]

    def book(node, name):
        return node.Histo1D((name, "", args.nbins, lo, hi),
                            "HiggsCandCorrMass", "wgt")

    def book_score(node, mva, name):
        return node.Histo1D((name, "", args.score_bins, edges[0], edges[-1]),
                            mva, "wgt")

    b_incl, b_slice, b_cls, b_x, b_m = [], [], [], [], []
    b_sbkg, b_ssig, b_sslice = [], [], []
    b_sbsb, b_bsbslice, b_sdata, b_dslice = [], [], [], []
    blo, bhi = args.blind
    sideband = (f"!(HiggsCandCorrMass > {blo} && HiggsCandCorrMass < {bhi})")
    persample, used = [], []
    for sid, path in paths:
        opened = open_sample(sid, path, args,
                             ("HiggsCandCorrMass", "w_allSF", "nGoodMuons"))
        if opened is None:
            continue
        df, inr, mva = opened
        # per-sample totals inside the mass range, so the contribution of
        # each sample to the sum is visible rather than assumed
        persample.append((
            sid, df.Count(), inr.Count(), inr.Sum("wgt"),
            inr.Filter("nGoodMuons == 2").Count(),
            inr.Filter("nGoodMuons >= 3").Count(),
            inr.Filter("nGoodMuons >= 3").Sum("wgt")))
        b_incl.append((sid, book(inr, f"x_{sid}")))
        b_sbkg.append((sid, book_score(inr, mva, f"sb_{sid}")))
        if data_paths:
            side = inr.Filter(sideband, "sidebands")
            b_sbsb.append((sid, book_score(side, mva, f"sbsb_{sid}")))
            for k in range(len(edges) - 1):
                a, b = edges[k], edges[k + 1]
                b_bsbslice.append((k, side.Filter(f"{mva} >= {a} && {mva} < {b}")
                                          .Sum("wgt")))
        for k in range(len(edges) - 1):
            a, b = edges[k], edges[k + 1]
            inbin = inr.Filter(f"{mva} >= {a} && {mva} < {b}", f"bdt{k}")
            b_slice.append((sid, k, book(inbin, f"t_{sid}_{k}")))
            if args.mult:
                for cname, cut in classes:
                    b_cls.append((k, cname, book(inbin.Filter(cut),
                                                 f"t_{sid}_{k}_{cname}")))
        if args.mult:
            # BDT-inclusive: the pairing acting on its own, with no BDT
            # selection to confound it. Kept per sample as well.
            for cname, cut in classes:
                b_x.append((sid, cname, book(inr.Filter(cut),
                                             f"x_{sid}_{cname}")))
            for mname, cut in MULT:
                b_m.append((sid, mname, book(inr.Filter(cut),
                                             f"m_{sid}_{mname}")))
        used.append((sid, mva))

    for sid, path in sig_paths:
        opened = open_sample(sid, path, args)
        if opened is None:
            continue
        _, inr, mva = opened
        b_ssig.append((sid, book_score(inr, mva, f"ss_{sid}")))
        for k in range(len(edges) - 1):
            a, b = edges[k], edges[k + 1]
            b_sslice.append((k, inr.Filter(f"{mva} >= {a} && {mva} < {b}")
                                   .Sum("wgt")))

    for sid, path in data_paths:
        opened = open_sample(sid, path, args, ("HiggsCandCorrMass",))
        if opened is None:
            continue
        _, inr, mva = opened
        # blinded: the sideband filter comes before anything is booked
        side = inr.Filter(sideband, "sidebands")
        b_sdata.append((sid, side.Histo1D(
            (f"sd_{sid}", "", args.score_bins, edges[0], edges[-1]), mva)))
        for k in range(len(edges) - 1):
            a, b = edges[k], edges[k + 1]
            b_dslice.append((k, side.Filter(f"{mva} >= {a} && {mva} < {b}")
                                    .Count()))

    if not b_incl:
        return None

    res = dict(used=used)
    res["incl"] = {sid: hh.GetPtr().Clone(f"incl_{sid}")   # triggers the loops
                   for sid, hh in b_incl}
    res["slice"] = {(sid, k): hh.GetPtr().Clone(f"sl_{sid}_{k}")
                    for sid, k, hh in b_slice}
    res["score_bkg"] = {sid: hh.GetPtr().Clone(f"scb_{sid}")
                        for sid, hh in b_sbkg}
    res["score_sig"] = {sid: hh.GetPtr().Clone(f"scs_{sid}")
                        for sid, hh in b_ssig}
    res["sig_slice"] = {k: 0.0 for k in range(len(edges) - 1)}
    for k, v in b_sslice:
        res["sig_slice"][k] += v.GetValue()
    res["score_bkg_sb"] = {sid: hh.GetPtr().Clone(f"scbsb_{sid}")
                           for sid, hh in b_sbsb}
    res["score_data"] = {sid: hh.GetPtr().Clone(f"scd_{sid}")
                         for sid, hh in b_sdata}
    res["bkg_sb_slice"] = {k: 0.0 for k in range(len(edges) - 1)}
    for k, v in b_bsbslice:
        res["bkg_sb_slice"][k] += v.GetValue()
    res["data_slice"] = {k: 0 for k in range(len(edges) - 1)}
    for k, v in b_dslice:
        res["data_slice"][k] += v.GetValue()
    res["rows"] = [dict(sid=sid, n_tot=a.GetValue(), n=b.GetValue(),
                        y=c.GetValue(), n2=d.GetValue(), n3=e.GetValue(),
                        y3=f.GetValue())
                   for sid, a, b, c, d, e, f in persample]
    if not args.mult:
        return res

    hists = {}
    for k in range(len(edges) - 1):
        for cname, _ in classes:
            hists[(k, cname)] = summed({}, [], f"h_{k}_{cname}", args)
    for k, cname, hh in b_cls:
        hists[(k, cname)].Add(hh.GetPtr())
    cls_incl = {cname: summed({}, [], f"hincl_{cname}", args)
                for cname, _ in classes}
    hsamp = {}
    for sid, cname, hh in b_x:
        h = hh.GetPtr().Clone(f"s_{sid}_{cname}")
        hsamp[(sid, cname)] = h
        cls_incl[cname].Add(h)
    hmult = {(sid, mname): hh.GetPtr().Clone(f"mm_{sid}_{mname}")
             for sid, mname, hh in b_m}

    table = {}
    for (k, cname), h in hists.items():
        y = h.Integral()
        b1 = h.FindFixBin(PEAK[0] + 1e-6)
        b2 = h.FindFixBin(PEAK[1] - 1e-6)
        ypk = h.Integral(b1, b2)
        table[(k, cname)] = dict(entries=h.GetEntries(), neff=n_eff(h),
                                 yield_=y, peak=ypk,
                                 frac=(ypk / y if y else 0.0),
                                 err=peak_frac(h)[1])
    res.update(hists=hists, table=table, cls_incl=cls_incl, hsamp=hsamp,
               hmult=hmult)
    return res


def summed(hmap, keys, name, args):
    """Sum hmap[key] over keys into a fresh histogram (empty for no keys)."""
    lo, hi = args.mass_range
    h = ROOT.TH1D(name, "", args.nbins, lo, hi)
    h.Sumw2()
    for key in keys:
        if key in hmap:
            h.Add(hmap[key])
    return h


def build_groups(res, args):
    """Group-level and total histograms, inclusive and per BDT bin."""
    hincl, hslice = res["incl"], res["slice"]
    nk = len(args.bdt_bins) - 1
    sids = sorted(hincl)
    names = [g for g, _ in GROUPS] + ["other"]
    members = {g: [s for s in sids if group_of(s) == g] for g in names}
    members = {g: m for g, m in members.items() if m}

    gincl = {g: summed(hincl, m, f"gi_{tag_of(g)}", args)
             for g, m in members.items()}
    gslice = {(g, k): summed(hslice, [(s, k) for s in m],
                             f"gs_{tag_of(g)}_{k}", args)
              for g, m in members.items() for k in range(nk)}
    tincl = summed(hincl, sids, "tot_incl", args)
    tslice = {k: summed(hslice, [(s, k) for s in sids], f"tot_{k}", args)
              for k in range(nk)}
    return members, gincl, gslice, tincl, tslice


def report_samples(rows, args):
    lo, hi = args.mass_range
    print(f"\nper-sample contribution to the sum   "
          f"(mass {lo:.0f}-{hi:.0f} GeV, no BDT cut)")
    print(f"  {'sample':<17}{'group':<12}{'snapshot':>11}{'in range':>10}"
          f"{'yield':>11}{'share':>8}{'n=2':>9}{'n>=3':>8}{'yield n>=3':>12}")
    print("  " + "-" * 98)
    tot = dict(n_tot=0, n=0, y=0.0, n2=0, n3=0, y3=0.0)
    for r in rows:
        for k in tot:
            tot[k] += r[k]
    for r in sorted(rows, key=lambda r: -r["y"]):
        share = r["y"] / tot["y"] if tot["y"] else 0.0
        print(f"  {str(SAMPLE_LABEL.get(r['sid'], r['sid'])):<17}"
              f"{group_of(r['sid']):<12}{r['n_tot']:>11,}{r['n']:>10,}"
              f"{r['y']:>11.4g}{100*share:>7.1f}%"
              f"{r['n2']:>9,}{r['n3']:>8,}{r['y3']:>12.4g}")
    print("  " + "-" * 98)
    print(f"  {'SUM':<29}{tot['n_tot']:>11,}{tot['n']:>10,}{tot['y']:>11.4g}"
          f"{'':>8}{tot['n2']:>9,}{tot['n3']:>8,}{tot['y3']:>12.4g}")
    print("  all samples above are added into one histogram, weighted by "
          "w_allSF")
    if args.mult and tot["n3"] == 0:
        print("  WARNING: no event has three or more good muons, so the "
              "pairing control\n           cannot be formed at all")
    if not args.mult and tot["n3"] > 0:
        print(f"  WARNING: {tot['n3']} events with nGoodMuons >= 3; the muon "
              f"veto of {args.cat}\n           should leave none (--mult "
              "runs the pairing control on them)")
    return tot


def report_composition(members, gslice, tslice, args):
    edges = args.bdt_bins
    groups = list(members)
    print("\ncomposition per BDT bin   (share of the summed yield)")
    wd = max(12, max(len(g) for g in groups) + 2)
    print(f"  {'BDT bin':<12}{'yield':>11}"
          + "".join(f"{g:>{wd}}" for g in groups))
    print("  " + "-" * (23 + wd * len(groups)))
    for k in range(len(edges) - 1):
        y = tslice[k].Integral()
        if tslice[k].GetEntries() == 0:
            continue
        cells = "".join(
            f"{100 * gslice[(g, k)].Integral() / y if y else 0:>{wd - 1}.1f}%"
            for g in groups)
        print(f"  {edges[k]:.2f}-{edges[k+1]:<7.2f}{y:>11.4g}{cells}")


def report_peak(members, gincl, gslice, tincl, tslice, args):
    edges = args.bdt_bins
    nk = len(edges) - 1
    print(f"\nfraction of the {args.mass_range[0]:.0f}-"
          f"{args.mass_range[1]:.0f} range inside {PEAK[0]:.0f}-{PEAK[1]:.0f} GeV,"
          " per BDT bin,\ncompared with the same sample with no BDT cut")
    rows = [("all background", tincl, tslice)]
    if len(members) > 1:
        rows += [(g, gincl[g], {k: gslice[(g, k)] for k in range(nk)})
                 for g in members]
    for name, hi, hs in rows:
        f0, e0 = peak_frac(hi)
        print(f"\n  {name}: no BDT cut {100*f0:.1f}% +-{100*e0:.1f}  "
              f"({hi.GetEntries():.0f} entries, {n_eff(hi):.0f} effective)")
        print(f"    {'BDT bin':<12}{'entries':>10}{'eff.':>9}{'yield':>12}"
              f"{'fraction':>14}{'minus no-BDT':>18}")
        print("    " + "-" * 75)
        for k in range(nk):
            h = hs[k]
            if n_eff(h) < args.min_entries:
                continue
            f, ef = peak_frac(h)
            # the slice is a subset of the inclusive sample, so the two errors
            # are correlated; the slice error alone is the honest one to quote
            d = f - f0
            sig = abs(d) / ef if ef else 0.0
            print(f"    {edges[k]:.2f}-{edges[k+1]:<7.2f}{h.GetEntries():>10.0f}"
                  f"{n_eff(h):>9.0f}{h.Integral():>12.4g}{100*f:>9.1f}% +-{100*ef:.1f}"
                  f"{100*d:>+10.1f}  ({sig:.1f} sigma)")
    print(f"\n  (BDT bins with fewer than {args.min_entries} effective entries "
          "omitted; eff. = (sum w)^2 / sum w^2)")
    print("  A process whose fraction is flat in BDT is not being sculpted;")
    print("  if only the summed fraction moves, the BDT is changing the mix.")


def report(table, args):
    edges = args.bdt_bins
    print(f"\n{'BDT bin':<14}{'class':>7}{'entries':>10}{'eff.':>9}{'yield':>12}"
          f"{'in 120-130':>12}{'fraction':>11}")
    print("-" * 75)
    for k in range(len(edges) - 1):
        for cname in ("n=2", "n>=3"):
            t = table[(k, cname)]
            if t["entries"] == 0:
                continue
            f, ef = t["frac"], t["err"]
            print(f"{edges[k]:.1f}-{edges[k+1]:<9.1f}{cname:>7}"
                  f"{t['entries']:>10.0f}{t['neff']:>9.0f}{t['yield_']:>12.4g}"
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
        ea, eb = a["err"], b["err"]
        d, ed = fb - fa, np.sqrt(ea*ea + eb*eb)
        sig = abs(d) / ed if ed else 0.0
        print(f"{edges[k]:.1f}-{edges[k+1]:<9.1f}{100*fa:>11.1f}%"
              f"{100*fb:>11.1f}%{100*d:>+9.1f} +-{100*ed:.1f}"
              f"  ({sig:.1f} sigma)")
    print("\na difference compatible with zero means the pairing is not "
          "building\nstructure under the peak in this sample")


def peak_frac(h):
    """(fraction of the range inside the peak window, its error).

    Each event's uncertainty is scaled by its own weight, so the error is
    propagated from the sums of squared weights inside and outside the
    window (kept by Sumw2):

        var(f) = [(1-f)^2 sum_in w^2 + f^2 sum_out w^2] / (sum w)^2

    With equal weights this is the binomial sqrt(f(1-f)/N); with a few
    heavy events it is larger, as it should be.
    """
    y = h.Integral()
    if y <= 0:
        return 0.0, 0.0
    b1 = h.FindFixBin(PEAK[0] + 1e-6)
    b2 = h.FindFixBin(PEAK[1] - 1e-6)
    f = h.Integral(b1, b2) / y
    w2 = [h.GetBinError(b) ** 2 for b in range(1, h.GetNbinsX() + 1)]
    w2_in = sum(w2[b1 - 1:b2])
    w2_out = sum(w2) - w2_in
    var = ((1 - f) ** 2 * w2_in + f ** 2 * w2_out) / (y * y)
    return f, np.sqrt(max(var, 0.0))


def n_eff(h):
    """Effective number of MC events, (sum w)^2 / sum w^2, over the range.

    Equal to the raw count for equal weights; smaller when a few heavy
    events carry the yield, which is what sets the statistical power.
    """
    w2 = sum(h.GetBinError(b) ** 2 for b in range(1, h.GetNbinsX() + 1))
    y = h.Integral()
    return y * y / w2 if w2 > 0 else 0.0


def overlay(entries, xtitle, header, stem, args, legend_right=True,
            normalize=True, logy=False, legend_x0=None, colors=None):
    """entries: list of (TH1D, legend label).

    colors overrides the default categorical colours, one per entry.

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
            ROOT.kViolet + 6, ROOT.kCyan + 2, ROOT.kMagenta + 2,
            ROOT.kSpring - 1, ROOT.kGray + 2, ROOT.kPink + 7]
    x0 = legend_x0 if legend_x0 is not None else (0.60 if legend_right else 0.17)
    leg = ROOT.TLegend(x0, 0.88 - 0.055 * len(live), x0 + 0.34, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(43)
    leg.SetTextSize(LAB)
    hs = []
    if colors:
        cols = [colors[k] for k, (h, _) in enumerate(entries) if h.Integral() > 0]
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
        leg.AddEntry(h, f"{lab}  ({n_eff(h0):.0f} eff. ev, "
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
        # headroom for the legend, which grows by one row per entry: the
        # curves sit below it as long as the axis spans enough decades
        rows = 0.055 * len(live) + 0.05
        ytop = ymax * (ymax / ylo) ** (rows / (1.0 - rows)) if ymax > ylo else 20 * ymax
        ytop = max(ytop, 20.0 * ymax)
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
            if n_eff(h0) < args.min_entries:
                skipped.append(f"{label} {m} ({n_eff(h0):.0f} eff.)")
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
            if n_eff(h) < args.min_entries:
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
        if n_eff(h3) < args.min_entries:
            continue
        f2, e2 = peak_frac(h2)
        f3, e3 = peak_frac(h3)
        dd, edd = f3 - f2, np.sqrt(e2*e2 + e3*e3)
        print(f"    {SAMPLE_LABEL.get(sid, sid):<12}{100*f2:>10.1f}%"
              f"{100*f3:>11.1f}%{100*dd:>+11.1f} +-{100*edd:.1f}")
    print(f"    (samples with fewer than {args.min_entries} effective n>=3 "
          "entries omitted)")


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

    panels = [(f"BDT {edges[k]:.1f}-{edges[k+1]:.1f}",
               [(hists[(k, "n=2")],
                 f"n=2 ({hists[(k, 'n=2')].GetEntries():.0f})", ROOT.kAzure + 2),
                (hists[(k, "n>=3")],
                 f"n#geq3 ({hists[(k, 'n>=3')].GetEntries():.0f})", ROOT.kRed + 1)])
              for k in live]
    keep.append(draw_grid(panels, os.path.join(outdir, "bkgMass_pairingControl"),
                          "c_ctrl"))
    return keep


def draw_grid(panels, stem, cname):
    """One pad per BDT bin, the shapes in it overlaid, each scaled to unit area.

    panels: list of (pad title, [(TH1D, legend label, colour), ...]); the
    first histogram of a pad sets its axes.
    """
    nx = min(4, len(panels))
    ny = int(np.ceil(len(panels) / nx))

    # Divided pads are small, and TitleOffset with pixel fonts (font 43) puts
    # the axis titles a fixed number of pixels from the axis regardless of pad
    # size, which is what pushed them outside the margin. The titles are drawn
    # by hand in pad-local NDC instead, and the margins are set from the pixel
    # widths that have to fit inside them.
    PW, PH = 430, 370                  # pad size in pixels
    lab_px, tit_px, leg_px = 15, 17, 14
    mL, mR, mB, mT = 0.20, 0.05, 0.19, 0.13

    keep = []
    c = ROOT.TCanvas(cname, "", PW * nx, PH * ny)
    c.Divide(nx, ny, 0.004, 0.004)
    for p, (title, entries) in enumerate(panels, 1):
        c.cd(p)
        ROOT.gPad.SetLeftMargin(mL)
        ROOT.gPad.SetRightMargin(mR)
        ROOT.gPad.SetBottomMargin(mB)
        ROOT.gPad.SetTopMargin(mT)
        hs = []
        for i, (h0, _, col) in enumerate(entries):
            h = h0.Clone(f"{cname}_{p}_{i}")
            if h.Integral() > 0:
                h.Scale(1.0 / h.Integral())
            h.SetLineColor(col)
            h.SetLineWidth(2)
            h.SetMarkerColor(col)
            h.SetMarkerStyle(20)
            h.SetMarkerSize(0.6)
            hs.append(h)
        ymax = max(h.GetMaximum() + h.GetBinError(h.GetMaximumBin())
                   for h in hs)
        h1 = hs[0]
        h1.SetMaximum(1.45 * ymax)
        h1.SetMinimum(0.0)
        # titles are hand-drawn below, so the axes carry labels only
        h1.GetXaxis().SetTitle("")
        h1.GetYaxis().SetTitle("")
        for ax in (h1.GetXaxis(), h1.GetYaxis()):
            ax.SetLabelFont(43)
            ax.SetLabelSize(lab_px)
        h1.GetXaxis().SetNdivisions(505)   # fewer labels, they crowd at this size
        h1.GetYaxis().SetNdivisions(505)
        for i, h in enumerate(hs):
            h.Draw("HIST E" if i == 0 else "HIST E SAME")

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

        lg = ROOT.TLegend(mL + 0.42 * span_x,
                          1.0 - mT - 0.10 * len(hs) * span_y,
                          1.0 - mR - 0.01, 1.0 - mT - 0.01)
        lg.SetBorderSize(0)
        lg.SetFillStyle(0)
        lg.SetTextFont(43)
        lg.SetTextSize(leg_px)
        for h, (_, label, _) in zip(hs, entries):
            lg.AddEntry(h, label, "lep")
        lg.Draw()

        tx = ROOT.TLatex()
        tx.SetNDC()
        tx.SetTextFont(43)
        tx.SetTextSize(tit_px)
        tx.SetTextAlign(12)
        tx.DrawLatex(mL, 1.0 - 0.5 * mT, title)
        keep.append((hs, lg, tx, xt, yt))
    c.SaveAs(f"{stem}.png")
    print(f"    wrote {stem}.png")
    keep.append(c)
    return keep


def score_sum(hmap, name, args):
    edges = args.bdt_bins
    h = ROOT.TH1D(name, "", args.score_bins, edges[0], edges[-1])
    h.Sumw2()
    for hh in hmap.values():
        h.Add(hh)
    return h


def report_score(res, tslice, args):
    """Signal and background yield per BDT bin, inside the mass range."""
    edges = args.bdt_bins
    nk = len(edges) - 1
    lo, hi = args.mass_range
    print(f"\nsignal and background per BDT bin   (MC, mass {lo:.0f}-{hi:.0f} GeV)")
    print(f"  {'BDT bin':<12}{'S':>11}{'B':>12}{'S/B':>10}{'S/sqrt(B)':>11}"
          f"{'above: S':>11}{'B':>12}{'S/sqrt(B)':>11}")
    print("  " + "-" * 88)
    for k in range(nk):
        S, B = res["sig_slice"][k], tslice[k].Integral()
        # cumulative from the top: what a cut at this bin's lower edge keeps
        Sa = sum(res["sig_slice"][j] for j in range(k, nk))
        Ba = sum(tslice[j].Integral() for j in range(k, nk))
        print(f"  {edges[k]:.2f}-{edges[k+1]:<7.2f}{S:>11.4g}{B:>12.4g}"
              f"{(S / B if B > 0 else 0):>10.3g}"
              f"{(S / np.sqrt(B) if B > 0 else 0):>11.3g}"
              f"{Sa:>11.4g}{Ba:>12.4g}"
              f"{(Sa / np.sqrt(Ba) if Ba > 0 else 0):>11.3g}")
    print("  'above' columns: everything with a score at or above the bin's "
          "lower edge")


def report_sidebands(res, args):
    """Data against MC background in the sidebands, per BDT bin."""
    edges = args.bdt_bins
    nk = len(edges) - 1
    nd = sum(res["data_slice"].values())
    nb = sum(res["bkg_sb_slice"].values())
    print(f"\ndata and MC background in the sidebands per BDT bin   "
          f"({sideband_text(args)} GeV)")
    print(f"  {'BDT bin':<12}{'data':>10}{'MC bkg':>12}{'data/MC':>10}"
          f"{'data/MC, shape':>16}")
    print("  " + "-" * 60)
    for k in range(nk):
        d, b = res["data_slice"][k], res["bkg_sb_slice"][k]
        r = d / b if b > 0 else 0.0
        # shape: both normalised to their totals, so only the BDT-score
        # distribution is compared and the overall normalisation drops out
        rs = (d / nd) / (b / nb) if nd > 0 and nb > 0 and b > 0 else 0.0
        ers = rs / np.sqrt(d) if d > 0 else 0.0
        print(f"  {edges[k]:.2f}-{edges[k+1]:<7.2f}{d:>10,}{b:>12.4g}"
              f"{r:>10.3f}{rs:>10.3f} +-{ers:.3f}")
    print(f"  {'total':<12}{nd:>10,}{nb:>12.4g}{(nd / nb if nb else 0):>10.3f}")
    print("  data/MC is the overall normalisation; the shape column divides it "
          "out.\n  Errors are data statistics only.")


def sideband_text(args):
    """'110-120, 130-150' for the default mass range and blinded window."""
    lo, hi = args.mass_range
    blo, bhi = args.blind
    parts = []
    if lo < blo:
        parts.append(f"{lo:.0f}-{min(blo, hi):.0f}")
    if bhi < hi:
        parts.append(f"{max(bhi, lo):.0f}-{hi:.0f}")
    return ", ".join(parts) or "none"


def draw_score(res, args, outdir):
    """Signal against background BDT score, all scaled to unit area.

    With data: the data sidebands as points, and the MC background in the
    same sidebands as a dashed line, which is what the data compares to.
    """
    edges = args.bdt_bins
    lo, hi = args.mass_range
    hb0 = score_sum(res["score_bkg"], "score_bkg", args)
    hs0 = score_sum(res["score_sig"], "score_sig", args)
    hbs0 = score_sum(res["score_bkg_sb"], "score_bkg_sb", args)
    hd0 = score_sum(res["score_data"], "score_data", args)
    have_sig, have_data = hs0.Integral() > 0, hd0.Integral() > 0
    if hb0.Integral() <= 0 or not (have_sig or have_data):
        print("    no signal or data next to the background in the mass range, "
              "so the score plot is skipped")
        return None
    sbtext = sideband_text(args)
    keep = []
    for logy in (False, True):
        tag = "_log" if logy else ""
        c = ROOT.TCanvas(f"c_score{tag}", "", 880, 650)
        c.SetLeftMargin(0.14)
        c.SetRightMargin(0.04)
        c.SetBottomMargin(0.14)
        c.SetTopMargin(0.09)
        c.SetLogy(logy)
        hb = hb0.Clone(f"score_bkg_n{tag}")
        hb.Scale(1.0 / hb.Integral())
        hb.SetFillColorAlpha(ROOT.kAzure - 9, 0.7)
        hb.SetLineColor(ROOT.kAzure + 2)
        hb.SetLineWidth(2)
        shown = [hb]
        hs = hbs = hd = None
        if have_sig:
            hs = hs0.Clone(f"score_sig_n{tag}")
            hs.Scale(1.0 / hs.Integral())
            hs.SetLineColor(ROOT.kRed + 1)
            hs.SetLineWidth(3)
            hs.SetMarkerSize(0)
            shown.append(hs)
        if have_data:
            hd = hd0.Clone(f"score_data_n{tag}")
            hd.Scale(1.0 / hd.Integral())
            hd.SetMarkerStyle(20)
            hd.SetMarkerSize(0.9)
            hd.SetMarkerColor(ROOT.kBlack)
            hd.SetLineColor(ROOT.kBlack)
            shown.append(hd)
            if hbs0.Integral() > 0:
                hbs = hbs0.Clone(f"score_bkg_sb_n{tag}")
                hbs.Scale(1.0 / hbs.Integral())
                hbs.SetLineColor(ROOT.kAzure + 3)
                hbs.SetLineStyle(2)
                hbs.SetLineWidth(3)
                hbs.SetMarkerSize(0)
                shown.append(hbs)
        ymax = max(h.GetMaximum() for h in shown)
        # the share of the frame, from the top, that the legend takes (one
        # row per curve, plus the sideband header with data), with a small
        # gap; the axis is stretched so every curve stays below it
        nrows = len(shown) + bool(hd)
        frame_lo, frame_hi = 0.14, 0.91           # bottom and 1 - top margin
        leg_lo = 0.88 - 0.06 * nrows
        f = min((frame_hi - leg_lo) / (frame_hi - frame_lo) + 0.04, 0.8)
        if logy:
            pos = [h.GetBinContent(b) for h in shown
                   for b in range(1, h.GetNbinsX() + 1) if h.GetBinContent(b) > 0]
            ylo = 0.5 * min(pos)
            ytop = ylo * (ymax / ylo) ** (1.0 / (1.0 - f))
        else:
            ylo, ytop = 0.0, ymax / (1.0 - f)
        hb.SetMinimum(ylo)
        hb.SetMaximum(ytop)
        for ax, t in ((hb.GetXaxis(), "BDT score"),
                      (hb.GetYaxis(), "Normalized Events")):
            ax.SetTitle(t)
            ax.SetTitleFont(43)
            ax.SetTitleSize(TIT)
            ax.SetLabelFont(43)
            ax.SetLabelSize(LAB)
        hb.GetXaxis().SetTitleOffset(1.4)
        hb.GetYaxis().SetTitleOffset(1.6)
        hb.Draw("HIST")
        if hbs:
            hbs.Draw("HIST SAME")
        if hs:
            hs.Draw("HIST E SAME")
        if hd:
            hd.Draw("E1 P SAME")
        # the slice edges used by every other plot, so the two can be read
        # against each other
        lines = []
        for x in edges[1:-1]:
            ln = ROOT.TLine(x, ylo, x, ytop if not logy else ymax * 3)
            ln.SetLineColor(ROOT.kGray + 1)
            ln.SetLineStyle(3)
            ln.Draw()
            lines.append(ln)
        leg = ROOT.TLegend(0.36, leg_lo, 0.96, 0.88)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextFont(43)
        leg.SetTextSize(LAB)
        if hd:
            leg.SetHeader(f"sidebands: {sbtext} GeV")
        if hs:
            leg.AddEntry(hs, f"signal H#rightarrow#mu#mu  ({hs0.Integral():.3g} exp, "
                             f"{hs0.GetEntries():.0f} MC)", "l")
        leg.AddEntry(hb, f"background  ({hb0.Integral():.4g} exp, "
                         f"{hb0.GetEntries():.0f} MC)", "f")
        if hbs:
            leg.AddEntry(hbs, f"background, sidebands  ({hbs0.Integral():.4g} exp)",
                         "l")
        if hd:
            leg.AddEntry(hd, f"data, sidebands  ({hd0.GetEntries():.0f} events)",
                         "pe")
        leg.Draw()
        c.RedrawAxis()
        lat = ROOT.TLatex()
        lat.SetNDC()
        lat.SetTextFont(43)
        lat.SetTextSize(HDR)
        lat.SetTextAlign(12)
        what = "MC and data" if hd else "MC"
        lat.DrawLatex(0.14, 0.955, f"{args.cat}, {what}, "
                                   f"{lo:.0f} < m_{{#mu#mu}} < {hi:.0f} GeV")
        stem = os.path.join(outdir, f"bkgMass_BDTscore{tag}")
        c.SaveAs(f"{stem}.png")
        print(f"    wrote {stem}.png")
        keep.append((c, hb, hs, hbs, hd, leg, lat, lines))
    return keep


def draw_composition(members, gslice, tslice, args, outdir):
    """Stacked share of each process group per BDT bin."""
    edges = args.bdt_bins
    nk = len(edges) - 1
    cols = [ROOT.kAzure + 1, ROOT.kOrange - 2, ROOT.kGreen - 6,
            ROOT.kGreen + 3, ROOT.kViolet - 4, ROOT.kOrange + 7,
            ROOT.kYellow - 7, ROOT.kRed - 7, ROOT.kGray + 1]
    c = ROOT.TCanvas("c_comp", "", 880, 650)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.24)
    c.SetBottomMargin(0.16)
    c.SetTopMargin(0.09)
    stack = ROOT.THStack("st_comp", "")
    leg = ROOT.TLegend(0.77, 0.90 - 0.06 * len(members), 0.99, 0.90)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(43)
    leg.SetTextSize(LAB)
    hs = []
    # colour by the group's place in GROUPS, so a group keeps its colour
    # whichever others are present
    order = [g for g, _ in GROUPS] + ["other"]
    for g in members:
        i = order.index(g)
        h = ROOT.TH1D(f"comp_{tag_of(g)}", "", nk, 0, nk)
        for k in range(nk):
            y = tslice[k].Integral()
            h.SetBinContent(k + 1, gslice[(g, k)].Integral() / y if y > 0 else 0)
        h.SetFillColor(cols[i % len(cols)])
        h.SetLineColor(ROOT.kBlack)
        h.SetLineWidth(1)
        stack.Add(h)
        hs.append(h)
    for h, g in reversed(list(zip(hs, members))):
        leg.AddEntry(h, g, "f")
    stack.Draw("HIST")
    stack.SetMaximum(1.0)
    stack.SetMinimum(0.0)
    ax = stack.GetXaxis()
    for k in range(nk):
        ax.SetBinLabel(k + 1, f"{edges[k]:.1f}-{edges[k+1]:.1f}")
    ax.SetTitle("BDT bin")
    stack.GetYaxis().SetTitle("Fraction of background")
    for a in (ax, stack.GetYaxis()):
        a.SetTitleFont(43)
        a.SetTitleSize(TIT)
        a.SetLabelFont(43)
        a.SetLabelSize(LAB)
    ax.SetTitleOffset(1.7)
    stack.GetYaxis().SetTitleOffset(1.3)
    c.Modified()
    leg.Draw()
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(43)
    lat.SetTextSize(HDR)
    lat.SetTextAlign(12)
    lat.DrawLatex(0.12, 0.955, f"{args.cat} background composition, "
                               f"{args.mass_range[0]:.0f}-"
                               f"{args.mass_range[1]:.0f} GeV")
    stem = os.path.join(outdir, "bkgMass_composition")
    c.SaveAs(f"{stem}.png")
    print(f"    wrote {stem}.png")
    return c, stack, hs, leg, lat


def slice_colors(n):
    """n colours in order, dark blue to green to yellow-green (viridis without
    its pale yellow end, which does not show on white)."""
    ROOT.gStyle.SetPalette(ROOT.kViridis)
    pal = ROOT.TColor.GetPalette()
    top = int(0.85 * (pal.GetSize() - 1))
    return [pal.At(int(top * i / max(n - 1, 1))) for i in range(n)]


def draw_processes(members, gincl, gslice, tincl, tslice, args, outdir):
    """Control 1: the mass in BDT slices, summed and per process group."""
    edges = args.bdt_bins
    nk = len(edges) - 1
    keep = []

    scol = slice_colors(nk)

    def slices(hs):
        """(entries, colours) for one overlay, in BDT order."""
        ks = [k for k in range(nk) if n_eff(hs[k]) >= args.min_entries]
        return ([(hs[k], f"{edges[k]:.1f}-{edges[k+1]:.1f}") for k in ks],
                [scol[k] for k in ks])

    # 1. the processes with no BDT cut, as yields, so their size is visible
    keep.append(overlay(
        [(gincl[g], g) for g in members], "m_{#mu#mu} [GeV]",
        f"{args.cat} background, NO BDT cut, by process (yields, log scale)",
        os.path.join(outdir, "bkgMass_noBDT_byProcess"), args,
        normalize=False, logy=True, legend_x0=0.48))

    # 2. all background summed, shape per BDT slice
    ents, cols = slices(tslice)
    keep.append(overlay(
        ents, "m_{#mu#mu} [GeV]",
        f"{args.cat} background, {nk} BDT slices, all processes",
        os.path.join(outdir, "bkgMass_BDTslices_all"), args,
        legend_x0=0.52, colors=cols))

    # 2b. the same, one pad per BDT bin against the shape with no BDT cut,
    #     the pairing-control layout without the multiplicity split
    live = [k for k in range(nk) if n_eff(tslice[k]) >= args.min_entries]
    if live:
        panels = [(f"BDT {edges[k]:.1f}-{edges[k+1]:.1f}",
                   [(tincl, "no BDT cut", ROOT.kAzure + 2),
                    (tslice[k], f"this bin ({n_eff(tslice[k]):.0f} eff.)",
                     ROOT.kRed + 1)])
                  for k in live]
        keep.append(draw_grid(panels, os.path.join(outdir, "bkgMass_BDTbins"),
                              "c_bins"))

    # 3. the same, one process group at a time: the composition is fixed
    #    inside each plot, so any change left is the BDT shaping that process
    if len(members) > 1:
        for g in members:
            ents, cols = slices({k: gslice[(g, k)] for k in range(nk)})
            if len(ents) < 2:
                print(f"    {g}: fewer than two populated BDT slices, skipped")
                continue
            keep.append(overlay(
                ents, "m_{#mu#mu} [GeV]",
                f"{args.cat} {g}, in {nk} BDT slices",
                os.path.join(outdir, f"bkgMass_BDTslices_{tag_of(g)}"), args,
                legend_x0=0.52, colors=cols))

        keep.append(draw_composition(members, gslice, tslice, args, outdir))
    return keep


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("year")
    p.add_argument("-d", "--indir",
                   default=f"/work/submit/{getpass.getuser()}/HmumuRun3/ROOTFILES",
                   help="snapshot directory, <dir>/<cat>/ (default: %(default)s, "
                        "the Hmm.py output)")
    p.add_argument("--cat", default="VLcat", choices=CATS,
                   help="category dir (default: %(default)s)")
    p.add_argument("-s", "--samples", default=None,
                   help="comma-separated background ids (default: every MC "
                        "snapshot in <indir>/<cat>/ except signal and the "
                        "alternative samples)")
    p.add_argument("--sig", default=None,
                   help="signal ids for the score plot, '' for none (default: "
                        "every signal_hmm snapshot in <sig-indir>/<cat>/)")
    p.add_argument("--sig-indir", default=None,
                   help="signal snapshot directory, <dir>/<cat>/, when it is "
                        "not -d (default: the same as -d)")
    p.add_argument("--data", default=None,
                   help="data ids for the sidebands in the score plot, written "
                        "with '=' since they are negative (--data=-41,-42), '' "
                        "for none (default: every data snapshot in -d)")
    p.add_argument("--blind", nargs=2, type=float, default=[120.0, 130.0],
                   metavar=("LO", "HI"),
                   help="mass window where data is never looked at "
                        "(default: %(default)s, the same as the peak window)")
    p.add_argument("--score-bins", type=int, default=50,
                   help="bins of the BDT score plot (default: %(default)s)")
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
                   help="filter on the category column, applied only where the "
                        "snapshot has one (default: %(default)s)")
    p.add_argument("--no-cat", action="store_true",
                   help="drop the category filter entirely, to recover statistics: "
                        "the pairing does not care about the lepton counting")
    p.add_argument("--filter", default=None,
                   help="extra RDataFrame filter applied to every sample")
    p.add_argument("--mult", action=argparse.BooleanOptionalAction, default=None,
                   help="run the good-muon multiplicity control (default: on "
                        f"for {', '.join(sorted(MULT_CATS))}, off elsewhere)")
    p.add_argument("--no-ptweight", action="store_true",
                   help="do not apply boson_ptWeight to the DY samples")
    p.add_argument("-o", "--plotdir",
                   default=f"/home/submit/{getpass.getuser()}/public_html/"
                           f"HmumuRun3/bkgshape",
                   help="plots go to <plotdir>/<cat> (default: %(default)s)")
    p.add_argument("--min-entries", type=int, default=20,
                   help="omit a BDT slice, or a sample from the per-sample "
                        "multiplicity plots, when its effective number of MC "
                        "events, (sum w)^2 / sum w^2, is below this "
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
    if args.mult is None:
        args.mult = args.cat in MULT_CATS
    if args.ncores != 1:
        ROOT.ROOT.EnableImplicitMT(args.ncores)

    if args.samples:
        paths = explicit([int(x) for x in args.samples.split(",")],
                         args.indir, args.cat, args.year, "background")
    else:
        found = discover(args.indir, args.cat, args.year)
        paths = [(sid, p) for sid, p in found
                 if sid not in SIGNAL and sid not in ALTERNATIVES]
        left = [sid for sid, _ in found if sid in ALTERNATIVES]
        if left:
            print("  left out as alternatives to a nominal sample (-s to use): "
                  + ", ".join(str(SAMPLE_LABEL.get(x, x)) for x in left))
    if not paths:
        sys.exit(f"no background snapshots in {os.path.join(args.indir, args.cat)}")

    if args.sig == "":
        sig_paths = []
    elif args.sig:
        sig_paths = explicit([int(x) for x in args.sig.split(",")],
                             args.sig_indir or args.indir, args.cat, args.year,
                             "signal")
    else:
        sig_paths = [(sid, p) for sid, p in
                     discover(args.sig_indir or args.indir, args.cat, args.year)
                     if sid in SIGNAL]

    print(f"{args.cat}, {len(paths)} sample(s), "
          f"mass {args.mass_range[0]:.0f}-{args.mass_range[1]:.0f} GeV, "
          f"{len(args.bdt_bins)-1} BDT slices")
    print("samples: " + ", ".join(str(SAMPLE_LABEL.get(s, s)) for s, _ in paths))
    if args.data == "":
        data_paths = []
    elif args.data:
        data_paths = explicit([int(x) for x in args.data.split(",")],
                              args.indir, args.cat, args.year, "data")
    else:
        data_paths = discover(args.indir, args.cat, args.year, data=True)

    print("signal:  " + (", ".join(str(SAMPLE_LABEL.get(s, s)) for s, _ in sig_paths)
                         + f"  ({os.path.dirname(sig_paths[0][1])})"
                         if sig_paths else
                         "off (--sig '')" if args.sig == "" else
                         "none found, so no score plot"))
    if args.cat_filter:
        print(f"category filter: {args.cat_filter}")
    else:
        print("category filter: none (--no-cat)")
    if args.filter:
        print(f"extra filter: {args.filter}")
    print("multiplicity control: " + ("on" if args.mult else "off"))

    print(f"data:    {len(data_paths)} file(s), sidebands {sideband_text(args)} GeV "
          f"only (blinded {args.blind[0]:.0f}-{args.blind[1]:.0f})" if data_paths
          else "data:    " + ("off (--data '')" if args.data == "" else "none found"))
    res = collect(paths, args, sig_paths, data_paths)
    if res is None:
        sys.exit("nothing to plot")
    mvas = sorted({m for _, m in res["used"]})
    print(f"discriminant: {', '.join(mvas)}")

    members, gincl, gslice, tincl, tslice = build_groups(res, args)
    report_samples(res["rows"], args)
    if len(members) > 1:
        report_composition(members, gslice, tslice, args)
    report_peak(members, gincl, gslice, tincl, tslice, args)
    if res["score_sig"]:
        report_score(res, tslice, args)
    if res["score_data"]:
        report_sidebands(res, args)
    if args.mult:
        report_nobdt(res["cls_incl"], res["hsamp"], args)
        report(res["table"], args)
    print()
    outdir = os.path.join(args.plotdir, args.cat)
    os.makedirs(outdir, exist_ok=True)
    _keeps = (draw_score(res, args, outdir)
              if res["score_sig"] or res["score_data"] else None)
    _keep = draw_processes(members, gincl, gslice, tincl, tslice, args,
                           outdir)
    if args.mult:
        _keep1 = draw_nobdt(res["cls_incl"], res["hsamp"], args, outdir)
        _keepd = draw_dominant(res["hmult"], args, outdir)
        _keep2 = draw(res["hists"], args, outdir, len(args.bdt_bins) - 1)


if __name__ == "__main__":
    main()
