"""
triggerEff.py -- signal trigger efficiency of HLT_IsoMu24 || HLT_Mu50, and what
adding the Mu17_Mu8 trigger would gain.

Reads the snapshots made with `Hmm.py --no-trigger` (written to
<indir>/<category>_notrig/), so the denominator is the full offline selection
of the category with no HLT requirement, and the numerators are the HLT bits.

    python triggerEff.py ggHcat 2024
    python triggerEff.py VLcat 2024 --samples 12,13,14 --region isH

Prints the integrated efficiency per sample and for all signal combined, a
non-overlapping breakdown (IsoMu24 only / both / Mu17_Mu8 only / neither),
plots the leading-muon pT per class, and writes efficiency-vs-variable histograms to <outdir> (ROOT file) and
<plotdir>/<category><year>/ (PNGs, public_html by default).
"""

import argparse
import getpass
import os
import sys

import ROOT

ROOT.gROOT.SetBatch(True)

DEFAULT_INDIR = f"/work/submit/{getpass.getuser()}/HmumuRun3/ROOTFILES/"
DEFAULT_OUTDIR = f"/work/submit/{getpass.getuser()}/HmumuRun3/TRIGEFF/"
DEFAULT_PLOTDIR = f"/home/submit/{getpass.getuser()}/public_html/HmumuRun3/triggerEff/"

SIGNAL = {"10": "VBF", "11": "ggH", "12": "WmH", "13": "WpH",
          "14": "ZH", "15": "ttH", "17": "bbH"}

MU17_MU8 = "HLT_Mu17_TrkIsoVVL_Mu8_TrkIsoVVL_DZ_Mass3p8"

TRIGGERS = {
    "IsoMu24":    "HLT_IsoMu24",
    "Mu50":       "HLT_Mu50",
    "OR":         "HLT_IsoMu24 || HLT_Mu50",
    "Mu50only":   "HLT_Mu50 && !HLT_IsoMu24",   # what Mu50 adds on top of IsoMu24
    "Mu17_Mu8":     MU17_MU8,
    "OR_Mu17_Mu8":  f"HLT_IsoMu24 || HLT_Mu50 || {MU17_MU8}",
    "Mu17_Mu8only": f"{MU17_MU8} && !(HLT_IsoMu24 || HLT_Mu50)",  # what Mu17_Mu8 adds on top of the OR
}

# curves drawn on the efficiency plots: (key in TRIGGERS, color, legend text)
CURVES = [
    ("IsoMu24", ROOT.kBlue + 1,    "IsoMu24"),
    ("OR",      ROOT.kRed + 1,     "IsoMu24 || Mu50"),
    ("Mu17_Mu8",    ROOT.kGreen + 2, "Mu17_Mu8"),
    ("OR_Mu17_Mu8", ROOT.kBlack,     "IsoMu24 || Mu50 || Mu17_Mu8"),
]

# Non-overlapping split of the denominator: IsoMu24 alone vs Mu17_Mu8 (no Mu50).
SINGLE = "HLT_IsoMu24"
CLASSES = {
    "IsoMu24only": f"{SINGLE} && !{MU17_MU8}",
    "both":        f"{SINGLE} && {MU17_MU8}",
    "Mu17_Mu8only": f"!{SINGLE} && {MU17_MU8}",
    "neither":     f"!{SINGLE} && !{MU17_MU8}",
}

# leading-muon pT distributions, overlaid: (selection, color, legend text)
DISTS = {
    "IsoMu24_or_Mu17_Mu8": (f"{SINGLE} || {MU17_MU8}", ROOT.kBlack, "IsoMu24 || Mu17_Mu8"),
    "IsoMu24only":  (CLASSES["IsoMu24only"],  ROOT.kRed + 1,   "IsoMu24 only"),
    "Mu17_Mu8only": (CLASSES["Mu17_Mu8only"], ROOT.kGreen + 2, "Mu17_Mu8 only"),
}
DIST_BINNING = (54, 15., 150.)   # 2.5 GeV bins

# variable -> (expression, nbins, xmin, xmax, axis label)
VARS = {
    "leadRawPt":  ("leadRawPt", 37, 15., 200., "leading muon p_{T} (raw) [GeV]"),
    "subRawPt":   ("subRawPt", 37, 15., 200., "subleading muon p_{T} (raw) [GeV]"),
    "mass":       ("HiggsCandCorrMass", 26, 70., 200., "m_{#mu#mu} [GeV]"),
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("category", help="e.g. ggHcat, VBFcat, VLcat, ...")
    p.add_argument("year", help="e.g. 2024, 12022, 22023")
    p.add_argument("-i", "--indir", default=DEFAULT_INDIR)
    p.add_argument("-o", "--outdir", default=DEFAULT_OUTDIR,
                   help="ROOT file with the histograms (default: %(default)s)")
    p.add_argument("-p", "--plotdir", default=DEFAULT_PLOTDIR,
                   help="PNGs go to <plotdir>/<category><year>/ (default: %(default)s)")
    p.add_argument("-s", "--samples", default=",".join(SIGNAL),
                   help="comma-separated MC ids (default: all H->mumu signal)")
    p.add_argument("--lead-pt-cut", type=float, default=0., metavar="GEV",
                   help="require leading-muon raw pT > GEV in the denominator. The "
                        "--no-trigger snapshots have no such cut; 26 reproduces the "
                        "nominal selection (default: %(default)s)")
    p.add_argument("--from-higgs", action="store_true",
                   help="keep only events where both selected muons are gen-matched "
                        "to the Higgs (pairFromHiggs)")
    p.add_argument("--region", default=None,
                   help="optional extra filter on the denominator, e.g. isH")
    p.add_argument("--weight", default="w_allSF",
                   help="event weight (default: %(default)s; no trigger SF applied)")
    p.add_argument("-j", "--ncores", type=int, default=0)
    return p.parse_args()


def load(indir, category, year, ids):
    chains = {}
    for i in ids:
        path = f"{indir}{category}_notrig/snapshot_mc_{i}_{year}_{category}.root"
        if not os.path.exists(path):
            print(f"  missing: {path}")
            continue
        c = ROOT.TChain("events")
        c.Add(path)
        chains[i] = c
    return chains


def define(df, region, lead_pt_cut, from_higgs):
    df = (df.Define("leadIs1", "Muon1_rawpt >= Muon2_rawpt")
            .Define("leadRawPt", "leadIs1 ? Muon1_rawpt : Muon2_rawpt")
            .Define("subRawPt", "leadIs1 ? Muon2_rawpt : Muon1_rawpt"))
    if from_higgs:
        df = df.Filter("pairFromHiggs", "both muons from H")
    if lead_pt_cut > 0:
        df = df.Filter(f"leadRawPt > {lead_pt_cut}", f"leadRawPt>{lead_pt_cut}")
    if region:
        df = df.Filter(region, region)
    return df


def book(df, weight, tag):
    """Integrated sums and binned histograms for one sample."""
    out = {"den": df.Sum(weight), "nraw": df.Count(), "num": {}, "h": {}}
    for t, sel in TRIGGERS.items():
        out["num"][t] = df.Filter(sel).Sum(weight)
    out["cls_n"] = {k: df.Filter(sel).Count() for k, sel in CLASSES.items()}
    out["cls_w"] = {k: df.Filter(sel).Sum(weight) for k, sel in CLASSES.items()}
    out["dist"] = {k: df.Filter(sel).Histo1D((f"leadRawPt_dist_{k}_{tag}", "", *DIST_BINNING),
                                             "leadRawPt", weight)
                   for k, (sel, _, _) in DISTS.items()}
    out["dist_n"] = {k: df.Filter(sel).Count() for k, (sel, _, _) in DISTS.items()}
    for v, (expr, nb, lo, hi, _) in VARS.items():
        model = (f"{v}_den_{tag}", "", nb, lo, hi)
        out["h"][(v, "den")] = df.Histo1D(model, expr, weight)
        for t, _, _ in CURVES:
            model = (f"{v}_{t}_{tag}", "", nb, lo, hi)
            out["h"][(v, t)] = df.Filter(TRIGGERS[t]).Histo1D(model, expr, weight)
    return out


def ratio(num, den, name):
    h = num.Clone(name)
    h.Divide(num, den, 1., 1., "B")
    return h


def main():
    args = parse_args()
    ROOT.ROOT.EnableImplicitMT(args.ncores)
    ids = [s.strip() for s in args.samples.split(",") if s.strip()]

    chains = load(args.indir, args.category, args.year, ids)
    if not chains:
        sys.exit(f"ERROR: no *_notrig snapshots under {args.indir}{args.category}_notrig/ "
                 f"-- run `python Hmm.py <year> <mode> --no-trigger` first")

    booked = {i: book(define(ROOT.RDataFrame(c), args.region, args.lead_pt_cut, args.from_higgs), args.weight, i)
              for i, c in chains.items()}
    ROOT.RDF.RunGraphs([r["den"] for r in booked.values()])

    # ---- integrated table --------------------------------------------
    # recovered = Mu17_Mu8only / (1 - OR): of the events IsoMu24||Mu50 misses,
    # the fraction Mu17_Mu8 brings back
    cols = list(TRIGGERS)

    def row(den, nums):
        effs = [nums[c] / den if den else float("nan") for c in cols]
        missed = den - nums["OR"]
        effs.append(nums["Mu17_Mu8only"] / missed if missed else float("nan"))
        return "".join(f"{e:>14.4f}" for e in effs)

    print(f"\nTrigger efficiency  category={args.category}  year={args.year}"
          f"  region={args.region or 'all'}  leadRawPt>{args.lead_pt_cut:g}"
          f"  fromHiggs={args.from_higgs}  weight={args.weight}")
    print(f"{'sample':<8}{'Nraw':>9}" + "".join(f"{c:>14}" for c in cols + ["recovered"]))
    tot_den, tot_num = 0., {c: 0. for c in cols}
    for i, r in booked.items():
        den = r["den"].GetValue()
        nums = {c: r["num"][c].GetValue() for c in cols}
        tot_den += den
        for c in cols:
            tot_num[c] += nums[c]
        print(f"{SIGNAL.get(i, i):<8}{r['nraw'].GetValue():>9}" + row(den, nums))
    print(f"{'all':<8}{'':>9}" + row(tot_den, tot_num))

    # ---- non-overlapping breakdown, all samples summed ----------------
    n_tot = sum(r["nraw"].GetValue() for r in booked.values())
    cls_n = {k: sum(r["cls_n"][k].GetValue() for r in booked.values()) for k in CLASSES}
    cls_w = {k: sum(r["cls_w"][k].GetValue() for r in booked.values()) for k in CLASSES}
    w_any = tot_den - cls_w["neither"]
    print(f"\nBreakdown (all samples summed)   IsoMu24 vs Mu17_Mu8, Mu50 not used")
    print(f"{'class':<14}{'Nraw':>10}{'weighted':>14}{'frac of all':>14}{'frac of any':>14}")
    for k in CLASSES:
        f_all = cls_w[k] / tot_den if tot_den else float("nan")
        f_any = cls_w[k] / w_any if (w_any and k != "neither") else float("nan")
        print(f"{k:<14}{cls_n[k]:>10}{cls_w[k]:>14.4g}{f_all:>14.4f}{f_any:>14.4f}")
    print(f"{'total':<14}{n_tot:>10}{tot_den:>14.4g}{1.:>14.4f}")
    missed = cls_w["Mu17_Mu8only"] + cls_w["neither"]
    if missed:
        print(f"Mu17_Mu8 recovers {cls_w['Mu17_Mu8only'] / missed:.4f} of the events IsoMu24 misses "
              f"({cls_n['Mu17_Mu8only']} of {cls_n['Mu17_Mu8only'] + cls_n['neither']} raw)")

    # ---- binned efficiencies, all samples summed ---------------------
    os.makedirs(args.outdir, exist_ok=True)
    plotdir = os.path.join(args.plotdir, f"{args.category}{args.year}")
    os.makedirs(plotdir, exist_ok=True)
    suffix = f"{args.category}_{args.year}" + (f"_{args.region}" if args.region else "")
    if args.lead_pt_cut > 0:
        suffix += f"_lead{args.lead_pt_cut:g}"
    if args.from_higgs:
        suffix += "_fromH"
    fout = ROOT.TFile(os.path.join(args.outdir, f"trigEff_{suffix}.root"), "RECREATE")
    ROOT.gStyle.SetOptStat(0)
    for v, (_, _, _, _, label) in VARS.items():
        summed = {}
        for key in ["den"] + [t for t, _, _ in CURVES]:
            hs = [r["h"][(v, key)].GetValue() for r in booked.values()]
            h = hs[0].Clone(f"{v}_{key}")
            for other in hs[1:]:
                h.Add(other)
            summed[key] = h
            h.Write()

        c = ROOT.TCanvas(f"c_{v}", "", 700, 600)
        frame = c.DrawFrame(summed["den"].GetXaxis().GetXmin(), 0.5,
                            summed["den"].GetXaxis().GetXmax(), 1.05)
        frame.GetXaxis().SetTitle(label)
        frame.GetYaxis().SetTitle("trigger efficiency")
        keep = []
        leg = ROOT.TLegend(0.45, 0.15, 0.88, 0.38)
        leg.SetBorderSize(0)
        for key, color, text in CURVES:
            e = ratio(summed[key], summed["den"], f"eff_{v}_{key}")
            e.SetLineColor(color)
            e.SetMarkerColor(color)
            e.SetMarkerStyle(20)
            e.Draw("E1 same")
            e.Write()
            keep.append(e)
            leg.AddEntry(e, text, "lp")
        leg.Draw()
        ROOT.SetOwnership(leg, False)
        c.SaveAs(os.path.join(plotdir, f"trigEff_{v}_{suffix}.png"))
    # ---- leading-pT distributions per trigger class ------------------
    c = ROOT.TCanvas("c_dist", "", 700, 600)
    c.SetLogy()
    keep = []
    leg = ROOT.TLegend(0.38, 0.72, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetTextSize(0.03)
    for k, (_, color, text) in DISTS.items():
        hs = [r["dist"][k].GetValue() for r in booked.values()]
        h = hs[0].Clone(f"leadRawPt_dist_{k}")
        for other in hs[1:]:
            h.Add(other)
        h.SetLineColor(color)
        h.SetLineWidth(2)
        h.GetXaxis().SetTitle("leading muon p_{T} (raw) [GeV]")
        h.GetYaxis().SetTitle("weighted events / 2.5 GeV")
        h.Write()
        keep.append(h)
        # N = raw events; fraction = weighted share of all selected events (incl. under/overflow)
        n_raw = sum(r["dist_n"][k].GetValue() for r in booked.values())
        frac = h.Integral(0, h.GetNbinsX() + 1) / tot_den if tot_den else float("nan")
        leg.AddEntry(h, f"{text}: N = {n_raw}, {100 * frac:.2f}%", "l")
    ymax = max(h.GetMaximum() for h in keep)
    positive = [h.GetBinContent(b) for h in keep for b in range(1, h.GetNbinsX() + 1)
                if h.GetBinContent(b) > 0]
    ymin = min(positive, default=ymax * 1e-4)
    for n, h in enumerate(keep):
        h.SetMinimum(ymin * 0.5)
        h.SetMaximum(ymax * 5)
        h.Draw("hist" if n == 0 else "hist same")
    leg.Draw()
    ROOT.SetOwnership(leg, False)
    c.SaveAs(os.path.join(plotdir, f"leadRawPt_byTrigger_{suffix}.png"))

    fout.Close()
    print(f"\nwrote {os.path.join(args.outdir, f'trigEff_{suffix}.root')}")
    print(f"PNGs in {plotdir}")


if __name__ == "__main__":
    main()
