"""
triggerEff.py -- signal trigger efficiency of HLT_IsoMu24 || HLT_Mu50, and what
adding the Mu17_Mu8 trigger would gain.

Reads the snapshots made with `Hmm.py --no-trigger` (written to
<indir>/<category>_notrig/), so the denominator is the full offline selection
of the category with no HLT requirement, and the numerators are the HLT bits.

    python triggerEff.py ggHcat 2024
    python triggerEff.py VLcat 2024 --samples 12,13,14 --region isH

Prints the integrated efficiency per sample and for all signal combined, and
writes efficiency-vs-variable histograms to <outdir> (ROOT file) and
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

# variable -> (expression, nbins, xmin, xmax, axis label)
VARS = {
    "leadRawPt":  ("leadRawPt", 37, 15., 200., "leading muon p_{T} (raw) [GeV]"),
    "leadAbsEta": ("leadAbsEta", 12, 0., 2.4, "leading muon |#eta|"),
    "subRawPt":   ("subRawPt", 37, 15., 200., "subleading muon p_{T} (raw) [GeV]"),
    "nPV":        ("PV_npvsGood", 20, 0., 80., "N_{PV}"),
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
            .Define("subRawPt", "leadIs1 ? Muon2_rawpt : Muon1_rawpt")
            .Define("leadAbsEta", "std::abs(leadIs1 ? Muon1_eta : Muon2_eta)"))
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
    fout.Close()
    print(f"\nwrote {os.path.join(args.outdir, f'trigEff_{suffix}.root')}")
    print(f"PNGs in {plotdir}")


if __name__ == "__main__":
    main()
