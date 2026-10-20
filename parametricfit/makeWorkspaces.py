"""
makeWorkspaces.py -- build the Combine input workspaces from the fit JSONs.

Step two of two. The fits run in conda (ROOT >= 6.28): sigFit.py and
bkgFit.py write every fitted number to JSON, and bkgFit.py also writes the
data as a plain TH1D. This script runs inside the Combine container, rebuilds
the PDFs from those numbers with the classes the container has, and writes
the workspaces bwsHrare.py reads, with tag = "<cat>_<bin>_<year>":

    Signal_<tag>_workspace.root   crystal_ball_<tag>_<mode>  (+ _norm)
    Bkg_<tag>_workspace.root      multipdf_<tag>_bkg, pdfindex_<tag>,
                                  multipdf_<tag>_bkg_norm, datahist_<tag>

Why the split: the container's ROOT (6.22) has no RooCrystalBall, and RooFit
objects written by a newer ROOT are not guaranteed to read back in it; and
before 6.28, Range("sbLow,sbHigh") normalises each sideband separately,
which breaks the blinded background fit. Numbers and a TH1D cross between
versions safely, and nothing is refitted here.

The signal Signal_<tag>_workspace.root written by sigFit.py in conda is
OVERWRITTEN with the container-readable version.

Usage (inside the Combine container, from parametricfit/)
-----
    python makeWorkspaces.py -c ggHcat -b incl -y 2024 --wsdir WS_2024
    python makeWorkspaces.py -c ggHcat -b incl -y 2024 --only bkg
"""

import argparse
import json
import os

import ROOT

import pdfDefinitions as pdefs

ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)


def read_json(path):
    if not os.path.isfile(path):
        raise SystemExit(f"missing {path}; run the conda fit step first")
    with open(path) as fh:
        return json.load(fh)


def observable(category, window):
    """The mass observable, named as sigFit.py, bkgFit.py and bwsHrare.py
    expect."""
    x = ROOT.RooRealVar(f"mh{category}", "m_{#mu#mu}", *window)
    x.setRange("full", *window)
    return x


def pin(var, value, const=True):
    """Set a parameter to its fitted value. setVal clamps silently to the
    range, and the container's bounds can be tighter (the RooDoubleCBFast
    n floor), so widen the range rather than move the value."""
    lo, hi = var.getMin(), var.getMax()
    if not lo <= value <= hi:
        print(f"  note: {var.GetName()} = {value:.4g} outside "
              f"[{lo:g}, {hi:g}]; widening the range")
        var.setRange(min(lo, value), max(hi, value))
    var.setVal(value)
    var.setConstant(const)


def build_signal(category, tag, sigjson, wsfile):
    """One frozen DSCB per production mode, with the shape nuisances left
    floating for Combine -- the same content sigFit.py writes."""
    rec = read_json(sigjson)
    x = observable(category, rec["provenance"]["fit_range"])
    w = ROOT.RooWorkspace("w", "workspace")
    keep = []
    for mode, r in sorted(rec["fits"].items()):
        if "parameters" not in r:
            print(f"  {mode}: {r.get('skipped', 'no fit')}, skipping")
            continue
        cfg = r.get("shape_model", {})
        pdf, bundle = pdefs.create_signal_pdf(
            x, f"{tag}_{mode}", scale_unc=cfg.get("scale_unc"),
            res_unc=cfg.get("res_unc"),
            correlated=cfg.get("correlated", True))
        for k, v in bundle["nom"].items():
            pin(v, r["parameters"][k]["value"])
        for v in bundle["nuis"].values():
            v.setVal(0.0)
            v.setConstant(False)
        norm = ROOT.RooRealVar(f"{pdf.GetName()}_norm",
                               f"{pdf.GetName()}_norm",
                               r["normalization"]["integral"])
        norm.setConstant(True)
        getattr(w, "import")(pdf, ROOT.RooFit.Silence())
        getattr(w, "import")(norm, ROOT.RooFit.Silence())
        keep.append((pdf, bundle, norm))
        print(f"  {pdf.GetName():<40} norm {norm.getVal():10.3f}  "
              f"({bundle['cfg']['pdf_class']})")
        if not r.get("fit_quality", {}).get("ok", True):
            print("    !! sigFit flagged this fit; check it before use")
    if not keep:
        raise SystemExit(f"no fitted signal modes in {sigjson}")
    w.writeToFile(wsfile)
    print(f"wrote {wsfile}")


def build_bkg(category, tag, bkgjson, wsdir, wsfile):
    """The RooMultiPdf envelope, each member at its sideband-fit values
    with its parameters floating, plus the floating norm and the data."""
    rec = read_json(bkgjson)
    prov = rec["provenance"]
    x = observable(category, prov["window"])

    keep, store = [], ROOT.RooArgList()
    for nm in rec["envelope"]:
        pdf, k = pdefs.create_bkg_pdf(x, tag, nm)
        fitted = rec["candidates"][nm]["params"]
        for p in k["params"]:
            if p.GetName() not in fitted:
                raise SystemExit(f"{bkgjson}: no fitted value for "
                                 f"{p.GetName()} ({nm})")
            pin(p, fitted[p.GetName()]["value"], const=False)
        store.add(pdf)
        keep.append((pdf, k))

    if not pdefs.load_combine():
        raise SystemExit("needs libHiggsAnalysisCombinedLimit: run inside "
                         "the Combine container")
    cat = ROOT.RooCategory(f"pdfindex_{tag}", f"pdfindex_{tag}")
    multi = ROOT.RooMultiPdf(f"multipdf_{tag}_bkg", "multipdf", cat, store)

    nrm_cfg = rec["normalization"]
    norm_val = nrm_cfg["value"]
    lo, hi = nrm_cfg["norm_range"]
    nrm = ROOT.RooRealVar(f"multipdf_{tag}_bkg_norm",
                          f"multipdf_{tag}_bkg_norm", norm_val,
                          lo * norm_val, hi * norm_val)
    nrm.setConstant(False)              # Combine profiles it

    datafile = os.path.join(wsdir, prov["data_file"])
    fin = ROOT.TFile.Open(datafile)
    if not fin or fin.IsZombie():
        raise SystemExit(f"cannot open {datafile}")
    h = fin.Get(f"data_{tag}")
    if not h:
        raise SystemExit(f"{datafile} has no data_{tag}")
    h.SetDirectory(0)
    fin.Close()
    data = ROOT.RooDataHist(f"datahist_{tag}", "data", ROOT.RooArgList(x), h)

    w = ROOT.RooWorkspace("w", "w")
    for obj in (cat, multi, nrm, data):  # data is FULL range, for data_obs
        getattr(w, "import")(obj, ROOT.RooFit.Silence())
    w.writeToFile(wsfile)
    print(f"  multipdf_{tag}_bkg        {', '.join(rec['envelope'])}")
    print(f"  multipdf_{tag}_bkg_norm   = {norm_val:.1f}, floating")
    print(f"  datahist_{tag}            {data.sumEntries():.0f} events")
    print(f"wrote {wsfile}")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--category", default="ggHcat")
    ap.add_argument("-b", "--bin", default="incl",
                    help="BDT bin label, as given to sigFit/bkgFit")
    ap.add_argument("-y", "--year", default="Run3")
    ap.add_argument("--wsdir", default="WS_LOCAL",
                    help="where the JSONs are, and the workspaces go")
    ap.add_argument("--only", choices=("signal", "bkg"), default=None)
    args = ap.parse_args()

    tag = f"{args.category}_{args.bin}_{args.year}"
    print(f"tag      : {tag}\nROOT     : {ROOT.gROOT.GetVersion()}")

    if args.only != "bkg":
        print("\n--- signal ---")
        build_signal(args.category, tag,
                     os.path.join(args.wsdir, f"Signal_{tag}_params.json"),
                     os.path.join(args.wsdir, f"Signal_{tag}_workspace.root"))
    if args.only != "signal":
        print("\n--- background ---")
        build_bkg(args.category, tag,
                  os.path.join(args.wsdir, f"Bkg_{tag}_params.json"),
                  args.wsdir,
                  os.path.join(args.wsdir, f"Bkg_{tag}_workspace.root"))


if __name__ == "__main__":
    main()
