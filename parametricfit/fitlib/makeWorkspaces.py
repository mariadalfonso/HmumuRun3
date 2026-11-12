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

With --members families the multipdf holds every family's F-test pick
(bkgFit.py's family_choice), the ones the goodness-of-fit cut rejected
included, and goes to Bkg_<tag>_families_workspace.root instead. That is
the input of biasStudy.py --mode single (the Run 2 note's truth x fit
matrix).

With --members core the multipdf is the core-PDF of corePdfFit.py
(Bkg_<tag>_core_params.json): core x this bin's SMF per core, the core
parameters and the index pdfindex_<cat>_core_<year> named for the
category, so the bins of the combined card share them. Written to
Bkg_<tag>_core_workspace.root.

Why the split: the container's ROOT (6.22) has no RooCrystalBall, and RooFit
objects written by a newer ROOT are not guaranteed to read back in it; and
before 6.28, Range("sbLow,sbHigh") normalises each sideband separately,
which breaks the blinded background fit. Numbers and a TH1D cross between
versions safely, and nothing is refitted here.

The signal Signal_<tag>_workspace.root written by sigFit.py in conda is
OVERWRITTEN with the container-readable version.

Usage (inside the Combine container, from parametricfit/)
-----
    python -m fitlib.makeWorkspaces -c ggHcat -b incl -y 2024 --wsdir WS_2024
    python -m fitlib.makeWorkspaces -c ggHcat -b incl -y 2024 --only bkg
    python -m fitlib.makeWorkspaces -c ggHcat -b bdt0 -y Run3 --members families

buildWorkspace.py runs it for every bin; call it directly only to debug.
"""

import argparse
import json
import os

import ROOT

from . import config, corePdf
from . import pdfDefinitions as pdefs

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
    range, so if a value sits outside the rebuilt bounds (they should
    match the fit's), widen the range rather than move the value."""
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


# file-name suffix of each --members choice
MEMBERS_SUFFIX = {"envelope": "", "families": "_families", "core": "_core"}


def bkg_members(rec, members):
    """The multipdf members, in pdfindex order. "envelope": bkgFit.py's
    envelope, best first. "families": the envelope, then every other
    family's F-test pick (rejected by the goodness-of-fit cut), lowest
    penalised NLL first; index 0 stays the envelope's best, which the
    norm comes from."""
    env = list(rec["envelope"])
    if members == "envelope":
        return env
    cand = rec["candidates"]
    rest = sorted({nm for nm in rec["family_choice"].values()
                   if nm not in env},
                  key=lambda n: cand[n].get("nll_pen", cand[n]["nll"]))
    return env + rest


def build_bkg(category, tag, bkgjson, wsdir, wsfile, members="envelope"):
    """The RooMultiPdf, each member at its sideband-fit values with its
    parameters floating, plus the floating norm and the data."""
    rec = read_json(bkgjson)
    prov = rec["provenance"]
    x = observable(category, prov["window"])
    names = bkg_members(rec, members)

    mc_core = (rec.get("mc_core") or {}).get("cores")   # mc* members
    keep, store = [], ROOT.RooArgList()
    for nm in names:
        pdf, k = pdefs.create_bkg_pdf(x, tag, nm, mc_core)
        pin_all(k["params"], rec["candidates"][nm]["params"], bkgjson, nm)
        store.add(pdf)
        keep.append((pdf, k))

    cat = ROOT.RooCategory(f"pdfindex_{tag}", f"pdfindex_{tag}")
    write_bkg(tag, x, cat, store, rec, wsdir, wsfile, names,
              "" if members == "envelope" else
              f"   (all families; envelope = {', '.join(rec['envelope'])})")


def build_core_bkg(category, tag, corejson, wsdir, wsfile):
    """The core-PDF of corePdfFit.py: one member per core, core x this
    bin's SMF, at the fitted values with the parameters floating. The core
    parameters and the index carry the category's name, not the bin's, so
    text2workspace merges them across the bins of the combined card."""
    rec = read_json(corejson)
    prov = rec["provenance"]
    x = observable(category, prov["window"])
    year = prov["year"]
    keep, store = [], ROOT.RooArgList()
    for core in rec["members"]:
        core_pdf, ck = corePdf.create_core(x, category, year, core)
        pin_all(ck["params"], rec["core_params"][core], corejson, core)
        pdf, smf, k = corePdf.create_member(x, core_pdf, core, tag,
                                            prov["smf_order"])
        pin_all(smf, rec["smf_params"][core], corejson, f"{core} SMF")
        store.add(pdf)
        keep.append((core_pdf, ck, pdf, smf, k))
    index = prov["index"]
    cat = ROOT.RooCategory(index, index)
    write_bkg(tag, x, cat, store, rec, wsdir, wsfile, rec["members"],
              f"   (core-PDF: cores shared as {prov['core_tag']}, index "
              f"{index}, SMF order {prov['smf_order']})")


def pin_all(pars, fitted, source, what):
    """Every floating parameter at its fitted value (still floating)."""
    for p in pars:
        if p.GetName() not in fitted:
            raise SystemExit(f"{source}: no fitted value for "
                             f"{p.GetName()} ({what})")
        pin(p, fitted[p.GetName()]["value"], const=False)


def write_bkg(tag, x, cat, store, rec, wsdir, wsfile, names, note):
    """multipdf_<tag>_bkg over `store`, its floating norm and the data."""
    prov = rec["provenance"]
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
    print(f"  multipdf_{tag}_bkg        {', '.join(names)}{note}")
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
    ap.add_argument("--wsdir", default=None,
                    help="where the JSONs are, and the workspaces go "
                         "(default WS_<year>)")
    ap.add_argument("--only", choices=("signal", "bkg"), default=None)
    ap.add_argument("--members", choices=tuple(MEMBERS_SUFFIX),
                    default="envelope",
                    help="multipdf members: the envelope (default), every "
                         "family's F-test pick (families: "
                         "Bkg_<tag>_families_workspace.root, for "
                         "biasStudy.py --mode single), or the core-PDF of "
                         "corePdfFit.py (core: Bkg_<tag>_core_workspace.root)")
    args = ap.parse_args()
    args.wsdir = args.wsdir or config.default_wsdir(args.year)

    # before writing anything: the workspaces must be built with Combine's classes
    if not pdefs.load_combine(quiet=True):
        raise SystemExit("needs libHiggsAnalysisCombinedLimit: run inside "
                         "the Combine container")

    tag = f"{args.category}_{args.bin}_{args.year}"
    print(f"tag      : {tag}\nROOT     : {ROOT.gROOT.GetVersion()}")

    if args.only != "bkg":
        print("\n--- signal ---")
        build_signal(args.category, tag,
                     os.path.join(args.wsdir, f"Signal_{tag}_params.json"),
                     os.path.join(args.wsdir, f"Signal_{tag}_workspace.root"))
    if args.only != "signal":
        print("\n--- background ---")
        sfx = MEMBERS_SUFFIX[args.members]
        wsfile = os.path.join(args.wsdir, f"Bkg_{tag}{sfx}_workspace.root")
        if args.members == "core":
            build_core_bkg(args.category, tag, os.path.join(
                args.wsdir, f"Bkg_{tag}_core_params.json"), args.wsdir,
                wsfile)
        else:
            build_bkg(args.category, tag,
                      os.path.join(args.wsdir, f"Bkg_{tag}_params.json"),
                      args.wsdir, wsfile, args.members)


if __name__ == "__main__":
    main()
