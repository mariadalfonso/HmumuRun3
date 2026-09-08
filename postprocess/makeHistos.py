"""
makeHistos.py -- histogramming stage for the HmumuRun3 analysis.

Runs one RDataFrame event loop over the snapshot chain and writes every
histogram the plotting scripts need into a single ROOT file:

    histos_<category>_<year>.root
        <region>/<var>_<process>

Usage:
    python makeHistos.py ggHcat 2024
    python makeHistos.py ggHcat 2024 --unblind
    python makeHistos.py ggHcat 2024 -i <snapshotdir> -o <histodir>
"""

import ROOT
import os
import sys
import argparse
import getpass
from datetime import datetime

from utils.LoadTree import loadTree
from utils import plot_vars
from utils import prepareHisto
from utils import histo_config as cfg
from utils.plot_style import lumis

ROOT.gROOT.SetBatch()
ROOT.ROOT.EnableImplicitMT()

DEFAULT_INDIR = f"/work/submit/{getpass.getuser()}/HmumuRun3/ROOTFILES/"
DEFAULT_OUTDIR = f"/work/submit/{getpass.getuser()}/HmumuRun3/HISTOS/"


def parse_args():
    p = argparse.ArgumentParser(
        description="Book and fill all histograms for one (category, year) "
                    "in a single event loop."
    )
    p.add_argument("category", choices=cfg.CATEGORIES,
                   help="analysis category (snapshot subfolder)")
    p.add_argument("year",
                   help="data-taking year, e.g. 2024, 2025, 12022 ...")
    p.add_argument("-i", "--indir", default=DEFAULT_INDIR,
                   help="input directory with snapshot ROOT files "
                        "(default: %(default)s)")
    p.add_argument("-o", "--outdir", default=DEFAULT_OUTDIR,
                   help="output directory for the histogram file "
                        "(default: %(default)s)")
    p.add_argument("--blind", dest="blind", action="store_true", default=True,
                   help="blind data in the signal region (default: on; "
                        "applies to blindable regions only)")
    p.add_argument("--unblind", dest="blind", action="store_false",
                   help="disable blinding (use with care)")
    p.add_argument("--regions", nargs="+", default=None,
                   choices=list(cfg.REGIONS),
                   help="only book these regions (default: all)")
    p.add_argument("--vars", nargs="+", default=None,
                   help="only book these variables (default: the category's "
                        "active set from histo_config.get_active_vars)")
    return p.parse_args()


def main():
    args = parse_args()

    category = args.category
    year = "_" + args.year          # internal convention: leading underscore

    indir = args.indir if args.indir.endswith("/") else args.indir + "/"
    outdir = args.outdir if args.outdir.endswith("/") else args.outdir + "/"

    # Guard: year must be known to the lumi table (used later by the CMS label).
    if year not in lumis:
        sys.exit(f"ERROR: year '{args.year}' not found in lumis table "
                 f"(known: {[k.lstrip('_') for k in lumis if k.startswith('_')]})")

    t0 = datetime.now()
    print(f"[makeHistos] category={category} year={args.year}")
    print(f"[makeHistos] threads: {ROOT.GetThreadPoolSize()}")
    print(f"[makeHistos] blinding: "
          f"{'ON' if args.blind else 'OFF (UNBLINDED)'}")

    # ---- input chain --------------------------------------------------
    chain = loadTree(ROOT.TChain("events"), indir, category, year)
    nfiles = chain.GetListOfFiles().GetEntries() if chain.GetListOfFiles() else 0
    print(f"[makeHistos] chained {nfiles} snapshot files")
    if nfiles == 0:
        sys.exit("ERROR: no input files -- check --indir, category and year")

    # ---- base node: built ONCE ----------------------------------------
    base = prepareHisto.make_base_node(chain, year)

    # ---- BOOK: region -> variable -> process --------------------------
    # Mirrors FastFrames' nesting (region outside variable, variable outside
    # process). NOTHING here is dereferenced.
    regions = args.regions or list(cfg.REGIONS)
    active_vars = args.vars or cfg.get_active_vars(category)

    booked = {}     # (region, varname, process) -> RResultPtr[TH1D]
    skipped = []

    for region_name in regions:
        region = cfg.REGIONS[region_name]

        node = base
        if region["filter"] is not None:
            node = node.Filter(region["filter"], f"region:{region_name}")

        # Blinding applies only where the region doesn't already exclude the
        # signal core (see histo_config.REGIONS).
        blind = args.blind and region["blindable"]

        for varname in active_vars:
            binning = plot_vars.get_binning(varname)
            if binning is None:
                skipped.append((region_name, varname, "no binning defined"))
                continue

            ptrs = prepareHisto.book_variable(
                node, category, varname, binning, blind=blind
            )
            if ptrs is None:
                skipped.append((region_name, varname, "branch not available"))
                continue

            for proc, ptr in ptrs.items():
                booked[(region_name, varname, proc)] = ptr

    if not booked:
        sys.exit("ERROR: nothing booked -- check --vars / snapshot contents")

    print(f"[makeHistos] booked {len(booked)} histograms "
          f"({len(regions)} regions x {len(active_vars)} vars "
          f"x {len(cfg.ALL_PROCESSES)} processes)")
    for region_name, varname, why in skipped:
        print(f"   -> skipped {region_name}/{varname}: {why}")

    # ---- TRIGGER: one event loop for everything -----------------------
    print("[makeHistos] triggering the event loop ...")
    ROOT.RDF.RunGraphs(list(booked.values()))
    t1 = datetime.now()
    print(f"[makeHistos] event loop done in {t1 - t0}")

    # ---- finalize + write ---------------------------------------------
    # Overflow folding and negative-bin removal happen HERE, after the loop.
    # In the old code these lived inside getHisto() and were part of what
    # forced a separate loop per plot.
    os.makedirs(outdir, exist_ok=True)
    outpath = outdir + cfg.histo_filename(category, year)
    fout = ROOT.TFile(outpath, "RECREATE")

    # Provenance, so a plot can always be traced back to how it was made.
    meta = ROOT.TNamed(
        "provenance",
        f"category={category} year={year} blind={args.blind} "
        f"indir={indir} nfiles={nfiles} created={t0:%Y-%m-%d %H:%M:%S}"
    )
    fout.WriteTObject(meta, "provenance")

    written = 0
    for region_name in regions:
        d = fout.mkdir(region_name)
        for (reg, varname, proc), ptr in booked.items():
            if reg != region_name:
                continue
            hist = ptr.GetValue()       # already computed; no loop here
            prepareHisto.finalize(hist)
            d.WriteTObject(hist, f"{varname}_{proc}")
            written += 1

    fout.Close()

    print(f"[makeHistos] wrote {written} histograms -> {outpath}")
    print(f"[makeHistos] total {datetime.now() - t0}")


if __name__ == "__main__":
    main()
