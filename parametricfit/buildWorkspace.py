"""
buildWorkspace.py -- from the fit JSONs to the Combine workspace, in one go.

    python buildWorkspace.py -c ggHcat -b bdt0 bdt1 bdt2 -y Run3 --wsdir WS_Run3

Steps 3-5 of the chain (README), inside the Combine container, after sigFit.py
and bkgFit.py have run in conda:

  1. makeWorkspaces.py for each bin
  2. bwsHrare.py for each bin: <wsdir>/datacard_<tag>.txt and
     <wsdir>/workspace_<tag>_workspace.root
  3. combineCards.py over the bins (more than one bin only):
     <wsdir>/datacard_<name>.txt
  4. text2workspace.py: <wsdir>/<name>.root

then writes the yields of every bin to <wsdir>/Yields_<name>.json, prints them,
and prints the combine commands for the expected limit and significance.

<name> is <category>_<bins>_<year>, e.g. ggHcat_bdt_Run3 for bdt0 bdt1 bdt2
(--name overrides); for one bin it is that bin's tag. Everything runs from
parametricfit/: the shapes paths in the cards are relative to it, so run
combine from here too.
"""

import argparse
import json
import os
import re
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))

MH = "125"
# discrete profiling, as in biasStudy.py and the Combine documentation
COMBINE_OPTS = [
    "--cminDefaultMinimizerStrategy", "0",
    "--X-rtd", "MINIMIZER_freezeDisassociatedParams",
    "--X-rtd", "MINIMIZER_multiMin_hideConstants",
    "--X-rtd", "MINIMIZER_multiMin_maskConstraints",
    "--X-rtd", "MINIMIZER_multiMin_maskChannels=2",
]


def in_container():
    return bool(os.environ.get("APPTAINER_CONTAINER")
                or os.environ.get("SINGULARITY_CONTAINER"))


def chain_name(args):
    """ggHcat_bdt_Run3 for bdt0 bdt1 bdt2, ggHcat_incl_2024 for one bin."""
    if args.name:
        return args.name
    if len(args.bins) == 1:
        label = args.bins[0]
    elif all(re.fullmatch(r"bdt\d+", b) for b in args.bins):
        label = "bdt"
    else:
        label = "_".join(args.bins)
    return f"{args.category}_{label}_{args.year}"


def tags(args):
    return [f"{args.category}_{b}_{args.year}" for b in args.bins]


def rel(path):
    """path relative to parametricfit/, as the cards and combine see it"""
    return os.path.relpath(path, HERE)


def run(cmd, args, cwd=HERE, stdout=None):
    print(f"\n$ (cd {os.path.relpath(cwd, HERE)}) {' '.join(cmd)}", flush=True)
    if args.dry_run:
        return
    subprocess.run(cmd, cwd=cwd, stdout=stdout, check=True)


# ---------------------------------------------------------------------------
# steps
# ---------------------------------------------------------------------------

def make_workspaces(args):
    for b in args.bins:
        run([sys.executable, "makeWorkspaces.py", "-c", args.category,
             "-b", b, "-y", args.year, "--wsdir", args.wsdir], args)


def make_cards(args):
    ws = rel(args.wsdir)
    for b, tag in zip(args.bins, tags(args)):
        # --output goes into the shapes lines as given: relative to here
        run([sys.executable, "bwsHrare.py",
             f"--whichCat={args.category}", f"--whichBIN={b}",
             f"--whichYear={args.year}",
             f"--inputFileSig={ws}/Signal_{tag}_workspace.root",
             f"--inputFileBKG={ws}/Bkg_{tag}_workspace.root",
             f"--output={ws}/workspace_{tag}_workspace.root",
             f"--datCardName={ws}/datacard_{tag}.txt"], args)


def combine_cards(args):
    """combineCards.py prefixes each shapes path with its card's directory,
    so it runs next to the cards (bare file names, no prefix)."""
    card = os.path.join(args.wsdir, f"datacard_{chain_name(args)}.txt")
    cmd = ["combineCards.py"] + [f"{b}=datacard_{t}.txt"
                                 for b, t in zip(args.bins, tags(args))]
    print(f"\n   > {rel(card)}")
    if args.dry_run:
        run(cmd, args, cwd=args.wsdir)
        return
    with open(card, "w") as fh:
        run(cmd, args, cwd=args.wsdir, stdout=fh)
    # guard: a doubled prefix if combineCards ever adds the card directory
    ws = rel(args.wsdir)
    with open(card) as fh:
        text = fh.read()
    if f"{ws}/{ws}/" in text:
        with open(card, "w") as fh:
            fh.write(text.replace(f"{ws}/{ws}/", f"{ws}/"))


def text2workspace(args):
    name, ws = chain_name(args), rel(args.wsdir)
    run(["text2workspace.py", "-m", MH, f"{ws}/datacard_{name}.txt",
         "-o", f"{ws}/{name}.root"], args)


# ---------------------------------------------------------------------------
# yields
# ---------------------------------------------------------------------------

def read_json(path):
    with open(path) as fh:
        return json.load(fh)


def signal_yields(rec):
    """{mode: yield} and the window, from a Signal_<tag>_params.json. JSONs
    written before sigFit.py saved 'yields' are summed from the fits."""
    if "yields" in rec:
        return rec["yields"]["by_mode"], rec["yields"]["window"]
    by_mode = {s: r["normalization"]["integral"]
               for s, r in rec["fits"].items() if "normalization" in r}
    return by_mode, None


def write_yields(args):
    name = chain_name(args)
    out = {"name": name, "category": args.category, "year": args.year,
           "bins": {}}
    tot_sig, tot_bkg = {}, 0.0
    for b, tag in zip(args.bins, tags(args)):
        sj = os.path.join(args.wsdir, f"Signal_{tag}_params.json")
        bj = os.path.join(args.wsdir, f"Bkg_{tag}_params.json")
        entry = {}
        if os.path.isfile(sj):
            srec = read_json(sj)
            by_mode, window = signal_yields(srec)
            out.setdefault("lumi_fb", sum(srec["provenance"]["lumi"].values()))
            entry["signal"] = {"by_mode": by_mode,
                               "total": sum(by_mode.values())}
            if window:
                out["window"] = window
            for s, v in by_mode.items():
                tot_sig[s] = tot_sig.get(s, 0.0) + v
        if os.path.isfile(bj):
            brec = read_json(bj)
            norm = brec["normalization"]
            entry["background"] = {
                # blinded: the sideband count extrapolated over the window
                # with the best envelope member
                "total": norm["value"],
                "sideband_count": norm["sideband_count"],
                "best": brec["best"], "envelope": brec["envelope"],
                "envelope_spread": norm.get("envelope_spread"),
            }
            out["window"] = brec["provenance"]["window"]
            tot_bkg += norm["value"]
        entry["inputs"] = [os.path.basename(p) for p in (sj, bj)
                           if os.path.isfile(p)]
        out["bins"][b] = entry
    out["total"] = {"signal": {"by_mode": tot_sig,
                               "total": sum(tot_sig.values())},
                    "background": tot_bkg}

    print_yields(out)
    path = os.path.join(args.wsdir, f"Yields_{name}.json")
    if args.dry_run:
        print(f"\n(dry run, {path} not written)")
        return
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2, sort_keys=True)
    print(f"\nwrote {path}")


def print_yields(out):
    w = out.get("window", ["?", "?"])
    modes = list(out["total"]["signal"]["by_mode"])
    head = ["bin"] + modes + ["signal", "background"]
    rows = []
    for b, e in list(out["bins"].items()) + [("total", None)]:
        if e is None:
            s, bkg = out["total"]["signal"], out["total"]["background"]
        else:
            s, bkg = e.get("signal", {}), e.get("background", {}).get("total")
        by = s.get("by_mode", {})
        rows.append([b] + ["%.2f" % by[m] if m in by else "-" for m in modes]
                    + ["%.2f" % s["total"] if "total" in s else "-",
                       "%.0f" % bkg if bkg else "-"])
    widths = [max(len(r[i]) for r in rows + [head]) for i in range(len(head))]
    print(f"\n=== {out['name']}: yields over {w[0]}-{w[1]} GeV, "
          f"{out.get('lumi_fb', 0):.2f} fb^-1 ===")
    for r in [head] + rows:
        print("  " + "  ".join(c.rjust(n) for c, n in zip(r, widths)))


def print_next(args):
    name, ws = chain_name(args), rel(args.wsdir)
    opts = " ".join(COMBINE_OPTS)
    print(f"\nnext, from parametricfit/ (blinded, -t -1):\n"
          f"  combine -M AsymptoticLimits -m {MH} -t -1 --run expected "
          f"--rMin 0 -d {ws}/{name}.root -n _{name} {opts}\n"
          f"  combine -M Significance -m {MH} -t -1 --expectSignal 1 "
          f"-d {ws}/{name}.root -n _{name}_sig {opts}\n"
          f"  python -u biasStudy.py -c {args.category} "
          f"-b {' '.join(args.bins)} -y {args.year} --carddir {ws} "
          f"--ntoys 1000 -j 10")


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--category", default="ggHcat")
    ap.add_argument("-b", "--bins", nargs="+", default=["incl"])
    ap.add_argument("-y", "--year", default="Run3")
    ap.add_argument("--wsdir", default=None,
                    help="fit JSONs and workspaces (default WS_<year>)")
    ap.add_argument("--name", default=None,
                    help="name of the combined card and workspace "
                         "(default <category>_<bins>_<year>)")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the commands, run nothing")
    args = ap.parse_args()

    args.wsdir = os.path.abspath(os.path.join(
        HERE, args.wsdir or f"WS_{args.year}"))
    if not in_container() and not args.dry_run:
        raise SystemExit(
            "run this inside the Combine container (leave conda first):\n"
            "  singularity shell --cleanenv -B /work -B /ceph --nv "
            "/work/submit/jaeyserm/software/docker/"
            "combine-standalone_v9.2.1.sif")

    make_workspaces(args)
    make_cards(args)
    if len(args.bins) > 1:
        combine_cards(args)
    text2workspace(args)
    write_yields(args)
    print_next(args)


if __name__ == "__main__":
    main()
