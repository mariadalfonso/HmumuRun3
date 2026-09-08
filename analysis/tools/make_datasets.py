#!/usr/bin/env python3
"""
make_datasets.py — generate the committed sample lists in datasets/.

Reads config/samples.yaml (via tools/datasets.py) and writes one text file
per category and year: datasets/<tag>_<year>.txt. Those files are committed
to git, so this only needs running when the datasets change.

    python tools/make_datasets.py --write-all           # every category and year
    python tools/make_datasets.py --write 2024 isVBF    # just one
    python tools/make_datasets.py --show 2024 isVBF     # print, write nothing

Nothing in the analysis imports this module; Hmm.py and run_all.py read the
generated files through datasets.read_list.
"""

import argparse
import os
from datetime import datetime
from pathlib import Path

from tools.datasets import (CATEGORY, LISTDIR, SAMPLES_YAML, VALID_MODES,
                            VALID_YEARS, getDataList, getMCList, list_path,
                            read_list, sample_label)


def write_list(year, mode, path=None, ids=None, stream="prompt"):
    """Generate the committed sample list for one (year, mode).

    Each line is a sample id with its dataset name as a comment, so the file
    is self-documenting, greppable, and diffable in git. Comment a line out
    to skip that sample without touching samples.yaml.
    """
    year = str(year)
    if ids is None:
        ids = getMCList(year, mode) + getDataList(year, mode, stream)

    path = Path(path or list_path(year, mode))
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        fh.write(f"# {CATEGORY[mode]} {year} -- sample list for {mode}\n")
        fh.write(f"# generated {datetime.now().strftime('%Y-%m-%d %H:%M')} "
                 f"from config/{os.path.basename(SAMPLES_YAML)} "
                 f"by tools/make_datasets.py ({stream} data stream)\n")
        fh.write("#\n")
        fh.write("# Committed to git: regenerate only when the datasets change.\n")
        fh.write("# Comment a line out to skip that sample.\n")
        mc = [i for i in ids if i > 0]
        data = [i for i in ids if i < 0]
        if mc:
            fh.write(f"\n# --- MC ({len(mc)}) ---\n")
            for sid in mc:
                fh.write(f"{sid:<6} # {sample_label(sid, year)}\n")
        if data:
            fh.write(f"\n# --- data ({len(data)}) ---\n")
            for sid in data:
                fh.write(f"{sid:<6} # {sample_label(sid, year)}\n")
    return path


def write_all_lists(listdir=None, stream="prompt"):
    """Regenerate datasets/<tag>_<year>.txt for every (year, mode)."""
    written = []
    for year in VALID_YEARS:
        for mode in VALID_MODES:
            path = write_list(year, mode,
                              list_path(year, mode, listdir), stream=stream)
            n = len(read_list(path, year, mode, stream))
            written.append((path, n))
    return written


# ===========================================================================
#  CLI: regenerate the committed sample lists
# ===========================================================================
def _main():
    ap = argparse.ArgumentParser(
        prog="python tools/make_datasets.py",
        description="Generate the committed sample lists in datasets/ from "
                    "config/samples.yaml. Run this only when the datasets change; "
                    "the .txt files are in git so normal running needs no regeneration.")
    ap.add_argument("--write", nargs=2, metavar=("YEAR", "MODE"),
                    help="regenerate one list")
    ap.add_argument("--write-all", action="store_true",
                    help="regenerate every (year, mode) list")
    ap.add_argument("--listdir", default=None,
                    help=f"output directory (default: {LISTDIR})")
    ap.add_argument("--stream", default="prompt", choices=["prompt", "parking"])
    ap.add_argument("--show", nargs=2, metavar=("YEAR", "MODE"),
                    help="print a list's resolved ids without writing anything")
    args = ap.parse_args()

    if args.show:
        year, mode = args.show
        ids = getMCList(year, mode) + getDataList(year, mode, args.stream)
        for sid in ids:
            print(f"{sid:<6} # {sample_label(sid, year)}")
        print(f"# {len(ids)} samples for {year} {mode}")
        return

    if args.write:
        year, mode = args.write
        if year not in VALID_YEARS:
            raise SystemExit(f"unknown year {year!r}; expected one of {', '.join(VALID_YEARS)}")
        if mode not in VALID_MODES:
            raise SystemExit(f"unknown mode {mode!r}; expected one of {', '.join(VALID_MODES)}")
        path = write_list(year, mode, list_path(year, mode, args.listdir),
                          stream=args.stream)
        print(f"wrote {path} ({len(read_list(path, year, mode, args.stream))} samples)")
        return

    if args.write_all:
        for path, n in write_all_lists(args.listdir, args.stream):
            print(f"  {os.path.basename(path):<20} {n:>4} samples")
        print("\ncommit datasets/ so nobody needs to regenerate these")
        return

    ap.print_help()


if __name__ == "__main__":
    _main()
