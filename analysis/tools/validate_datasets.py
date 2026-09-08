#!/usr/bin/env python3
"""
validate_datasets.py — prove the YAML-backed datasets.py is equivalent to the
legacy hard-coded one, without touching the filesystem or needing ROOT.

Compares, for every period and every sample ID:
  * the resolved path pattern
  * the cross-section
  * the per-mode MC and data selections

    python tools/validate_datasets.py --legacy /path/to/old/analysis
"""

import argparse
import importlib.util
import sys
import types


def stub_root():
    stub = types.ModuleType("ROOT")
    stub.vector = lambda t: list
    sys.modules.setdefault("ROOT", stub)


def load_module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def norm(p):
    while "//" in p:
        p = p.replace("//", "/")
    return p.rstrip("/")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--legacy", default="/mnt/project",
                    help="directory holding the legacy datasets.py")
    ap.add_argument("--new", default="tools/datasets.py",
                    help="path to the YAML-backed datasets.py")
    args = ap.parse_args()

    stub_root()

    legacy = load_module("legacy_datasets", f"{args.legacy}/datasets.py")
    legacy.findDIR = lambda directory, useXROOTD=False: directory

    new = load_module("new_datasets", args.new)

    n_path = n_xsec = n_sel = 0
    fails = []

    for year in new.VALID_YEARS:
        old_d = legacy.BuildDict(year)
        new_d = new.BuildDict(year)

        missing = set(old_d) - set(new_d)
        extra = set(new_d) - set(old_d)
        if missing:
            fails.append(f"{year}: ids missing from YAML: {sorted(missing)}")
        if extra:
            fails.append(f"{year}: ids only in YAML: {sorted(extra)}")

        for sid in sorted(set(old_d) & set(new_d)):
            old_pattern, old_xsec = old_d[sid]
            new_pattern, new_xsec = new_d[sid]

            n_path += 1
            if norm(old_pattern) != norm(new_pattern):
                fails.append(f"{year} {sid}: path\n    old {norm(old_pattern)}\n"
                             f"    new {norm(new_pattern)}")

            if sid > 0:
                n_xsec += 1
                if old_xsec == 0:
                    if new_xsec != 0:
                        fails.append(f"{year} {sid}: xsec {old_xsec} -> {new_xsec}")
                elif abs(new_xsec - old_xsec) / abs(old_xsec) > 1e-9:
                    fails.append(f"{year} {sid}: xsec {old_xsec} -> {new_xsec}")

        for mode in new.VALID_MODES:
            n_sel += 1
            if sorted(legacy.getMCList(year, mode)) != sorted(new.getMCList(year, mode)):
                a, b = set(legacy.getMCList(year, mode)), set(new.getMCList(year, mode))
                fails.append(f"{year} {mode}: MC selection differs "
                             f"(legacy only {sorted(a - b)}, YAML only {sorted(b - a)})")
            if sorted(legacy.getDataList(year, mode)) != sorted(new.getDataList(year, mode)):
                fails.append(f"{year} {mode}: data selection differs "
                             f"{sorted(legacy.getDataList(year, mode))} vs "
                             f"{sorted(new.getDataList(year, mode))}")

    print(f"compared {n_path} path patterns, {n_xsec} cross-sections, "
          f"{n_sel} (period, mode) selections")
    if fails:
        print(f"\n{len(fails)} MISMATCHES:\n")
        for f in fails:
            print("  " + f)
        sys.exit(1)
    print("\nOK — YAML-backed datasets.py is equivalent to the legacy version")


if __name__ == "__main__":
    main()
