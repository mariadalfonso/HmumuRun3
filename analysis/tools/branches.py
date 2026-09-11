"""
branches.py — the snapshot output schema.

Reads config/branches.yaml and builds the branch list that Hmm.py hands to
Snapshot. Separate from datasets.py, which is about input samples: this is
about what the analysis writes out.

    branch_list(mode, is_mc)        -> list of branch names for one category
    check_branches(df, names, mode) -> the same list, after checking every
                                       name is a column of df
"""

import difflib
import os

import yaml

from tools.datasets import VALID_MODES

HERE = os.path.dirname(os.path.abspath(__file__))
ANALYSIS = os.path.dirname(HERE)                       # tools/ -> analysis/
BRANCHES_YAML = os.path.join(ANALYSIS, "config", "branches.yaml")

REQUIRED_BLOCKS = ["base", "per_mode", "per_mode_mc_only"]


# ===========================================================================
#  config/branches.yaml, parsed once at import
# ===========================================================================
def _load_branches_yaml(path):
    """Read and sanity-check config/branches.yaml."""
    if not os.path.isfile(path):
        raise SystemExit(f"ERROR: branch definitions not found: {path}")

    try:
        with open(path) as f:
            cfg = yaml.safe_load(f)
    except yaml.YAMLError as exc:
        raise SystemExit(f"ERROR: {path} is not valid YAML\n{exc}")

    if not isinstance(cfg, dict):
        raise SystemExit(f"ERROR: {path} should be a mapping of top-level blocks, "
                         f"got {type(cfg).__name__}")

    missing = [k for k in REQUIRED_BLOCKS if k not in cfg]
    if missing:
        raise SystemExit(f"ERROR: {path} is missing required block(s): "
                         f"{', '.join(missing)}\n"
                         f"       expected: {', '.join(REQUIRED_BLOCKS)}")

    absent = [m for m in VALID_MODES if m not in cfg["per_mode"]]
    if absent:
        raise SystemExit(f"ERROR: {path} has no per_mode entry for: "
                         f"{', '.join(absent)}")

    return cfg


CONFIG = _load_branches_yaml(BRANCHES_YAML)


# ===========================================================================
#  Building and checking the branch list
# ===========================================================================
def branch_list(mode, is_mc):
    """Snapshot output branches for one category.

    base + per_mode[mode], plus per_mode_mc_only[mode] for MC. Duplicates are
    dropped while order is preserved, so a name may appear in more than one
    block without upsetting Snapshot.
    """
    if mode not in CONFIG["per_mode"]:
        raise KeyError(f"unknown mode {mode!r} in branches.yaml "
                       f"(have: {', '.join(sorted(CONFIG['per_mode']))})")

    out = list(CONFIG["base"]) + list(CONFIG["per_mode"][mode] or [])
    if is_mc:
        out += list(CONFIG["per_mode_mc_only"].get(mode) or [])
    return list(dict.fromkeys(out))


def check_branches(df, branches, mode):
    """Fail early if a requested branch is not a column of df.

    Called before Snapshot, so a typo in branches.yaml is reported in the
    first seconds of a job rather than after the event loop. Returns the
    list unchanged when everything is present.
    """
    have = set(str(c) for c in df.GetColumnNames())
    absent = [b for b in branches if b not in have]
    if not absent:
        return branches

    print(f"ERROR: {len(absent)} branch(es) in config/branches.yaml are not "
          f"defined for mode {mode}:")
    for b in absent:
        near = difflib.get_close_matches(b, have, 3)
        print(f"   {b}" + (f"   did you mean: {', '.join(near)}" if near else ""))
    raise SystemExit("fix config/branches.yaml, or define the missing columns")
