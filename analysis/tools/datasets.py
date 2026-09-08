"""
datasets.py — dataset selection and definitions for the HmumuRun3 analysis.

Sample definitions are in config/samples.yaml, parsed into CONFIG at import.
This module only builds file lists from config/samples.yaml.
utilsAna.SwitchSample calls findDIR when a sample is picked.
Execution lives in Hmm.loopOnDataset.

  Selection  (which IDs to run for a year/mode):
    - getMCList(year, mode)   -> list of MC sample IDs
    - getDataList(year, mode) -> list of data (negative) sample IDs

  Sample lists (datasets/<tag>_<year>.txt):
    - list_path(year, mode)   -> canonical path of the committed list
    - read_list(path, ...)    -> ids from a list file
    (writing them lives in make_datasets.py)

  Definition (what each ID is):
    - BuildDict(year)     : {ID: (path pattern, xsec)}   [no I/O]
    - findDIR(pattern)    : ROOT.vector<string> of files [glob the files]
    - getXsec(name, run)  : named cross-section in fb ('run3' default, 'run2')
    - getBR(name)         : branching ratio
"""

import glob
import os
import sys
from subprocess import check_output

import ROOT
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
ANALYSIS = os.path.dirname(HERE)                       # tools/ -> analysis/
SAMPLES_YAML = os.path.join(ANALYSIS, "config", "samples.yaml")
LISTDIR = os.path.join(ANALYSIS, "datasets")           # committed sample lists

VALID_YEARS = ["12022", "22022", "12023", "22023", "2024", "2025", "2026"]

# mode -> snapshot category suffix. Single source of truth: Hmm.py and
# run_all.py both import this instead of keeping their own copies.
MODE_MAP = {
    "isVBF":   "VBFcat",
    "isGGH":   "ggHcat",
    "isZinv":  "Zinvcat",
    "isVlep":  "VLcat",
    "isVhad":  "VHcat",
    "isTTlep": "TTLcat",
    "isTThad": "TTHcat",
}
VALID_MODES = list(MODE_MAP)

# mode -> short tag used in datasets/<tag>_<year>.txt  (VBFcat -> VBF)
CATEGORY = {m: c[:-3] for m, c in MODE_MAP.items()}


# ===========================================================================
#  config/samples.yaml, parsed once at import
# ===========================================================================
REQUIRED_BLOCKS = ["periods", "xsecs", "branching_ratios", "groups", "samples"]


def _load_samples_yaml(path):
    """Read and sanity-check config/samples.yaml.

    Anything wrong with the file is fatal: this module is useless without it,
    and a clear message now beats a confusing KeyError later.
    """
    if not os.path.isfile(path):
        raise SystemExit(f"ERROR: sample definitions not found: {path}")

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

    for run in ("run3",):
        if run not in cfg["xsecs"]:
            raise SystemExit(f"ERROR: {path} has no 'xsecs.{run}' table")

    return cfg


CONFIG = _load_samples_yaml(SAMPLES_YAML)


DEFAULT_RUN = "run3"


def getXsec(name=None, run=DEFAULT_RUN):
    """Named cross-sections from samples.yaml, in fb.

    getXsec()                    -> the whole Run 3 table as a plain dict
    getXsec('ggH')               -> 52230.0   (13.6 TeV)
    getXsec('ggH', run='run2')   -> 48580.0   (13 TeV)
    """
    tables = CONFIG["xsecs"]
    if run not in tables:
        raise KeyError(f"unknown run {run!r} in samples.yaml "
                       f"(have: {sorted(tables)})")
    table = tables[run]
    if name is None:
        return dict(table)
    if name not in table:
        raise KeyError(f"unknown cross-section {name!r} for {run} in samples.yaml "
                       f"(have: {sorted(table)})")
    return float(table[name])


def getBR(name):
    """A branching ratio (or decay fraction) from samples.yaml."""
    table = CONFIG["branching_ratios"]
    if name not in table:
        raise KeyError(f"unknown branching ratio {name!r} in samples.yaml "
                       f"(have: {sorted(table)})")
    return float(table[name])


def resolve_xsec(spec):
    """Turn a sample's 'xsec' entry into a number in fb.

    Accepted forms:
        2219000                             a plain value
        {ref: W}                            named process from xsecs.run3
        {ref: VBFH, br: H_to_mumu}          times a branching ratio
        {ref: Wm, br: [H_to_WW, W_to_qq]}   times several
        {value: 21650, factor: 0.6}         times an empirical scaling
    'br' entries come from branching_ratios and multiply together.
    'factor' is for anything that is not a physics constant.
    'run' may be given alongside 'ref' to indicate Run2 vs Run3.
    """
    if not isinstance(spec, dict):
        return float(spec)

    if "ref" in spec:
        xsec = getXsec(spec["ref"], spec.get("run", DEFAULT_RUN))
    elif "value" in spec:
        xsec = float(spec["value"])
    else:
        raise KeyError(f"xsec entry needs 'ref' or 'value': {spec!r}")

    br = spec.get("br", [])
    for name in ([br] if isinstance(br, str) else br):
        xsec *= getBR(name)

    return xsec * float(spec.get("factor", 1.0))


# ===========================================================================
#  File discovery
# ===========================================================================
def findDIR(directory, useXROOTD=False):
    """Expand one path pattern into a sorted list of ROOT files."""
    print("resolving: %s" % directory)
    rootFiles = ROOT.vector("string")()

    if useXROOTD and "/data/submit/cms" in directory:
        xrd = "root://submit50.mit.edu/"
        xrdpath = directory.replace("/data/submit/cms", "")
        out = check_output(["xrdfs", xrd, "ls", xrdpath]).decode(sys.stdout.encoding)
        for entry in sorted(out.split()):
            if any(skip in entry for skip in ("failed/", "log/", ".txt")):
                continue
            rootFiles.push_back(xrd + entry)
    else:
        for f in sorted(glob.glob("{}/*.root".format(directory))):
            rootFiles.push_back(f)

    print("  -> %d files" % len(rootFiles))
    return rootFiles


# ===========================================================================
#  Dataset definition: ID -> (path pattern, xsec)
# ===========================================================================
def BuildDict(year):
    """Map every sample ID for this period to (pattern, xsec).

    Patterns are strings, not file lists: no filesystem access happens here.
    SwitchSample resolves the one sample it is asked for.
    """
    year = str(year)

    if year not in CONFIG["periods"]:
        raise KeyError(f"period {year!r} not in samples.yaml "
                       f"(have: {sorted(CONFIG['periods'])})")
    period = CONFIG["periods"][year]

    def fill(raw, prefix=""):
        return raw.format(ceph=period["ceph"],
                          scratch=period["scratch"],
                          year=year,
                          campaign=prefix + period["campaign"])

    thisdict = {}

    for sid, spec in CONFIG["samples"].items():
        cp = spec.get("campaign_prefix")
        prefix = cp["value"] if cp and year in cp.get("years", []) else ""
        thisdict[int(sid)] = (fill(spec["path"], prefix), resolve_xsec(spec["xsec"]))

    for sid, raw in CONFIG.get("data", {}).get(year, {}).items():
        thisdict[int(sid)] = (fill(raw), -1)

    return thisdict


# ===========================================================================
#  Sample selection: which IDs to run for a given (year, mode)
# ===========================================================================
def getMCList(year, mode):
    """MC sample IDs to process for this (year, mode).

    The ID lists come from the 'groups' table in samples.yaml.
    """
    year = str(year)
    groups = CONFIG["groups"]

    def g(name):
        if name not in groups:
            raise KeyError(f"group {name!r} not in samples.yaml "
                           f"(have: {sorted(groups)})")
        return list(groups[name])

    is_v12 = year in ("12022", "22022", "12023", "22023")
    is_2024 = year == "2024"

    mc = g("signal_hmm")

    if is_v12:
        mc += g("zgamma_v12")
    elif is_2024:
        mc += g("zgamma_v15")

    if mode == "isVhad":
        if is_2024:
            mc += g("dy_ptbinned_v15")
        elif is_v12:
            mc += g("dy_ptbinned_v12")
    else:
        mc += g("dy_mu_tau") if is_2024 else g("dy_inclusive")

    mc += g("dy_ewk")
    if mode == "isVBF":
        mc += g("dy_ewk_mass")

    mc += g("ttbar_2l")
    mc += g("vv")
    mc += g("vvv")
    mc += g("ttv")
    mc += g("ttz_qq_v15") if is_2024 else g("ttz_qq_v12")
    mc += g("ttgamma_v15") if is_2024 else g("ttgamma_v12")
    mc += g("top_other")

    return list(dict.fromkeys(mc))   # de-duplicate, preserve order


def getDataList(year, mode, stream="prompt"):
    """Data (negative) sample IDs for this period.

    'prompt' is the Muon primary datasets
    'parking' is the ParkingDoubleMuonLowMass streams.
    """
    sel = CONFIG.get("data_selection", {}).get(str(year), {})
    return list(sel.get(stream, []))


def sample_label(sid, year, short=True):
    """Human-readable dataset name for a sample id, from its path pattern."""
    if sid < 0:
        raw = CONFIG.get("data", {}).get(str(year), {}).get(sid, "")
    else:
        raw = CONFIG["samples"].get(sid, {}).get("path", "")
    if not raw:
        return "?"

    parts = raw.split("/")
    if "{year}" not in parts:
        return "?"

    # everything between {year} and the campaign/NANOAODSIM level is the name
    name_parts = []
    for p in parts[parts.index("{year}") + 1:]:
        if not p or p == "*" or p == "NANOAODSIM" or p.startswith("{"):
            break
        name_parts.append(p)
    name = "/".join(name_parts) or "?"

    if short:
        for cut in ("_TuneCP5", "_TuneCH3", "_Tune", "-Tune"):
            if cut in name:
                name = name.split(cut)[0]
                break
    return name


def resolve_ids(year, mode, spec, stream="prompt", _seen=None):
    """Turn a comma-separated selection string into a list of sample ids.

    Each token may be:
        103, -41        a numeric id
        signal_hmm      a group name from 'groups' in samples.yaml
        mc / data       everything of that kind for this (year, mode)
        'DYto2Mu*'      a glob matched against the dataset name
        @lists/vbf.txt  a file of tokens, one per line, '#' comments allowed

    Tokens are resolved against every sample in samples.yaml, not just this
    mode's default selection, so an id outside it can still be run deliberately.
    Anything unmatched raises, so a typo fails loudly instead of silently
    processing nothing.
    """
    year = str(year)
    all_mc = sorted(int(k) for k in CONFIG["samples"])
    all_data = sorted((int(k) for k in CONFIG.get("data", {}).get(year, {})), reverse=True)
    _seen = _seen or set()

    out = []
    for token in (t.strip() for t in spec.split(",")):
        if not token:
            continue

        if token.startswith("@"):
            path = os.path.abspath(token[1:])
            if path in _seen:
                raise ValueError(f"sample list {token} includes itself")
            if not os.path.isfile(path):
                raise FileNotFoundError(f"sample list not found: {path}")
            _seen.add(path)
            with open(path) as fh:
                for line in fh:
                    line = line.split("#", 1)[0].strip()
                    if line:
                        out += resolve_ids(year, mode, line, stream, _seen)
            continue

        if token.lstrip("-").isdigit():
            out.append(int(token))
            continue

        if token == "mc":
            out += getMCList(year, mode)
            continue
        if token == "data":
            out += getDataList(year, mode, stream)
            continue

        if token in CONFIG["groups"]:
            out += list(CONFIG["groups"][token])
            continue

        if any(c in token for c in "*?["):
            import fnmatch
            hits = [s for s in all_mc + all_data
                    if fnmatch.fnmatch(sample_label(s, year, short=False), token)
                    or fnmatch.fnmatch(sample_label(s, year), token)]
            if not hits:
                raise ValueError(f"--samples pattern {token!r} matched nothing")
            out += hits
            continue

        raise ValueError(
            f"--samples token {token!r} is not an id, a group name, 'mc', 'data', "
            f"a glob, or @file. Known groups: {', '.join(sorted(CONFIG['groups']))}")

    return list(dict.fromkeys(out))


def list_path(year, mode, listdir=None):
    """Canonical path of the committed sample list: datasets/<tag>_<year>.txt"""
    if mode not in CATEGORY:
        raise KeyError(f"unknown mode {mode!r} (have: {', '.join(VALID_MODES)})")
    return os.path.join(listdir or LISTDIR, f"{CATEGORY[mode]}_{year}.txt")


def read_list(path, year, mode, stream="prompt"):
    """Read a sample list file into a list of ids."""
    if not os.path.isfile(path):
        raise FileNotFoundError(
            f"sample list not found: {path}\n"
            f"  generate it with:  python make_datasets.py --write {year} {mode}")
    return resolve_ids(year, mode, f"@{path}", stream)
