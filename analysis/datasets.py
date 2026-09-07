"""
datasets.py — dataset selection and definitions for the HmumuRun3 analysis.

Sample definitions are inconfig/samples.yaml.
This module only builds file lists from config/samples.yaml.
utilsAna.SwitchSample calls findDIR when a sample is picked.
Execution lives in Hmm.loopOnDataset.

  Selection  (which IDs to run for a year/mode):
    - getMCList(year, mode)   -> list of MC sample IDs
    - getDataList(year, mode) -> list of data (negative) sample IDs

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
REGISTRY = os.path.join(HERE, "config", "samples.yaml")

VALID_YEARS = ["12022", "22022", "12023", "22023", "2024", "2025", "2026"]
VALID_MODES = ["isVBF", "isGGH", "isZinv", "isVlep", "isVhad", "isTTlep", "isTThad"]


# ===========================================================================
#  Registry
# ===========================================================================
_REG = None


def registry(path=None):
    """Parse config/samples.yaml once per process."""
    global _REG
    if _REG is None:
        path = path or REGISTRY
        if not os.path.isfile(path):
            raise FileNotFoundError(f"sample registry not found: {path}")
        with open(path) as f:
            _REG = yaml.safe_load(f)
    return _REG


DEFAULT_RUN = "run3"


def getXsec(name=None, run=DEFAULT_RUN):
    """Named cross-sections from samples.yaml, in fb.

    getXsec()                    -> the whole Run 3 table as a plain dict
    getXsec('ggH')               -> 52230.0   (13.6 TeV)
    getXsec('ggH', run='run2')   -> 48580.0   (13 TeV)
    """
    tables = registry()["xsecs"]
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
    table = registry()["branching_ratios"]
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
    reg = registry()

    if year not in reg["periods"]:
        raise KeyError(f"period {year!r} not in samples.yaml "
                       f"(have: {sorted(reg['periods'])})")
    period = reg["periods"][year]

    def fill(raw, prefix=""):
        return raw.format(ceph=period["ceph"],
                          scratch=period["scratch"],
                          year=year,
                          campaign=prefix + period["campaign"])

    thisdict = {}

    for sid, spec in reg["samples"].items():
        cp = spec.get("campaign_prefix")
        prefix = cp["value"] if cp and year in cp.get("years", []) else ""
        thisdict[int(sid)] = (fill(spec["path"], prefix), resolve_xsec(spec["xsec"]))

    for sid, raw in reg.get("data", {}).get(year, {}).items():
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
    groups = registry()["groups"]

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
    sel = registry().get("data_selection", {}).get(str(year), {})
    return list(sel.get(stream, []))
