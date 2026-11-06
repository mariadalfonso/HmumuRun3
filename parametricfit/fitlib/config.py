"""
config.py -- settings shared by every script of the chain: the fit window,
its binning, the production modes per category, the luminosity lookup and
the git commit recorded in the outputs.

BINS_PER_GEV is read at call time, so sigFit.py --bins-per-gev (or a study)
can override it by assigning config.BINS_PER_GEV.
"""

import subprocess

XLOW, XHIGH = 110.0, 150.0
BINS_PER_GEV = 10         # 0.1 GeV bins, matching SIGfits.py and fitBkg.
                          # Checked with a closure-test scan (see README):
                          # unbiased at -0.11% +/- 0.25% on sigma. Anything
                          # from 0.1 to 0.5 GeV is fine; 1.0 GeV biases sigma
                          # by +2%. Override with --bins-per-gev.

MC_LIST = {
    "ggHcat":  ["ggH", "qqH", "VH", "ttH"],
    "VBFcat":  ["ggH", "qqH", "VH", "ttH"],
    "Zinvcat": ["qqH", "VH", "ttH"],
    "VLcat":   ["VH", "ttH"],
    "VHcat":   ["VH", "ttH"],
    "TTLcat":  ["VH", "ttH"],
    "TTHcat":  ["VH", "ttH"],
}


def get_lumi(year, full_year=False):
    """Integrated luminosity (fb^-1) for a year.

    By default that of the signal MC, so Run3 -> 2022 + 2023 + 2024 only.
    full_year=True gives that of the data (Run3 -> 2022 ... 2026), which the
    background fit reads and the signal yield is extrapolated to.
    Imported lazily so --selftest does not pull in prepareFits.
    """
    from .prepareFits import (get_lumi as _get_lumi, SIGNAL_YEAR_GROUPS,
                              DATA_YEAR_GROUPS)
    return _get_lumi(year, DATA_YEAR_GROUPS if full_year
                     else SIGNAL_YEAR_GROUPS)


def default_wsdir(year):
    """WS_<year>: where every step of the chain reads and writes for a year
    when --wsdir is not given, e.g. WS_Run3, WS_2024."""
    return f"WS_{str(year).lstrip('_')}"


def git_commit():
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "--short", "HEAD"],
            stderr=subprocess.DEVNULL).decode().strip()
    except Exception:
        return "unknown"
