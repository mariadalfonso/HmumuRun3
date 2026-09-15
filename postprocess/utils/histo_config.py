"""
histo_config.py -- shared configuration for the histogramming + plotting stages.

"""

from .LoadTree import (hmumu, hzgamma, hww, vv, tt2l, ttV2223, ttV24, top,
                      dyewk, dy_2223, dy_24, dy_pt2223, dy_pt24, dy_minllo,
                      dy_j)


# ===========================================================================
#  Categories
# ===========================================================================
CATEGORIES = ["VBFcat", "ggHcat", "VLcat", "TTLcat", "TTHcat",
              "VHcat", "Zinvcat"]

# Display names for the on-plot category line.
CATEGORY_LABELS = {
    "VBFcat":  "VBF cat.",
    "ggHcat":  "ggF cat.",
    "VLcat":   "VH-lep cat.",
    "TTLcat":  "ttH-lep cat.",
    "TTHcat":  "ttH-had cat.",
    "VHcat":   "VH-had cat.",
    "Zinvcat": "Zinv cat.",
}

# Categories whose signal legend shows the leptonic modes (ttH/WH/ZH) rather
# than the hadronic ones (VBF/ggH).
LEP_CATEGORIES = ["VLcat", "TTLcat", "TTHcat", "VHcat", "Zinvcat"]


# ===========================================================================
#  Processes
# ===========================================================================
# mc IDs per process. Data (mc<0) is separate: it is the only process with
# blinding logic attached.
#   hDY  = all DY QCD flavours/binnings
#   hTop = single-top + tt+V + 4-top (everything top-like except tt->2l)
#   hZg  = H->Zgamma and H->WW
MC_PROCESSES = {
    "hDY":   dy_2223 + dy_24 + dy_pt2223 + dy_pt24 + dy_minllo + dy_j,
    "hEWK":  dyewk,
    "hTT2L": tt2l,
    "hTop":  top + ttV2223 + ttV24,
    "hVV":   vv,
    "hVBFH": ["10"],
    "hggH":  ["11"],
    "hWH":   ["12", "13"],
    "hZH":   ["14"],
    "hTTH":  ["15"],
    "hZg":   hzgamma + hww,
}

DATA_PROCESS = "hData"

# Stacking order for the background stack (bottom -> top), and the signal
# processes drawn unstacked on top of it.
BKG_PROCS = ["hTT2L", "hTop", "hZg", "hVV", "hEWK", "hDY"]
SIG_PROCS = ["hWH", "hTTH", "hZH", "hVBFH", "hggH"]

# All processes written to / read from the histogram file.
ALL_PROCESSES = list(MC_PROCESSES.keys()) + [DATA_PROCESS]

# Human-readable labels, used by both plotting scripts' legends.
PROCESS_LABELS = {
    "hData": "Data",
    "hDY":   "DY+jets(QCD)",
    "hEWK":  "DY+jets(EWK)",
    "hVV":   "VV+VVV",
    "hTT2L": "t#bar{t} 2l",
    "hTop":  "Top(1l,tW/tZq,ttV/4t)",
    "hZg":   "H#rightarrowZ#gamma+jets",
    "hTTH":  "ttH",
    "hWH":   "WH",
    "hZH":   "ZH",
    "hVBFH": "VBF H",
    "hggH":  "ggH",
}

# Legend order for the stacked plots (process, draw-style).
LEGEND_BKG_ORDER = ["hData", "hDY", "hEWK", "hVV", "hTT2L", "hTop", "hZg"]
LEGEND_SIG_LEP_ORDER = ["hTTH", "hWH", "hZH"]
LEGEND_SIG_HAD_ORDER = ["hVBFH", "hggH"]


# ---------------------------------------------------------------------------
# mc -> process-group index
# ---------------------------------------------------------------------------
# Integer process index, defined once on the base node so each histogram
# filters on "procGroup==k" instead of a jitted "mc==x || mc==y ..." chain.
# Index 0 is unassigned MC: an ID in the chain but in no group is excluded
# rather than landing in a real process.
PROC_INDEX = {name: i for i, name in enumerate(MC_PROCESSES, start=1)}
DATA_INDEX = -1


def proc_group_expr():
    """C++ expression defining the procGroup column (see PROC_INDEX)."""
    clauses = []
    for name, ids in MC_PROCESSES.items():
        if not ids:
            continue
        cond = " || ".join(f"mc=={i}" for i in ids)
        clauses.append(f"({cond}) ? {PROC_INDEX[name]}")
    clauses.append(f"(mc<0) ? {DATA_INDEX}")
    return " : ".join(clauses) + " : 0"


# ===========================================================================
#  Regions
# ===========================================================================
# Region windows are boolean columns defined in Hmm.py, written to every
# snapshot, so "filter" is just a branch name:
#
#   isZ    (76, 106)     isH    (115, 135)
#   isHSB  [110,115] U [135,150]   (disjoint from isH)
#
# "filter"    : region branch name, or None for the full snapshot window
# "blindable" : whether --blind applies
# "data"      : book a data histogram at all
# "signal"    : book/draw the signal processes
# "label"     : on-plot mass-range line
#
# "Unrestricted" supplies the full signal prediction for SR_sideband plots,
# where the resonance sits in the window isHSB excludes.
REGIONS = {
    "Inclusive": {
        "filter": None,
        "blindable": True,
        "data": True,
        "signal": True,
        "label": "m_{#mu#mu}#in[70,200]GeV",
    },
    "SR_sideband": {
        "filter": "isHSB",
        "blindable": False,
        "data": True,
        "signal": True,
        "label": "m_{#mu#mu}#in[110,115]#cup[135,150]GeV",
    },
    "Zboson_CR": {
        "filter": "isZ",
        "blindable": False,
        "data": True,
        "signal": False,
        "label": "m_{#mu#mu}#in[76,106]GeV",
    },
    "SR_plus_sideband": {
        # full fit range (isH and isHSB are disjoint); contains the SR
        "filter": "isH || isHSB",
        "blindable": False,
        "data": False,
        "signal": True,
        "label": "m_{#mu#mu}#in[110,150]GeV",
    },
    "Unrestricted": {
        "filter": None,
        "blindable": False,
        "data": True,
        "signal": True,
        "label": "m_{#mu#mu}#in[70,200]GeV",
    },
}

# Where a region's signal curves come from. None = no signal drawn.
SIGNAL_REGION_FOR = {
    "Inclusive":         "Inclusive",
    "SR_sideband":       "Unrestricted",
    "Zboson_CR":         None,
    "SR_plus_sideband":  "SR_plus_sideband",
    "Unrestricted":      "Unrestricted",
}


# ===========================================================================
#  Variables
# ===========================================================================
def get_active_vars(category):
    """The variables actually plotted for `category`.

    Was maintained in two places (SummaryPlots.py's GROUP_FUNCS dispatch and
    SummaryPlotsShapes.py's own copy). Now one list drives the histogramming
    stage AND both plotting scripts, so they cannot fall out of sync.

    Grouped so SummaryPlots.py can still organise output into per-group
    subdirectories; see VAR_GROUPS.
    """
    names = ["dimu_mass", "mva"]

    if category == "VLcat":
        names.append("category_vlcat")
    elif category == "TTLcat":
        names.append("category_ttlcat")
    elif category == "TTHcat":
        names.append("category_tthcat")

    names += [
        "muon1_pt", "muon2_pt", "muon1_eta", "muon2_eta",
        "dimuon_pt", "dimuon_eta", "dimuon_rapidity",
        "muon1_norm_pt", "muon2_norm_pt", "costhetacs", "phistarcs",
    ]
    if category in ("Zinvcat", "VLcat", "TTHcat"):
        names += ["muon1_sip3d", "muon2_sip3d"]
    if category in ("ggHcat", "TTHcat"):
        names.append("njets")
    names.append("deta_muons")

    # jets and MET.
    names.append("met_pt")
    if category in ("ggHcat", "TTHcat", "TTLcat"):
        names += ["jet1_pt", "jet1_eta"]

    return names


# Plot groups, by how data is treated: "mass" has the blind window cut out,
# "mva" is blinded above the per-category score cut, "objects" is not blinded.
# Anything unlisted falls into "objects".
VAR_GROUPS = {
    "dimu_mass":       "mass",
    "mva":             "mva",
    "category_vlcat":  "mva",
    "category_ttlcat": "mva",
    "category_tthcat": "mva",
}
DEFAULT_VAR_GROUP = "objects"

# Which regions each group is DRAWN in; everything in REGIONS is always booked.
# Output tree is <cat>_<year>/<group>/<region>/, or <cat>_<year>/<group>/ for a
# single-region group.
PLOT_REGIONS = {
    "mass":    ["Inclusive"],
    "mva":     ["Inclusive", "SR_sideband"],
    "objects": ["Inclusive", "SR_sideband", "Zboson_CR"],
}
DEFAULT_PLOT_REGIONS = ["Inclusive"]


def plot_regions_for(group):
    return PLOT_REGIONS.get(group, DEFAULT_PLOT_REGIONS)


def group_of(varname):
    return VAR_GROUPS.get(varname, DEFAULT_VAR_GROUP)


def vars_in_group(category, group):
    return [v for v in get_active_vars(category) if group_of(v) == group]


def groups_for(category):
    """Groups present for this category, in a stable display order."""
    order = ["mass", "mva", "objects"]
    present = {group_of(v) for v in get_active_vars(category)}
    return [g for g in order if g in present]


# ===========================================================================
#  Histogram file naming / layout
# ===========================================================================
# Layout inside the file, one directory per region:
#     <region>/<varname>_<process>
def histo_filename(category, year):
    """year is the internal '_<year>' convention (leading underscore)."""
    return f"histos_{category}{year}.root"


def histo_path(region, varname, process):
    return f"{region}/{varname}_{process}"
