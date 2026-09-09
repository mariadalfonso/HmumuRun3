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
# Each MC process is a list of mc IDs. Data is handled separately (mc<0), since
# it is the only process with blinding logic attached.
#
# NOTE the grouping is exactly the one prepareHisto.getHisto() used before:
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
# Booking one histogram per process used to mean one jitted filter per process,
# each a long "mc==x || mc==y || ..." chain (make_filter). Instead we Define a
# single integer column once and every histogram filters on "procGroup==k".
# Far less for Cling to compile, and one integer compare per event.
#
# Index 0 is reserved for "unassigned MC" so an ID that is in the chain but in
# no group is silently excluded rather than silently landing in a real process.
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
# SR-sideband mass restriction: m_mumu in [110,120] or [130,150] GeV,
# deliberately excluding the innermost [120,130] core. Independent of
# plot_vars.get_blind_range("dimu_mass") == (110,150), which is the FULL blind
# window used for the Inclusive region's mass-window blinding logic; this is a
# narrower, two-piece window.
SR_SIDEBAND_FILTER = ("((HiggsCandCorrMass>=110 && HiggsCandCorrMass<=120) || "
                      "(HiggsCandCorrMass>=130 && HiggsCandCorrMass<=150))")

# "filter"      : extra mass restriction, or None to use the full snapshot window
# "blindable"   : whether the --blind flag applies here. False means data is
#                 shown unblinded because the region restriction is ITSELF the
#                 safety mechanism (it already excludes the signal core).
# "label"       : on-plot mass-range line
#
# "Unrestricted" exists so SR_sideband plots can show the FULL predicted signal:
# a narrow resonance concentrates almost entirely in [120,130], the exact core
# SR_sideband excludes, so a region-restricted signal would be a misleading
# sliver. It used to be a second getHisto() call per variable (doubling the
# event loops); now it is just another branch of the same computation graph.
REGIONS = {
    "Inclusive": {
        "filter": None,
        "blindable": True,
        "label": "m_{#mu#mu}#in[70,200]GeV",
    },
    "SR_sideband": {
        "filter": SR_SIDEBAND_FILTER,
        "blindable": False,
        "label": "m_{#mu#mu}#in[110,120]#cup[130,150]GeV",
    },
    "Unrestricted": {
        "filter": None,
        "blindable": False,
        "label": "m_{#mu#mu}#in[70,200]GeV",
    },
}

# Which region supplies the SIGNAL curves when drawing a given region (see the
# "Unrestricted" note above).
SIGNAL_REGION_FOR = {
    "Inclusive":    "Inclusive",
    "SR_sideband":  "Unrestricted",
    "Unrestricted": "Unrestricted",
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
    if category in ("VLcat", "TTLcat", "TTHcat"):
        names += ["muon1_sip3d", "muon2_sip3d"]
    if category in ("ggHcat", "TTHcat"):
        # nGoodJetsAll is only written to the snapshot for isGGH/isTThad
        # (see Hmm.py's mode_branches).
        names.append("njets")
    names.append("deta_muons")

    return names


# Output-directory grouping for SummaryPlots.py. A variable not listed here
# falls into "muons".
VAR_GROUPS = {
    "dimu_mass":       "mass",
    "mva":             "mva",
    "category_vlcat":  "category",
    "category_ttlcat": "category",
    "category_tthcat": "category",
}
DEFAULT_VAR_GROUP = "muons"

# The "mass" group deliberately skips the region split: the SR-sideband window
# is just a truncated/gapped view of the same dimu_mass spectrum, so a second
# copy adds nothing the way it does for other variables.
GROUPS_WITHOUT_REGION_SPLIT = ["mass"]


def group_of(varname):
    return VAR_GROUPS.get(varname, DEFAULT_VAR_GROUP)


def vars_in_group(category, group):
    return [v for v in get_active_vars(category) if group_of(v) == group]


def groups_for(category):
    """Groups present for this category, in a stable display order."""
    order = ["mass", "mva", "category", "muons"]
    present = {group_of(v) for v in get_active_vars(category)}
    return [g for g in order if g in present]


# ===========================================================================
#  Histogram file naming / layout
# ===========================================================================
# Layout inside the file (one directory per region):
#     <region>/<varname>_<process>
# Modelled on FastFrames' <syst>/<var>_<region> convention. A <syst> level can
# be inserted above <region> later without touching the readers, provided they
# go through histo_path().
def histo_filename(category, year):
    """year is the internal '_<year>' convention (leading underscore)."""
    return f"histos_{category}{year}.root"


def histo_path(region, varname, process):
    return f"{region}/{varname}_{process}"
