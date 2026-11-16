"""
bdtBins.py -- the BDT-bin cuts of every category: SELMVA, and the check that
a fit output was made with the cuts in force now.

No ROOT here, so the Combine-container side (buildWorkspace.py) can import
it. prepareFits re-exports SELMVA and get_selection.

Where the cuts come from, per category:
  1. fitlib/bdt_bins.json, if it has the category: the split that
     optimizeBins.py --apply wrote, with how it was chosen. It replaces the
     category's whole entry (a 3-bin split never keeps an old bdt3).
  2. otherwise SELMVA_DEFAULT below.
SELMVA_SOURCE says which, per category. The bins file is tracked in git:
commit it with the change of cuts, and rerun the chain.

Every fit writes provenance["selection"], the cut string of its bin (and
"selection_source"). check_selection() compares such a record with the cut
in force and stops on a mismatch, so files made with the old cuts are not
combined with new ones by mistake.
"""

import json
import os

BINS_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "bdt_bins.json")

# BDT-bin selection thresholds, keyed by category then bin label; the
# defaults when fitlib/bdt_bins.json does not list the category.
# VLcat: 1 W->e; 2 Z->ee; 3 W->mu; 4 Z->mumu
SELMVA_DEFAULT = {
    "ggHcat":  {"bdt2":"discrMVA>=0.78", "bdt1":"discrMVA>=0.5 && discrMVA<0.78",
                "bdt0":"discrMVA<0.5",   "incl":"true", "":"true"},
    "VBFcat":  {"bdt3":"discrMVA>=0.96", "bdt2":"discrMVA>=0.90 && discrMVA<0.96",
                "bdt1":"discrMVA>=0.72 && discrMVA<0.90",
                "bdt0":"discrMVA<0.72",  "incl":"true", "":"true"},
    "VLcat":   {"bdt3":"discrMVA>=0.80", "bdt2":"discrMVA>=0.68 && discrMVA<0.80",
                "bdt1":"discrMVA>=0.34 && discrMVA<0.68",
                "bdt0":"discrMVA<0.34",  "incl":"true", "":"true"},
    "Zinvcat": {"bdt2":"discrMVA>=0.98", "bdt1":"discrMVA>=0.72 && discrMVA<0.98",
                "bdt0":"discrMVA<0.72",  "incl":"true", "":"true"},
    "VHcat":   {"bdt2":"discrMVA>=0.90", "bdt1":"discrMVA>=0.50 && discrMVA<0.90",
                "bdt0":"discrMVA<0.50",  "incl":"true", "":"true"},
    "TTHcat":  {"bdt2":"discrMVA>=0.94", "bdt1":"discrMVA>=0.78 && discrMVA<0.94",
                "bdt0":"discrMVA<0.78",  "incl":"true", "":"true"},
    # TODO: thresholds NOT synced. analysis/FITS has bdt2 >= 0.58 and
    # bdt1 "discrMVA>=0.58 && discrMVA<0.36", which is empty -- keep the old
    # two-bin split until those are fixed there.
    "TTLcat":  {"bdt1":"discrMVA>=0.54", "bdt0":"discrMVA<0.54",
                "incl":"true", "":"true"},
}
INCLUSIVE = {"incl": "true", "": "true"}


def read_bins_file(path=BINS_FILE):
    """{category: entry} of the bins file, {} if there is none."""
    if not os.path.isfile(path):
        return {}
    with open(path) as fh:
        return json.load(fh)


def _load():
    selmva = {c: dict(v) for c, v in SELMVA_DEFAULT.items()}
    source = {c: "SELMVA_DEFAULT (fitlib/bdtBins.py)" for c in selmva}
    for cat, entry in read_bins_file().items():
        selmva[cat] = {**entry["bins"], **INCLUSIVE}
        source[cat] = (f"fitlib/bdt_bins.json ({len(entry['bins'])} bins, "
                       f"applied {entry.get('applied_utc', '?')})")
    return selmva, source


SELMVA, SELMVA_SOURCE = _load()


def get_selection(category, binMVA):
    """Return the RDataFrame filter string for a (category, BDT bin)."""
    if category not in SELMVA:
        raise ValueError(f"Unknown category {category}")
    cuts = SELMVA[category]
    if binMVA not in cuts:
        raise ValueError(f"Unknown BDT bin '{binMVA}' for category {category} "
                         f"(have: {sorted(cuts)}; source "
                         f"{SELMVA_SOURCE[category]})")
    return cuts[binMVA]


def selection_record(category, binMVA):
    """The provenance fields a fit output records."""
    return {"selection": get_selection(category, binMVA or ""),
            "selection_source": SELMVA_SOURCE.get(category)}


def check_selection(prov, category, binMVA, what):
    """Stop if `what` (a fit output's provenance) was made with another cut
    than the one in force. Files written before the cuts were recorded give
    a warning only, since they cannot be checked."""
    cur = get_selection(category, binMVA or "")
    rec = prov.get("selection")
    if rec is None:
        print(f"warning: {what} does not record its BDT cut (made before "
              f"the cuts were recorded); cannot check it against {cur!r}")
        return
    if rec != cur:
        raise SystemExit(
            f"{what} was made with the cut {rec!r}, but {category} "
            f"{binMVA} is now {cur!r} ({SELMVA_SOURCE.get(category)}).\n"
            "Rerun the chain for the new cuts (in a new --wsdir), or restore "
            "the cuts it was made with.")
