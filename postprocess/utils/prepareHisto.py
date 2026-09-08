"""
prepareHisto.py -- RDataFrame booking for the HmumuRun3 histogramming stage.

Here, book_variable() only books. The single event loop is triggered by
makeHistos.py via RDF.RunGraphs, and the post-processing (overflow folding,
negative-bin removal) happens afterwards on plain TH1Ds.
"""

import ROOT

from . import plot_vars
from . import plot_style
from . import histo_config as cfg

plot_style.setup_style()

# lumis kept importable from here for back-compat; canonical copy is in plot_style
lumis = plot_style.lumis

RDataFrame = ROOT.RDataFrame

ROOT.gInterpreter.Declare("""
float deltaPhi(float phi1, float phi2) {
    float dphi = phi1 - phi2;
    while (dphi >  M_PI) dphi -= 2*M_PI;
    while (dphi <= -M_PI) dphi += 2*M_PI;
    return abs(dphi);
}
""")


# ===========================================================================
#  Post-loop histogram fixes (operate on plain TH1, after RunGraphs)
# ===========================================================================
def removeNegBin(hist):
    """Zero out negative bins (can arise from negative MC weights)."""
    for ibin in range(1, hist.GetNbinsX() + 1):
        if hist.GetBinContent(ibin) < 0.0:
            hist.SetBinContent(ibin, 0.0)
            hist.SetBinError(ibin, 0.0)


def addOverflow(h, addUnderflow=False):
    """Fold the overflow bin into the last visible bin (and optionally the
    underflow into the first), then clear them."""
    nb = h.GetNbinsX()

    h.SetBinContent(nb, h.GetBinContent(nb) + h.GetBinContent(nb + 1))
    h.SetBinError(nb, (h.GetBinError(nb)**2 + h.GetBinError(nb + 1)**2)**0.5)

    if addUnderflow:
        h.SetBinContent(1, h.GetBinContent(1) + h.GetBinContent(0))
        h.SetBinError(1, (h.GetBinError(1)**2 + h.GetBinError(0)**2)**0.5)

    h.SetBinContent(nb + 1, 0.)
    h.SetBinError(nb + 1, 0.)

    if addUnderflow:
        h.SetBinContent(0, 0.)
        h.SetBinError(0, 0.)


def finalize(hist):
    """Standard post-loop treatment applied to every histogram before it is
    written: drop negative bins, then fold in the overflow."""
    removeNegBin(hist)
    addOverflow(hist)
    return hist


# ===========================================================================
#  Base node: everything region- and variable-independent
# ===========================================================================
def make_base_node(chain, year):
    """Build the shared RDataFrame node.

    Defines, ONCE for the whole run:
      procGroup -- integer process-group index (replaces the long
                   "mc==x || mc==y || ..." filter chains; see
                   histo_config.proc_group_expr)
      weight    -- w_allSF, with the DY pT reweighting and the year
                   normalisation folded in

    Doing this once instead of per-plot is most of the JIT saving: Cling used
    to recompile ~12 filter chains plus the Defines on every getHisto() call.
    """
    df = RDataFrame(chain)

    # DY pT reweighting (DYturbo), applied only to the DY samples that carry it
    ggHcorr = "(mc==100 || mc==103 || mc==104 || mc==109) ? boson_ptWeight : 1."

    # 2025/2026 snapshots are made with 2024 MC
    if year in ("_2025", "_2026"):
        ratio = plot_style.lumis[year] / plot_style.lumis["_2024"]
        norm = f"(mc>0) ? {ratio} : 1."
    else:
        norm = "1."

    return (df.Define("procGroup", cfg.proc_group_expr())
              .Define("weight", f"w_allSF * ({ggHcorr}) * ({norm})"))


# ===========================================================================
#  Per-variable booking
# ===========================================================================
def book_variable(node, category, varname, binning, blind=True):
    """Book one histogram per process for `varname` on `node`.

    Returns {process_name: RResultPtr[TH1D]}, or None if the variable's
    expression cannot be built (typically a branch absent from this
    category's snapshot -- e.g. discrMVA0 when the MVA classifier was not
    run, or nGoodJetsAll outside isGGH/isTThad).

    Does NOT trigger the event loop.
    """
    nbin, low, high = binning
    var = plot_vars.get_expr(varname)

    # Let ROOT validate the expression here. If it references a missing
    # branch the JIT fails at Define time; catch it and skip this variable
    # cleanly rather than aborting the whole run.
    #
    # NOTE: catch BaseException -- PyROOT JIT failures don't always surface as
    # a plain Exception subclass, so a narrower `except Exception` can miss
    # them.
    try:
        node = node.Define(f"var_{varname}", var)
    except BaseException as e:
        print(f"⚠️  '{varname}': cannot build variable '{var}' "
              f"(likely a missing branch: {type(e).__name__})")
        return None

    col = f"var_{varname}"
    ptrs = {}

    # --- MC ---
    for proc in cfg.MC_PROCESSES:
        idx = cfg.PROC_INDEX[proc]
        ptrs[proc] = (node
                      .Filter(f"procGroup=={idx}")
                      .Histo1D((f"{varname}_{proc}", "", nbin, low, high),
                               col, "weight"))

    # --- data, with optional blinding ---
    # Blinding only ever affects data; MC is never blinded. Three cases,
    # preserving the original getHisto() logic:
    #   1. this variable has a blind mass window  -> cut that window out
    #   2. varname == "mva"                       -> cut above the per-category
    #                                                MVA score threshold
    #   3. otherwise                              -> nothing to blind
    data_filter = f"procGroup=={cfg.DATA_INDEX}"
    if blind:
        blind_range = plot_vars.get_blind_range(varname)
        if blind_range is not None:
            lo, hi = blind_range
            data_filter += f" && ({col}<{lo} || {col}>{hi})"
        elif varname == "mva":
            cut = _mva_blind_cut(category)
            if cut is not None:
                data_filter += f" && {col}<{cut}"
    else:
        if plot_vars.get_blind_range(varname) is not None:
            print(f"[book_variable] UNBLINDED data for '{varname}'")

    ptrs[cfg.DATA_PROCESS] = (node
                              .Filter(data_filter)
                              .Histo1D((f"{varname}_{cfg.DATA_PROCESS}", "",
                                        nbin, low, high),
                                       col, "weight"))

    return ptrs


# Per-category MVA score above which data is blinded (the signal-rich tail).
_MVA_BLIND_CUT = {
    "ggHcat":  0.50,
    "VBFcat":  0.64,
    "VLcat":   0.32,
    "VHcat":   0.86,
    "TTHcat":  0.80,
    "Zinvcat": 0.78,
}


def _mva_blind_cut(category):
    return _MVA_BLIND_CUT.get(category)
