import ROOT
from array import array
import math
import os
import getpass

ROOTFILES_DIR = f'/work/submit/{getpass.getuser()}/HmumuRun3/ROOTFILES'

ROOT.ROOT.EnableImplicitMT()
RDataFrame = ROOT.RDataFrame
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()

# Integrated luminosity in fb^-1, keyed by concrete era (bare, no underscore).
# Composite years (2022, 2023, Run3) are NOT listed: get_lumi() sums their
# eras via the same year groups used for the file lookup.
# TODO: check the correct Run 3 lumi
LUMIS = {
    # Run 2
    '12016': 19.52,  # APV (B-F for 2016 pre)
    '22016': 16.80,  # postVFP
    '2016': 35.9,
    '2017': 41.5,
    '12017': 7.7,    # (F for 2017) for VBF
    '2018': 59.70,
    '12018': 39.54,
    'all': 86.92,    # 19.52 + 7.7 + 59.70
    'Run2': 138.,
    # Run 3
    '12022': 7.99,   # C-D
    '22022': 26.68,  # E, F, G
    '12023': 17.96,  # C
    '22023': 9.68,   # D
    '2024': 109.82,  # C-I
    '2025': 110.59,  # C-G
    '2026': 25.31,   # C, B
}

# ---------------------------------------------------------------------------
# shared configuration and helpers
# ---------------------------------------------------------------------------

# BDT-bin selection thresholds, keyed by category then bin label.
# Single source of truth: both getHisto and getHistoSignal read from here.
# VLcat: 1 W->e; 2 Z->ee; 3 W->mu; 4 Z->mumu
SELMVA = {
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

# production mode -> file tag(s)
SIGNAL_TAGS = {
    "ggH": ["11"],              # ggH
    "qqH": ["10"],              # VBF
    "VH":  ["12", "13", "14"],  # VH
    "ttH": ["15"],              # TTL
}

# data file tags, keyed by year
DATA_TAGS = {
    "12022": ["-11", "-13", "-14"],
    "22022": ["-15", "-16", "-17"],
    "12023": ["-23", "-24"],
    "22023": ["-31", "-32"],
    "2024" : [str(i) for i in range(-41, -55, -1)],  # -41 ... -54
    "2025" : [str(i) for i in range(-61, -73, -1)],  # -61 ... -72
    "2026" : [str(i) for i in range(-81, -89, -1)],  # -81 ... -88
}

# background MC for the frozen-core fit (mcCoreFit.py --all-bkg-mc), file
# tags per MC era ("*" = every era without its own entry); datasets as in
# analysis/config/samples.yaml. Non-overlapping samples only, so the
# weighted sum is the background. The Drell-Yan samples carry a generator
# filter 105 < m_mumu < 160 GeV, around the 110-150 window, so nearly all
# their events are in it (the MLL-50 samples spend most on the Z peak);
# the 103 vs 109 shapes agree across 110-150, the 160 edge only shows above
# 150. NOT used, because they overlap the samples below:
#   DY QCD  100, 103 (MLL-50 inclusive), 108, 110 (powheg MLL-50to120 /
#           120to200), 111-113, 119-121 (jet-binned), 114-117, 122-125
#           (pT-binned), 126-128 (MiNNLO: never used for DY)
#   DY EWK  VBF: 101 (inclusive 2L2J); ggH: 98 (not in the ggH lists);
#           99 (M2Mu-105to160 with M2J > 120, pythia8)
#   tt      141-149 (other generators, tune / hdamp variations)
#   VV      207-212 (amcatnlo, overlap the powheg set)
_DY_QCD = {
    # DYto2Mu-2Jets MLL-105to160 + DYto2Tau-2Jets MLL-50, amcatnloFXFX
    "2024": ["109", "104"],
    # DYto2Mu-2Jets MLL-105to160, amcatnloFXFX; no tautau sample for
    # 2022-23 in the lists -- negligible at 110-150
    "*":    ["109"],
}
_DY_EWK = {
    "ggHcat": ["101"],          # EWK Z 2L2J, madgraph-herwig7
    "VBFcat": ["98"],           # EWK Z 2Mu2J MLL-105to160, madgraph-herwig7
}
_BKG_OTHER = [
    "140",                      # TTto2L2Nu, powheg
    "107",                      # TTtoLNu2Q, powheg
    "105", "106",               # TWminus / TbarWplus to2L2Nu, powheg
    "201", "202", "203",        # WZto2L2Q, ZZto2L2Nu, ZZto2L2Q, powheg
    "204", "205", "206",        # WZto3LNu, ZZto4L, WWto2L2Nu, powheg
]
BKG_MC_TAGS = {
    cat: {era: dy + ewk + _BKG_OTHER for era, dy in _DY_QCD.items()}
    for cat, ewk in _DY_EWK.items()
}

# the Drell-Yan subset, mcCoreFit.py's default. VBF keeps the EWK Z
# (vector-boson fusion Z, the irreducible VBF-like DY); in ggH it is
# negligible and left with the rest.
BKG_DY_TAGS = {
    "ggHcat": dict(_DY_QCD),
    "VBFcat": {era: dy + _DY_EWK["VBFcat"] for era, dy in _DY_QCD.items()},
}

# Composite years expand into the concrete eras whose snapshots exist on disk.
# Any year not listed here is used as-is (e.g. '2024' -> ['2024']).
RUN3_SIGNAL_ERAS = ["12022", "22022", "12023", "22023", "2024"]
RUN3_DATA_ERAS   = RUN3_SIGNAL_ERAS + ["2025", "2026"]

_COMMON_YEAR_GROUPS = {
    "2022": ["12022", "22022"],
    "2023": ["12023", "22023"],
}
SIGNAL_YEAR_GROUPS = {**_COMMON_YEAR_GROUPS, "Run3": RUN3_SIGNAL_ERAS}
DATA_YEAR_GROUPS   = {**_COMMON_YEAR_GROUPS, "Run3": RUN3_DATA_ERAS}

# Eras with no signal MC of their own, mapped to the MC era whose signal yield
# per fb^-1 they borrow (same detector configuration and NanoAOD version).
# Only the YIELD is extrapolated: the signal shape is fitted on the MC eras.
SIGNAL_PROXY = {"2025": "2024", "2026": "2024"}


def normalize_year(year):
    """'Run3', 2022, '_2022' -> 'Run3', '2022', '2022' (bare, idempotent)."""
    return str(year).lstrip('_')


def expand_year(year, groups):
    """Return the list of concrete eras for a requested year."""
    year = normalize_year(year)
    return list(groups.get(year, [year]))


def signal_yield_scales(year):
    """{MC era: factor} that turns each MC era's signal yield into its share
    of the yield for the FULL requested year, i.e. the same eras as the data.

    An era with its own MC contributes 1 to itself; an era in SIGNAL_PROXY
    adds L_era / L_proxy to its proxy. So Run3 -> the 2022 and 2023 eras x 1,
    and 2024 x (L2024 + L2025 + L2026) / L2024; '2025' -> 2024 x L2025/L2024.
    """
    scales = {}
    for era in expand_year(year, DATA_YEAR_GROUPS):
        mc = SIGNAL_PROXY.get(era, era)
        missing = [e for e in (era, mc) if e not in LUMIS]
        if missing:
            raise ValueError(f"no luminosity for era(s) {missing} "
                             f"(requested {normalize_year(year)})")
        scales[mc] = scales.get(mc, 0.0) + LUMIS[era] / LUMIS[mc]
    return scales


def get_lumi(year, groups=SIGNAL_YEAR_GROUPS):
    """Integrated luminosity (fb^-1) for a requested year.

    Summed over the same eras the file lookup uses, so the label always
    matches the files that were read. Pass DATA_YEAR_GROUPS for data.
    Returns 0.0 if any era is missing from LUMIS.
    """
    year = normalize_year(year)
    if year in LUMIS and year not in groups:
        return LUMIS[year]
    eras = expand_year(year, groups)
    missing = [e for e in eras if e not in LUMIS]
    if missing:
        print(f"⚠️ no luminosity for era(s) {missing} (requested {year})")
        return 0.0
    return sum(LUMIS[e] for e in eras)


def get_selection(category, binMVA):
    """Return the RDataFrame filter string for a (category, BDT bin)."""
    if category not in SELMVA:
        raise ValueError(f"Unknown category {category}")
    cuts = SELMVA[category]
    if binMVA not in cuts:
        raise ValueError(f"Unknown BDT bin '{binMVA}' for category {category} "
                         f"(have: {sorted(cuts)})")
    return cuts[binMVA]


def _snapshot_path(dirLOCAL_, tag, era, category):
    """Build the snapshot path for one file tag in one concrete era."""
    return f"{dirLOCAL_}snapshot_mc_{tag}_{era}_{category}.root"


def _collect_files(tags_per_era, category, rootfiles_dir):
    """tags_per_era: iterable of (era, [tags]).

    Strict: every (era, tag) must yield at least one readable file with the
    tree. Otherwise raise, so a partially loaded year can never be fit and
    labelled with the full luminosity.
    """
    dirLOCAL_ = f'{rootfiles_dir}/{category}/'
    files, missing = [], []
    for era, tags in tags_per_era:
        for tag in tags:
            path = _snapshot_path(dirLOCAL_, tag, era, category)
            n_before = len(files)
            safe_add_tree(files, path)
            if len(files) == n_before:
                missing.append(path)
    if missing:
        raise FileNotFoundError(
            f"{len(missing)} required snapshot(s) missing or without tree:\n  "
            + "\n  ".join(missing))
    return files


def get_signal_files(sig, category, year, rootfiles_dir=ROOTFILES_DIR):
    """Return the signal snapshot files for one production mode.

    year may be a concrete era ('12022', '2024') or a group
    ('2022', '2023', 'Run3'); see SIGNAL_YEAR_GROUPS.
    """
    if sig not in SIGNAL_TAGS:
        raise ValueError(f"Unknown signal '{sig}' (have: {sorted(SIGNAL_TAGS)})")
    eras = list(signal_yield_scales(year))   # MC eras; proxied eras -> 2024
    return _collect_files([(era, SIGNAL_TAGS[sig]) for era in eras],
                          category, rootfiles_dir)


def bkg_mc_tags(category, era, dy_only=False):
    """Background MC sample ids of ONE concrete MC era; see BKG_MC_TAGS.

    dy_only: the Drell-Yan subset instead, BKG_DY_TAGS.
    """
    table, name = ((BKG_DY_TAGS, "BKG_DY_TAGS") if dy_only
                   else (BKG_MC_TAGS, "BKG_MC_TAGS"))
    if category not in table:
        raise ValueError(f"no background MC list for {category} "
                         f"(have {sorted(table)}); add it to "
                         f"{name} in prepareFits.py")
    per_era = table[category]
    return list(per_era.get(era, per_era["*"]))


def get_bkg_mc_files(category, era, rootfiles_dir=ROOTFILES_DIR,
                     dy_only=False):
    """Background MC snapshot files of ONE concrete MC era."""
    return _collect_files([(era, bkg_mc_tags(category, era, dy_only))],
                          category, rootfiles_dir)


def get_data_files(year, category, rootfiles_dir=ROOTFILES_DIR):
    """Return the data snapshot files for a year.

    year may be a concrete era or a group; see DATA_YEAR_GROUPS.
    """
    eras = expand_year(year, DATA_YEAR_GROUPS)
    missing = [era for era in eras if era not in DATA_TAGS]
    if missing:
        raise ValueError(f"No DATA_TAGS for era(s) {missing} "
                         f"(requested {normalize_year(year)})")
    return _collect_files([(era, DATA_TAGS[era]) for era in eras],
                          category, rootfiles_dir)


def fill_mass_histo(df, name, title, nbin, low, high):
    """Define the weight and fill the HiggsCandCorrMass histogram.

    Returns a detached TH1 (SetDirectory(0)) so it survives past the RDF.
    """
    df = df.Define("weight", "w_allSF")  # lumiIntegrated is in the weight already
    h = df.Histo1D((name, title, nbin, low, high), "HiggsCandCorrMass", "weight")
    h.GetValue().SetDirectory(0)
    return h.GetValue()


def safe_add_tree(file_list, filepath_pattern, treename="events"):
    """Safely add ROOT files to a list if they contain the TTree."""
    import os
    import glob

    #print('to search filepath = ',filepath)
    files = glob.glob(filepath_pattern)

    if not files:
        print(f"⚠️ File not found (glob): {filepath_pattern}")
        return False

    n_added = 0

    for filepath in sorted(files):
        f = ROOT.TFile.Open(filepath)
        if not f or f.IsZombie():
            print(f"⚠️ Could not open: {filepath}")
            continue

        if not f.GetListOfKeys().Contains(treename):
            print(f"⚠️ No TTree '{treename}' in {filepath}")
            f.Close()
            continue

        f.Close()
        file_list.append(filepath)
        n_added += 1
        print(f"✅ Added {filepath}")

    return True

def getFWHM(h):

      max_bin = h.GetMaximumBin()
      maxP = h.GetBinContent(max_bin)
      x_maximum = h.GetBinCenter(max_bin)

      print('x_maximum = ',x_maximum)
      
      half_max = maxP / 2

      # Search left for half max
      fwhm_left_bin = max_bin
      while fwhm_left_bin > 1 and h.GetBinContent(fwhm_left_bin) > half_max:
         fwhm_left_bin -= 1
         
      # Search right for half max
      fwhm_right_bin = max_bin
      while fwhm_right_bin < h.GetNbinsX() and h.GetBinContent(fwhm_right_bin) > half_max:
         fwhm_right_bin += 1
            
      # Convert bin to x value
      fwhm_left_x = h.GetBinCenter(fwhm_left_bin)
      fwhm_right_x = h.GetBinCenter(fwhm_right_bin)
      
      # Calculate width and integral
      width = fwhm_right_x - fwhm_left_x
      myRange = h.Integral(fwhm_left_bin, fwhm_right_bin)

      print("FWHM width =", width, "Integral in range =", myRange)

      return fwhm_left_x,fwhm_right_x


def getHistoSignal(nbin, low, high, category, year, binMVA, sig,
                   rootfiles_dir=ROOTFILES_DIR):
    """Signal-only mass histogram for one (category, BDT bin, production mode).

    signal is never blinded and never region-restricted. 
    Returns a detached TH1, or None if no signal files were found.
    """
    print("getHistoSignal getting called for year: ", year)

    year = normalize_year(year)
    files = get_signal_files(sig, category, year, rootfiles_dir)
    if not files:
        print(f"⚠️ no signal files for {sig} {category} {year}")
        return None

    df = ROOT.RDataFrame("events", files)
    print(f"✅ Loaded {df.Count().GetValue()} entries from {len(files)} files")

    df = df.Filter(get_selection(category, binMVA), "selection cut")
    df = df.Filter("!isnan(HiggsCandCorrMass)", "Valid mass")
    df = df.Filter("mc >= 10 && mc <= 15", "Signal (ggH, VBF, VH, ttH)")

    return fill_mass_histo(df, f"h_{category}_{year}", f"{category} {year}",
                           nbin, low, high)


def getHistoSignalByEra(nbin, low, high, category, year, binMVA, sig,
                        rootfiles_dir=ROOTFILES_DIR):
    """Signal mass histograms for one (category, BDT bin, production mode),
    one per MC era of signal_yield_scales(year), UNSCALED.

    Their sum is the MC the shape is fitted on; each era's integral times
    its scale gives the signal yield of the full requested year.
    Returns {MC era: detached TH1}.
    """
    year = normalize_year(year)
    if sig not in SIGNAL_TAGS:
        raise ValueError(f"Unknown signal '{sig}' (have: {sorted(SIGNAL_TAGS)})")
    out = {}
    for era in signal_yield_scales(year):
        files = _collect_files([(era, SIGNAL_TAGS[sig])], category,
                               rootfiles_dir)
        df = (ROOT.RDataFrame("events", files)
              .Filter(get_selection(category, binMVA), "selection cut")
              .Filter("!isnan(HiggsCandCorrMass)", "Valid mass")
              .Filter("mc >= 10 && mc <= 15", "Signal (ggH, VBF, VH, ttH)"))
        out[era] = fill_mass_histo(df, f"h_{category}_{binMVA}_{sig}_{era}",
                                   f"{category} {era}", nbin, low, high)
        print(f"✅ {sig} {era}: {len(files)} file(s)")
    return out


def getHisto(nbin, low, high, doLog, category, year, doSignal, binMVA, sig='',
             rootfiles_dir=ROOTFILES_DIR):
    """Mass histogram for signal MC (doSignal=True) or data (doSignal=False).

    Only the file group that will be kept is loaded: signal fits do not read
    the data snapshots and vice versa.
    """

    print("getHisto getting called for year: ", year)

    year = normalize_year(year)

    # Build only the file group we will actually keep.
    files = []
    if doSignal:
        if sig != '':
            files += get_signal_files(sig, category, year, rootfiles_dir)
    else:
        files += get_data_files(year, category, rootfiles_dir)

    if not files:
        kind = f"signal {sig}" if doSignal else "data"
        print(f"⚠️ no {kind} files for {category} {year}")
        return None

    selection_cut = get_selection(category, binMVA)

    # --- Build the RDataFrame ---
    df = ROOT.RDataFrame("events", files)
    print(f"✅ Loaded {df.Count().GetValue()} entries from {len(files)} files")

    df = df.Filter("{}".format(selection_cut), "selection cut")

    # --- Filters ---
    # Sanity check: skip NaN masses
    df = df.Filter("!isnan(HiggsCandCorrMass)", "Valid mass")

    # Select signal or background
    if doSignal:
        df_sel = df.Filter("mc == 10 || mc == 11 || mc == 12 || mc == 13 || mc == 14 || mc == 15", "Signal (ggH, VBF, VH, ttH)")
    else:
        df_sel = df.Filter("mc < 0 && mc > -100 ", "data ")

    # --- Histogram creation ---
    return fill_mass_histo(df_sel, f"h_{category}_{year}", f"{category} {year}",
                           nbin, low, high)


def _require_branch(files, branch, category, treename="events"):
    """Raise, naming every file whose tree lacks `branch`.

    RDataFrame reads the column list off the first file only, so without
    this a stale snapshot surfaces as a cling 'undeclared identifier'.
    """
    stale = []
    for path in files:
        f = ROOT.TFile.Open(path)
        t = f.Get(treename)
        if not t.GetBranch(branch):
            stale.append(path)
        f.Close()
    if stale:
        tags = sorted({os.path.basename(p).split("_")[2] for p in stale})
        raise RuntimeError(
            f"{len(stale)} of {len(files)} background-MC snapshot(s) have no "
            f"'{branch}' branch -- written before the {category} classifier; "
            f"rerun Hmm.py on them (sample ids {','.join(tags)}):\n  "
            + "\n  ".join(stale))


def getHistoBkgMC(nbin, low, high, category, year, binMVA,
                  rootfiles_dir=ROOTFILES_DIR, dy_only=False):
    """Weighted background-MC mass histogram, for the frozen-core fit.

    Summed over the MC eras of `year` with the same factors as the signal
    (signal_yield_scales): an era without MC of its own (2025, 2026) is
    represented by its proxy, so the era mix follows the data. Only the
    SHAPE is used downstream. dy_only: the Drell-Yan samples alone.
    Returns (TH1, {era: (factor, n_files)}).
    """
    year = normalize_year(year)
    sel = get_selection(category, binMVA)
    per_era = {era: get_bkg_mc_files(category, era, rootfiles_dir, dy_only)
               for era in signal_yield_scales(year)}
    if "discrMVA" in sel:
        _require_branch([f for fs in per_era.values() for f in fs],
                        "discrMVA", category)
    total, info = None, {}
    for era, factor in signal_yield_scales(year).items():
        files = per_era[era]
        df = ROOT.RDataFrame("events", files)
        df = df.Filter(sel, "selection cut")
        df = df.Filter("!isnan(HiggsCandCorrMass)", "Valid mass")
        df = df.Filter("mc > 15", "background MC")
        h = fill_mass_histo(df, f"h_bkgmc_{category}_{binMVA}_{era}",
                            f"{category} {era} background MC",
                            nbin, low, high)
        h.Scale(factor)
        print(f"✅ background MC {era}: {len(files)} file(s), "
              f"x {factor:.3f}, sum w = {h.Integral():.1f}")
        info[era] = (factor, len(files))
        if total is None:
            total = h.Clone(f"h_bkgmc_{category}_{binMVA}_{year}")
            total.SetDirectory(0)
        else:
            total.Add(h)
    return total, info
