import ROOT
from array import array
import math
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
    "ggHcat":  {"bdt2":"discrMVA0>=0.76", "bdt1":"discrMVA0>=0.5 && discrMVA0<0.76",
                "bdt0":"discrMVA0<0.5",   "incl":"true", "":"true"},
    "VBFcat":  {"bdt2":"discrMVA0>=0.92", "bdt1":"discrMVA0>=0.64 && discrMVA0<0.92",
                "bdt0":"discrMVA0<0.64",  "incl":"true", "":"true"},
    "VLcat":   {"bdt2":"discrMVA0>=0.78", "bdt1":"discrMVA0>=0.32 && discrMVA0<0.78",
                "bdt0":"discrMVA0<0.32",  "incl":"true", "":"true"},
    "Zinvcat": {"bdt2":"discrMVA0>=0.98", "bdt1":"discrMVA0>=0.78 && discrMVA0<0.98",
                "bdt0":"discrMVA0<0.78",  "incl":"true", "":"true"},
    "VHcat":   {"bdt2":"discrMVA0>=0.94", "bdt1":"discrMVA0>=0.86 && discrMVA0<0.94",
                "bdt0":"discrMVA0<0.86",  "incl":"true", "":"true"},
    "TTHcat":  {"bdt2":"discrMVA0>=0.98", "bdt1":"discrMVA0>=0.8 && discrMVA0<0.98",
                "bdt0":"discrMVA0<0.8",   "incl":"true", "":"true"},
    "TTLcat":  {"bdt1":"discrMVA0>=0.54", "bdt0":"discrMVA0<0.54",
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


def normalize_year(year):
    """'Run3', 2022, '_2022' -> 'Run3', '2022', '2022' (bare, idempotent)."""
    return str(year).lstrip('_')


def expand_year(year, groups):
    """Return the list of concrete eras for a requested year."""
    year = normalize_year(year)
    return list(groups.get(year, [year]))


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
    eras = expand_year(year, SIGNAL_YEAR_GROUPS)
    return _collect_files([(era, SIGNAL_TAGS[sig]) for era in eras],
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
