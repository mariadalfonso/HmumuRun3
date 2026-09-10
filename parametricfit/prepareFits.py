import ROOT
from array import array
import math
import getpass

ROOTFILES_DIR = f'/work/submit/{getpass.getuser()}/HmumuRun3/ROOTFILES'

ROOT.ROOT.EnableImplicitMT()
RDataFrame = ROOT.RDataFrame
ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()

lumis={
    '_12016': 19.52, #APV #(B-F for 2016 pre)
    '_22016': 16.80, #postVFP
    '_2016': 35.9,
    '_2017': 41.5,
    '_12017': 7.7, #(F for 2017) for VBF
    '_2018': 59.70,
    '_12018': 39.54,
    '_all': 86.92,      #19.52 + 7.7 + 59.70
    '_Run2': 138.,      #19.52 + 7.7 + 59.70
      #RUN3
    '_12022':7.98, # C-D
    '_22022':26.67, # E, F, G
    '_12023':17.794, #C
    '_22023':9.451, #D
    '_2024':108.95, #C-I
    '_2025':115.65, #C-G
    '_2026':25.8, #from Gluillelmo
    '_Run3':312, #C-G
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
    "_12022": ["-11", "-13", "-14"],
    "_22022": ["-15", "-16", "-17"],
    "_12023": ["-23", "-24"],
    "_22023": ["-31", "-32"],
    "_2024" : [str(i) for i in range(-41, -55, -1)],  # -41 ... -54
    "_2025" : [str(i) for i in range(-61, -73, -1)],  # -61 ... -72
    "_2026" : [str(i) for i in range(-81, -89, -1)],  # -81 ... -88
    "_Run3" : ["-11", "-13", "-14", "-15", "-16", "-17", "-23", "-24", "-31", "-32"]
              + [str(i) for i in range(-41, -55, -1)]
              + [str(i) for i in range(-61, -73, -1)]
              + [str(i) for i in range(-81, -89, -1)],
}


def get_selection(category, binMVA):
    """Return the RDataFrame filter string for a (category, BDT bin)."""
    if category not in SELMVA:
        raise ValueError(f"Unknown category {category}")
    cuts = SELMVA[category]
    if binMVA not in cuts:
        raise ValueError(f"Unknown BDT bin '{binMVA}' for category {category} "
                         f"(have: {sorted(cuts)})")
    return cuts[binMVA]


def _snapshot_path(dirLOCAL_, tag, year, category):
    """Build one snapshot glob/path for a given file tag."""
    if year != '_Run3':
        return f"{dirLOCAL_}snapshot_mc_{tag}{year}_{category}.root"
    return f"{dirLOCAL_}snapshot_mc_{tag}_*_{category}.root"


def get_signal_files(sig, category, year, rootfiles_dir=ROOTFILES_DIR):
    """Return the list of signal snapshot files for one production mode."""
    dirLOCAL_ = f'{rootfiles_dir}/{category}/'
    files = []
    for tag in SIGNAL_TAGS.get(sig, []):
        safe_add_tree(files, _snapshot_path(dirLOCAL_, tag, year, category))
    return files


def get_data_files(year, category, rootfiles_dir=ROOTFILES_DIR):
    """Return the list of data snapshot files for a year."""
    dirLOCAL_ = f'{rootfiles_dir}/{category}/'
    files = []
    for tag in DATA_TAGS.get(year, []):
        safe_add_tree(files, _snapshot_path(dirLOCAL_, tag, year, category))
    return files


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

    year = '_' + str(year)
    files = get_signal_files(sig, category, year, rootfiles_dir)
    if not files:
        print(f"⚠️ no signal files for {sig} {category} {year}")
        return None

    df = ROOT.RDataFrame("events", files)
    print(f"✅ Loaded {df.Count().GetValue()} entries from {len(files)} files")

    df = df.Filter(get_selection(category, binMVA), "selection cut")
    df = df.Filter("!isnan(HiggsCandCorrMass)", "Valid mass")
    df = df.Filter("mc >= 10 && mc <= 15", "Signal (ggH, VBF, VH, ttH)")

    return fill_mass_histo(df, f"h_{category}{year}", f"{category} {year}",
                           nbin, low, high)


def getHisto(nbin, low, high, doLog, category, year, doSignal, binMVA, sig='',
             rootfiles_dir=ROOTFILES_DIR):
    """Mass histogram for signal MC (doSignal=True) or data (doSignal=False).

    Only the file group that will be kept is loaded: signal fits do not read
    the data snapshots and vice versa.
    """

    print("getHisto getting called for year: ", year)

    year = '_'+str(year)

    # Build only the file group we will actually keep.
    files = []
    if doSignal:
        if sig != '':
            files += get_signal_files(sig, category, year, rootfiles_dir)
    else:
        files += get_data_files(year, category, rootfiles_dir)

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
