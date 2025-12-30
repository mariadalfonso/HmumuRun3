import ROOT
from array import array
import math

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
    '_Run3':286.5, #C-G
}

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

def getHisto(nbin, low, high, doLog, category, year, doSignal, binMVA, sig=''):

   print("getHisto getting called for year: ", year)

   year = '_'+str(year)
#   dirLOCAL_='/work/submit/mariadlf/HmumuRun3/ROOTFILES/'+category+'/WS/'
   dirLOCAL_='/work/submit/mariadlf/HmumuRun3/ROOTFILES/'+category+'/'
   files = []

   mytree = ROOT.TChain('events')   
   # Map category → file tags
   signal_files = {
         "ggH":  ["11"],             # ggH
         "qqH":  ["10"],             # VBF
         "VH":   ["12", "13", "14"], # VH
         "ttH":  ["15"],             # TTL
   }

   # Add files safely
   if sig != '':
       for tag in signal_files.get(sig, []):
           if year != '_Run3': path = f"{dirLOCAL_}snapshot_mc_{tag}{year}_{category}.root"
           else: path = f"{dirLOCAL_}snapshot_mc_{tag}_*_{category}.root"
           safe_add_tree(files, path)

   # Year-specific samples
   data_samples = {
       "_12022": ["-11", "-13", "-14"],
       "_22022": ["-15", "-16", "-17"],
       "_12023": ["-23", "-24"],
       "_22023": ["-31", "-32"],
       "_2024" : [str(i) for i in range(-41, -55, -1)],  # -41 … -54
       "_2025" : [str(i) for i in range(-61, -73, -1)],  # -61 … -72
       "_Run3" : ["-11", "-13", "-14", "-15", "-16", "-17", "-23", "-24", "-31", "-32"]+ [str(i) for i in range(-41, -55, -1)] + [str(i) for i in range(-61, -73, -1)],
   }

   # Add files safely
   for tag in data_samples.get(year, []):
         if year != '_Run3': path = f"{dirLOCAL_}snapshot_mc_{tag}{year}_{category}.root"
         else: path = f"{dirLOCAL_}snapshot_mc_{tag}_*_{category}.root"
         safe_add_tree(files, path)

   #-------------- selection

   selMVAvbf = {
       "bdt0" : "discrMVA0>=0.75",
       "bdt1" : "discrMVA0>=0.4 && discrMVA0<0.75",
       "bdt2" : "discrMVA0<0.4",
       "" : "true"
   }

   selMVAggh = {
       "bdt0" : "discrMVA0>=0.65",
       "bdt1" : "discrMVA0>=0.35 && discrMVA0<0.65",
       "bdt2" : "discrMVA0<0.35",
       "" : "true"
   }

   selection = {
       "VLcat":   "true",
       "ggHcat":  ""+selMVAggh[binMVA],
       "VBFcat":  "Mjj>200 && "+selMVAvbf[binMVA],
       "VHcat":   "goodWjj_discr>0.9",
       "Zinvcat": "true",
       "TTHcat":  "true",
       "TTLcat":  "HT>100",
   }
   print('sel=',selection[category])

   # --- Build the RDataFrame ---
   df = ROOT.RDataFrame("events", files)
   print(f"✅ Loaded {df.Count().GetValue()} entries from {len(files)} files")

   df = df.Filter("{}".format(selection[category]),"selection cut")

   # --- Filters ---
   # Sanity check: skip NaN masses
   df = df.Filter("!isnan(HiggsCandCorrMass)", "Valid mass")

   # Select signal or background
   if doSignal:
         df_sel = df.Filter("mc == 10 || mc == 11 || mc == 12 || mc == 13 || mc == 14 || mc == 15", "Signal (ggH, VBF, VH, ttH)")
         df_sel = df_sel.Define("weight", "w_allSF * lumiIntegrated")
   else:
         df_sel = df.Filter("mc < 0 && mc > -100 ", "data ")
         df_sel = df_sel.Define("weight", "w_allSF")

   # --- Histogram creation ---
   h = df_sel.Histo1D(
         (f"h_{category}_{year}", f"{category} {year}", nbin, low, high),
         "HiggsCandCorrMass",
         "weight"
   )

   # Convert to TH1 for later ROOT use
   h.GetValue().SetDirectory(0)
   return h.GetValue()
