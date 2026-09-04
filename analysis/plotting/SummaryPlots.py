import ROOT
import os
import sys
import argparse
import getpass
from array import array
import math

from LoadTree import loadTree
from prepareHisto import getHisto
import plot_style
from plot_style import lumis
import plot_vars

plot_style.setup_style()

# ---------------------------------------------------------------------------
# Configuration / CLI
# ---------------------------------------------------------------------------

CATEGORIES = ["VBFcat", "ggHcat", "VLcat", "TTLcat", "TTHcat",
              "VHcat", "Zinvcat"]

# Which plot groups each category can draw. "mass" always runs.
# Names map to the plot*() functions / plotMuons() below.
DEFAULT_INDIR  = f"/work/submit/{getpass.getuser()}/HmumuRun3/ROOTFILES/"
DEFAULT_OUTDIR = f"/home/submit/{getpass.getuser()}/public_html/HmumuRun3/"

# Plot groups the user can request via --plots. "all" expands to the
# category-appropriate set (mass + the category's own function + muons).
PLOT_GROUPS = ["mass", "mva", "category", "muons", "all"]


def parse_args():
    p = argparse.ArgumentParser(
        description="Make stacked data/MC comparison plots from Hmm.py snapshots."
    )
    p.add_argument("category", choices=CATEGORIES,
                   help="analysis category (snapshot subfolder)")
    p.add_argument("year",
                   help="data-taking year, e.g. 2024, 2025, 12022 ...")
    p.add_argument("-i", "--indir", default=DEFAULT_INDIR,
                   help="input directory with snapshot ROOT files "
                        "(default: %(default)s)")
    p.add_argument("-o", "--outdir", default=DEFAULT_OUTDIR,
                   help="output directory for PNG plots (default: %(default)s)")
    p.add_argument("--plots", nargs="+", default=["all"], choices=PLOT_GROUPS,
                   help="which plot groups to draw (default: all)")
    p.add_argument("--blind", dest="blind", action="store_true", default=True,
                   help="blind data in the signal region (default: on)")
    p.add_argument("--unblind", dest="blind", action="store_false",
                   help="disable blinding (use with care)")
    return p.parse_args()


args = parse_args()

category = args.category
year     = "_" + args.year          # internal convention: leading underscore
dirLOCAL_ = args.indir
myOutDir  = args.outdir
if not dirLOCAL_.endswith("/"): dirLOCAL_ += "/"
if not myOutDir.endswith("/"):  myOutDir  += "/"

# baseOutDir is the root the user gave; run_group() redirects myOutDir into a
# per-group subdirectory under it. Individual folders are created lazily in
# plot() right before each SaveAs.
baseOutDir = myOutDir
os.makedirs(baseOutDir, exist_ok=True)

# Guard: year must be known to the lumi table used in the CMS label.
if year not in lumis:
    sys.exit(f"ERROR: year '{args.year}' not found in lumis table "
             f"(known: {[k.lstrip('_') for k in lumis if k.startswith('_')]})")

mytree = ROOT.TChain('events')
mytree = loadTree(mytree, dirLOCAL_, category, year)


# Process groups (stacking order matters; legend order matches original).
BKG_PROCS = ["hTT2L", "hTop", "hZg", "hVV", "hEWK", "hDY"]
SIG_PROCS = ["hWH", "hTTH", "hZH", "hVBFH", "hggH"]

LEGEND_BKG = [
    ("hData", "Data",                        "lep"),
    ("hDY",   "DY+jets (QCD)",                "f"),
    ("hEWK",  "DY+jets (EWK)",                "f"),
    ("hVV",   "VV + VVV",                     "f"),
    ("hTT2L", "t#bar{t} 2l",                  "f"),
    ("hTop",  "Top (1l, tW/tZq, ttV/4t)",     "f"),
    ("hZg",   "H#rightarrowZ#gamma + jets",   "f"),
]
LEGEND_SIG_LEP = [("hTTH", "ttH", "l"), ("hWH", "WH", "l"), ("hZH", "ZH", "l")]
LEGEND_SIG_HAD = [("hVBFH", "VBF H", "l"), ("hggH", "ggH", "l")]
LEP_CATEGORIES = ["VLcat", "TTLcat", "TTHcat", "VHcat", "Zinvcat"]


def _style_ratio_axis(axis, offset):
   """Apply the shared absolute-pixel ratio-pad text style to one axis."""
   axis.SetTitleFont(plot_style.RATIO_FONT_ABS)
   axis.SetLabelFont(plot_style.RATIO_FONT_ABS)
   axis.SetTitleSize(plot_style.RATIO_TITLE_SIZE_PX)
   axis.SetLabelSize(plot_style.RATIO_LABEL_SIZE_PX)
   axis.SetTitleOffset(offset)


def plot(varname):

   nbin, low, high = plot_vars.get_binning(varname)
   doLog = varname in plot_vars.get_logy_vars()
   titleX = plot_vars.get_xlabel(varname)

   listHisto = getHisto(mytree, category, varname, year, nbin, low, high, blind=args.blind)

   if not listHisto:                       # variable missing from the snapshot
      print(f"   -> skipped '{varname}': input not available")
      return

   hists = {obj.name: obj.hOBJ for obj in listHisto}
   hData = hists.get('hData')
   hDY = hists.get('hDY')

   c, pad1, pad2 = plot_style.make_canvas_pads(doLog)

   # Draw stackXS: backgrounds stacked+filled; signal drawn unstacked as
   # solid outline lines on top (true scale, not summed into the background).
   BKGstack = ROOT.THStack()
   SIGstack = ROOT.THStack()

   bin_labels = plot_vars.get_bin_labels(varname)   # None unless categorical

   def _clean_and_label(h):
      """Remove negative bins and apply category-axis bin labels; returns the TH1."""
      hist = h.GetValue()
      for ibin in range(1, hist.GetNbinsX() + 1):
         if hist.GetBinContent(ibin) < 0.0:
            hist.SetBinContent(ibin, 0.0)
            hist.SetBinError(ibin, 0.0)
      if bin_labels:
         for i, lab in enumerate(bin_labels, start=1):
            hist.GetXaxis().SetBinLabel(i, lab)
      return hist

   for name in BKG_PROCS:
      h = hists.get(name)
      if not h:
         continue
      print('Integral ',h.GetName(), " = ", h.Integral())
      BKGstack.Add(_clean_and_label(h))

   for name in SIG_PROCS:
      h = hists.get(name)
      if not h:
         continue
      print('Integral ',h.GetName(), " = ", h.Integral())
      hist = _clean_and_label(h)
      hist.SetFillStyle(0)          # outline only, not filled/stacked
      SIGstack.Add(hist)

   stack = BKGstack

   rangeYax = 10
   if not doLog: rangeYax = 2
   if hData and hDY: stack.SetMaximum(rangeYax*max(hData.GetValue().GetMaximum(),hDY.GetValue().GetMaximum()))
   if hDY: stack.SetMinimum(hDY.GetValue().GetMaximum()/1000000);

   if varname == "mva":
      stack.SetMinimum(hDY.GetValue().GetMaximum()/1000000000);

   pad1.cd()
   # Draw data first
   if hData:  print('Integral ',hData.GetName(), " = ", hData.Integral())      
   if hData: hData.SetMarkerStyle(20)
   if hData: hData.SetMarkerSize(1.2)
   if hData: hData.SetLineWidth(2)
   if hData: hData.SetLineColor(ROOT.kBlack)
   if hData: hData.Draw("E ")

   pad1.cd()
   stack.Draw("HIST")
   
   stack.GetXaxis().SetLabelSize(0)     # glued layout: x-axis shown only on ratio pad
   stack.GetXaxis().SetTitleSize(0)

   if varname == "mass":
      stack.GetYaxis().SetTitle("Events/ 1 [GeV]")
   stack.GetYaxis().SetTitleOffset(1.1)
   stack.GetYaxis().SetLabelSize(0.04)
   stack.GetYaxis().SetTitleSize(0.045)
   stack.GetYaxis().ChangeLabel(1, -1, 0)

   SIGstack.Draw("HIST NOSTACK SAME")   # signal: unstacked solid outlines, true scale
   
   if hData: hData.Draw("E SAME") #comment for the plots now
   
   pad2.cd()

   ratio = hData.Clone("dataratio")
   mcTOT = BKGstack.GetStack().Last()
   print("ALL mcTOT integral(): ",mcTOT.Integral())
   print("ALL data integral(): ",hData.Integral())
    
   ratio.Divide(mcTOT)
   ratio.GetYaxis().SetTitle("data/MC")
   ratio.GetYaxis().SetRangeUser(*plot_vars.get_ratio_range(varname))
   if 'CR' in dirLOCAL_: ratio.GetYaxis().SetRangeUser(0.,2.5)
   ratio.GetXaxis().SetTitle(titleX)    # single x-axis title on the bottom pad
   _style_ratio_axis(ratio.GetXaxis(), plot_style.RATIO_X_TITLE_OFFSET)
   _style_ratio_axis(ratio.GetYaxis(), plot_style.RATIO_Y_TITLE_OFFSET)
   
   ratio.Draw("pe")
   lineZero = ROOT.TLine(mcTOT.GetXaxis().GetXmin(), 1.,  mcTOT.GetXaxis().GetXmax(), 1.)
   lineZero.SetLineColor(11)
   lineZero.Draw("same")
   
   pad1.cd()
   legend = ROOT.TLegend(*plot_style.LEGEND_POS)

   legend.SetNColumns(plot_style.LEGEND_NCOLUMNS)
   legend.SetColumnSeparation(plot_style.LEGEND_COLUMN_SEP)
   legend.SetFillStyle(0)
   legend.SetBorderSize(0)
   legend.SetTextSize(plot_style.LEGEND_TEXT_SIZE)
#   legend.SetTextAlign(32)
   legend.SetTextAlign(12)  # left align for readability

   sig_entries = LEGEND_SIG_LEP if category in LEP_CATEGORIES else LEGEND_SIG_HAD
   for name, label, style in LEGEND_BKG + sig_entries:
      h = hists.get(name)
      if h and h.Integral() > 0:
         legend.AddEntry(h.GetValue(), label, style)
   legend.Draw();

   
   # CMS label + lumi/energy (official style via plot_style / cmsstyle)
   _cmslabel = plot_style.cms_label(pad1, year)

   # Add TLine(s) marking the blinded region, if this variable is blinded
   blind_range = plot_vars.get_blind_range(varname)
   if blind_range is not None:
      lo, hi = blind_range
      line1 = ROOT.TLine(lo, 0, lo, 500000.)
      line1.SetLineColor(11)
      line1.Draw()
      line2 = ROOT.TLine(hi, 0, hi, 500000.)
      line2.SetLineColor(11)
      line2.Draw()

   os.makedirs(myOutDir, exist_ok=True)   # create the (per-group) output dir on demand
   c.SaveAs(f"{myOutDir}{varname}_{category}{year}_Stack.png")
   print(varname+".png")


def plotVBF():
   
   plot("mva")
   '''
   plot("jetvbf1_pt")
   plot("jetvbf2_pt")
   plot("jetvbf1_eta")
   plot("jetvbf2_eta")

   plot("mjj")
   plot("detajj")
   plot("rpt")
   plot("zepvar")
   plot("puppimet_pt_vbf")
   '''

def plotVHlep():

   plot("mva")
   plot("category_vlcat")

   '''
   plot("rpt")

   plot("lepton_pt")
   plot("puppimet_pt_vhlep")
   plot("mt_vhlep")
   '''

def plotTTHlep():

   plot("mva")
   plot("category_ttlcat")

   '''
   plot("lepton_pt")
   plot("lepton_eta")
   plot("lepton2_pt")
   plot("puppimet_pt_tthlep_tthad")
   plot("mt_tthlep")
   plot("ht")
   plot("jet1_pt_tthlep")
   plot("jet1_eta")
   '''

def plotVHhad():

   plot("mva")

   '''
   plot("goodwjj_mass")
   plot("goodwjj_discr")
   plot("goodwjj_pt")
   plot("goodwjj_eta")
   plot("goodwjjptoverhpt")
   plot("detawjjh")
   plot("dphiwjjh")
   plot("rpt")
   '''

def plotTThad():

   plot("category_tthcat")
   plot("mva")
   '''
   plot("jet1_pt_tthad")
   plot("jet1_eta")
   plot("wtopjetmass")
   plot("wtopjetdiscr")
   plot("ht")
   plot("njets")
#   plot("nbmjets")

   plot("puppimet_pt_tthlep_tthad")
   plot("topmassreco")
   plot("deta_j1j2")

   plot("mindr_h_bjet")
   plot("mindr_h_anyjet")
   '''

def plotZinvH():

   plot("mva")
   '''
   plot("dphimeth")
   plot("rpt")

   plot("puppimet_ptoverhpt")

   plot("puppimet_pt_zinv")
   plot("deltaphimetmu1")
   plot("deltaphimetmu2")
   '''

def plotMuons():

   plot("muon1_pt")
   plot("muon2_pt")
   plot("muon1_eta")
   plot("muon2_eta")
   plot("dimuon_pt")
   plot("dimuon_eta")
   plot("dimuon_rapidity")
   plot("muon1_norm_pt")
   plot("muon2_norm_pt")
   plot("costhetacs")
   plot("phistarcs")
   if category == "VLcat" or category == "TTLcat" or category == "TTHcat":
      plot("muon1_sip3d")
      plot("muon2_sip3d")
#   plot("fsrph_pt")
#   plot("muon1_phi")
#   plot("muon2_phi")
   plot("deta_muons")


# ---------------------------------------------------------------------------
# Dispatch: which plot groups run, driven by --plots
# ---------------------------------------------------------------------------
#
# Groups are non-overlapping and each maps to one small draw_* function:
#   mass     -> dimuon Higgs candidate mass          (all categories)
#   mva      -> BDT/MVA discriminant                 (all categories)
#   category -> category-index plot (category_vlcat/ttlcat/tthcat)  (VLcat/TTLcat/TTHcat only)
#   muons    -> muon kinematics                      (all categories)
#
# The per-category plot*() functions above (plotVBF, plotVHlep, ...) are kept
# as a reference menu of detailed kinematic plots (mostly commented out). They
# are NOT part of the default dispatch; call them manually or uncomment lines
# inside them if you want those extra distributions.


def draw_mass():
    plot("mass")


def draw_mva():
    plot("mva")


def draw_category():
    # category-index plot; binning differs per category, only defined for some
    if category == "VLcat":
        plot("category_vlcat")
    elif category == "TTLcat":
        plot("category_ttlcat")
    elif category == "TTHcat":
        plot("category_tthcat")
    else:
        print(f"[plots] 'category' index plot not defined for {category} — skipping")


# group name -> the function that draws it
GROUP_FUNCS = {
    "mass":     draw_mass,
    "mva":      draw_mva,
    "category": draw_category,
    "muons":    plotMuons,
}

# order used when --plots all (or default) is requested
ALL_GROUPS = ["mass", "mva", "category", "muons"]


def run_group(name):
    """Point the output at a clear per-group subdirectory, then draw it.

    Output layout:  <outdir>/<category>_<year>/<group>/Stack<var>_<cat>_<year>.png
    The directory is created lazily inside plot() right before each SaveAs, so
    no empty folders are left for groups that draw nothing (e.g. 'category' in
    a category that has no index plot).
    """
    global myOutDir
    myOutDir = os.path.join(baseOutDir, f"{category}{year}", name) + "/"
    print(f"[plots] group '{name}' -> {myOutDir}")
    GROUP_FUNCS[name]()


def main():
    requested = args.plots
    if "all" in requested:
        groups = ALL_GROUPS
    else:
        # preserve a sensible order, drop duplicates, keep only known groups
        groups = [g for g in ALL_GROUPS if g in set(requested)]

    print(f"[plots] category={category} year={year.lstrip('_')} "
          f"groups={groups}")
    print(f"[plots] blinding: {'ON (data blinded in signal region)' if args.blind else 'OFF (UNBLINDED)'}")
    print(f"[plots] base output dir: {baseOutDir}")

    for name in groups:
        run_group(name)


if __name__ == "__main__":
    main()
