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

plot_style.setup_style()

# ---------------------------------------------------------------------------
# Configuration / CLI
# ---------------------------------------------------------------------------

CATEGORIES = ["VBFcat", "ggHcat", "VLcat", "TTLcat", "TTHcat",
              "VHcat", "Zinvcat"]

# Which plot groups each category can draw. "mass" always runs.
# Names map to the plot*() functions / plotMuons() below.
DEFAULT_INDIR  = f"/work/submit/{getpass.getuser()}/HmumuRun3/ROOTFILES/"
DEFAULT_OUTDIR = f"/home/submit/{getpass.getuser()}/HmumuRun3/output_plots/"

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


def plot(item, nbin, low, high, doLog, plotString, titleX):

   listHisto = getHisto(mytree, category, item, year, nbin, low, high, blind=args.blind)

   if not listHisto:                       # variable missing from the snapshot
      print(f"   -> skipped '{plotString}' (item {item}): input not available")
      return

   for obj in listHisto:
      if obj.name == 'hData': hData = obj.hOBJ
      #
      if obj.name == 'hZH': hZH = obj.hOBJ
      if obj.name == 'hWH': hWH = obj.hOBJ
      if obj.name == 'hTTH': hTTH = obj.hOBJ
      if obj.name == 'hVBFH': hVBFH = obj.hOBJ
      if obj.name == 'hggH': hggH = obj.hOBJ
      if obj.name == 'hZg': hZg = obj.hOBJ
      #
      if obj.name == 'hDY': hDY = obj.hOBJ
      if obj.name == 'hTT2L': hTT2L = obj.hOBJ
      if obj.name == 'hTop': hTop = obj.hOBJ
      if obj.name == 'hVV': hVV = obj.hOBJ
      if obj.name == 'hEWK': hEWK = obj.hOBJ
   
   c, pad1, pad2 = plot_style.make_canvas_pads(doLog)

   # Draw stackXS: backgrounds stacked+filled; signal drawn unstacked as
   # solid outline lines on top (true scale, not summed into the background).
   BKGstack = ROOT.THStack()
   SIGstack = ROOT.THStack()

   labelsH = ["", "n_top", "n_W", "resolved (5jets)"]
   labelsL = ["", "H_{#mu#mu}+e", "H_{#mu#mu}+ee", "H_{#mu#mu}+#mu", "H_{#mu#mu}+#mu#mu", "H_{#mu#mu}+e#mu"]

   def _clean_and_label(h):
      """Remove negative bins and apply category-axis bin labels; returns the TH1."""
      hist = h.GetValue()
      for ibin in range(1, hist.GetNbinsX() + 1):
         if hist.GetBinContent(ibin) < 0.0:
            hist.SetBinContent(ibin, 0.0)
            hist.SetBinError(ibin, 0.0)
      if item==210:
         for i, lab in enumerate(labelsL, start=1):
            hist.GetXaxis().SetBinLabel(i, lab)
      if item==211:
         for i, lab in enumerate(labelsH, start=1):
            hist.GetXaxis().SetBinLabel(i, lab)
      return hist

   for h in [hTT2L, hTop, hZg, hVV, hEWK, hDY]:
      if not h:
         continue
      print('Integral ',h.GetName(), " = ", h.Integral())
      BKGstack.Add(_clean_and_label(h))

   for h in [hWH, hTTH, hZH, hVBFH, hggH]:
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

   if item==99:
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

   if item==4:
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
   ratio.GetYaxis().SetRangeUser(0.5,1.5)
#   if item==4: ratio.GetYaxis().SetRangeUser(0.75,1.25)
   if item==4: ratio.GetYaxis().SetRangeUser(0.90,1.10) # mass
   if item==5: ratio.GetYaxis().SetRangeUser(0.75,1.25) # PT
   if 'CR' in dirLOCAL_: ratio.GetYaxis().SetRangeUser(0.,2.5)
   if (item==205): ratio.GetYaxis().SetRangeUser(0.,2.)
   ratio.GetXaxis().SetTitle(titleX)    # single x-axis title on the bottom pad
   ratio.GetXaxis().SetTitleFont(plot_style.RATIO_FONT_ABS)
   ratio.GetXaxis().SetLabelFont(plot_style.RATIO_FONT_ABS)
   ratio.GetXaxis().SetTitleSize(plot_style.RATIO_TITLE_SIZE_PX)
   ratio.GetXaxis().SetLabelSize(plot_style.RATIO_LABEL_SIZE_PX)
   ratio.GetXaxis().SetTitleOffset(plot_style.RATIO_X_TITLE_OFFSET)

   ratio.GetYaxis().SetTitleFont(plot_style.RATIO_FONT_ABS)
   ratio.GetYaxis().SetLabelFont(plot_style.RATIO_FONT_ABS)
   ratio.GetYaxis().SetTitleSize(plot_style.RATIO_TITLE_SIZE_PX)
   ratio.GetYaxis().SetLabelSize(plot_style.RATIO_LABEL_SIZE_PX)
   ratio.GetYaxis().SetTitleOffset(plot_style.RATIO_Y_TITLE_OFFSET)
   
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

   if hData and hData.Integral()>0: legend.AddEntry(hData.GetValue(), "Data" ,"lep")
   if hDY and hDY.Integral()>0: legend.AddEntry(hDY.GetValue(), "DY+jets (QCD)", "f")
   if hEWK and hEWK.Integral()>0: legend.AddEntry(hEWK.GetValue(), "DY+jets (EWK)", "f")
   if hVV and hVV.Integral()>0: legend.AddEntry(hVV.GetValue(), "VV + VVV", "f")
   if hTT2L and hTT2L.Integral()>0: legend.AddEntry(hTT2L.GetValue(), "t#bar{t} 2l", "f")
   if hTop and hTop.Integral()>0: legend.AddEntry(hTop.GetValue(), "Top (1l, tW/tZq, ttV/4t)", "f")
   if hZg and hZg.Integral()>0: legend.AddEntry(hZg.GetValue(), "H#rightarrowZ#gamma + jets", "f")
   if category in ["VLcat", "TTLcat", "TTHcat", "VHcat", "Zinvcat"]:
      if hTTH and hTTH.Integral()>0: legend.AddEntry(hTTH.GetValue(), "ttH", "l")
      if hWH and hWH.Integral()>0: legend.AddEntry(hWH.GetValue(), "WH", "l")
      if hZH and hZH.Integral()>0: legend.AddEntry(hZH.GetValue(), "ZH", "l")
   else:
      if hVBFH and hVBFH.Integral()>0: legend.AddEntry(hVBFH.GetValue(), "VBF H", "l")
      if hggH and hggH.Integral()>0: legend.AddEntry(hggH.GetValue(), "ggH", "l")
   legend.Draw();

   
   # CMS label + lumi/energy (official style via plot_style / cmsstyle)
   _cmslabel = plot_style.cms_label(pad1, year)

   # Add TLine blind
   line1 = ROOT.TLine( 110, 0, 110, 500000.)
   line1.SetLineColor(11);
   if item==4: line1.Draw()
   line2 = ROOT.TLine( 150, 0, 150, 500000.)
   line2.SetLineColor(11);
   if item==4: line2.Draw()

   string = category+year
   os.makedirs(myOutDir, exist_ok=True)   # create the (per-group) output dir on demand
   c.SaveAs(myOutDir+"Stack"+plotString+"_"+string+".png")
#   c.SaveAs(myOutDir+"Stack"+plotString+"_"+string+"_HSB.png")
#   c.SaveAs(myOutDir+"Stack"+plotString+"_"+string+"_ZCR.png")
   print(plotString+".png")


def plotVBF():
   
   plot(99, 100, 0. , 1., True, "discrMVA", "MVA discr")
   '''
   plot(104, 200, 0. , 200., True, "jetVBF1_Pt","jetVBF1_Pt")
   plot(105, 200, 0. , 200., True, "jetVBF2_Pt","jetVBF2_Pt")   
   plot(106, 100, -5. , 5., True, "jetVBF1_Eta","#eta jetVBF1")
   plot(107, 100, -5. , 5., True, "jetVBF2_Eta","#eta jetVBF2")

   plot(100, 100, 0. , 1000., True, "Mjj","Mjj")
   plot(101, 100, 0. , 10., True, "dEtaJJ","dEtaJJ")
   plot(102, 100, 0. , 1., True, "Rpt","Rpt")
   plot(103, 100, 0. , 2., True, "ZepVar","ZepVar")
   plot(301, 100, 0. , 500., False, "PuppiMET_Pt","PuppiMET_Pt")
   '''

def plotVHlep():

   plot(99, 100, 0. , 1., True, "discrMVA", "MVA discr")
   plot(210, 5, 0. , 5., False, "category","category")

   '''
   plot(102, 100, 0. , 1., True, "Rpt","Rpt")

   plot(201, 100, 10. , 100., False, "Lepton_Pt","Lepton_Pt")
   plot(301, 300, 0. , 300., True, "PuppiMET_Pt","PuppiMET_Pt")
   plot(212, 150, 0. , 150., False, "mt","mt")
   '''

def plotTTHlep():

   plot(99, 100, 0. , 1., True, "discrMVA", "MVA discr")
   plot(210, 6, 0. , 6., False, "category","category")

   '''
   plot(201, 100, 10. , 100., False, "Lepton_Pt","Lepton_Pt")
   plot(203, 60, -3. , 3., False, "Lepton_Eta","Lepton_Eta")
   plot(202, 100, 10. , 100., False, "Lepton2_Pt","Lepton2_Pt")
   plot(301, 100, 0. , 500., True, "PuppiMET_Pt","PuppiMET_Pt")
   plot(212, 30, 0. , 150., False, "mt","mt")
   plot(263, 200, 0. , 1000., True, "HT","HT")
   plot(260, 150, 10. , 160., False, "Jet1_Pt","Jet1_Pt")
   plot(265, 100, -5. , 5., False, "Jet1_Eta","Jet1_Eta")
   '''

def plotVHhad():

   plot(99, 100, 0. , 1., True, "discrMVA", "MVA discr")

   '''
   plot(251, 50, 60. , 110., True, "goodWjj_mass","goodWjj_mass")
   plot(252, 50, 0.75, 1., True, "goodWjj_discr","goodWjj_discr")
   plot(253, 40, 150. , 550., True, "goodWjj_pt","goodWjj_pt")
   plot(255, 60, -3. , 3., True, "goodWjj_eta","goodWjj_eta")
   plot(254, 50, 0. , 5., True, "goodWjjPtOverHpt","goodWjjPtOverHpt")
   plot(255, 60, 0. , 2.5, True, "dEtaWjjH","dEtaWjjH")
   plot(256, 100, 0. , 6.28, True, "dPhiWjjH","dPhiWjjH")
   plot(102, 100, 0. , 1., True, "Rpt","Rpt")
   '''

def plotTThad():

   plot(211, 4, 0. , 4., False, "category","category")
   plot(99, 100, 0. , 1., True, "discrMVA", "MVA discr")
   '''
   plot(260, 100, 10. , 200., False, "Jet1_Pt","Jet1_Pt")
   plot(265, 100, -5. , 5., False, "Jet1_Eta","Jet1_Eta")
   plot(261, 150, 50. , 200., True, "WTopJetMass","WTopJetMass")
   plot(262, 50, 0.5 , 1., False, "WTopJetDiscr","WTopJetDiscr")
   plot(263, 200, 0. , 1000., True, "HT","HT")
   plot(264, 10, 0. , 10., True, "Njets","Njets")
#   plot(267, 10, 0. , 10., False, "nBMjets","nBMjets")

   plot(301, 100, 0. , 500., True, "PuppiMET_Pt","PuppiMET_Pt")
   plot(266, 35, 0. , 350., True, "TopMassReco", "TopMassReco")
   plot(268, 30, 0. , 3., True, "dEta_j1j2", "dEta_j1j2")

   plot(269, 30, 0. , 6.5, True, "mindR_H_BJet", "mindR_H_BJet")
   plot(270, 30, 0. , 6.5, True, "mindR_H_AnyJet", "mindR_H_AnyJet")
   '''

def plotZinvH():

   plot(99, 100, 0. , 1., True, "discrMVA", "MVA discr")
   '''
   plot(305, 50, 0. , 3.5, True, "dPhiMETH","dPhiMETH")
   plot(102, 100, 0. , 1., True, "Rpt","Rpt")

   plot(302, 50, 0. , 5., True, "PuppiMET_PtOverHpt","PuppiMET_PtOverHpt")

   plot(301, 400, 0. , 400., False, "PuppiMET_Pt","PuppiMET_Pt")
   plot(303, 50, 0. , 3.5, True, "deltaPhiMETMu1","deltaPhiMETMu1")
   plot(304, 50, 0. , 3.5, True, "deltaPhiMETMu2","deltaPhiMETMu2")
   '''

def plotMuons():

   plot(10, 100, 0. , 200., True, "Muon1_pt", "p_{T}^{#mu_{1}} [GeV]")
   plot(11, 100, 0. , 200., True, "Muon2_pt", "p_{T}^{#mu_{2}} [GeV]")
   plot(12, 60, -3. , 3., True, "Muon1_eta", "#eta^{#mu_{1}}")
   plot(13, 60, -3. , 3., True, "Muon2_eta", "#eta^{#mu_{2}}")
   if category == "VLcat" or category == "TTLcat" or category == "TTHcat":
      plot(14, 100, -0. , 20., True, "Muon1_sip3d", "Muon1_sip3d")
      plot(15, 100, -0. , 20., True, "Muon2_sip3d", "Muon2_sip3d")
#   plot(16, 200, 0.1 , 20.1, False, "FsrPH_pt", "p^{T}_{#gammaFSR} [GeV]")
#   plot(18, 60, -3.14 , 3.14, False, "Muon1_phi", "#phi_{#mu^{1}}")
#   plot(19, 60, -3.14 , 3.14, False, "Muon2_phi", "#phi_{#mu^{2}}")
   plot(20, 60, 0. , 6., True, "dEtaMuons", "|#Delta#eta(#mu^{1}, #mu^{2})|")


# ---------------------------------------------------------------------------
# Dispatch: which plot groups run, driven by --plots
# ---------------------------------------------------------------------------
#
# Groups are non-overlapping and each maps to one small draw_* function:
#   mass     -> dimuon Higgs candidate mass          (all categories)
#   mva      -> BDT/MVA discriminant                 (all categories)
#   category -> category-index plot (item 210/211)   (VLcat/TTLcat/TTHcat only)
#   muons    -> muon kinematics                      (all categories)
#
# The per-category plot*() functions above (plotVBF, plotVHlep, ...) are kept
# as a reference menu of detailed kinematic plots (mostly commented out). They
# are NOT part of the default dispatch; call them manually or uncomment lines
# inside them if you want those extra distributions.


def draw_mass():
    plot(4, 130, 70., 200., True, "HCandCorrMass", "m_{#mu#mu} [GeV]")


def draw_mva():
    plot(99, 100, 0., 1., True, "discrMVA", "MVA discr")


def draw_category():
    # category-index plot; binning differs per category, only defined for some
    if category == "VLcat":
        plot(210, 5, 0., 5., False, "category", "category")
    elif category == "TTLcat":
        plot(210, 6, 0., 6., False, "category", "category")
    elif category == "TTHcat":
        plot(211, 4, 0., 4., False, "category", "category")
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
