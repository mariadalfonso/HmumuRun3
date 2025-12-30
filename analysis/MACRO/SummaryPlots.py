import ROOT
import os
from array import array
import math
import sys

from LoadTree import loadTree
from prepareHisto import getHisto, createCanvasPads, lumis

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()

category = sys.argv[1]
year = '_'+sys.argv[2]

myOutDir  = "/home/submit/mariadlf/public_html/HMUMU_NOV/"
dirLOCAL_ = '/work/submit/mariadlf/HmumuRun3/ROOTFILES/'

mytree = ROOT.TChain('events')
mytree = loadTree(mytree, dirLOCAL_, category, year )


def plot(item, nbin, low, high, doLog, plotString, titleX):

   listHisto = getHisto(mytree, category, item, nbin, low, high)

   print(listHisto)
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
   
   c, pad1, pad2 = createCanvasPads(doLog)

   # Draw stackXS with MC contributions
   allstack = ROOT.THStack()
   BKGstack = ROOT.THStack()
   SIGstack = ROOT.THStack()

   for h in [hTT2L, hTop, hZg, hVV, hEWK, hDY]:

      print('Integral ',h.GetName(), " = ", h.Integral())
      BKGstack.Add(h.GetValue())

   for h in [hWH, hTTH, hZH, hVBFH, hggH]:
      print('Integral ',h.GetName(), " = ", h.Integral())      
      SIGstack.Add(h.GetValue())


   labelsH = ["", "n_top", "n_W", "resolved (3jets)"]
   labelsL = ["", "H_{#mu#mu}+e", "H_{#mu#mu}+ee", "H_{#mu#mu}+#mu", "H_{#mu#mu}+#mu#mu"]

   for h in [hWH, hTTH, hZH, hVBFH, hggH, hTT2L, hTop, hZg, hVV, hEWK, hDY]:

      if not h:
         continue  # skip null histograms

      if item==210:
         # Set labels for the 4 bins
         for i, lab in enumerate(labelsL, start=1):
            h.GetXaxis().SetBinLabel(i, lab)

      if item==211:
         # Set labels for the 4 bins
         for i, lab in enumerate(labelsH, start=1):
            h.GetXaxis().SetBinLabel(i, lab)

      allstack.Add(h.GetValue())

   stack = allstack

   rangeYax = 10
   if not doLog: rangeYax = 2
   if hData and hDY: stack.SetMaximum(rangeYax*max(hData.GetValue().GetMaximum(),hDY.GetValue().GetMaximum()))
   if hDY: stack.SetMinimum(hDY.GetValue().GetMaximum()/1000000);

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
   
   stack.GetXaxis().SetLabelSize(0.04)
   stack.GetXaxis().SetTitleSize(0.045)
   stack.GetXaxis().SetTitle(titleX)

   if item==4:
      stack.GetYaxis().SetTitle("Events/ 1 [GeV]")
   stack.GetYaxis().SetTitleOffset(1.1)
   stack.GetYaxis().SetLabelSize(0.04)
   stack.GetYaxis().SetTitleSize(0.045)
   stack.GetYaxis().ChangeLabel(1, -1, 0)
   
   if hData: hData.Draw("E SAME")
   
   pad2.cd()

   ratio = hData.Clone("dataratio")
   mcTOT = BKGstack.GetStack().Last()
   print("ALL mcTOT integral(): ",mcTOT.Integral())
   print("ALL data integral(): ",hData.Integral())
    
   ratio.Divide(mcTOT)
   ratio.GetYaxis().SetTitle("data/MC")
   ratio.GetYaxis().SetRangeUser(0.5,1.5)
   if item==4: ratio.GetYaxis().SetRangeUser(0.75,1.25)
#   if item==4: ratio.GetYaxis().SetRangeUser(0.90,1.10)
   if 'CR' in dirLOCAL_: ratio.GetYaxis().SetRangeUser(0.,2.5)
   if (item==205): ratio.GetYaxis().SetRangeUser(0.,2.)
   ratio.GetXaxis().SetTitleOffset(4.)
   ratio.GetXaxis().SetTitleSize(0.15)
   ratio.GetXaxis().SetLabelSize(0.12)
   
   ratio.GetYaxis().SetTitleOffset(0.3)
   ratio.GetYaxis().SetTitleSize(0.15)
   ratio.GetYaxis().SetLabelSize(0.12)
   
   ratio.Draw("pe")
   lineZero = ROOT.TLine(mcTOT.GetXaxis().GetXmin(), 1.,  mcTOT.GetXaxis().GetXmax(), 1.)
   lineZero.SetLineColor(11)
   lineZero.Draw("same")
   
   pad1.cd()
   legend = ROOT.TLegend(0.60, 0.60, 0.95, 0.9)  # adjust as needed

   legend.SetNColumns(2)
   legend.SetFillStyle(0)
   legend.SetBorderSize(0)
   legend.SetTextSize(0.035)
#   legend.SetTextAlign(32)
   legend.SetTextAlign(12)  # left align for readability

   if hData and hData.Integral()>0: legend.AddEntry(hData.GetValue(), "Data" ,"lep")
   if hDY and hDY.Integral()>0: legend.AddEntry(hDY.GetValue(), "DY+jets (QCD)", "f")
   if hEWK and hEWK.Integral()>0: legend.AddEntry(hEWK.GetValue(), "DY+jets (EWK)", "f")
   if hVV and hVV.Integral()>0: legend.AddEntry(hVV.GetValue(), "VV + VVV", "f")
   if hTT2L and hTT2L.Integral()>0: legend.AddEntry(hTT2L.GetValue(), "ttbar 2l+ jets", "f")
   if hTop and hTop.Integral()>0: legend.AddEntry(hTop.GetValue(), "tt1l, tW, ttV", "f")
   if hZg and hZg.Integral()>0: legend.AddEntry(hZg.GetValue(), "H#rightarrowZ#gamma + jets", "f")
   if category in ["VLcat", "TTLcat", "TTHcat", "VHcat", "Zinvcat"]:
      if hTTH and hTTH.Integral()>0: legend.AddEntry(hTTH.GetValue(), "ttH", "f")
      if hWH and hWH.Integral()>0: legend.AddEntry(hWH.GetValue(), "WH", "f")
      if hZH and hZH.Integral()>0: legend.AddEntry(hZH.GetValue(), "ZH", "f")
   else:
      if hVBFH and hVBFH.Integral()>0: legend.AddEntry(hVBFH.GetValue(), "VBF H", "f")
      if hggH and hggH.Integral()>0: legend.AddEntry(hggH.GetValue(), "ggH", "f")
   legend.Draw();

   
   # Add label
   text = ROOT.TLatex()
   text.SetNDC()
   text.SetTextFont(72)
   text.SetTextSize(0.045)
   text.DrawLatex(0.15, 0.93, "CMS")
   text.SetTextFont(42)
#   text.DrawLatex(0.15 + 0.10, 0.93, "Simulation")
   text.DrawLatex(0.15 + 0.10, 0.93, "Internal")
   text.SetTextSize(0.04)
   text.DrawLatex(0.65, 0.93, "#sqrt{s} = 13.6 TeV,%0.2f fb^{-1}"% (lumis[year]))

   # Add TLine blind
   line1 = ROOT.TLine( 110, 0, 110, 500000.)
   line1.SetLineColor(11);
   if item==4: line1.Draw()
   line2 = ROOT.TLine( 150, 0, 150, 500000.)
   line2.SetLineColor(11);
   if item==4: line2.Draw()

   string = category+year
#   c.SaveAs(myOutDir+"Stack"+plotString+string+"_bothBB.png")
   c.SaveAs(myOutDir+"Stack"+plotString+"_"+string+".png")
   print(plotString+".png")

def plotVBF():
   
   plot(100, 100, 0. , 1000., True, "Mjj","Mjj")
   plot(101, 100, 0. , 10., True, "dEtaJJ","dEtaJJ")
   plot(102, 100, 0. , 1., True, "Rpt","Rpt")
   plot(103, 100, 0. , 2., True, "ZepVar","ZepVar")
   plot(104, 200, 0. , 200., True, "jetVBF1_Pt","jetVBF1_Pt")
   plot(105, 200, 0. , 200., True, "jetVBF2_Pt","jetVBF2_Pt")   
   plot(106, 100, -5. , 5., True, "jetVBF1_Eta","#eta jetVBF1")
   plot(107, 100, -5. , 5., True, "jetVBF2_Eta","#eta jetVBF2")
   plot(301, 100, 0. , 500., False, "PuppiMET_Pt","PuppiMET_Pt")

def plotVHlep():

   plot(201, 100, 10. , 100., False, "Lepton_Pt","Lepton_Pt")
   plot(210, 5, 0. , 5., False, "category","category")
   plot(301, 300, 0. , 300., True, "PuppiMET_Pt","PuppiMET_Pt")
   plot(212, 150, 0. , 150., False, "mt","mt")

def plotTTHlep():

   plot(201, 100, 10. , 100., False, "Lepton_Pt","Lepton_Pt")
   plot(210, 5, 0. , 5., False, "category","category")
   plot(301, 100, 0. , 500., True, "PuppiMET_Pt","PuppiMET_Pt")
   plot(212, 150, 0. , 150., False, "mt","mt")
   plot(260, 150, 10. , 160., False, "Jet1_Pt","Jet1_Pt")
   plot(263, 200, 0. , 1000., True, "HT","HT")

def plotVHhad():

   plot(251, 50, 60. , 110., True, "goodWjj_mass","goodWjj_mass")
   plot(252, 50, 0.75, 1., True, "goodWjj_discr","goodWjj_discr")
   plot(253, 40, 150. , 550., True, "goodWjj_pt","goodWjj_pt")
   plot(254, 50, 0. , 5., True, "goodWjjPtOverHpt","goodWjjPtOverHpt")

def plotTThad():

   plot(211, 4, 0. , 4., False, "category","category")
   plot(260, 150, 10. , 160., False, "Jet1_Pt","Jet1_Pt")
   plot(261, 150, 50. , 200., True, "WTopJetMass","WTopJetMass")
   plot(262, 50, 0.5 , 1., False, "WTopJetDiscr","WTopJetDiscr")
   plot(263, 200, 0. , 1000., True, "HT","HT")
   plot(264, 10, 0. , 10., True, "Njets","Njets")

def plotZinvH():

   plot(301, 400, 0. , 400., False, "PuppiMET_Pt","PuppiMET_Pt")
   plot(302, 50, 0. , 5., True, "PuppiMET_PtOverHpt","PuppiMET_PtOverHpt")

def plotMuons():

   plot(10, 100, 0. , 200., True, "Muon1_pt", "p^{T}_{#mu^{1}} [GeV]")
   plot(11, 100, 0. , 200., True, "Muon2_pt", "p^{T}_{#mu^{2}} [GeV]")
   plot(12, 60, -3. , 3., False, "Muon1_eta", "#eta_{#mu^{1}}")
   plot(13, 60, -3. , 3., False, "Muon2_eta", "#eta_{#mu^{2}}")
   if category == "VLcat" or category == "TTLcat" or category == "TTHcat":
      plot(14, 100, -0. , 20., True, "Muon1_sip3d", "Muon1_sip3d")
      plot(15, 100, -0. , 20., True, "Muon2_sip3d", "Muon2_sip3d")
#   plot(16, 200, 0.1 , 20.1, False, "FsrPH_pt", "p^{T}_{#gammaFSR} [GeV]")
#   plot(18, 60, -3.14 , 3.14, False, "Muon1_phi", "#phi_{#mu^{1}}")
#   plot(19, 60, -3.14 , 3.14, False, "Muon2_phi", "#phi_{#mu^{2}}")
   plot(20, 60, 0. , 6., False, "dEtaMuons", "#eta_{#mu^{1}}-#eta_{#mu^{2}}")

if __name__ == "__main__":

   #   plot(4, 130, 70. , 200., True, "HCandCorrMass", "m_{#mu^{+},#mu^{-}}^{H} [GeV]")

   if category == "VLcat" or category == "TTLcat" or category == "TTHcat":
      plot(4, 130, 70. , 200., True, "HCandCorrMass", "m_{#mu^{+},#mu^{-}}^{H} [GeV]")
   else:
      plot(4, 130, 70. , 200., True, "HCandCorrMass", "m_{#mu^{+},#mu^{-}}^{H} [GeV]")

   if category == "VBFcat" or category == "ggHcat": plot(99, 100, 0. , 1., True, "discrMVA", "MVA discr")

   plot(7, 100, 0. , 0.001, True, "HCandCorrMassErr", "#sigma(m_{#mu^{+},#mu^{-}}^{H})/m_{#mu^{+},#mu^{-}}^{H}")

   plotMuons()

   plot(95, 100, 0. , 100., False, "nVTX","nVTX")
   plot(6, 60, -3 , 3., True, "HCandCorrRapidity", "y_{#mu^{+},#mu^{-}}^{H}")
   if category == "VBFcat" or category == "TTLcat":
      plot(5, 200, 0. , 200., False, "HCandCorrPt", "p^{T}_{#mu^{+},#mu^{-}}^{H} [GeV]")
   else:
      plot(5, 100, 0. , 100., True, "HCandCorrPt", "p^{T}_{#mu^{+},#mu^{-}}^{H} [GeV]")

   if category == "VBFcat": plotVBF()
   if category == "VLcat": plotVHlep()
   if category == "TTLcat": plotTTHlep()
   if category == "VHcat": plotVHhad()
   if category == "TTHcat": plotTThad()
   if category == "Zinvcat" : plotZinvH()
