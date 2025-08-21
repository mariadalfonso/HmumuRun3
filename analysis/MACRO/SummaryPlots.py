import ROOT
import os
from array import array
import math
import sys

from LoadTree import loadTree
from prepareHisto import getHisto, createCanvasPads, lumis

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()

category = '_'+sys.argv[1]
year = '_'+sys.argv[2]

myOutDir  = "/home/submit/mariadlf/public_html/HMUMU_July/"
#dirLOCAL_ = '/home/submit/mariadlf/Hrare/CMSSW_10_6_27_new/src/HmumuRun3/analysis/ROOTFILES/'
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
      if obj.name == 'hTT12L': hTT12L = obj.hOBJ
      if obj.name == 'hVV': hVV = obj.hOBJ
      if obj.name == 'hEWK': hEWK = obj.hOBJ
   
   c, pad1, pad2 = createCanvasPads(doLog)

   # Draw stackXS with MC contributions
   allstack = ROOT.THStack()
   BKGstack = ROOT.THStack()
   SIGstack = ROOT.THStack()
   for h in [hTT12L, hZg, hVV, hEWK, hDY]:
      print('Integral ',h.GetName(), " = ", h.Integral())
      BKGstack.Add(h.GetValue())

   for h in [hZH, hWH, hTTH, hVBFH, hggH]:
      print('Integral ',h.GetName(), " = ", h.Integral())      
      SIGstack.Add(h.GetValue())

   for h in [hZH, hWH, hTTH, hVBFH, hggH, hTT12L, hZg, hVV, hEWK, hDY]:
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
   legend = ROOT.TLegend()
   legend.SetX1NDC(0.85)  # left
   legend.SetY1NDC(0.80)  # bottom
   legend.SetX2NDC(0.95)  # right
   legend.SetY2NDC(0.98)  # top
   
   legend.SetFillStyle(0)
   legend.SetBorderSize(0)
   legend.SetTextSize(0.04)
   legend.SetTextAlign(32)

   if hDY and hDY.Integral()>0: legend.AddEntry(hDY.GetValue(), "DY + jets", "f")
   if hTT12L and hDY.Integral()>0: legend.AddEntry(hTT12L.GetValue(), "ttbar 2l + jets", "f")
   if hVV and hVV.Integral()>0: legend.AddEntry(hVV.GetValue(), "VV + VVV", "f")
   if hZg and hZg.Integral()>0: legend.AddEntry(hZg.GetValue(), "Z#gamma + jets", "f")
   if hVBFH and hVBFH.Integral()>0: legend.AddEntry(hVBFH.GetValue(), "VBF H", "f")
   if hggH and hggH.Integral()>0: legend.AddEntry(hggH.GetValue(), "ggH", "f")
   if hData and hData.Integral()>0: legend.AddEntry(hData.GetValue(), "Data" ,"lep")
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
   text.DrawLatex(0.65, 0.93, "#sqrt{s} = 13 TeV,%0.2f fb^{-1}"% (lumis[year]))

   # Add TLine blind
   line1 = ROOT.TLine( 110, 0, 110, 500000.)
   line1.SetLineColor(11);
   if item==4: line1.Draw()
   line2 = ROOT.TLine( 150, 0, 150, 500000.)
   line2.SetLineColor(11);
   if item==4: line2.Draw()

   string = category+year
   c.SaveAs(myOutDir+"Stack"+plotString+string+".png")
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

def plotVH():

   plot(201, 100, 0. , 100., False, "Lepton_Pt","Lepton_Pt")

def plotMuons():

   plot(10, 200, 0. , 200., False, "Muon1_pt", "p^{T}_{#mu^{1}} [GeV]")
   plot(11, 200, 0. , 200., False, "Muon2_pt", "p^{T}_{#mu^{2}} [GeV]")
   plot(12, 200, 0.1 , 20.1, False, "FsrPH_pt", "p^{T}_{#gammaFSR} [GeV]")


if __name__ == "__main__":

   if category == "_WHcat" :
      plot(4, 150, 50. , 200., False, "HCandCorrMass", "m_{#mu^{+},#mu^{-}}^{H} [GeV]")
   else:
      plot(4, 150, 50. , 200., True, "HCandCorrMass", "m_{#mu^{+},#mu^{-}}^{H} [GeV]")

   plot(95, 100, 0. , 100., False, "nVTX","nVTX")
   plot(6, 60, -3 , 3., True, "HCandCorrRapidity", "y_{#mu^{+},#mu^{-}}^{H}")
   if category == "_VBFcat" :
      plot(5, 200, 0. , 200., False, "HCandCorrPt", "p^{T}_{#mu^{+},#mu^{-}}^{H} [GeV]")
   else:
      plot(5, 100, 0. , 100., False, "HCandCorrPt", "p^{T}_{#mu^{+},#mu^{-}}^{H} [GeV]")

   plotMuons()

   if category == "_VBFcat" : plotVBF()
   if category == "_WHcat" : plotVH()
