import ROOT
import os
from array import array
import math
import sys

from LoadTree import loadTree

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch()

class MyHisto():

    def __init__(self, name, hOBJ):
        self.name = name
        self.hOBJ = hOBJ

    def __repr__(self):
        return self.name

ROOT.ROOT.EnableImplicitMT()
RDataFrame = ROOT.RDataFrame

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
    '_12022':7.98, # C-D
    '_22022':26.67, # E, F, G
    '_12023':17.794, #C
    '_22023':9.451, #D
    '_2024':107.3, #C-I
}

def deltaPhi(phi1,phi2):
   result = phi1 - phi2
   if result > float(M_PI): result -= float(2*M_PI)
   if result <= -float(M_PI): result += float(2*M_PI)
   return result


# Create the plot
def createCanvasPads(doLog):

   # Create canvas with pad
   c = ROOT.TCanvas("c", "", 600, 600)
   pad = ROOT.TPad("upper_pad", "", 0, 0, 1, 1)
#   gStyle.SetOptStat(0)                                                                                                                                                                                                                                                                                             
   if doLog : pad.SetLogy(1)
   pad.SetTickx(False)
   pad.SetTicky(False)
   pad.Draw()
   pad.cd()

   ydiv = 0.2
   pad1 = ROOT.TPad("upper_pad", "", 0.,ydiv,1.,1.)
   if doLog : pad1.SetLogy(1)
   pad2 = ROOT.TPad("lower_pad", "", 0.,0.,1.,ydiv)

   pad1.Draw()
   pad2.Draw()

   return c, pad1, pad2


def getHisto(mytree, category, item, nbin, low, high):

   darkBlue = (0, 0, 255)
   redDark = (191, 34, 41)
#   redMed = (255, 40, 0)
   redMed = (237, 41, 57)
#   redMed = (220,0,5)
#   redLight = (255,61,65)
   redLight = (255,82,82)
   orange = (255, 204, 153)
   orangeDark = (255,150,79)
   gray = (136,139,141)
   azure = (100, 192, 232)
   azureDark = (96, 147, 172)
   green = (144, 238, 144)
   greenDark = (98, 172, 141)
   gold = (212 ,175, 55)
   violet = (181 ,100, 227)

   ####

   df = RDataFrame(mytree)

   ## -------------------   
   ## PUT here some PRESELECTION
   ## -------------------
   
   ## -------------------   
   ## PUT here some Variable
   ## -------------------   
   if item == 4 : var = "HiggsCandCorrMass"
   if item == 5 : var = "HiggsCandCorrPt"
   if item == 6 : var = "HiggsCandCorrRapidity"
   if item == 7 : var = "HiggsCandMassErr/HiggsCandCorrMass"
   ##
   if item == 10 : var = "Muon1_pt"
   if item == 11 : var = "Muon2_pt"
   if item == 12 : var = "Muon1_eta"
   if item == 13 : var = "Muon2_eta"
   if item == 14 : var = "Muon1_sip3d"
   if item == 15 : var = "Muon2_sip3d"
   if item == 16 : var = "FsrPH_pt"
   if item == 17 : var = "FsrPH_eta"
   if item == 95 : var = "PV_npvsGood"

   # for VBF
   if item == 100 : var = "Mjj"   
   if item == 101 : var = "dEtaJJ"
   if item == 102 : var = "RPt"
   if item == 103 : var = "ZepVar"
   if item == 104 : var = "jetVBF1_Pt"
   if item == 105 : var = "jetVBF2_Pt"
   if item == 106 : var = "jetVBF1_Eta"
   if item == 107 : var = "jetVBF2_Eta"
   if item == 108 : var = "jetVBF1_Phi"
   if item == 109 : var = "jetVBF2_Phi"

   # for VH lep
   if item == 201: var = "Lepton_Pt"
   if item == 210: var = "category"

   # for VH had
   if item == 251: var = "goodWjj_mass[0]"
   if item == 252: var = "FatJet_pNet_WvsQCD_discr[0]"
   if item == 253: var = "FatJet_pNet_mass[0]"

   # for ZinvH
   if item == 301: var = "PuppiMET_pt"
   if item == 302: var = "PuppiMET_pt/HiggsCandCorrPt"

   ## -------------------
   ## FILL the histograms
   ## -------------------

   df_common = df.Define("var","{}".format(var)).Define("weight","w_allSF")
   #.Filter("int(category)==3")
   #.Filter("abs(HiggsCandCorrMass-125)<15")
   hDY = df_common.Filter("mc==100 or mc==103 or mc==104").Histo1D(("hDY","h",nbin, low, high),"var","weight")
   hEWK = df_common.Filter("mc==101").Histo1D(("hEWK","h",nbin, low, high),"var","weight")
   hTT2L = df_common.Filter("mc==102").Histo1D(("hTT2L","h",nbin, low, high),"var","weight")
   hTop = df_common.Filter("mc==105 or mc==106 or mc==107 or mc==221 or mc==222 or mc==223 or mc==224 or mc==225").Histo1D(("hTop","h",nbin, low, high),"var","weight")
   hVV = df_common.Filter("mc==201 or mc==202 or mc==203 or mc==204 or mc==205 or mc==211 or mc==212 or mc==213 or mc==214").Histo1D(("hVV","h",nbin, low, high),"var","weight")

   hVBFH = df_common.Filter("(mc==10)").Histo1D(("hVBFH","h",nbin, low, high),"var","weight")
   hggH = df_common.Filter("(mc==11)").Histo1D(("hggH ","h",nbin, low, high),"var","weight")
   hWH = df_common.Filter("mc==12 or mc==13").Histo1D(("hWH","h",nbin, low, high),"var","weight")
   hZH = df_common.Filter("mc==14").Histo1D(("hZH","h",nbin, low, high),"var","weight")
   hTTH = df_common.Filter("mc==15").Histo1D(("hTTH","h",nbin, low, high),"var","weight")
   hZg = df_common.Filter("mc==20 or mc==21 or mc==22 or mc==23 or mc==24 or mc==25 or mc==26").Histo1D(("hZg","h",nbin, low, high),"var","weight")
   
   if False: hData = df_common.Filter("mc<0").Histo1D(("hData","h",nbin, low, high),"var","weight")
   else:
       print('HELLO')
       if (item==4): hData = df_common.Filter("mc<0 and (var<110 or var>150)").Histo1D(("hData","h",nbin, low, high),"var","weight")
       else: hData = df_common.Filter("mc<0").Histo1D(("hData","h",nbin, low, high),"var","weight")

   if hData: hData.SetMarkerStyle(20)
   if hData: hData.SetMarkerSize(1.2)
   if hData: hData.SetLineWidth(2)      
   if hData: hData.SetLineColor(ROOT.kBlack)

   for h, color in zip([hDY, hTT2L, hTop, hVV, hEWK, hVBFH, hggH, hWH, hZH, hTTH, hZg], [azure, green, greenDark, orange, gold,redMed, redDark, redLight, redLight, redDark, orangeDark]):
       if h:
           h.SetLineWidth(3)
           h.SetLineColor(ROOT.TColor.GetColor(*color))
           h.SetFillColor(ROOT.TColor.GetColor(*color))

   if hData: hData_ = MyHisto('hData', hData)
   #
   hDY_ = MyHisto('hDY', hDY)
   hTT2L_ = MyHisto('hTT2L',hTT2L)
   hTop_ = MyHisto('hTop',hTop)
   hVV_ = MyHisto('hVV', hVV)
   hEWK_ = MyHisto('hEWK', hEWK)
   #
   hVBFH_ = MyHisto('hVBFH', hVBFH)
   hggH_ = MyHisto('hggH', hggH)
   hZH_ = MyHisto('hZH', hZH)
   hWH_ = MyHisto('hWH', hWH)
   hTTH_ = MyHisto('hTTH', hTTH)
   #
   hZg_ = MyHisto('hZg', hZg)

   listHisto = [hDY_, hTT2L_, hTop_, hVV_, hEWK_, hData_, hVBFH_, hggH_, hWH_, hZH_, hTTH_, hZg_]
   print(item)

   return listHisto
