import numpy as np
#import pandas as pd
import ROOT
ROOT.gROOT.SetBatch(1)
from ROOT import TCanvas, TH1F, TF1, TLegend, gPad, THStack, TColor
import tdrstyle
import os.path
import os

from ROOT import TCanvas, TGraph
from ROOT import gROOT
from math import sin
from array import array

BR = 0.00022
xsecRun3={
    'ggH':52230, # 0.4 × σ(13 TeV) + 0.6 × σ(14 TeV)                                                                                                                                           $
    'VBFH':4078,
    'Wm':567.7,
    'Wp':888.9,
    'ZH':943.9,
    'TTH':570,
}

print('theo(ggH) (fb)',xsecRun3['ggH']*BR)
print('theo(qqH) (fb)',xsecRun3['VBFH']*BR)

dir_ = 'workspaces/'

ax = range(1,8) #
#     ax = range(1,6) #
labels = [ 'VBF'
           ,'ggH'
           ,'VL'
           ,'VH'
           ,'ZinvH'
           ,'TTL'
           ,'TTH'
          ]

limitfiles = ["workspaces_GGHVBF_bdt_DEC16/higgsCombine_DEC16_VBF.AsymptoticLimits.mH125.root",
              "workspaces_GGHVBF_bdt_DEC16/higgsCombine_DEC16_ggH.AsymptoticLimits.mH125.root",
              "workspaces_MINOR_DEC16/higgsCombine_DEC16_VL.AsymptoticLimits.mH125.root",
              "workspaces_MINOR_DEC16/higgsCombine_DEC16_VH.AsymptoticLimits.mH125.root",
              "workspaces_MINOR_DEC16/higgsCombine_DEC16_Zinv.AsymptoticLimits.mH125.root",
              "workspaces_MINOR_DEC16/higgsCombine_DEC16_TTL.AsymptoticLimits.mH125.root",
              "workspaces_MINOR_DEC16/higgsCombine_DEC16_TTH.AsymptoticLimits.mH125.root",
              ]

c44 = ROOT.TCanvas("c44","c44",1200,800)
aexl = [0.5 for a in ax]
aexh = [0.5 for a in ax]
ay      = [  ]
aeyl    = [  ]
aeyh    = [  ]
aeyl_y  = [  ]
aeyh_y  = [  ]
obs = [ ]

for f in limitfiles:
     fo= ROOT.TFile.Open(f)
     tr = fo.Get("limit")
     tr.GetEntry(0)
     aeyl_y.append(tr.limit)
     tr.GetEntry(1)
     aeyl.append(tr.limit)
     tr.GetEntry(2)
     ay.append(tr.limit)
     tr.GetEntry(3)
     aeyh.append(tr.limit)
     tr.GetEntry(4)
     aeyh_y.append(tr.limit)
#     tr.GetEntry(5)
#     obs.append(tr.limit)
    
n = len(ax)

col = ROOT.TColor()
green = col.GetColor('#19C405')
yellow = col.GetColor('#FEC30A')

legend_args = (0.64, 0.7, 0.92, 0.85, '', 'NDC')

legend = TLegend(*legend_args)
legend.SetFillStyle(0)
legend.SetTextFont(42)

x, y, exl, exh, eyl, eyh, eyl_y, eyh_y, eyl_n, eyh_n = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ) , array( 'd' ), array( 'd' )
#x, y, exl, exh, eyl, eyh = array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' ), array( 'd' )
print(ay[0])

yo = array( 'd' )
yow = array( 'd' )
for i in range( n ): 
     
    print('i=',i,' ;limit=',ay[i])
    x.append( ax[i])
    y.append( ay[i])
    
    exh.append( aexh[i] )
    exl.append( aexl[i] )
    
    eyh.append( (aeyh[i]-ay[i]) )
    eyl.append( (ay[i]-aeyl[i]) )
    
    eyh_y.append( (aeyh_y[i]-ay[i]) )
    eyl_y.append( (ay[i]-aeyl_y[i]) )
    
    eyh_n.append( 0 )
    eyl_n.append( 0 )
    
#    yow.append(0.5)
#    yo.append(obs[i])
    
gae = ROOT.TGraphAsymmErrors(n, x, y, exl, exh, eyl, eyh)
yae = ROOT.TGraphAsymmErrors(n, x, y, exl, exh, eyl_y, eyh_y)
nae = ROOT.TGraphAsymmErrors(n, x, y, exl, exh, eyl_n, eyh_n)
#gro = ROOT.TGraphErrors(n, x, yo,yow,eyh_n)
#gro.SetMarkerStyle(20)
#gro.SetMarkerColor(1)
#gro.SetLineColor(1)
#gro.SetLineWidth(2)
#gro.SetMarkerSize(1.5)

gae.SetFillColor(green)
#gae.SetFillStyle(3001)
gae.SetMarkerStyle(0)

nae.SetLineStyle(2)
nae.SetLineWidth(2)
nae.SetMarkerSize(0)

yae.SetFillColor(yellow)
#yae.SetFillStyle(3001)
yae.SetMarkerStyle(0)


#legend.AddEntry(gro, "Observed", 'lp')
legend.AddEntry(nae, "Median expected", 'lp')
legend.AddEntry(gae, "68% expected", 'f')
legend.AddEntry(yae, "95% expected", 'f')

dummyHistogram = ROOT.TH1F("dummy","",n,0.5,n+.5)
dummyHistogram.GetXaxis().SetTitleOffset(1.8)
dummyHistogram.GetXaxis().SetLabelSize(0.045)
dummyHistogram.GetXaxis().SetRangeUser(0.5, ax[-1]+0.5)
dummyHistogram.GetYaxis().SetRangeUser(0, 15) # all
#if doChekSync: dummyHistogram.GetYaxis().SetRangeUser(0, 0.5) # GF and VBF
dummyHistogram.GetYaxis().SetTitle("signal strength (H#rightarrow #mu #mu)")
dummyHistogram.GetYaxis().SetTitleSize(0.04)
dummyHistogram.GetYaxis().SetTitleOffset(1.8)

xax = dummyHistogram.GetXaxis()
for i in range(1, n+1):
#for i in range(0, n):     
    binIndex = xax.FindBin(i)
    print(binIndex, " ", str(labels[i-1]))
    xax.SetBinLabel(binIndex, str(labels[i-1]))
#xax.LabelsOption("v")
xax.SetTickLength(0)
dummyHistogram.Draw("AXIS")

#.LabelsOption("v", "X")
gae.Draw("p e2")
yae.Draw("e2 same")
gae.Draw("e2 same")
nae.Draw("p same")
#gro.Draw("p same")  # comment observed

legend.Draw()

gPad.RedrawAxis()
#if doChekSync: tdrstyle.cmsPrel(39540, energy= 13,simOnly=False)
tdrstyle.cmsPrel(286500, energy= 13.6,simOnly=False)
c44.Draw()

latex = ROOT.TLatex()
latex.SetTextSize(0.04)
#if mesonCat=='Rho': latex.DrawLatex(2, 0.85*dummyHistogram.GetMaximum(), "Rho")
#if mesonCat=='Phi': latex.DrawLatex(2, 0.85*dummyHistogram.GetMaximum(), "Phi")

#mySTR=''
#if doChekSync: mySTR="_synchJune22"
#else: mySTR='_final'
mySTR='_15dec'

c44.SaveAs("~/public_html/HMUMU_FITS/limit"+mySTR+".pdf")
c44.SaveAs("~/public_html/HMUMU_FITS/limit"+mySTR+".png")
