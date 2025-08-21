import ROOT
import os
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
    '_2024':107.3, #C-I
}

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

def getHisto(item, nbin, low, high, doLog, category, year, doSignal):

   print("getHisto getting called")

   year = '_'+str(year)
   dirLOCAL_='/work/submit/mariadlf/HmumuRun3/ROOTFILES/'
   
   mytree = ROOT.TChain('events')   
   if (category == '_ggHcat'): mytree.Add(dirLOCAL_+'snapshot_mc10'+year+category+'.root')  #ggH
   if (category == '_VBFcat'): mytree.Add(dirLOCAL_+'snapshot_mc11'+year+category+'.root')  #VBF

   mytree.SetBranchStatus("*", 0)
   mytree.SetBranchStatus("HiggsCandCorrMass", 1)
   mytree.SetBranchStatus("mc", 1)
   mytree.SetBranchStatus("w_allSF", 1)
   mytree.SetBranchStatus("lumiIntegrated", 1)

   ROOT.gROOT.cd() 
   h = ROOT.TH1F("hSig", "hSig", nbin, low, high)
   h.SetDirectory(0)
   h1 = ROOT.TH1F("hBkg", "hBkg", nbin, low, high)
   h1.SetDirectory(0)   

   #-------------- BELOW RDF

   nEntries = mytree.GetEntries()
   print('nEntries=',nEntries)
   
   for i in range(nEntries):
         mytree.GetEntry(i)
         var = mytree.HiggsCandCorrMass
         wei = mytree.w_allSF * mytree.lumiIntegrated
         if doSignal:
               if mytree.mc==10 or mytree.mc==11  and not math.isnan(var):
                     h.Fill( var, wei )
         else:
               if (mytree.mc==100) and not math.isnan(var): # DY mc=100
                     h.Fill( var, wei )               

   return h
