import ROOT
import os
from array import array
import math
import sys

from LoadTree import loadTree
import plot_style

plot_style.setup_style()

# lumis kept importable from here for back-compat; canonical copy is in plot_style
lumis = plot_style.lumis

class MyHisto():

    def __init__(self, name, hOBJ):
        self.name = name
        self.hOBJ = hOBJ

    def __repr__(self):
        return self.name

ROOT.ROOT.EnableImplicitMT()
RDataFrame = ROOT.RDataFrame

def deltaPhi(phi1,phi2):
   result = phi1 - phi2
   if result > float(M_PI): result -= float(2*M_PI)
   if result <= -float(M_PI): result += float(2*M_PI)
   return result

ROOT.gInterpreter.Declare("""
float deltaPhi(float phi1, float phi2) {
    float dphi = phi1 - phi2;
    while (dphi >  M_PI) dphi -= 2*M_PI;
    while (dphi <= -M_PI) dphi += 2*M_PI;
    return abs(dphi);
}
""")

def removeNegBin(hist):

    # Remove negative bins
    for ibin in range(1, hist.GetNbinsX() + 1):
        if hist.GetBinContent(ibin) < 0.0:
            hist.SetBinContent(ibin, 0.0)
            hist.SetBinError(ibin, 0.0)   # or keep the original error, see below

def addOverflow(h, addUnderflow=False):
    nb = h.GetNbinsX()

    # Overflow -> last visible bin
    h.SetBinContent(
        nb,
        h.GetBinContent(nb) + h.GetBinContent(nb + 1)
    )
    h.SetBinError(
        nb,
        (h.GetBinError(nb)**2 + h.GetBinError(nb + 1)**2)**0.5
    )

    if addUnderflow:
        h.SetBinContent(
            1,
            h.GetBinContent(1) + h.GetBinContent(0)
        )
        h.SetBinError(
            1,
            (h.GetBinError(1)**2 + h.GetBinError(0)**2)**0.5
        )

    # Optional: clear overflow/underflow
    h.SetBinContent(nb + 1, 0.)
    h.SetBinError(nb + 1, 0.)

    if addUnderflow:
        h.SetBinContent(0, 0.)
        h.SetBinError(0, 0.)


# createCanvasPads now lives in plot_style.make_canvas_pads (glued ratio pads).
# Colours now live in plot_style (PROCESS_COLORS / _RGB).



from LoadTree import hmumu, hzgamma, hww
from LoadTree import dy_2223,dy_24,dy_pt2223,dy_pt24,dy_minllo,dy_j
from LoadTree import dyewk
from LoadTree import vv,tt2l,ttV2223,ttV24,top
import plot_vars

def make_filter(ids):
    return " || ".join([f"mc=={x}" for x in ids])


def getHisto(mytree, category, varname, year, nbin, low, high, blind=True, region_filter=None):

   ####

   df = RDataFrame(mytree)

   ## -------------------   
   ## PRESELECTION
   ## -------------------

   var = plot_vars.get_expr(varname)

   ## -------------------
   ## FILL the histograms
   ## -------------------

   selectionReg = "true"

   #ZCR
   #selectionReg = "HiggsCandCorrMass>70 and HiggsCandCorrMass<110"
   #HiggsSideband
   #selectionReg = "(HiggsCandCorrMass>110 and HiggsCandCorrMass<115) or ((HiggsCandCorrMass>135 and HiggsCandCorrMass<150))"

   # fit range
   #selectionMVA = "HiggsCandCorrMass>110 and HiggsCandCorrMass<150"

   selectionMVA = "true"
   if varname == "mva":
       print(category)
       if category in ['ggHcat']: selectionMVA = '(var<0.5)'
       if category in ['VBFcat']: selectionMVA = '(var<0.64)'
       if category in ['VLcat']: selectionMVA = '(var<0.32)'
       if category in ['VHcat']: selectionMVA = '(var<0.86)'
       if category in ['TTHcat']: selectionMVA = '(var<0.80)'
       if category in ['Zinvcat']: selectionMVA = '(var<0.78)'

   # DY pT reweighting
   ggHcorr = "(mc==100 || mc==103 || mc==104 || mc==109) ? boson_ptWeight : 1."

   # year-dependent normalization
   if year == "_2025":
       norm = "(mc>0) ? (110.59/109.82) : 1."
   elif year == "_2026":
       norm = "(mc>0) ? (25.31/109.82) : 1."
   else:
       norm = "1."

   weightSTD = "w_allSF"
   weightExpr = f"{weightSTD} * ({ggHcorr}) * ({norm})"

   # Let ROOT validate the variable expression. If it references a branch that
   # isn't in the snapshot (e.g. discrMVA0 when the MVA classifier wasn't run),
   # the JIT compilation fails here; catch it and skip this plot cleanly.
   # NOTE: catch BaseException — PyROOT JIT failures don't always surface as a
   # plain Exception subclass, so a narrower `except Exception` can miss them.
   try:
       df_common = df.Define("var","{}".format(var)).Define("weight","{}".format(weightExpr)).Filter(selectionReg)
       if region_filter:
           # applied identically to every process (backgrounds, signals, data)
           # so it's a shared, single-source-of-truth restriction, e.g. a
           # Higgs-mass sideband/window split for MC shape comparisons.
           df_common = df_common.Filter(region_filter)
   except BaseException as e:
       print(f"⚠️  '{varname}': cannot build variable '{var}' "
             f"(likely a missing branch: {type(e).__name__}) — skipping this plot")
       return None

   #.Filter(selection)
   #.Filter("int(category)==3")
   #.Filter("abs(HiggsCandCorrMass-125)<15")
   hDY = df_common.Filter(make_filter(dy_2223+dy_24+dy_pt2223+dy_pt24+dy_minllo+dy_j)).Histo1D(("hDY","h",nbin, low, high),"var","weight")
   hEWK = df_common.Filter(make_filter(dyewk)).Histo1D(("hEWK","h",nbin, low, high),"var","weight")
   hTT2L = df_common.Filter(make_filter(tt2l)).Histo1D(("hTT2L","h",nbin, low, high),"var","weight")
   hTop = df_common.Filter(make_filter(top + ttV2223 + ttV24)).Histo1D(("hTop","h",nbin, low, high),"var","weight")
   hVV = df_common.Filter(make_filter(vv)).Histo1D(("hVV","h",nbin, low, high),"var","weight")

   hVBFH = df_common.Filter("(mc==10)").Histo1D(("hVBFH","h",nbin, low, high),"var","weight")
   hggH = df_common.Filter("(mc==11)").Histo1D(("hggH ","h",nbin, low, high),"var","weight")
   hWH = df_common.Filter("mc==12 or mc==13").Histo1D(("hWH","h",nbin, low, high),"var","weight")
   hZH = df_common.Filter("mc==14").Histo1D(("hZH","h",nbin, low, high),"var","weight")
   hTTH = df_common.Filter("mc==15").Histo1D(("hTTH","h",nbin, low, high),"var","weight")
   hZg = df_common.Filter(make_filter(hzgamma+hww)).Histo1D(("hZg","h",nbin, low, high),"var","weight")

   # data histogram, with optional blinding of the signal-sensitive region.
   # blinding only affects data (hData); MC histograms are never blinded.
   blind_range = plot_vars.get_blind_range(varname)   # e.g. (110,150) for "mass"
   if not blind:
       if blind_range is not None:
           print('[getHisto] UNBLINDED data')
       hData = df_common.Filter("mc<0").Histo1D(("hData","h",nbin, low, high),"var","weight")
   elif blind_range is not None:
       lo, hi = blind_range
       hData = df_common.Filter(f"mc<0 and (var<{lo} or var>{hi})").Histo1D(("hData","h",nbin, low, high),"var","weight")
   elif varname == "mva":
       # blind above the per-category MVA cut
       hData = df_common.Filter("mc<0 and {}".format(selectionMVA)).Histo1D(("hData","h",nbin, low, high),"var","weight")
   else:
       hData = df_common.Filter("mc<0").Histo1D(("hData","h",nbin, low, high),"var","weight")

   if hData: hData.SetMarkerStyle(20)
   if hData: hData.SetMarkerSize(1.2)
   if hData: hData.SetLineWidth(2)      
   if hData: hData.SetLineColor(ROOT.kBlack)
   if hData: addOverflow(hData)

   _proc_order = ["hDY", "hTT2L", "hTop", "hVV", "hEWK",
                  "hVBFH", "hggH", "hWH", "hZH", "hTTH", "hZg"]
   for h, pname in zip([hDY, hTT2L, hTop, hVV, hEWK, hVBFH, hggH, hWH, hZH, hTTH, hZg], _proc_order):
       if h:
           h.SetLineWidth(3)
           col = plot_style.process_color(pname)
           h.SetLineColor(col)
           h.SetFillColor(col)
           addOverflow(h)
           removeNegBin(h)

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
   print(varname)

   return listHisto
