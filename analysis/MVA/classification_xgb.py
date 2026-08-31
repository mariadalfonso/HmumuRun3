import ROOT
import json

import xgboost as xgb
import os
import sys

ROOT.gROOT.SetBatch(True)

from settings import *

import matplotlib
matplotlib.use("Agg")     # IMPORTANT for batch mode / no display
import matplotlib.pyplot as plt
import scipy.stats

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score, roc_curve, auc
from scipy.ndimage import gaussian_filter1d

from sklearn.metrics import confusion_matrix
from sklearn.metrics import ConfusionMatrixDisplay

from LoadTree import loadTree
from LoadTree import vv, tt2l, ttV

#myDir='/work/submit/mariadlf/HmumuRun3/ROOTFILES/newID_'
myDir='/work/submit/mariadlf/HmumuRun3/ROOTFILES/'

years = ['_12022', '_22022', '_12023', '_22023', '_2024']

#category = sys.argv[1]

category="ggHcat"
#category="VBFcat"
#category="TTLcat"
#category="VLcat"
#
#category="VHcat"
#category="TTHcat"
#category="Zinvcat"
mytree = ROOT.TChain('events')
for year in years:
    mytree = loadTree(mytree, myDir, category, year )

doMultiClass=True
if category=="Zinvcat": doMultiClass=False
if category=="TTHcat": doMultiClass=False
if category=="VHcat": doMultiClass=False

if category in ["Zinvcat","TTLcat"]:
    params = paramsClass2.copy()
else:
    params = paramsClass1.copy()

signal_map = {
    "VBFcat": ["10"],
    "ggHcat": ["11","10"],
    "VHcat": ["12","13","14"],
    "VLcat": ["12","13","14"],
    "Zinvcat": ["14"],
    "TTHcat": ["15"],
    "TTLcat": ["15"]
}

class_map = {
    "ggHcat": {
        0: ["11"],                         # ggH  # -- main signal is 0
        1: ["10"],                         # VBF
        2: ["100","103","104","109"],      # BKGA
    },
    "VBFcat": {
        0: ["10"],                         # VBF  # -- main signal is 0
        1: ["11"],                         # ggH
        2: ["100","103","104","109"],      # BKGA
        3: ["101","99","98"]               # BKGB
    },
    "VLcat": {
	0: ["12","13"],                    # 3l
        1: ["14"],                         # 4l
        2: ["201","202","203","206","207","208","209","210"] + ["204","211","215"],  # BKGA 2l-3l
        3: ["205","212","213","214","216"] # BKGA 4l
    },
    "TTLcat": {
	0: ["15"],                    # 3l
        1: ["12","13","14"],
        2: tt2l+ ["221", "222", "224","226","227","228","229"] + ["107","105","106"] + ["242","243","244","245"] + ["225","227","228","230"],
        3: ["223","225"] + ["230","231","232","233","234","235","236"],
    }
}

labels_map = {
    "ggHcat": ["ggH", "qqH", "DY-allJets"],
    "VBFcat": ["qqH", "ggH", "DY-QCD", "DY-EWK"],
    "VLcat": ["WH", "ZH", "2l3l", "4l"],
    "TTLcat": ["ttH","VH","tt2l+ttW","ttZ"]
}

if doMultiClass:

    sizeClasses = len(labels_map[category])
    signal_classes = [0, 1]
    background_classes = list(range(2, sizeClasses))

else:

    signal_classes = [1]
    background_classes = [0]


# note inclusive 100, 103, 104 might disturb
bkg_map = {
    "VBFcat": ["100","103","104","109"] + ["101","99","98"],# DY QCD + EWK
    "ggHcat": ["100","103","104","109"], # DY (incl and mass binned)
    "VHcat": ["109"] + ["114","115","116","117","122","123","124","125"] + tt2l + vv,
    "VLcat": vv, # diboson
    "Zinvcat": tt2l, # top 2l
    "TTHcat": tt2l + ["109"], # ttbar2l both powheg nominal and altern ; DY as well
    "TTLcat": tt2l + ttV , # ttbar2l + ttV + tt1l + singleTop
}

#labelForPNG = "_DYQCD"
#labelForPNG = "_DYEWK"
#if doMultiClass: labelForPNG = "_multiclass"
#else: labelForPNG = ""
#labelForPNG = "_RelSigma1"

labelForPNG = ""

extraFolder = ""
if doMultiClass: extraFolder = "multiclass/"

dir_map = {
    "VBFcat": f"/home/submit/mariadlf/public_html/HMUMU_MVA/VBF/{extraFolder}",
    "ggHcat": f"/home/submit/mariadlf/public_html/HMUMU_MVA/ggH/{extraFolder}",
    "VHcat": f"/home/submit/mariadlf/public_html/HMUMU_MVA/VH/{extraFolder}",
    "VLcat": f"/home/submit/mariadlf/public_html/HMUMU_MVA/VL/{extraFolder}",
    "Zinvcat": f"/home/submit/mariadlf/public_html/HMUMU_MVA/Zinv/{extraFolder}",
    "TTHcat": f"/home/submit/mariadlf/public_html/HMUMU_MVA/TTH/{extraFolder}",
    "TTLcat": f"/home/submit/mariadlf/public_html/HMUMU_MVA/TTL/{extraFolder}"
}

variables_map = {
    "VBFcat": ["HiggsCandCorrPt", "RPt", "Mjj", "dEtaJJ", "ZepVar", "minDetaDiMuVBF", "dPhiJJ", "Muon1_norm_pt", "Muon2_norm_pt","jetVBF2_Pt","jetVBF1_Pt","jetVBF1_Eta", "jetVBF2_Eta","CenEta","CenPt"],#,"ptAsy","log_Mjj"],#"minDphiDiMuVBF"],
    "ggHcat": ["HiggsCandCorrPt", "Muon1_norm_pt", "Muon2_norm_pt", "nGoodJetsAll","Jet1_Pt","Jet1_Eta","deltaRJet1H"],
    "Zinvcat": ["HiggsCandCorrPt", "Muon1_norm_pt","Muon2_norm_pt","PuppiMET_pt","dPhiMETH","RPt"],
    "VLcat": ["HiggsCandCorrPt", "category","Muon1_norm_pt","Muon2_norm_pt","Lepton_Pt","Lepton2_Pt","RPt","PuppiMET_pt","dPhiVH","dEtaVH","VMass"],  ##,"dPhiWH","dEtaWH","dPhiZH","dEtaZH","ZMassPull","WMassPull"]
    #    ,"dPhiWH","dEtaWH","dPhiZH","dEtaZH","ZMassPull","WMassPull"
    "TTLcat": ["HiggsCandCorrPt", "category","Muon1_norm_pt","Muon2_norm_pt"]+["Lepton_Pt","Lepton2_Pt","PuppiMET_pt","dEtaLepH","mt","dPhiMETH","MetBisectorProj"] + ["HT","dEta_j1j2","mbb","Centrality"] + ["Jet1_Pt","Lepton_Eta","Lepton2_Eta"],
    "VHcat": ["HiggsCandCorrPt", "goodWjj_discr", "goodWjj_mass", "dEtaWjjH","dPhiWjjH","Muon1_norm_pt","Muon2_norm_pt","RPt"],
    "TTHcat": ["HiggsCandCorrPt", "HT", "nGoodJetsAll","category","Centrality","Jet1_Eta"] + ["PuppiMET_pt", "MetBisectorProj","dPhiMETH"] + [ "WTopJetDiscr","TopMassReco","TopPairChi2","dEta_j1j2","mindR_H_BJet"],
}

# here the variable that can help with the resolution
variables_resolution = {
    "VBFcat": ["HiggsCandCorrRapidity", "cosThetaCS", "phiStarCS"],
    "ggHcat": ["HiggsCandCorrRapidity", "cosThetaCS", "phiStarCS", "Muon1_eta","Muon2_eta","PuppiMET_pt"],
    "VHcat": ["HiggsCandCorrRapidity","cosThetaCS", "phiStarCS", "PuppiMET_pt"],
    "VLcat": ["HiggsCandCorrRapidity","cosThetaCS","phiStarCS"],
    "Zinvcat": [ "HiggsCandCorrRapidity","cosThetaCS","phiStarCS"],
    "TTLcat": ["HiggsCandCorrRapidity","phiStarCS"],
    "TTHcat": ["HiggsCandCorrRapidity","cosThetaCS","phiStarCS","WTopJetMass","Jet1_Pt"],
}

# here the no discrimination and very weak
variables_notUseful_map = {
    "ggHcat": ["HiggsCandCorrRapidity", "cosThetaCS", "phiStarCS", "Muon1_eta","Muon2_eta"], #"Jet1_Eta" has a peak at 0
    "VBFcat": ["HiggsCandCorrRapidity", "cosThetaCS","phiStarCS","Muon1_eta","Muon2_eta"],
    "Zinvcat": [ "HiggsCandCorrRapidity","cosThetaCS","phiStarCS","Muon1_eta","Muon2_eta"],
    "TTLcat": ["HiggsCandCorrRapidity", "cosThetaCS", "phiStarCS","Muon1_eta","Muon2_eta","ST","Lepton_MVAid","Lepton2_MVAid","Lepton_sip3d","Lepton2_sip3d","Jet1_Pt","Lepton_Eta","Lepton2_Eta","mindR_H_BJet","dR_H_LeadB","Lepton_charge","Lepton2_charge"],
    "TTHcat": ["HiggsCandCorrRapidity", "cosThetaCS", "phiStarCS", "Muon1_eta","Muon2_eta","Jet1_Pt","Muon1_sip3d","Muon2_sip3d","WTopJetMass","Muon1_norm_pt","Muon2_norm_pt"]+["mindR_H_AnyJet","mbb","nBMJets","LeadBJetPt","dR_H_LeadB"],
    "VLcat": ["Muon1_eta","Muon2_eta","Lepton_sip3d","Lepton2_sip3d","Lepton_MVAid","Lepton2_MVAid"],
    "VHcat": ["Muon1_eta","Muon2_eta","goodWjj_pt","goodWjj_eta","Muon2_pt"],
}

def get_var_range(df, var, pad=0.05):
    """Return (xmin, xmax) for a variable using RDataFrame Min/Max."""
    min_val = df.Min(var).GetValue()
    max_val = df.Max(var).GetValue()

    # Add small padding around edges (5% by default)
    span = max_val - min_val
    xmin = min_val - span * pad
    xmax = max_val + span * pad

    return xmin, xmax

def make_plots(df):

    # Define weight
    df = df.Define("weight", "w_allSF")

    # Define signal & background filters
    sig_ids = signal_map[category]
    bkg_ids = bkg_map[category]
    sig_filter = " || ".join([f"mc == {mid}" for mid in sig_ids])
    bkg_filter = " || ".join([f"mc == {mid}" for mid in bkg_ids])
    df_sig = df.Filter(sig_filter, f"Signal ({sig_filter})")
    df_bkg = df.Filter(bkg_filter, f"Signal ({bkg_filter})")
    df_sig = df_sig.Define("absweight", "fabs(weight)")
    df_bkg = df_bkg.Define("absweight", "fabs(weight)")

    for var in variables_map[category] + variables_resolution[category]:
        print(f"Plotting {var}...")

        # Get automatic min/max from **both** samples together
        xmin_s, xmax_s = get_var_range(df_sig, var)
        xmin_b, xmax_b = get_var_range(df_bkg, var)
        xmin = min(xmin_s, xmin_b)
        xmax = max(xmax_s, xmax_b)

        nbins = 40  # if you want different binning per variable, I can automate this too

        # Create histograms
        hsig = df_sig.Histo1D(
            (f"h_{var}_sig", f"{var} (Signal);{var};Normalized Events", nbins, xmin, xmax),
            var, "weight"
        )
        hbkg = df_bkg.Histo1D(
            (f"h_{var}_bkg", f"{var} (Background);{var};Normalized Events", nbins, xmin, xmax),
            var, "weight"
        )
        hsigABS = df_sig.Histo1D(
            (f"h_{var}_abssig", f"{var} (Signal) ABSw;{var};Normalized Events", nbins, xmin, xmax),
            var, "absweight"
        )
        hbkgABS = df_bkg.Histo1D(
            (f"h_{var}_absbkg", f"{var} (Background) ABSw;{var};Normalized Events", nbins, xmin, xmax),
            var, "absweight"
        )

        # Draw
        c = ROOT.TCanvas(f"c_{var}", "", 800, 600)
        hsig.SetLineColor(ROOT.kOrange)
        hsigABS.SetLineColor(ROOT.kRed)
        hbkg.SetLineColor(ROOT.kBlue)
        hbkgABS.SetLineColor(ROOT.kGreen+1)

        if hsig.Integral() > 0: hsig.Scale(1.0 / hsig.Integral())
        if hsigABS.Integral() > 0: hsigABS.Scale(1.0 / hsigABS.Integral())
        if hbkg.Integral() > 0: hbkg.Scale(1.0 / hbkg.Integral())
        if hbkgABS.Integral() > 0: hbkgABS.Scale(1.0 / hbkgABS.Integral())
        hsig.SetMaximum(1.5*max(hsig.GetMaximum(),hbkg.GetMaximum()))
        hsig.Draw("HIST")
        hsigABS.Draw("HIST SAME")
        hbkg.Draw("HIST SAME")
        hbkgABS.Draw("HIST SAME")

        # Explicit legend (middle right)
        leg = ROOT.TLegend(0.60, 0.40, 0.88, 0.65)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.04)

        leg.AddEntry(hsigABS.GetValue(), "sig", "l")
        leg.AddEntry(hbkg.GetValue(), "bkg", "l")
        leg.Draw()
        #c.BuildLegend()

        c.SaveAs(dir_map[category]+f"VARS/{var}_sig_bkg_{category}{labelForPNG}.png")

        print(f" → Saved {var}_sig_bkg.png")

def plotCovMatrix(matrixTXT,sample_name):

    df = pd.read_csv(matrixTXT, sep="\t", index_col=0)

    plt.figure(figsize=(8, 7))
    ax = plt.gca()

    im = plt.imshow(df, vmin=-1, vmax=1)

    plt.xticks(range(len(df.columns)), df.columns, rotation=45, ha="right")
    plt.yticks(range(len(df.index)), df.index)

    plt.colorbar(im, label="Correlation")

    # -----------------------------
    # ADD SEPARATION LINES
    # -----------------------------
    split = len(df) - 2

    ax.axhline(split - 0.5, color="white", linewidth=2)
    ax.axvline(split - 0.5, color="white", linewidth=2)

    # -----------------------------
    # values in cells
    # -----------------------------
    for i in range(len(df.index)):
        for j in range(len(df.columns)):
            plt.text(j, i, f"{df.iloc[i, j]:.2f}",
                     ha="center", va="center", fontsize=8)

    plt.title("Feature Correlation Matrix")
    plt.tight_layout()

    plt.savefig(f"{dir_map[category]}/correlation_matrix_{category}_{sample_name}{labelForPNG}.png", dpi=200)
    plt.close()

    print("✅ Saved correlation_matrix.png")

def load_process_class(class_id, variables, drawPlot=False):

    # Heavy-lifting in C++ and remote access of data
    df = ROOT.RDataFrame(mytree)
    df = df.Filter("(HiggsCandCorrMass>(110) and HiggsCandCorrMass<(150))","HiggsMass within reasonable range 125+-20") # Run2 selection
#    df = df.Define("HiggsCandCorrPt_norm", "HiggsCandCorrMass>0 ? HiggsCandCorrPt/HiggsCandCorrMass: 0.f")
    if category == "TTLcat":
        df = df.Define("Muon1_norm_pt", "HiggsCandCorrPt>0 ? Muon1_pt/HiggsCandCorrPt: 0.f")
        df = df.Define("Muon2_norm_pt", "HiggsCandCorrPt>0 ? Muon2_pt/HiggsCandCorrPt: 0.f")
    else:
        df = df.Define("Muon1_norm_pt", "HiggsCandCorrMass>0 ? Muon1_pt/HiggsCandCorrMass: 0.f")
        df = df.Define("Muon2_norm_pt", "HiggsCandCorrMass>0 ? Muon2_pt/HiggsCandCorrMass: 0.f")

#    df = df.Define("log_Mjj", "log(1.+Mjj)")
    df = df.Define("ptAsy", "(Muon1_pt-Muon2_pt)/(Muon1_pt+Muon2_pt)")
#    df = df.Define("log_HT", "log(1.+HT)")
#    df = df.Define("log_MET", "log(1.+PuppiMET_pt)")
#    df = df.Define("log_HiggsPt", "log(1.+HiggsCandCorrPt)")


    conditionSig="true"
    if category == "VLcat" or category == "TTLcat":
        conditionSig="Muon1_genPartFlav==1 and Muon2_genPartFlav==1 and Lepton_genPartFlav==1"  # this does nothing

    # IMPORTANT: multiclass classID 0-1-signal and 2-3-BKG
    if doMultiClass and ( category == "VBFcat" or  category == "ggHcat"):
        ids = class_map[category][class_id]
        filt = " || ".join([f"mc == {x}" for x in ids])

        df = df.Filter(filt)

    elif doMultiClass and (category == "VLcat" or  category == "TTLcat"):  ## to evolve to subsplit in category
        ids = class_map[category][class_id]
        filt = " || ".join([f"mc == {x}" for x in ids])

        df = df.Filter(conditionSig)
        df = df.Filter(filt)

    else:

#        if class_id == 1 and drawPlot: make_plots(df) # only make the plot once

        sig_ids = signal_map[category]
        bkg_ids = bkg_map[category]

        # IMPORTANT: binary classID 1-signal and 0-BKG
        if class_id == 1:
            sig_filter = " || ".join([f"mc == {mid}" for mid in sig_ids])
            df = df.Filter(sig_filter)
        elif class_id == 0:
            bkg_filter = " || ".join([f"mc == {mid}" for mid in bkg_ids])
            df = df.Filter(bkg_filter)
        else:
            raise ValueError("Invalid class_id")

    # Define the weight column (if not already present)
    df = df.Define("weight", "w_allSF")
    nevts = df.Count().GetValue()
    print(f"class {class_id} -- evt counts {nevts} (1 signal and 0 for BKG)")

    nWevts = df.Sum("weight").GetValue()
    print(f"class {class_id} -- evt weights {nWevts} (1 signal and 0 for BKG)")

    cols = df.AsNumpy(variables + ["weight"])

    # Push data to scipy ecosystem
    pdf = pd.DataFrame(cols)
    pdf['target'] = class_id
    return pdf

def overtraining(bdt, train_data, train_labels, train_weights, test_data, test_labels, test_weights):

    if doMultiClass:
        # score = P(VBF)
        train_scores = bdt.predict_proba(train_data)[:,0]
        test_scores  = bdt.predict_proba(test_data)[:,0]

        sig_mask_train = (train_labels == 0)
        sig_mask_test  = (test_labels  == 0)

    else:
        train_scores = bdt.predict_proba(train_data)[:,1]
        test_scores  = bdt.predict_proba(test_data)[:,1]

        sig_mask_train = (train_labels == 1)
        sig_mask_test  = (test_labels  == 1)

    bkg_mask_train = ~sig_mask_train
    bkg_mask_test  = ~sig_mask_test

    train_sig = train_scores[sig_mask_train]
    train_bkg = train_scores[bkg_mask_train]

    test_sig  = test_scores[sig_mask_test]
    test_bkg  = test_scores[bkg_mask_test]

    train_w_sig = train_weights[sig_mask_train]
    train_w_bkg = train_weights[bkg_mask_train]

    test_w_sig  = test_weights[sig_mask_test]
    test_w_bkg  = test_weights[bkg_mask_test]

    bins = np.linspace(0, 1, 50)

    plt.figure(figsize=(8,6))

    # --- TRAIN (histograms) ---
    plt.hist(train_sig, bins=bins, weights=train_w_sig,
             density=True, histtype='step', linewidth=2,
             label='Train Signal')

    plt.hist(train_bkg, bins=bins, weights=train_w_bkg,
             density=True, histtype='step', linewidth=2,
             label='Train Background')

    # --- TEST (points with errors) ---
    def weighted_hist_with_err(values, weights, bins):
        hist, edges = np.histogram(values, bins=bins, weights=weights, density=True)

        # statistical uncertainty (approx)
        sumw2, _ = np.histogram(values, bins=bins, weights=weights**2)
#        err = np.sqrt(sumw2) / np.sum(weights)
        err = np.sqrt(sumw2) / np.diff(edges) / np.sum(weights)

        centers = (edges[:-1] + edges[1:]) / 2
        return centers, hist, err

    # Signal
    x_sig, y_sig, err_sig = weighted_hist_with_err(test_sig, test_w_sig, bins)
    plt.errorbar(x_sig, y_sig, yerr=err_sig, fmt='o', label='Test Signal')

    # Background
    x_bkg, y_bkg, err_bkg = weighted_hist_with_err(test_bkg, test_w_bkg, bins)
    plt.errorbar(x_bkg, y_bkg, yerr=err_bkg, fmt='o', label='Test Background')

    # Labels
    plt.xlabel('BDT Output Score')
    plt.ylabel('Probability Density')
    plt.title('XGBoost Overtraining Check')
    plt.legend()
    plt.grid(alpha=0.3)
    outfile = dir_map[category] + f"overtraining_{category}{labelForPNG}.png"
    plt.savefig(outfile, dpi=200)
    plt.close()

    print(f"📊 Saved overtraining to: {outfile}")

    # Plot Over Training curve
    results = bdt.evals_result()

    # Plot
    plt.plot(results['validation_0']['auc'], label='Train')
    plt.plot(results['validation_1']['auc'], label='Validation')
#    plt.ylim(0.6, 1.0)
    plt.xlabel('Boosting Rounds')
    plt.ylabel('AUC')
    plt.legend()
    plt.title('XGBoost Overfitting Plot')
    outfile = dir_map[category] + f"overtraining_AUCvsRound_{category}{labelForPNG}.png"
    plt.savefig(outfile, dpi=200)
    plt.close()

    print(f"📊 Saved overtraining to: {outfile}")

def plot_mass_correlations(data, variables, bkg_mask):

    corrs = []

    mass = data.loc[bkg_mask, "HiggsCandCorrMass"]

    for v in variables:

        x = data.loc[bkg_mask, v]

        corr = np.corrcoef(x, mass)[0,1]

        corrs.append((v, corr))

    # sort by absolute correlation
    corrs.sort(key=lambda x: abs(x[1]), reverse=True)

    names  = [x[0] for x in corrs]
    values = [x[1] for x in corrs]

    plt.figure(figsize=(10,6))

    plt.barh(names, values)

    plt.axvline(0, color="black")
    plt.xlabel(r"Corr(variable,$m_{\mu\mu}$)")
    plt.title("Background correlation with dimuon mass")

    plt.tight_layout()

    outfile = (
        dir_map[category]
        + f"mass_correlation_bkg_{category}{labelForPNG}.png"
    )

    plt.savefig(outfile, dpi=200)
    plt.close()

    print(f"Saved {outfile}")

def plot_score_vs_mass(score, sig_mask, bkg_mask, data):


    mass = data["HiggsCandCorrMass"].values

    # Check for NaN values
    valid_mask = ~(np.isnan(score) | np.isnan(mass))
    score = score[valid_mask]
    mass = mass[valid_mask]
    sig_mask = sig_mask[valid_mask]
    bkg_mask = bkg_mask[valid_mask]

    corr_sig = np.corrcoef(score[sig_mask], mass[sig_mask])[0,1]
    corr_bkg = np.corrcoef(score[bkg_mask], mass[bkg_mask])[0,1]

    print("corr(score,mass) signal =", corr_sig)
    print("corr(score,mass) bkg    =", corr_bkg)

    plot_mass_correlations(
        data,
        variables_map[category] + variables_resolution[category],
        bkg_mask
    )

    for name, mask in [("sig", sig_mask), ("bkg", bkg_mask)]:

        plt.figure(figsize=(8,6))

        plt.hist2d(
            mass[mask],
            score[mask],
            bins=[60,50],
            range=[[105,155],[0,1]]
        )

        plt.xlabel(r"$m_{\mu\mu}$ [GeV]")
        plt.ylabel("BDT score")
        plt.title(f"{name}: BDT score vs mass")
        plt.colorbar(label="Events")

        outfile = (
            dir_map[category]
            + f"score_vs_mass_{name}_{category}{labelForPNG}.png"
        )

        plt.tight_layout()
        plt.savefig(outfile, dpi=200)
        plt.close()

        print(f"Saved {outfile}")


def plot_mass_in_score_bins(score, sig_mask, bkg_mask, data, variables):

    mass = data["HiggsCandCorrMass"].values

    for v in variables:

        corr = np.corrcoef(
            data.loc[bkg_mask,v],
        data.loc[bkg_mask,"HiggsCandCorrMass"]
        )[0,1]

        print(f"{v:25s} {corr: .4f}")

    score_bins = np.arange(0, 1.01, 0.1)

    for name, mask in [("sig", sig_mask), ("bkg", bkg_mask)]:

        plt.figure(figsize=(8,6))

        for low, high in zip(score_bins[:-1], score_bins[1:]):

            sel = (
                mask
                & (score >= low)
                & (score < high)
            )

            if np.sum(sel) < 10:
                continue

            plt.hist(
                mass[sel],
                bins=60,
                range=(105,155),
                density=True,
                histtype="step",
                linewidth=1.5,
                label=f"{low:.1f}-{high:.1f}"
            )

        plt.xlabel(r"$m_{\mu\mu}$ [GeV]")
        plt.ylabel("Normalized events")
        plt.title(f"{category} {name}: mass in score bins")
        plt.legend(
            ncol=2,
            fontsize=8
        )

        outfile = (
            dir_map[category]
            + f"mass_in_score_bins_{name}_{category}{labelForPNG}.png"
        )

        plt.tight_layout()
        plt.savefig(outfile, dpi=200)
        plt.close()

        print(f"Saved {outfile}")

def plot_confusion_matrix(y_true, y_pred, test_weights):

    labels = labels_map[category]

    # weighted yields (if desired)
    cm = confusion_matrix(y_true, y_pred)

    # normalize rows
    cm_frac = cm.astype(float) / cm.sum(axis=1)[:, np.newaxis]

    plt.figure(figsize=(8,7))

    disp = ConfusionMatrixDisplay(
        confusion_matrix=cm_frac,
        display_labels=labels,
    )

    disp.plot(
        cmap="Blues",
        values_format=".2f"
    )

    plt.title("Confusion Matrix (row normalized)")

    outfile = dir_map[category] + f"confusion_matrix_fraction_{category}{labelForPNG}.png"

    plt.savefig(outfile, dpi=200)
    plt.close()

    print("saved:", outfile)

def diagnostic(bdt,proba,y_true_binary,variables):                                          

    P_sig = proba[:, signal_classes].sum(axis=1)
    P_bkg = proba[:, background_classes].sum(axis=1)

    y_score = P_sig / (P_sig + P_bkg + 1e-12)
    print("Score range:", y_score.min(), y_score.max())

    # -------------------------------------------------
    # ROC
    # -------------------------------------------------
    roc = roc_auc_score(y_true_binary, y_score)
    print(f"ROC AUC (sig vs all bkg): {roc:.4f}")

    fpr, tpr, thresholds = roc_curve(y_true_binary, y_score)
    roc_auc = auc(fpr, tpr)

    ##--

    # Plot ROC curve
    plt.figure(figsize=(8,6))
    plt.plot(fpr, tpr, color='blue', lw=2, label=f'ROC curve (AUC = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], color='gray', lw=1, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('XGBoost ROC Curve')
    plt.legend(loc="lower right")
    plt.grid(alpha=0.3)
    outfile = dir_map[category] + f"roc_{category}{labelForPNG}.png"
    plt.savefig(outfile, dpi=200)
    plt.close()
    print(f"📊 Saved roc to: {outfile}")

    # Set feature names for export
    booster = bdt.get_booster()
    booster.feature_names = variables
    print("Number of input vars:", len(variables))
    print("Booster feature names:", booster.feature_names)
    importance = booster.get_score(importance_type="gain")
    print("Non-zero features:", importance)

    importance_types = ["gain", "weight", "cover"]
    titles = {
        "gain": "XGBoost Feature Importance — Gain",
        "weight": "XGBoost Feature Importance — Weight (Split Count)",
        "cover": "XGBoost Feature Importance — Cover"
    }

    for imp in importance_types:
        plt.figure(figsize=(8, 6))
        xgb.plot_importance(bdt, importance_type=imp, max_num_features=20)
        plt.title(titles[imp])
#        plt.tight_layout()
        plt.subplots_adjust(left=0.35)

        outfile = dir_map[category] + f"feature_importance_{imp}_{category}{labelForPNG}.png"
        plt.savefig(outfile, dpi=200)
        plt.close()

        print(f"📊 Saved feature importance ({imp}) to: {outfile}")

def plot_score_signed_vs_abs(score, sig_mask, bkg_mask, data):

    # --------

    w_signed = data["weight"].values
    w_abs    = np.abs(w_signed)

    bins = np.linspace(0,1,50)

    plt.figure(figsize=(8,6))

    # signal
    plt.hist(score[sig_mask], bins=bins,
             weights=w_signed[sig_mask],
             density=True,
             histtype='step',
             linewidth=2,
             label="sig signed")

    plt.hist(score[sig_mask], bins=bins,
             weights=w_abs[sig_mask],
             density=True,
             histtype='step',
             linewidth=2,
             linestyle="--",
             label="sig abs")

    # background
    plt.hist(score[bkg_mask], bins=bins,
             weights=w_signed[bkg_mask],
             density=True,
             histtype='step',
             linewidth=2,
             label="bkg signed")

    plt.hist(score[bkg_mask], bins=bins,
             weights=w_abs[bkg_mask],
             density=True,
             histtype='step',
             linewidth=2,
             linestyle="--",
             label="bkg abs")

    plt.xlabel("BDT score")
    plt.ylabel("Normalized entries")
    plt.legend()
    plt.grid(alpha=0.3)

    outfile = dir_map[category]+f"score_signed_vs_abs_{category}{labelForPNG}.png"
    plt.savefig(outfile,dpi=200)
    plt.close()

    print("saved", outfile)

def compare_even_odd(score_even, score_odd,
                     y_even, y_odd,
                     w_even, w_odd):

    sig_mask_even = np.isin(y_even, signal_classes)
    sig_mask_odd  = np.isin(y_odd, signal_classes)

    bkg_mask_even = np.isin(y_even, background_classes)
    bkg_mask_odd  = np.isin(y_odd, background_classes)

    bins = np.linspace(0, 1, 50)

    plt.figure(figsize=(8,6))

    # -------------------------
    # SIGNAL
    # -------------------------
    plt.hist(score_even[sig_mask_even],
             bins=bins,
             weights=w_even[sig_mask_even],
             density=True,
             histtype='step',
             linewidth=2,
             label="Even (signal)")

    plt.hist(score_odd[sig_mask_odd],
             bins=bins,
             weights=w_odd[sig_mask_odd],
             density=True,
             histtype='step',
             linewidth=2,
             linestyle="--",
             label="Odd (signal)")

    # -------------------------
    # BACKGROUND
    # -------------------------
    plt.hist(score_even[bkg_mask_even],
             bins=bins,
             weights=w_even[bkg_mask_even],
             density=True,
             histtype='step',
             linewidth=2,
             label="Even (bkg)")

    plt.hist(score_odd[bkg_mask_odd],
             bins=bins,
             weights=w_odd[bkg_mask_odd],
             density=True,
             histtype='step',
             linewidth=2,
             linestyle="--",
             label="Odd model (bkg)")

    plt.xlabel("BDT score")
    plt.ylabel("Normalized entries")
    plt.title("Even/Odd training stability")
    plt.legend()
    plt.grid(alpha=0.3)

    outfile = dir_map[category]+f"score_even_vs_odd_{category}{labelForPNG}.png"
    plt.savefig(outfile,dpi=200)
    plt.close()

    print("saved", outfile)
    

def train_one_fold(data, train_mask, variables, verbose ): #, labelForPNG):

    # -------------------------------------------------
    # Split even odd
    # -------------------------------------------------

    variables = [str(v).strip().strip("'").strip('"') for v in variables]
    variables = list(dict.fromkeys(variables))  # also removes duplicates safely    
    
##    train_mask = (data["event"] % 2 == 0)
    test_mask  = ~train_mask

    print("Train even:", (data.loc[train_mask,"event"] % 2).unique())
    print("Test odd :", (data.loc[test_mask,"event"] % 2).unique())
    
    train_data   = data.loc[train_mask, variables]
    test_data    = data.loc[test_mask,  variables]

    train_labels = data.loc[train_mask, "target"]
    test_labels  = data.loc[test_mask,  "target"]

    train_weights = data.loc[train_mask, "weight_balanced"]
    test_weights  = data.loc[test_mask,  "weight_balanced"]

    # ✅ Keep test_data_df for later plotting with proper indexing
    test_data_df = data.loc[test_mask].copy()

    sig_mask = np.isin(test_labels, signal_classes)
    bkg_mask = np.isin(test_labels, background_classes)

    n_signal = sig_mask.sum()
    n_background = bkg_mask.sum()
    n_total = len(train_labels)

    print(f"Train set:")
    print(f"  Total:      {n_total}")
    print(f"  Signal (0): {n_signal} ({100*n_signal/n_total:.1f}%)")
    print(f"  Background: {n_background} ({100*n_background/n_total:.1f}%)")
    print(f"  Imbalance ratio (bkg/sig): {n_background/max(n_signal,1):.2f}")

    # -------------------------------------------------
    # TRAINING: BETTER PARAMETERS FOR CLASS IMBALANCE
    # -------------------------------------------------
    ## NOTE: scale_pos_weight is WRONG in multiclass

    params_local = params.copy()

    if not doMultiClass:
        scale_pos_weight = n_background / max(n_signal, 1)
        params_local["scale_pos_weight"] = scale_pos_weight
        print(f"✅ Scale pos weight: {scale_pos_weight:.2f}")
    else:
        # Ensure it's not set for multiclass
        params_local.pop("scale_pos_weight", None)

    print("\n" + "="*60)
    print("CLASS BALANCE DIAGNOSTICS")
    print("="*60)


    print(f"\nWeighted balance:")
    sig_w_sum = test_weights[sig_mask].sum()
    bkg_w_sum = test_weights[bkg_mask].sum()
    print(f"  Signal weight sum:     {sig_w_sum:.6f}")
    print(f"  Background weight sum: {bkg_w_sum:.6f}")
    print(f"  Weight ratio: {bkg_w_sum/max(sig_w_sum,1e-8):.2f}")


    print(f"\nLabel dtype: {train_labels.dtype}")
    print(f"Label unique: {np.unique(train_labels)}")
    print(f"Weights min/max: {train_weights.min():.2e} / {train_weights.max():.2e}")
    print("="*60 + "\n")

    # -------------------------------------------------
    # DIAGNOSTIC
    # -------------------------------------------------

    if verbose:
        print("Train/test sizes:", train_data.shape, test_data.shape)
        print("Train signal fraction:", train_labels.mean(),
              "Test signal fraction:", test_labels.mean())

        print("Train labels unique:", np.unique(train_labels))
        print("Test labels unique:", np.unique(test_labels))
        print("Train positive fraction:", train_labels.mean())
        print("Test positive fraction:", test_labels.mean())

        # Basic statistics
        for v in variables:
            col = train_data[v]  # Series
            print(
                f"{v:20s}",
                "min =", col.min(),
                "max =", col.max(),
                "std =", col.std(),
            )

        '''
        KS      Interpretation
        < 0.05  No separation
        0.05–0.15  Very weak
        0.15–0.30  Some discrimination
        > 0.30  Strong
        '''

        # KS test between signal and background
        for v in variables:
            sig = train_data.loc[train_labels == 1, v]
            bkg = train_data.loc[train_labels == 0, v]

            ks_stat = scipy.stats.ks_2samp(sig, bkg).statistic

            if ks_stat < 0.05: comment='No separation'
            elif ks_stat > 0.05 and ks_stat < 0.15: comment='very weak'
            elif ks_stat > 0.15 and ks_stat < 0.30: comment='some discrimination'
            elif ks_stat > 0.30: comment='STRONG'

            print(
                f"{v:20s}",
                "⟨sig⟩ =", sig.mean(),
                "⟨bkg⟩ =", bkg.mean(),
                "KS =", ks_stat,
                "==> ",comment
            )

        # Sanity checks
        print("Any NaNs in train_data:", train_data.isna().any().any())
        print("Any infs in train_data:", np.isinf(train_data.values).any())

        print("Mean train weight:", train_weights.mean())
        print("Min train weight :", train_weights.min())
        print("Max train weight :", train_weights.max())

    # -------------------------------------------------
    # TRAINING
    # -------------------------------------------------

    train_data = train_data.to_numpy()
    test_data = test_data.to_numpy()
    train_weights = train_weights.to_numpy()

    train_labels = train_labels.to_numpy()
    test_labels = test_labels.to_numpy()
    test_weights = test_weights.to_numpy()

    eval_set = [(train_data, train_labels), (test_data, test_labels)]
    bdt = xgb.XGBClassifier(**params_local)

    bdt.fit(train_data, train_labels,
            sample_weight=train_weights,
            verbose=True,
            eval_set=eval_set)

    print("\n" + "="*60)
    print("TRAINING IN PROGRESS")
    print("="*60)

    print("Training complete.")

    if verbose:
        if doMultiClass: fOutName = f"output/classification_model_multiclass_{category}.root"
        else: fOutName = f"output/classification_model_{category}.root"
        model_name = f"bdt_model_{category}"
        print("variables",variables)
        print("Export model ",model_name)

        ROOT.TMVA.Experimental.SaveXGBoost(bdt, model_name, fOutName, num_inputs=len(variables))
        print(f"output written to {fOutName} with name {model_name}")

        variables_ = ROOT.TList()
        for var in variables:
            print(var)
            variables_.Add(ROOT.TObjString(var))
        fOut = ROOT.TFile(fOutName, "UPDATE")
        fOut.WriteObject(variables_, "variables")
        fOut.Close()
        print('FILE SAVED')

    # -------------------------------------------------
    # DIAGNOSTIC
    # -------------------------------------------------

    if verbose:
        ## only the main signal
        overtraining(bdt, train_data, train_labels, train_weights, test_data, test_labels, test_weights)

    proba = bdt.predict_proba(test_data)

    if doMultiClass:
        y_true_binary = np.isin(test_labels, signal_classes).astype(int)
        y_pred = np.argmax(proba, axis=1)
        plot_confusion_matrix(test_labels, y_pred, test_weights)
    else:
        y_true_binary = (test_labels == 1).astype(int)

    diagnostic(bdt,proba,y_true_binary,variables)

    # -------------------------------------------------
    # FINAL SCORE
    # -------------------------------------------------

    P_sig = proba[:, signal_classes].sum(axis=1)
    P_bkg = proba[:, background_classes].sum(axis=1)

    score = P_sig

    # -------------------------------------------------
    # EXTRA FUNCTIONS diagnostic
    # -------------------------------------------------

    if verbose:
        plot_score_vs_mass(score, sig_mask, bkg_mask, test_data_df)
        plot_mass_in_score_bins(score, sig_mask, bkg_mask, test_data_df, variables)
        plot_score_signed_vs_abs(score, sig_mask, bkg_mask, test_data_df)

    return score,test_labels,test_weights,bdt


def correlation_diagnostics(train_variables):

    print("\n" + "="*80)
    print("CORRELATION DIAGNOSTICS")
    print("="*80)

    variables_hm = (
        train_variables
        + ["HiggsCandCorrMass"]
        + ["HiggsCandMassErr"]
    )

    # --------------------------
    # load samples
    # --------------------------

    if doMultiClass:
        # Load all classes defined in labels_map
        num_classes = len(labels_map[category])
        class_ids_to_load = list(range(num_classes))
    else:
        # For binary: load only signal (1) and background (0)
        class_ids_to_load = [0, 1]

    class_samples = {}
    for c in class_ids_to_load:
        print(f"  Loading class {c}...")
        class_samples[c] = load_process_class(c, variables_hm, False if not doMultiClass else None)

    data = pd.concat(
        [class_samples[c] for c in class_ids_to_load],
        ignore_index=True
    )

    # Create subsamples
    samples = {
        "all": data.copy(),
        "sig": pd.concat(
            [class_samples[c] for c in class_ids_to_load if c in signal_classes],
            ignore_index=True
        ),
        "bkg": pd.concat(
            [class_samples[c] for c in class_ids_to_load if c in background_classes],
            ignore_index=True
        ),
    }

    sig_mask = data["target"].isin(signal_classes)
    bkg_mask = data["target"].isin(background_classes)


    # --------------------------
    # save separately
    # --------------------------

    for sample_name, df in samples.items():

        print(f"\n--- {sample_name.upper()} ---")

        corr_df = df[variables_hm].corr()

        out_corr = (f"output/correlation_matrix{category}_{sample_name}.txt")
        corr_df.round(3).to_csv(out_corr,sep="\t")
        print(f"📐 saved {out_corr}")

        plotCovMatrix(out_corr,sample_name)

        # ----------------------
        # correlation wrt mass
        # ----------------------

        print("\nCorrelation with Higgs mass")

        mass_corr = (
            corr_df["HiggsCandCorrMass"]
            .drop("HiggsCandCorrMass")
            .sort_values(
                key=lambda x: np.abs(x),
                ascending=False
            )
        )

        for v,c in mass_corr.items():
            print(f"{v:25s} {c:8.4f}")

        print("\nCorrelation with mass error")

        err_corr = (
            corr_df["HiggsCandMassErr"]
            .drop("HiggsCandMassErr")
            .sort_values(
                key=lambda x: np.abs(x),
                ascending=False
            )
        )

        for v,c in err_corr.items():
            print(f"{v:25s} {c:8.4f}")

    # --------------------------
    # feature diagnostics
    # --------------------------

    print("\n=== FEATURE DIAGNOSTICS ===")

    for v in train_variables:

        x = samples["all"][v].values

        print(f"\n{v}")
        print("  min:", np.min(x))
        print("  max:", np.max(x))
        print("  std:", np.std(x))
        print("  unique:", len(np.unique(x)))



def build_training_weights(data,
                           doMultiClass=False,
                           flatten_background=True,
                           signal_sigma_power=1.0,
                           signal_relative_sigma=True):

    w = data["weight"].abs().to_numpy().clip(min=1e-10)

    weights = w.copy()

def build_training_weights(data,
                           doMultiClass=False,
                           flatten_background=True,
                           signal_sigma_power=1.0,
                           signal_relative_sigma=True):

    w = data["weight"].abs().to_numpy().clip(min=1e-10)
    weights = w.copy()

    if doMultiClass:
        nclass = len(np.unique(data["target"]))

        # ✅ CORRECT: Weight inversely by class frequency
        class_weights = {}
        total_weight = 0

        for k in range(nclass):
            mask = data["target"] == k
            class_sum = weights[mask].sum()
            class_weights[k] = 1.0 / max(class_sum, 1e-10)
            total_weight += class_weights[k]

        # Normalize so average weight across all classes = 1
        for k in range(nclass):
            class_weights[k] /= total_weight / nclass

        # Apply per-class weights
        for k in range(nclass):
            mask = data["target"] == k
            weights[mask] *= class_weights[k]

        # Final normalization to mean=1
        weights /= weights.mean()

        print("✅ Multiclass weights per class:")
        for k in range(nclass):
            mask = data["target"] == k
            print(f"   Class {k}: weight sum = {weights[mask].sum():.2e}")

        return weights

    ###############################################################
    # for Binary

    sig_mask = data["target"] == 1
    bkg_mask = data["target"] == 0

    ###############################################################
    # Signal weighting
    ###############################################################

    sigma = data["HiggsCandMassErr"].to_numpy()

    sigma_safe = np.where(
        (sigma <= 0) | (~np.isfinite(sigma)),
        1e-3,
        sigma,
    )

    if signal_relative_sigma:

        mass = data["HiggsCandCorrMass"].to_numpy()

        mass_safe = np.where(
            (mass <= 0) | (~np.isfinite(mass)),
            125.,
            mass,
        )

        sigma_used = sigma_safe / mass_safe

    else:

        sigma_used = sigma_safe

    weights[sig_mask] *= (
        1.0 / np.power(sigma_used[sig_mask], signal_sigma_power)
    )

    ###############################################################
    # Background flattening
    ###############################################################

    if flatten_background:

        mass = data["HiggsCandCorrMass"].to_numpy()

        bins = np.linspace(110,150,41)

        hist, edges = np.histogram(
            mass[bkg_mask],
            bins=bins,
            weights=w[bkg_mask]
        )

        #
        # smooth to avoid huge fluctuations
        #
        hist = gaussian_filter1d(hist.astype(float), sigma=1.5)

        hist = np.maximum(hist,1e-6)

        inv = 1./hist

        #
        # normalize so average weight = 1
        #
        inv *= hist.sum()/np.sum(hist*inv)

        ibin = np.clip(
            np.digitize(mass,edges)-1,
            0,
            len(inv)-1
        )

        weights[bkg_mask] *= inv[ibin][bkg_mask]

    ###############################################################
    # Final normalization
    ###############################################################

    weights[sig_mask] /= weights[sig_mask].mean()
    weights[bkg_mask] /= weights[bkg_mask].mean()

    return weights

def _test_XGB_class(label):

    # -------------------------------------------------
    # Main training dataset (no Higgs mass)
    # -------------------------------------------------

    train_variables = variables_map[category] + variables_resolution[category]

    aux_variables = ["HiggsCandCorrMass", "HiggsCandMassErr", "event"]

    variables = train_variables + aux_variables

    # --------------------------
    # load samples
    # --------------------------

    if doMultiClass:
        num_classes = len(labels_map[category])
        class_ids_to_load = list(range(num_classes))
    else:
        class_ids_to_load = [0, 1]  # binary: signal and background only

    class_samples = {}
    for c in class_ids_to_load:
        class_samples[c] = load_process_class(c, variables, False if not doMultiClass else None)

    data = pd.concat(
        [class_samples[c] for c in class_ids_to_load],
        ignore_index=True
    )

    ## full dataset
    sig_mask = data["target"].isin(signal_classes)
    bkg_mask = data["target"].isin(background_classes)

    # -------------------------------------------------
    # WEIGHTS: SIMPLER APPROACH
    # -------------------------------------------------

    w = data["weight"].abs().clip(lower=1e-10)

    data["weight_balanced"] = build_training_weights(
        data,
        doMultiClass=doMultiClass,
        flatten_background=False,
        signal_sigma_power=1,
        signal_relative_sigma=True,
    )

    print("\n" + "="*60)
    print("WEIGHT DIAGNOSTICS (AFTER SCALING)")
    print("="*60)

    train_weights_temp = data["weight_balanced"]
    print(f"Scaled weight range: {train_weights_temp.min():.2e} to {train_weights_temp.max():.2e}")
    print(f"Scaled weight mean: {train_weights_temp.mean():.2e}")
    print(f"Signal weight avg: {train_weights_temp[sig_mask].mean():.4f}")
    print(f"Bkg weight avg: {train_weights_temp[bkg_mask].mean():.4f}")
    print(f"Signal weight sum: {train_weights_temp[sig_mask].sum():.2e}")
    print(f"Bkg weight sum: {train_weights_temp[bkg_mask].sum():.2e}")
    print("="*60 + "\n")

    # -------------------------------------------------
    # Correlation check with HiggsCandCorrMass
    # -------------------------------------------------

    correlation_diagnostics(train_variables)

    ###
    # -------------------------

    print("Start training in 632")

    variables = [v.strip().strip("'").strip('"') for v in train_variables]
    score_odd, y_odd, w_odd, bdt_odd  = train_one_fold(data, data["event"] % 2 == 1, train_variables, False ) #, labelForPNG)
    score_even, y_even, w_even, bdt_even = train_one_fold(data, data["event"] % 2 == 0, train_variables, True ) #, labelForPNG)

    # below score is ok, but  not yet ok for multiclass
    compare_even_odd(score_even, score_odd, y_even, y_odd, w_even, w_odd)

if __name__ == "__main__":

    _test_XGB_class("default")
