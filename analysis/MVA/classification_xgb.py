import ROOT
import json

import xgboost as xgb
import os
import sys

ROOT.gROOT.SetBatch(True)

#from settings import paramsClass1, paramsClass2, paramsMulti
from settings import *

import matplotlib
matplotlib.use("Agg")     # IMPORTANT for batch mode / no display
import matplotlib.pyplot as plt
import scipy.stats

import pandas as pd
import numpy as np
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score, roc_curve, auc

from sklearn.metrics import confusion_matrix
from sklearn.metrics import ConfusionMatrixDisplay

from LoadTree import loadTree
from LoadTree import vv, tt2l, ttV

myDir='/work/submit/mariadlf/HmumuRun3/ROOTFILES/'


years = ['_12022', '_22022', '_12023', '_22023', '_2024']

category = sys.argv[1]

#category="VBFcat"
#category="ggHcat"
#category="VHcat"
#category="VLcat"
#category="TTHcat"
#category="TTLcat"
#category="Zinvcat"
mytree = ROOT.TChain('events')
for year in years:
    mytree = loadTree(mytree, myDir, category, year )

doMultiClass=False

if category in ["Zinvcat","TTLcat"]:
    params = paramsClass2.copy()
else:
    params = paramsClass1.copy()

signal_map = {
    "VBFcat": ["10"],
    "ggHcat": ["11"],
    "VHcat": ["12","13","14"],
    "VLcat": ["12","13","14"],
    "Zinvcat": ["14"],
    "TTHcat": ["15"],
    "TTLcat": ["15"]
}

class_map = {
    "VBFcat": {
        0: ["10"],                         # VBF
        1: ["11"],                         # ggH
        2: ["100","103","104","109","108","110"],   # BKGA
        3: ["101","99","98"]               # BKGB
    }
}

# note inclusive 100, 103, 104 might disturb
bkg_map = {
    "VBFcat": ["100","103","104","109","108","110"] + ["101","99","98"],# DY QCD + EWK
    "ggHcat": ["100","103","104","109","108","110"], # DY (incl and mass binned)
    "VHcat": ["109"] + ["114","115","116","117","122","123","124","125"] + tt2l,
    "VLcat": vv, # diboson
    "Zinvcat": tt2l, # top 2l
    "TTHcat": tt2l + ["109"], # ttbar2l both powheg nominal and altern ; DY as well
    "TTLcat": tt2l + ttV , # ttbar2l + ttV + tt1l + singleTop
}

#labelForPNG = "_DYQCD"
#labelForPNG = "_DYEWK"
#if doMultiClass: labelForPNG = "_multiclass"
#else: labelForPNG = ""
labelForPNG = "_Sigma1"

dir_map = {
    "VBFcat": "/home/submit/mariadlf/public_html/HMUMU_MVA/VBF/",
    "ggHcat": "/home/submit/mariadlf/public_html/HMUMU_MVA/ggH/",
    "VHcat": "/home/submit/mariadlf/public_html/HMUMU_MVA/VH/",
    "VLcat": "/home/submit/mariadlf/public_html/HMUMU_MVA/VL/",
    "Zinvcat": "/home/submit/mariadlf/public_html/HMUMU_MVA/Zinv/",
    "TTHcat": "/home/submit/mariadlf/public_html/HMUMU_MVA/TTH/",
    "TTLcat": "/home/submit/mariadlf/public_html/HMUMU_MVA/TTL/"
}

variables_map = {
    "VBFcat": ["HiggsCandCorrPt", "RPt", "Mjj", "dEtaJJ", "ZepVar", "minDetaDiMuVBF", "dPhiJJ", "Muon1_norm_pt", "Muon2_norm_pt","jetVBF2_Pt","jetVBF1_Pt","jetVBF1_Eta", "jetVBF2_Eta","CenEta","CenPt"],#"minDphiDiMuVBF"],
    "ggHcat": ["HiggsCandCorrPt", "Muon1_norm_pt", "Muon2_norm_pt", "nGoodJetsAll","Jet1_Pt"],
    "Zinvcat": ["HiggsCandCorrPt", "Muon1_norm_pt","Muon2_norm_pt","PuppiMET_pt","dPhiMETH","RPt"],
    "VLcat": ["HiggsCandCorrPt", "category","Muon1_norm_pt","Muon2_norm_pt","Lepton_Pt","Lepton2_Pt","VMass","PuppiMET_pt","dPhiVH","dEtaVH","RPt"],
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

def plotCovMatrix(matrixTXT):

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

    plt.savefig(f"{dir_map[category]}/correlation_matrix_{category}{labelForPNG}.png", dpi=200)
    plt.close()

    print("✅ Saved correlation_matrix.png")

def load_process_class(class_id, variables, drawPlot=False):

    # Heavy-lifting in C++ and remote access of data
    df = ROOT.RDataFrame(mytree)
    df = df.Filter("(HiggsCandCorrMass>(125-20) and HiggsCandCorrMass<(125+20))","HiggsMass within reasonable range 125+-20")
    # FOR VH - VL - TTH - TTL
    #df = df.Filter("(HiggsCandCorrMass>(125-25) and HiggsCandCorrMass<(125+25))","HiggsMass within reasonable range 125+-25")
#    df = df.Define("Muon1_norm_pt", "HiggsCandCorrPt>0 ? Muon1_pt/HiggsCandCorrMass: 0.f")
#    df = df.Define("Muon2_norm_pt", "HiggsCandCorrPt>0 ? Muon2_pt/HiggsCandCorrMass: 0.f")
    df = df.Define("Muon1_norm_pt", "HiggsCandCorrPt>0 ? Muon1_pt/HiggsCandCorrPt: 0.f")
    df = df.Define("Muon2_norm_pt", "HiggsCandCorrPt>0 ? Muon2_pt/HiggsCandCorrPt: 0.f")
#    df = df.Define("log_Mjj", "log(1.+Mjj)")
#    df = df.Define("log_HT", "log(1.+HT)")
#    df = df.Define("log_MET", "log(1.+PuppiMET_pt)")
#    df = df.Define("log_HiggsPt", "log(1.+HiggsCandCorrPt)")


    if doMultiClass:
        ids = class_map[category][class_id]
        filt = " || ".join([f"mc == {x}" for x in ids])

        df = df.Filter(filt)

    else:
        if class_id == 1 and drawPlot: make_plots(df) # only make the plot once

        sig_ids = signal_map[category]
        bkg_ids = bkg_map[category]

        if class_id == 1:
            sig_filter = " || ".join([f"mc == {mid}" for mid in sig_ids])
            df = df.Filter(sig_filter)
        elif class_id == 0 and not doMultiClass:
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

def plot_score_vs_mass(bdt, data, variables):

    # predict signal score
    if doMultiClass:
        score = bdt.predict_proba(data[variables].to_numpy())[:,0]
    else:
        score = bdt.predict_proba(data[variables].to_numpy())[:,1]

    mass = data["HiggsCandCorrMass"].values

    plt.figure(figsize=(8,6))

    plt.hist2d(
        mass,
        score,
        bins=[60,50],
        range=[[105,145],[0,1]]
    )

    plt.xlabel(r"$m_{\mu\mu}$ [GeV]")
    plt.ylabel("BDT score")
    plt.title("BDT score vs Higgs mass")
    plt.colorbar(label="Events")

    outfile = dir_map[category] + f"score_vs_mass_{category}{labelForPNG}.png"
    plt.tight_layout()
    plt.savefig(outfile,dpi=200)
    plt.close()

    print(f"📊 Saved: {outfile}")

def plot_confusion_matrix(bdt, X_test, y_test):

    y_pred = bdt.predict(X_test)

    labels = ["qqH", "ggH", "DY-QCD", "DY-EWK"]

    # weighted yields (if desired)
    cm = confusion_matrix(y_test, y_pred)

    # normalize rows
    cm_frac = cm.astype(float) / cm.sum(axis=1)[:, np.newaxis]

    plt.figure(figsize=(8,7))

    disp = ConfusionMatrixDisplay(
        confusion_matrix=cm_frac,
        display_labels=labels
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

#def diagnostic(bdt,test_data,test_labels,variables):
def diagnostic(bdt,proba,y_true_binary,variables):                                          

    # Probabilities (what you actually want)
#    proba = bdt.predict_proba(test_data)

    if doMultiClass:
        P_vbf  = proba[:,0]
        P_ggh  = proba[:,1]
        P_bkga = proba[:,2]
        P_bkgb = proba[:,3]

        y_score = (P_vbf + P_ggh) / (P_vbf + P_ggh + P_bkga + P_bkgb)
        y_score = P_vbf / (P_vbf + P_ggh + P_bkga + P_bkgb)
#        y_score = P_sig / (P_sig + P_bkg1 + P_bkg2 + 1e-12)

    else:
        P_sig  = proba[:, 1]
        P_bkg = proba[:, 0]
        y_score = P_sig / (P_sig + P_bkg + 1e-12)

    print("Score range:", y_score.min(), y_score.max())
    # Convert labels → binary (signal vs all)
#    y_true_binary = (test_labels == 1).astype(int)

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

def plot_score_signed_vs_abs(bdt, data, variables):

    #"this is FULL dataset (sig + bkg + even + odd mixed)"
    X = data[variables].to_numpy()
    score = bdt.predict_proba(X)[:,0] if doMultiClass else \
        bdt.predict_proba(X)[:,1]

    sig = data["target"] == 1
    bkg = data["target"] == 0

    w_signed = data["weight"].values
    w_abs    = np.abs(w_signed)

    bins = np.linspace(0,1,50)

    plt.figure(figsize=(8,6))

    # signal
    plt.hist(score[sig], bins=bins,
             weights=w_signed[sig],
             density=True,
             histtype='step',
             linewidth=2,
             label="sig signed")

    plt.hist(score[sig], bins=bins,
             weights=w_abs[sig],
             density=True,
             histtype='step',
             linewidth=2,
             linestyle="--",
             label="sig abs")

    # background
    plt.hist(score[bkg], bins=bins,
             weights=w_signed[bkg],
             density=True,
             histtype='step',
             linewidth=2,
             label="bkg signed")

    plt.hist(score[bkg], bins=bins,
             weights=w_abs[bkg],
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

    bins = np.linspace(0, 1, 50)

    plt.figure(figsize=(8,6))

    # -------------------------
    # SIGNAL
    # -------------------------
    plt.hist(score_even[y_even == 1],
             bins=bins,
             weights=w_even[y_even == 1],
             density=True,
             histtype='step',
             linewidth=2,
             label="Even (signal)")

    plt.hist(score_odd[y_odd == 1],
             bins=bins,
             weights=w_odd[y_odd == 1],
             density=True,
             histtype='step',
             linewidth=2,
             linestyle="--",
             label="Odd (signal)")

    # -------------------------
    # BACKGROUND
    # -------------------------
    plt.hist(score_even[y_even == 0],
             bins=bins,
             weights=w_even[y_even == 0],
             density=True,
             histtype='step',
             linewidth=2,
             label="Even (bkg)")

    plt.hist(score_odd[y_odd==0],
             bins=bins,
             weights=w_odd[y_odd == 0],             
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
    

def train_one_fold(data, train_mask, variables, verbose):

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

    print("Sum weights (train):", train_weights.sum())
    print("Signal weight sum  :", train_weights[train_labels==1].sum())
    print("Bkg weight sum     :", train_weights[train_labels==0].sum())

    # -------------------------------------------------
    # TRAINING: BETTER PARAMETERS FOR CLASS IMBALANCE
    # -------------------------------------------------
    ## NOTE: scale_pos_weight is WRONG in multiclass

    if not doMultiClass:

        n_signal = (train_labels == 1).sum()
        n_background = (train_labels == 0).sum()
        n_total = len(train_labels)

        print(f"Train: {n_signal} signal ({100*n_signal/n_total:.1f}%), {n_background} bkg ({100*n_background/n_total:.1f}%)")
        print(f"Imbalance: {n_background/n_signal:.3f}")

        scale_pos_weight = n_background / max(n_signal, 1)

        print(f"\n✅ Signal samples: {n_signal}, Background: {n_background}")
        print(f"✅ Scale pos weight: {scale_pos_weight:.2f}\n")

        params["scale_pos_weight"] = scale_pos_weight

    print("\n" + "="*60)
    print("CLASS BALANCE DIAGNOSTICS")
    print("="*60)

    n_signal = (train_labels == 1).sum()
    n_background = (train_labels == 0).sum()
    n_total = len(train_labels)

    print(f"Train set:")
    print(f"  Total:      {n_total}")
    print(f"  Signal (0): {n_signal} ({100*n_signal/n_total:.1f}%)")
    print(f"  Background: {n_background} ({100*n_background/n_total:.1f}%)")
    print(f"  Imbalance ratio (bkg/sig): {n_background/max(n_signal,1):.2f}")

    print(f"\nWeighted balance:")
    sig_w_sum = train_weights[train_labels == 1].sum()
    bkg_w_sum = train_weights[train_labels == 0].sum()
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
    # TRAIN
    # -------------------------------------------------

    train_data = train_data.to_numpy()
    test_data = test_data.to_numpy()
    train_weights = train_weights.to_numpy()

    train_labels = train_labels.to_numpy()
    test_labels = test_labels.to_numpy()
    test_weights = test_weights.to_numpy()

    eval_set = [(train_data, train_labels), (test_data, test_labels)]
    bdt = xgb.XGBClassifier(**params)

    bdt.fit(train_data, train_labels,
            sample_weight=train_weights,
            verbose=True,
            eval_set=eval_set)

    print("\n" + "="*60)
    print("TRAINING IN PROGRESS")
    print("="*60)

    print("Training complete.")

    if verbose:
        fOutName = f"output/classification_model_{category}.root"
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
        print('FILE SAVED')
        
    # -------------------------------------------------
    # DIAGNOSTIC
    # -------------------------------------------------

    if verbose:
    
        overtraining(bdt, train_data, train_labels, train_weights, test_data, test_labels, test_weights)

        proba = bdt.predict_proba(test_data)

        if doMultiClass:
            plot_confusion_matrix(bdt, test_data, test_labels)
            y_true_binary = (test_labels == 0).astype(int)

        else:
            y_true_binary = (test_labels == 1).astype(int)
        diagnostic(bdt,proba,y_true_binary,variables)

    if doMultiClass:
        proba = bdt.predict_proba(test_data)

        print("\nMean predicted probs")
        print("VBF :", proba[:,0].mean())
        print("ggH :", proba[:,1].mean())
        print("BKGA:", proba[:,2].mean())
        print("BKGB:", proba[:,3].mean())

        score = bdt.predict_proba(test_data)[:,0]

    else:
        score = bdt.predict_proba(test_data)[:,1]

    return score,test_labels,test_weights,bdt

def _test_XGB_class(label):

    # -------------------------------------------------
    # Main training dataset (no Higgs mass)
    # -------------------------------------------------

    train_variables = variables_map[category] + variables_resolution[category]

    aux_variables = ["HiggsCandCorrMass", "HiggsCandMassErr", "event"]

    variables = train_variables + aux_variables

    sig_df = load_process_class(1,  variables, False)
    if doMultiClass:
        c0 = load_process_class(0, variables)
        c1 = load_process_class(1, variables)
        c2 = load_process_class(2, variables)
        c3 = load_process_class(3, variables)
        data = pd.concat([c0,c1,c2,c3], ignore_index=True)

    else:
        bkg_df = load_process_class(0, variables, False)
        data = pd.concat([sig_df, bkg_df], ignore_index=True)

#    data = data.sample(frac=1, random_state=42)
#    data["event"] = np.arange(len(data))

    # -------------------------------------------------
    # WEIGHTS: SIMPLER APPROACH
    # -------------------------------------------------
    if doMultiClass:
        sig_mask = data["target"].values == 0
        bkg1_mask = data["target"].values == 1
        bkg2_mask = data["target"].values == 2
        bkg3_mask = data["target"].values == 3

    else:
        sig_mask = data["target"] == 1
        bkg_mask = data["target"] == 0

    w = data["weight"].abs().clip(lower=1e-10)

    if False:

        # Normalize per class to mean=1
        sig_avg = w[sig_mask].mean()
        bkg_avg = w[bkg_mask].mean()

        data.loc[sig_mask, "weight_balanced"] = w[sig_mask] / sig_avg
        data.loc[bkg_mask, "weight_balanced"] = w[bkg_mask] / bkg_avg

    else:

        if True:
            # Apply mass error weighting: 1 / sigma^2
            sigma = data["HiggsCandMassErr"].values
            sigma_safe = np.where((sigma <= 0) | (~np.isfinite(sigma)), 1e-3, sigma)
            w_with_sigma = w / (sigma_safe ** 2)
#            w_with_sigma = w / (sigma_safe)

        else:
            sigma = data["HiggsCandMassErr"].values
            mass  = data["HiggsCandCorrMass"].values

            sigma_safe = np.where(
                (sigma <= 0) | (~np.isfinite(sigma)),
                1e-3,
                sigma
            )

            mass_safe = np.where(
                (mass <= 0) | (~np.isfinite(mass)),
                125.0,
                mass
            )

            rel_sigma = sigma_safe / mass_safe
            w_with_sigma = w / (rel_sigma ** 2)

        # Weight by mass resolution: events with better resolution get higher weight
        w_with_sigma = w_with_sigma.clip(lower=1e-10)  # Safety

        if doMultiClass:

            w_balanced = np.zeros_like(w)

            for k in range(4):
                mask = data["target"] == k
                avg = w[mask].mean()
                w_balanced[mask] = w[mask] / avg

            data["weight_balanced"] = w_balanced

        else:
            # Normalize per class to mean=1
            sig_avg = w_with_sigma[sig_mask].mean()
            bkg_avg = w_with_sigma[bkg_mask].mean()

            data.loc[sig_mask, "weight_balanced"] = w_with_sigma[sig_mask] / sig_avg
            data.loc[bkg_mask, "weight_balanced"] = w_with_sigma[bkg_mask] / bkg_avg

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

    variables_hm = train_variables + ["HiggsCandCorrMass", "HiggsCandMassErr"]

    sig_df_hm = load_process_class(1,  variables_hm, True)
    bkg_df_hm = load_process_class(0, variables_hm, False) ## NOTE: need to be updated for the multiclasee

    data_hm = pd.concat([sig_df_hm, bkg_df_hm], ignore_index=True)
    data_hm["event"] = np.arange(len(data_hm))

    train_mask_hm = (data_hm["event"] % 2 == 0)

    train_hm = data_hm.loc[train_mask_hm, variables_hm].to_numpy()

    corr = np.corrcoef(train_hm, rowvar=False)
    corr_df = pd.DataFrame(corr, index=variables_hm, columns=variables_hm).round(3)

    out_corr = f"output/correlation_matrix{category}.txt"
    corr_df.to_csv(out_corr, sep="\t")
    # TO DO: this is mixed SIG-BKG correlation, split SIG and BKG

    print(f"📐 Correlation matrix saved to {out_corr}")

    plotCovMatrix(out_corr)

    print("=== FEATURE DIAGNOSTICS ===")

    for i, v in enumerate(variables_hm):
        col = train_hm[:, i]

        print(v)
        print("  min:", np.min(col))
        print("  max:", np.max(col))
        print("  std:", np.std(col))
        print("  unique (approx):", len(np.unique(col)))

    ###
    # -------------------------

    print("Start training in 632")

    variables = [v.strip().strip("'").strip('"') for v in train_variables]
    score_even, y_even, w_even, bdt_even = train_one_fold(data, data["event"] % 2 == 0, train_variables, True)

    score_odd, y_odd, w_odd, bdt_odd  = train_one_fold(data, data["event"] % 2 == 1, train_variables, False)

    # -------------------------------------------------
    # EXTRA FUNCTIONS
    # -------------------------------------------------

    mass_df = pd.concat([sig_df_hm, bkg_df_hm], ignore_index=True)
    plot_score_vs_mass(bdt_even, mass_df, train_variables)

    plot_score_signed_vs_abs(bdt_even, data, train_variables)

    compare_even_odd(score_even, score_odd, y_even, y_odd, w_even, w_odd)
    
if __name__ == "__main__":

    _test_XGB_class("default")
