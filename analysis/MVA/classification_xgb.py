import ROOT
import json

import matplotlib
matplotlib.use("Agg")     # IMPORTANT for batch mode / no display
import matplotlib.pyplot as plt

import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score

from LoadTree import loadTree
myDir='/work/submit/mariadlf/HmumuRun3/ROOTFILES/'

years = ['_12022', '_22022', '_12023', '_22023', '_2024']
category="VBFcat"
#category="ggHcat"
#category="VHcat"
#category="VLcat"
#category="TTHcat"
#category="TTLcat"
mytree = ROOT.TChain('events')
for year in years:
    mytree = loadTree(mytree, myDir, category, year )

ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT()

paramsClass = {
    "objective": "binary:logistic",
    "eval_metric": "auc",
    "max_depth": 4,
    "eta": 0.1,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "n_estimators": 400,
}

signal_map = {
    "VBFcat": ["10"],
    "ggHcat": ["11"],
    "VHcat": ["12","13","14"],
    "VLcat": ["12","13","14"],
    "TTHcat": ["15"],
    "TTLcat": ["15"]
}

bkg_map = {
    "VBFcat": ["100","103","109","108","110"], # DY
    "ggHcat": ["100","103","109","108","110"], # DY (incl and mass binned)
    "VHcat": ["100","103","109","108","110","111","112","113","114","115","116","117"],  # DY 11-13 jet binned; 14-17 ptbinned
    "VLcat": ["201","202","203","204","205"], # diboson
    "Zinvcat": ["201","202","203","204","205"], # diboson
    "TTHcat": ["102","118","100","103","109","108","110"], # ttbar2l both powheg nominal and spare
    "TTLcat": ["102","118","105","106","107"], # ttbar2l + tt1l + singleTop
}

dir_map = {
    "VBFcat": "/home/submit/mariadlf/public_html/HMUMU_MVA/VBF/",
    "ggHcat": "/home/submit/mariadlf/public_html/HMUMU_MVA/ggH/",
    "VHcat": "/home/submit/mariadlf/public_html/HMUMU_MVA/VH/",
    "VLcat": "/home/submit/mariadlf/public_html/HMUMU_MVA/VL/",
    "TTHcat": "/home/submit/mariadlf/public_html/HMUMU_MVA/TTH/",
    "TTLcat": "/home/submit/mariadlf/public_html/HMUMU_MVA/TTL/"
}

variables_map = {
    "VBFcat": ["HiggsCandCorrPt", "HiggsCandCorrRapidity", "HiggsCandMassErr", "cosThetaCS","RPt","Mjj","dEtaJJ","ZepVar","dPhiJJ"],
    "ggHcat": ["HiggsCandCorrPt", "HiggsCandCorrRapidity", "HiggsCandMassErr", "cosThetaCS", "phiStarCS", "Muon1_eta","Muon2_eta"],
#    "Zinvcat": []
    "VLcat": ["HiggsCandCorrPt", "HiggsCandMassErr", "cosThetaCS", "category","Lepton_Pt","mt"],
    "VHcat": ["HiggsCandCorrPt", "goodWjj_mass", "goodWjj_discr", "goodWjj_pt","dEtaWjjH"],
    "TTLcat": ["HiggsCandCorrPt", "HiggsCandCorrRapidity", "cosThetaCS", "phiStarCS", "Jet1_Pt","category","Lepton_Pt","mt","Muon1_eta","Muon2_eta"],
    "TTHcat": ["HiggsCandCorrPt", "cosThetaCS", "Jet1_Pt", "HT", "nGoodJetsAll"]
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

    for var in variables_map[category]:
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

        hsig.Scale(1.0 / hsig.Integral())
        hsigABS.Scale(1.0 / hsigABS.Integral())
        hbkg.Scale(1.0 / hbkg.Integral())
        hbkgABS.Scale(1.0 / hbkgABS.Integral())
        hsig.Draw("HIST")
        hsigABS.Draw("HIST SAME")
        hbkg.Draw("HIST SAME")
        hbkgABS.Draw("HIST SAME")

        c.BuildLegend()

        c.SaveAs(dir_map[category]+var+"_sig_bkg.png")

        print(f" → Saved {var}_sig_bkg.png")

def load_process_class(isSignal, variables, drawPlot=False):
#, target=0, weight=1.):
    
    # Heavy-lifting in C++ and remote access of data
    df = ROOT.RDataFrame(mytree)
    df = df.Filter("(HiggsCandCorrMass>(125-20) and HiggsCandCorrMass<(125+20))","HiggsMass within reasonable range 125+-20")
    # FOR VH - VL - TTH - TTL
    #df = df.Filter("(HiggsCandCorrMass>(125-50) and HiggsCandCorrMass<(125+50))","HiggsMass within reasonable range 125+-20")

    if isSignal and drawPlot: make_plots(df) # only make the plot once

    sig_ids = signal_map[category]
    sig_filter = " || ".join([f"mc == {mid}" for mid in sig_ids])
    bkg_ids = bkg_map[category]
    bkg_filter = " || ".join([f"mc == {mid}" for mid in bkg_ids])
    df = df.Filter(sig_filter if isSignal else bkg_filter)
    # Define the weight column (if not already present)
    df = df.Define("weight", "w_allSF")
    nevts = df.Count().GetValue()
    print('case isSignal=',isSignal,' -- nevts',nevts)

    cols = df.AsNumpy(variables + ["weight"])
    # Push data to scipy ecosystem
    pdf = pd.DataFrame(cols)
    pdf['target'] = 1 if isSignal else 0
    return pdf


def diagnostic(bdt,test_data,test_labels,variables):

    y_pred = bdt.predict(test_data)
    acc = accuracy_score(test_labels, y_pred)
    roc = roc_auc_score(test_labels, bdt.predict_proba(test_data)[:, 1])
    print(f"✅ Model trained successfully")
    print(f"Accuracy: {acc:.4f}")
    print(f"ROC AUC: {roc:.4f}")

    print("Number of input vars:", len(variables))
    print("Booster feature names:", bdt.get_booster().feature_names)
    print("Booster feature attributes:", bdt.get_booster().attributes())

    # Set feature names for export
    bdt.get_booster().feature_names = variables
    importance = bdt.get_booster().get_score(importance_type="gain")
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
        plt.tight_layout()

        outfile = dir_map[category] + f"feature_importance_{imp}.png"
        plt.savefig(outfile, dpi=200)
        plt.close()

        print(f"📊 Saved feature importance ({imp}) to: {outfile}")

def _test_XGB_class(label):

    # TO DO: 1) refine the variable

    ## $$$$$$$$$$$$$$$$$
    ## the first block is to check the correlation with the HiggsMass

    variables_plus_HM = variables_map[category] + ["HiggsCandCorrMass"]

    sig_df_hm = load_process_class(True, variables_plus_HM,True)
    bkg1_df_hm = load_process_class(False, variables_plus_HM,False)

    data = pd.concat([sig_df_hm, bkg1_df_hm], ignore_index=True)

    train_data_hm, test_data_hm, train_labels_hm, test_labels_hm, train_weights_hm, test_weights_hm  = train_test_split(
        data[variables_plus_HM], data['target'], data['weight'], test_size=0.2, random_state=42)

    train_data_hm = train_data_hm.to_numpy()

    correlation_matrix = np.corrcoef(train_data_hm, rowvar=False)
    correlation_df = pd.DataFrame(correlation_matrix, index=variables_plus_HM, columns=variables_plus_HM)
    correlation_df_rounded = correlation_df.round(3)
    fOutName = f"output/correlation_matrix{category}.txt"
    correlation_df_rounded.to_csv(fOutName, sep='\t')
    print("Correlation matrix saved to correlation_matrix.txt")

    ## $$$$$$$$$$$$$$$$$
    ## repeated block w/o the "HiggsCandCorrMass"

    variables = variables_map[category]

    sig_df = load_process_class(True, variables, False)
    bkg1_df = load_process_class(False, variables, False)

    data = pd.concat([sig_df, bkg1_df], ignore_index=True).sample(frac=1, random_state=42)
    data["event"] = data.index

    train_mask = (data["event"] % 2 == 0)
    test_mask  = ~train_mask

    train_data   = data.loc[train_mask, variables]
    test_data    = data.loc[test_mask, variables]

    train_labels = data.loc[train_mask, "target"]
    test_labels  = data.loc[test_mask, "target"]

    # -------------------------
    # adjust weights: sig and BKG the same and how to handle neg weights
    # -------------------------

    # compute total original weights (on the whole dataset or only training later; we do it on full dataset)
    Wsig = data.loc[data.target == 1, "weight"].sum()
    Wbkg = data.loc[data.target == 0, "weight"].sum()

    if Wsig == 0 or Wbkg == 0:
        raise RuntimeError("ERROR: Signal or background total weight = 0!")

    scale_bkg = Wsig / Wbkg
    print(f"  Total signal weight = {Wsig:.3f}")
    print(f"  Total background weight = {Wbkg:.3f}")
    print(f"  Scaling background weights by {scale_bkg:.5f}")

    # build balanced weights
    data["weight_balanced"] = data["weight"]
    data.loc[data.target == 0, "weight_balanced"] *= scale_bkg

    # clip negatives to zero (XGBoost requires non-negative weights)
    n_neg = (data["weight_balanced"] < 0).sum()
    #    if n_neg:
    #        print(f"  Warning: found {n_neg} negative balanced weights -> clipping to 0")
    #    data["weight_balanced"] = data["weight_balanced"].clip(lower=0.0)
    if n_neg:
        print(f"  Warning: found {n_neg} negative balanced weights -> taking absolute value")
    data["weight_balanced"] = data["weight_balanced"].abs()

    train_weights = data.loc[train_mask, "weight_balanced"]
    test_weights  = data.loc[test_mask, "weight_balanced"]

    # -------------------------
    # train
    # -------------------------

    train_data = train_data.to_numpy()
    test_data = test_data.to_numpy()
    train_weights = train_weights.to_numpy()

    train_labels = train_labels.to_numpy()
    test_labels = test_labels.to_numpy()
    test_weights = test_weights.to_numpy()

    print("Train/test sizes:", train_data.shape, test_data.shape)
    print("Train signal fraction:", train_labels.mean(),
          "Test signal fraction:", test_labels.mean())

    eval_set = [(train_data, train_labels), (test_data, test_labels)]

    ###

    print("Start training")
    bdt = xgb.XGBClassifier(**paramsClass)
    bdt.fit(train_data, train_labels, sample_weight=train_weights, verbose=True, eval_set=eval_set)

    print("Training complete.")

    fOutName = f"output/classification_model_{category}.root"
    model_name = f"bdt_model_{category}"
    print("variables",variables)
    print("Export model ",model_name)

    ROOT.TMVA.Experimental.SaveXGBoost(bdt, model_name, fOutName, num_inputs=len(variables))
    print(f"output written to {fOutName} with name {model_name}")

    # append the variables

    variables_ = ROOT.TList()
    for var in variables:
        print(var)
        variables_.Add(ROOT.TObjString(var))
    fOut = ROOT.TFile(fOutName, "UPDATE")
    fOut.WriteObject(variables_, "variables")

    ## call diagnostic
    diagnostic(bdt,test_data,test_labels,variables)

if __name__ == "__main__":

    _test_XGB_class("default")
