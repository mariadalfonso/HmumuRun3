import ROOT
import json

import matplotlib
matplotlib.use("Agg")     # IMPORTANT for batch mode / no display
import matplotlib.pyplot as plt
import scipy.stats

import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score, roc_curve, auc

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score

from LoadTree import loadTree
from LoadTree import vv, tt2l, ttV

myDir='/work/submit/mariadlf/HmumuRun3/ROOTFILES/'

years = ['_12022', '_22022', '_12023', '_22023', '_2024']
#category="VBFcat"
category="ggHcat"
#category="VHcat"
#category="VLcat"
#category="TTHcat"
#category="TTLcat"
#category="Zinvcat"
mytree = ROOT.TChain('events')
for year in years:
    mytree = loadTree(mytree, myDir, category, year )

ROOT.gROOT.SetBatch(True)
ROOT.ROOT.EnableImplicitMT()

paramsClass = {
    "lambda": 0, # to make the discriminator larger for VL
    "objective": "binary:logistic",
    "eval_metric": "auc",
    "max_depth": 4,
    "eta": 0.1,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "n_estimators": 400,
    "tree_method": "hist",
    "seed": 42,
    "min_child_weight": 1,
}

signal_map = {
    "VBFcat": ["10"],
    "ggHcat": ["11"],
    "VHcat": ["12","13","14"],
    "VLcat": ["12","13","14"],
    "Zinvcat": ["14"],
    "TTHcat": ["15"],
    "TTLcat": ["15"]
}

bkg_map = {
    "VBFcat": ["100","103","104","109","108","110"] + ["101","99","98"],# DY QCD + EWK
    "ggHcat": ["100","103","104","109","108","110"], # DY (incl and mass binned)
    "VHcat": ["109"] + ["114","115","116","117","122","123","124","125"] + tt2l,
    "VLcat": vv, # diboson
    "Zinvcat": tt2l, # top 2l
    "TTHcat": tt2l + ["109"], # ttbar2l both powheg nominal and altern ; DY as well
    "TTLcat": tt2l + ttV , # ttbar2l + ttV + tt1l + singleTop
}

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
    "VBFcat": ["HiggsCandCorrPt", "RPt", "Mjj", "dEtaJJ", "ZepVar", "minDetaDiMuVBF", "dPhiJJ", "Muon1_norm_pt", "Muon2_norm_pt"], ##"HiggsCandMassErr"],
    "ggHcat": ["HiggsCandCorrPt", "HiggsCandMassErr", "Muon1_pt", "Muon2_pt"],
    "Zinvcat": ["HiggsCandCorrPt", "Muon1_norm_pt","Muon2_norm_pt","PuppiMET_pt","dPhiMETH","RPt"], ##"HiggsCandMassErr"],
    "VLcat": ["HiggsCandCorrPt", "category","Muon1_norm_pt","Muon2_norm_pt","Lepton_Pt","Lepton2_Pt","VMass","PuppiMET_pt","dPhiVH","dEtaVH"], ##"HiggsCandMassErr"],
    "TTLcat": ["HiggsCandCorrPt", "category","Muon1_norm_pt","Muon2_norm_pt","Lepton_Pt","Lepton2_Pt","PuppiMET_pt","HT","dEtaLepH","mt","Lepton_MVAid","Lepton2_MVAid"], # very weak variables only
    "VHcat": ["HiggsCandCorrPt", "goodWjj_discr", "goodWjj_pt", "goodWjj_eta", "goodWjj_mass", "dEtaWjjH","dPhiWjjH","Muon1_norm_pt","Muon2_norm_pt","RPt"], ##"HiggsCandMassErr"]
    "TTHcat": ["HiggsCandCorrPt", "HT", "nGoodJetsAll", "PuppiMET_pt","Muon1_norm_pt","Muon2_norm_pt","WTopJetDiscr","category"], ##"HiggsCandMassErr"]
}

# here the variable that can help with the resolution
variables_resolution = {
    "VBFcat": ["HiggsCandCorrRapidity", "cosThetaCS", "phiStarCS"],
#    "ggHcat": ["HiggsCandMassErr","HiggsCandCorrRapidity", "cosThetaCS", "phiStarCS", "Muon1_eta","Muon2_eta"],
    "VHcat": ["HiggsCandCorrRapidity","cosThetaCS", "phiStarCS", "PuppiMET_pt"],
    "VLcat": ["HiggsCandCorrRapidity","cosThetaCS","phiStarCS"],
    "Zinvcat": [ "HiggsCandCorrRapidity","cosThetaCS","phiStarCS"],
    "TTLcat": ["HiggsCandCorrRapidity","cosThetaCS","phiStarCS"],
    "TTHcat": ["HiggsCandCorrRapidity","cosThetaCS","phiStarCS","WTopJetMass","Jet1_Pt","Jet1_Eta"],
}

# here the no discrimination and very weak
variables_notUseful_map = {
    "ggHcat": ["HiggsCandCorrRapidity", "cosThetaCS", "phiStarCS", "Muon1_eta","Muon2_eta", "PuppiMET_pt"],
    "VBFcat": ["HiggsCandCorrRapidity", "HiggsCandMassErr", "cosThetaCS","phiStarCS","Muon1_eta","Muon2_eta","dPhiJJ"],
    "Zinvcat": [ "HiggsCandCorrRapidity","cosThetaCS","phiStarCS","Muon1_eta","Muon2_eta"],
    "TTLcat": ["category","Lepton_Pt","Lepton_sip3d","Jet1_Pt","HiggsCandCorrRapidity", "cosThetaCS", "phiStarCS","Muon1_eta","Muon2_eta","PuppiMET_pt","HiggsCandMassErr"],
    "TTHcat": ["HiggsCandCorrRapidity", "cosThetaCS", "phiStarCS", "Muon1_eta","Muon2_eta","Jet1_Pt","Muon1_sip3d","Muon2_sip3d","WTopJetMass"],
    "VLcat": ["Lepton_sip3d","HiggsCandCorrRapidity","dPhiLepH","cosThetaCS","phiStarCS", "Muon1_eta","Muon2_eta"],
    "VHcat": ["cosThetaCS", "phiStarCS", "Muon1_eta","Muon2_eta","PuppiMET_pt","goodWjj_mass","Muon2_pt"],
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

        c.SaveAs(dir_map[category]+var+"_sig_bkg.png")

        print(f" → Saved {var}_sig_bkg.png")

def plotCovMatrix(matrixTXT):

    # Load the CSV
    df = pd.read_csv(matrixTXT, sep="\t", index_col=0)

    # Create figure
    plt.figure(figsize=(8, 7))

    # Plot heatmap
    im = plt.imshow(df, vmin=-1, vmax=1)

    # Axis labels
    plt.xticks(range(len(df.columns)), df.columns, rotation=45, ha="right")
    plt.yticks(range(len(df.index)), df.index)

    # Add colorbar
    plt.colorbar(im, label="Correlation")

    # Add values inside cells (optional but useful)
    for i in range(len(df.index)):
        for j in range(len(df.columns)):
            plt.text(j, i, f"{df.iloc[i, j]:.2f}",
                     ha="center", va="center", fontsize=8)

    plt.title("Feature Correlation Matrix")
    plt.tight_layout()

    # Save
    plt.savefig(dir_map[category]+"correlation_matrix_"+category+".png", dpi=200)
    plt.close()
    print("✅ Saved correlation_matrix.png")

def load_process_class(class_id, variables, drawPlot=False):
#, target=0, weight=1.):
    
    # Heavy-lifting in C++ and remote access of data
    df = ROOT.RDataFrame(mytree)
    df = df.Filter("(HiggsCandCorrMass>(125-20) and HiggsCandCorrMass<(125+20))","HiggsMass within reasonable range 125+-20")
    df = df.Define("Muon1_norm_pt", "HiggsCandCorrPt>0 ? Muon1_pt/HiggsCandCorrPt: 0.f")
    df = df.Define("Muon2_norm_pt", "HiggsCandCorrPt>0 ? Muon2_pt/HiggsCandCorrPt: 0.f")

    if class_id == 0 and drawPlot: make_plots(df) # only make the plot once

    sig_ids = signal_map[category]
    sig_filter = " || ".join([f"mc == {mid}" for mid in sig_ids])
    bkg_ids = bkg_map[category]
    bkg_filter = " || ".join([f"mc == {mid}" for mid in bkg_ids])

    if class_id == 0:
        df = df.Filter(sig_filter)
    elif class_id == 1:
        df = df.Filter(bkg_filter)
    else:
        raise ValueError("Invalid class_id")

    # Define the weight column (if not already present)
    df = df.Define("weight", "w_allSF")
    nevts = df.Count().GetValue()
    print(f"class {class_id} -- evt counts {nevts} (0 signal and 1,2 for BKG)")

    nWevts = df.Sum("weight").GetValue()
    print(f"class {class_id} -- evt weights {nWevts} (0 signal and 1,2 for BKG)")

    cols = df.AsNumpy(variables + ["weight"] + ["HiggsCandMassErr"])

    # Push data to scipy ecosystem
    pdf = pd.DataFrame(cols)
    pdf['target'] = class_id
    return pdf

def diagnostic(bdt,test_data,test_labels,variables):

    # Probabilities (what you actually want)
    y_score = bdt.predict_proba(test_data)[:, 1]

    print("Score range:", y_score.min(), y_score.max())

    # -------------------------------------------------
    # ROC
    # -------------------------------------------------
    roc = roc_auc_score(y_true_binary, y_score)
    print(f"ROC AUC (sig vs all bkg): {roc:.4f}")

    # Compute ROC curve and AUC
    fpr, tpr, thresholds = roc_curve(test_labels, y_score)
    roc_auc = auc(fpr, tpr)

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
    outfile = dir_map[category] + f"roc.png"
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

        outfile = dir_map[category] + f"feature_importance_{imp}.png"
        plt.savefig(outfile, dpi=200)
        plt.close()

        print(f"📊 Saved feature importance ({imp}) to: {outfile}")

def _test_XGB_class(label):

    # -------------------------------------------------
    # Main training dataset (no Higgs mass)
    # -------------------------------------------------

    variables = variables_map[category] + variables_resolution[category]

    sig_df = load_process_class(0,  variables, False)
    bkg_df = load_process_class(1, variables, False)

    data = data.sample(frac=1, random_state=42)
    data["event"] = np.arange(len(data))

    # -------------------------------------------------
    # Build balanced, positive-definite weights
    # -------------------------------------------------

    # Build balanced weight
    if False:
        data["weight_balanced"] = data["weight"]
    else:
        eps = 1e-8  # small positive number
        sigma = data["HiggsCandMassErr"].values

        # protect against bad values
        sigma_safe = np.where((sigma <= 0) | (~np.isfinite(sigma)), eps, sigma)

        data["weight_balanced"] = data["weight"] / (sigma_safe * sigma_safe)

    # Fix negative weights (MANDATORY)
    n_neg = (data["weight_balanced"] <= 0).sum()
    if n_neg:
        print(f"⚠️  Found {n_neg} non-positive weights → taking absolute value")

        data["weight_balanced"] = data["weight_balanced"].abs().clip(lower=1e-6)

    # -------------------------------------------------
    # Split even odd
    # -------------------------------------------------

    train_mask = (data["event"] % 2 == 0)
    test_mask  = ~train_mask

    train_data   = data.loc[train_mask, variables]
    test_data    = data.loc[test_mask,  variables]

    train_labels = data.loc[train_mask, "target"]
    test_labels  = data.loc[test_mask,  "target"]

    train_weights = data.loc[train_mask, "weight_balanced"]
    test_weights  = data.loc[test_mask,  "weight_balanced"]

    # -------------------------------------------------
    # Correlation check with HiggsCandCorrMass
    # -------------------------------------------------

    variables_hm = variables_map[category] + variables_resolution[category] + ["HiggsCandMassErr"] + ["HiggsCandCorrMass"]

    sig_df_hm = load_process_class(0,  variables_hm, True)
    bkg_df_hm = load_process_class(1, variables_hm, False)

    data_hm = pd.concat([sig_df_hm, bkg_df_hm], ignore_index=True)
    data_hm["event"] = np.arange(len(data_hm))

    train_mask_hm = (data_hm["event"] % 2 == 0)

    train_hm = data_hm.loc[train_mask_hm, variables_hm].to_numpy()

    corr = np.corrcoef(train_hm, rowvar=False)
    corr_df = pd.DataFrame(corr, index=variables_hm, columns=variables_hm).round(3)

    out_corr = f"output/correlation_matrix{category}.txt"
    corr_df.to_csv(out_corr, sep="\t")

    print(f"📐 Correlation matrix saved to {out_corr}")

    plotCovMatrix(out_corr)

    # -------------------------------------------------
    # Class imbalance (MANDATORY)
    # -------------------------------------------------

    print("Sum weights (train):", train_weights.sum())
    print("Signal weight sum  :", train_weights[train_labels==0].sum())
    print("Bkg weight sum     :", train_weights[train_labels==1].sum())

    # -------------------------
    # some diagnostic i.e. of the variables
    # -------------------------
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
        sig = train_data.loc[train_labels == 0, v]
        bkg = train_data.loc[train_labels == 1, v]

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

    train_data = train_data.to_numpy()
    test_data = test_data.to_numpy()
    train_weights = train_weights.to_numpy()

    train_labels = train_labels.to_numpy()
    test_labels = test_labels.to_numpy()
    test_weights = test_weights.to_numpy()

    ###
    # -------------------------

    print("Start training in 632")

    eval_set = [(train_data, train_labels), (test_data, test_labels)]
    bdt = xgb.XGBClassifier(**paramsClass)
    bdt.fit(train_data, train_labels, sample_weight=train_weights, verbose=True, eval_set=eval_set)

    print("Training complete.")

    fOutName = f"output/classification_model_{category}.root"
    model_name = f"bdt_model_{category}"
    print("variables",variables)
    print("Export model ",model_name)

    ROOT.TMVA.Experimental.SaveXGBoost(bdt, model_name, fOutName, num_inputs=len(variables))
    print(f"output written to {fOutName} with name {model_name}")

    ###
    # -------------------------

    # append the variables
    variables_ = ROOT.TList()
    for var in variables:
        print(var)
        variables_.Add(ROOT.TObjString(var))
    fOut = ROOT.TFile(fOutName, "UPDATE")
    fOut.WriteObject(variables_, "variables")
    print('FILE SAVED')

    ## call diagnostic
    diagnostic(bdt,test_data,test_labels,variables)

if __name__ == "__main__":

    _test_XGB_class("default")
