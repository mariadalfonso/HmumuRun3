import ROOT
import json

import pandas as pd
import numpy as np
import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, roc_auc_score

from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import cross_val_score

from LoadTree import loadTree
myDir='/work/submit/mariadlf/HmumuRun3/ROOTFILES/'

year = '_22023'
category="_VBFcat"
mytree = ROOT.TChain('events')
mytree = loadTree(mytree, myDir, category, year )

ROOT.gROOT.SetBatch(True)

weight_sf = 1e9

paramsClass = {
    'objective': 'binary:logistic',
    'eval_metric': 'logloss',
    'eta': 0.1,
    'max_depth': 4,
    'subsample': 0.5,
    'colsample_bytree': 0.5,
    'seed': 42,
    'n_estimators': 4000,
    'early_stopping_rounds': 25,
    'num_rounds': 20,
    'learning_rate': 0.1,
    'gamma': 3,
    'min_child_weight': 10,
    'max_delta_step': 0,
}

variables = [
    "RPt",
    "Mjj",
    "dEtaJJ",
    "ZepVar"
]

variables_map = {
    "_VBFcat": ["RPt","Mjj","dEtaJJ","ZepVar"]
#    "_ggHcat": []
#    "_Zinvcat": []
#    "_VLcat": []
#    "_VHcat": []
#    "_TTLcat": []
#    "_TTHcat": []
}


def load_process_class(mySel, variables, target=0, weight=1.):

    
    # Heavy-lifting in C++ and remote access of data
    df = ROOT.RDataFrame(mytree)
    df = df.Filter("(HiggsCandCorrMass>(125-15) and HiggsCandCorrMass<(125+15))","HiggsMass within reasonable range 125+-15")
    df = df.Filter('{} ? mc==10: mc>50'.format(mySel))
    # Define the weight column (if not already present)
    df = df.Define("weight", "w_allSF")
    nevts = df.Count().GetValue()
    print('case=',mySel,' -- nevts',nevts, 'target = ', target)
    
    cols = df.AsNumpy(variables + ["weight"])
    # Push data to scipy ecosystem
    pdf = pd.DataFrame(cols)
    pdf['target'] = target # add a target column to indicate signal (1) and background (0)
#    pdf['weight'] = weight
    return pdf


def _test_XGB_class(label):

    # TO DO: 1) refine the variable 2) split even and odd

    # ROC AUC: 0.5000; Accuracy: 0.1743; variables ['RPt', 'Mjj', 'dEtaJJ', 'ZepVar']
    variables_plus_HM = variables_map[category] + ["HiggsCandCorrMass"]

    #####
    sig_df_hm = load_process_class(1, variables_plus_HM, weight=weight_sf, target=1)
    bkg1_df_hm = load_process_class(0, variables_plus_HM, weight=weight_sf)
    
    data = pd.concat([sig_df_hm, bkg1_df_hm], ignore_index=True)

    train_data_hm, test_data_hm, train_labels_hm, test_labels_hm, train_weights_hm, test_weights_hm  = train_test_split(
        data[variables_plus_HM], data['target'], data['weight'], test_size=0.2, random_state=42)

    train_data_hm = train_data_hm.to_numpy()
    #####
    
    correlation_matrix = np.corrcoef(train_data_hm, rowvar=False)
    correlation_df = pd.DataFrame(correlation_matrix, index=variables_plus_HM, columns=variables_plus_HM)
    correlation_df_rounded = correlation_df.round(3)
    correlation_df_rounded.to_csv("output/correlation_matrix.txt", sep='\t')
    print("Correlation matrix saved to correlation_matrix.txt")
    
    ### repeated block
    sig_df_hm = load_process_class(1, variables_plus_HM, weight=weight_sf, target=1)
    bkg1_df_hm = load_process_class(0, variables_plus_HM, weight=weight_sf)
    
    data = pd.concat([sig_df_hm, bkg1_df_hm], ignore_index=True)
    
    train_data, test_data, train_labels, test_labels, train_weights, test_weights  = train_test_split(
        data[variables], data['target'], data['weight'], test_size=0.2, random_state=42)

    #XGBoost cannot handle negative sample weights
    train_weights = train_weights.abs()
    test_weights = test_weights.abs()

    train_data = train_data.to_numpy()
    test_data = test_data.to_numpy()
    train_labels = train_labels.to_numpy()
    test_labels = test_labels.to_numpy()
    train_weights = train_weights.to_numpy()
    test_weights = test_weights.to_numpy()

    eval_set = [(train_data, train_labels), (test_data, test_labels)]
    ###

    print("Start training")
    bdt = xgb.XGBClassifier(**paramsClass)
    bdt.fit(train_data, train_labels, verbose=True, eval_set=eval_set, sample_weight=train_weights)
    
    y_pred = bdt.predict(test_data)
#    y_pred = bdt.predict(train_labels)
    acc = accuracy_score(test_labels, y_pred)
    roc = roc_auc_score(test_labels, bdt.predict_proba(test_data)[:, 1])
    print(f"✅ Model trained successfully")
    print(f"Accuracy: {acc:.4f}")
    print(f"ROC AUC: {roc:.4f}")
    
    fOutName = f"output/classification_model{category}.root"
    model_name = f"bdt_model{category}"
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

if __name__ == "__main__":

    _test_XGB_class("default")
