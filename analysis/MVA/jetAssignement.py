import xgboost as xgb
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score
from sklearn.metrics import roc_auc_score

import ROOT
import pandas as pd
import numpy as np

myDir='/work/submit/mariadlf/HmumuRun3/ROOTFILES/'

ROOT.gROOT.SetBatch(True)

paramsClass = {
    'objective': 'binary:logistic',
    'eval_metric': 'logloss',
    'n_estimators':300,
    'learning_rate':0.1,
    'max_depth':6
}

variables = [
    "jetVBF1_Pt",
    "jetVBF2_Pt",
    "jetVBF1_Eta",    
    "jetVBF2_Eta"
]

def load_process_class(fIn, target=0, weight=1.):

    # Always make sure the truth variables are loaded
    truth_vars = ["jetVBF1_LHE", "jetVBF2_LHE"]
    all_vars = list(set(variables + truth_vars))  # avoid duplicates

    # Heavy-lifting in C++ and remote access of data
    df = ROOT.RDataFrame("events", fIn)
    cols = df.AsNumpy(all_vars)
    pdf = pd.DataFrame(cols)

    pdf["target"] = ((pdf["jetVBF1_LHE"] != -1) & (pdf["jetVBF2_LHE"] != -1)).astype(int)
#    pdf['target'] = target # add a target column to indicate signal (1) and background (0)                                                                                                
    pdf['weight'] = weight
    return pdf



def _test_XGB_assignement():

    sig_df = load_process_class(myDir+"snapshot_mc_10_12022_VBFcat.root")

    # Define X (features) and y (truth)
    X = sig_df[variables].copy()
    y = sig_df["target"]
    
    # X = feature matrix, y = true labels
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

    X_train = X_train.to_numpy()
    X_test = X_test.to_numpy()
    y_train = y_train.to_numpy()
    y_test = y_test.to_numpy()

    # Create model and train
    model = xgb.XGBClassifier(**paramsClass)
    model.fit(X_train, y_train)

    # ------------------------

    y_pred = model.predict(X_test)
    acc = accuracy_score(y_test, y_pred)
    roc = roc_auc_score(y_test, model.predict_proba(X_test)[:, 1])
    print(f"✅ Model trained successfully")
    print(f"Accuracy: {acc:.4f}")
    print(f"ROC AUC: {roc:.4f}")

    # ------------------------

    json_path = "output/truth_assignment_xgb.json"
    model.save_model(json_path)
    fOutName = "output/truth_assignment_xgb.root"
    model_name = "truth_assignment_vbfjets"
        

    ROOT.TMVA.Experimental.SaveXGBoost(model, model_name, fOutName, len(variables))
    print(f"📦 Model exported successfully:")
    print(f"   - JSON: {json_path}")
    print(f"   - ROOT: {fOutName} (TMVA model name: {model_name})")
    

#    variables_ = ROOT.TList()
##    for var in variables+[target_var]:                                                                                                                   
#    for var in variables:
#        print(var)
#        variables_.Add(ROOT.TObjString(var))
#    fOut = ROOT.TFile(fOutName+".root", "UPDATE")
#    fOut.WriteObject(variables_, "variables")
    
    
if __name__ == "__main__":

    _test_XGB_assignement()
