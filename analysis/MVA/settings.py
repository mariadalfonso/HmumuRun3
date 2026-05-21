import ROOT
import os

num_cores=64
ROOT.ROOT.EnableImplicitMT(num_cores)

# Set threading via environment variables (works in all versions)

os.environ['OMP_NUM_THREADS'] = str(num_cores)
os.environ['OPENBLAS_NUM_THREADS'] = str(num_cores)
os.environ['MKL_NUM_THREADS'] = str(num_cores)
os.environ['NUMEXPR_NUM_THREADS'] = str(num_cores)

print(f"XGBoost threading configured for {num_cores} cores")


paramsClass1 = {
    "objective": "binary:logistic",
    "eval_metric": "auc",
    ##
#    "objective": "multi:softprob",
#    "num_class": 3,   # signal + 2 backgrounds
#    "eval_metric": "mlogloss",
    ##
    # Model complexity
    "tree_method": "hist",
    "seed": 42,
    "n_estimators": 500,
    # Model complexity (reduce to fight overfitting)
    "max_depth": 4,
    "min_child_weight": 5,
    "eta": 0.1,
    # Sampling (reduce to fight overfitting)
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    # Regularization (increase to fight overfitting)
    "lambda": 1.0,            # L2 regularization
    "gamma": 1.0,             # Minimum loss reduction
    "alpha": 0.1,             # L1 regularization
    #
    "early_stopping_rounds": 50,
    "n_jobs": num_cores
}

paramsClass2 = {
    "objective": "binary:logistic",
    "eval_metric": "auc",
    ##
    # Model complexity
    "tree_method": "hist",
    "seed": 42,
    "n_estimators": 500,
    "max_bin": 512,
    # Model complexity (reduce to fight overfitting)
    "max_depth": 3,
    "min_child_weight": 15,
    "eta": 0.02,
    # Sampling (reduce to fight overfitting)
    "subsample": 0.9,
    "colsample_bytree": 0.9,
    # Regularization (increase to fight overfitting)
    "lambda": 5.0,            # L2 regularization
    "gamma": 2.0,             # Minimum loss reduction
    "alpha": 0.5,             # L1 regularization
    #
    "early_stopping_rounds": 50,
    "n_jobs": num_cores
}

paramsMulti = {
    # ====== MULTI-CLASS OBJECTIVE ======
    "objective": "multi:softprob",      # ← Multi-class (probabilities)
    "num_class": 4,                     # ← 4 classes!
    "eval_metric": "mlogloss",          # ← Multi-class loss (not AUC)
    # Anti-overfitting
    "max_depth": 4,
    "min_child_weight": 5,
    "eta": 0.01,
    "lambda": 5.0,
    "alpha": 0.5,
    "gamma": 2.0,
    # Sampling
    "subsample": 0.9,
    "colsample_bytree": 0.9,
    "colsample_bylevel": 0.9,
    # Tree method
    "tree_method": "hist",
    "grow_policy": "lossguide",
    "max_leaves": 15,
    "max_bin": 512,
    # Threading
    "n_jobs": num_cores,
    # Training
    "seed": 42,
    "n_estimators": 2000,
    "early_stopping_rounds": 100,
}
