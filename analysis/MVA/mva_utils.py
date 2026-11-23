import ROOT
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import xgboost as xgb
import scipy.stats
from sklearn.metrics import roc_auc_score, roc_curve, auc, accuracy_score
from sklearn.metrics import confusion_matrix, ConfusionMatrixDisplay

try:
    from scipy.ndimage import gaussian_filter1d
except ImportError:
    gaussian_filter1d = None

def get_var_range(df, var, pad=0.05):
    min_val = df.Min(var).GetValue()
    max_val = df.Max(var).GetValue()
    span = max_val - min_val
    return min_val - span * pad, max_val + span * pad

def make_input_plots(df, sig_ids, bkg_ids, variables, outdir, suffix, weight_col="w_allSF"):
    df = df.Define("weight", weight_col)
    sig_filter = " || ".join([f"mc == {mid}" for mid in sig_ids])
    bkg_filter = " || ".join([f"mc == {mid}" for mid in bkg_ids])
    df_sig = df.Filter(sig_filter).Define("absweight", "fabs(weight)")
    df_bkg = df.Filter(bkg_filter).Define("absweight", "fabs(weight)")

    for var in variables:
        print(f"Plotting {var}...")
        xmin_s, xmax_s = get_var_range(df_sig, var)
        xmin_b, xmax_b = get_var_range(df_bkg, var)
        xmin, xmax = min(xmin_s, xmin_b), max(xmax_s, xmax_b)
        nbins = 40

        hsig = df_sig.Histo1D((f"h_{var}_sig", f"{var};{var};Norm. Events", nbins, xmin, xmax), var, "weight")
        hbkg = df_bkg.Histo1D((f"h_{var}_bkg", f"{var};{var};Norm. Events", nbins, xmin, xmax), var, "weight")
        hsigABS = df_sig.Histo1D((f"h_{var}_abssig", f"{var};{var};Norm. Events", nbins, xmin, xmax), var, "absweight")
        hbkgABS = df_bkg.Histo1D((f"h_{var}_absbkg", f"{var};{var};Norm. Events", nbins, xmin, xmax), var, "absweight")

        c = ROOT.TCanvas(f"c_{var}", "", 800, 600)
        hsig.SetLineColor(ROOT.kOrange)
        hsigABS.SetLineColor(ROOT.kRed)
        hbkg.SetLineColor(ROOT.kBlue)
        hbkgABS.SetLineColor(ROOT.kGreen+1)

        for h in [hsig, hsigABS, hbkg, hbkgABS]:
            if h.Integral() > 0:
                h.Scale(1.0 / h.Integral())

        hsig.SetMaximum(1.5 * max(hsig.GetMaximum(), hbkg.GetMaximum()))
        hsig.Draw("HIST")
        hsigABS.Draw("HIST SAME")
        hbkg.Draw("HIST SAME")
        hbkgABS.Draw("HIST SAME")

        leg = ROOT.TLegend(0.60, 0.40, 0.88, 0.65)
        leg.SetBorderSize(0)
        leg.SetFillStyle(0)
        leg.SetTextSize(0.04)
        leg.AddEntry(hsigABS.GetValue(), "sig", "l")
        leg.AddEntry(hbkg.GetValue(), "bkg", "l")
        leg.Draw()

        c.SaveAs(f"{outdir}/VARS/{var}{suffix}_sig_bkg.png")
        print(f"  Saved {var}_sig_bkg.png")


def plot_roc(y_true, y_score, outfile):
    fpr, tpr, _ = roc_curve(y_true, y_score)
    roc_auc = auc(fpr, tpr)

    plt.figure(figsize=(8, 6))
    plt.plot(fpr, tpr, color='blue', lw=2, label=f'ROC curve (AUC = {roc_auc:.2f})')
    plt.plot([0, 1], [0, 1], color='gray', lw=1, linestyle='--')
    plt.xlim([0.0, 1.0])
    plt.ylim([0.0, 1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.title('XGBoost ROC Curve')
    plt.legend(loc="lower right")
    plt.grid(alpha=0.3)
    plt.savefig(outfile, dpi=200)
    plt.close()
    print(f"Saved ROC to: {outfile}")
    return roc_auc


def plot_feature_importance(bdt, variables, outdir, suffix=""):
    booster = bdt.get_booster()
    booster.feature_names = variables
    importance = booster.get_score(importance_type="gain")
    print("Non-zero features:", importance)

    titles = {
        "gain": "XGBoost Feature Importance - Gain",
        "weight": "XGBoost Feature Importance - Weight (Split Count)",
        "cover": "XGBoost Feature Importance - Cover",
    }

    for imp, title in titles.items():
        plt.figure(figsize=(8, 6))
        xgb.plot_importance(bdt, importance_type=imp, max_num_features=20)
        plt.title(title)
        plt.subplots_adjust(left=0.35)
        outfile = f"{outdir}feature_importance_{imp}{suffix}.png"
        plt.savefig(outfile, dpi=200)
        plt.close()
        print(f"Saved feature importance ({imp}) to: {outfile}")


def plot_correlation_matrix(corr_csv, outfile, split_from_end=2):
    df = pd.read_csv(corr_csv, sep="\t", index_col=0)

    plt.figure(figsize=(8, 7))
    ax = plt.gca()
    im = plt.imshow(df, vmin=-1, vmax=1)
    plt.xticks(range(len(df.columns)), df.columns, rotation=45, ha="right")
    plt.yticks(range(len(df.index)), df.index)
    plt.colorbar(im, label="Correlation")

    split = len(df) - split_from_end
    ax.axhline(split - 0.5, color="white", linewidth=2)
    ax.axvline(split - 0.5, color="white", linewidth=2)

    for i in range(len(df.index)):
        for j in range(len(df.columns)):
            plt.text(j, i, f"{df.iloc[i, j]:.2f}", ha="center", va="center", fontsize=8)

    plt.title("Feature Correlation Matrix")
    plt.tight_layout()
    plt.savefig(outfile, dpi=200)
    plt.close()
    print(f"Saved {outfile}")


def plot_overtraining(bdt, train_data, train_labels, train_weights,
                      test_data, test_labels, test_weights,
                      signal_classes, outfile, doMultiClass=False):
    if doMultiClass:
        train_scores = bdt.predict_proba(train_data)[:, 0]
        test_scores  = bdt.predict_proba(test_data)[:, 0]
        sig_mask_train = (train_labels == 0)
        sig_mask_test  = (test_labels  == 0)
    else:
        train_scores = bdt.predict_proba(train_data)[:, 1]
        test_scores  = bdt.predict_proba(test_data)[:, 1]
        sig_mask_train = (train_labels == 1)
        sig_mask_test  = (test_labels  == 1)

    bkg_mask_train = ~sig_mask_train
    bkg_mask_test  = ~sig_mask_test

    bins = np.linspace(0, 1, 50)

    def weighted_hist_with_err(values, weights, bins):
        hist, edges = np.histogram(values, bins=bins, weights=weights, density=True)
        sumw2, _ = np.histogram(values, bins=bins, weights=weights**2)
        err = np.sqrt(sumw2) / np.diff(edges) / np.sum(weights)
        centers = (edges[:-1] + edges[1:]) / 2
        return centers, hist, err

    plt.figure(figsize=(8, 6))

    plt.hist(train_scores[sig_mask_train], bins=bins, weights=train_weights[sig_mask_train],
             density=True, histtype='step', linewidth=2, label='Train Signal')
    plt.hist(train_scores[bkg_mask_train], bins=bins, weights=train_weights[bkg_mask_train],
             density=True, histtype='step', linewidth=2, label='Train Background')

    x_sig, y_sig, err_sig = weighted_hist_with_err(test_scores[sig_mask_test], test_weights[sig_mask_test], bins)
    plt.errorbar(x_sig, y_sig, yerr=err_sig, fmt='o', label='Test Signal')

    x_bkg, y_bkg, err_bkg = weighted_hist_with_err(test_scores[bkg_mask_test], test_weights[bkg_mask_test], bins)
    plt.errorbar(x_bkg, y_bkg, yerr=err_bkg, fmt='o', label='Test Background')

    plt.xlabel('BDT Output Score')
    plt.ylabel('Probability Density')
    plt.title('XGBoost Overtraining Check')
    plt.legend()
    plt.grid(alpha=0.3)
    plt.savefig(outfile, dpi=200)
    plt.close()
    print(f"Saved overtraining to: {outfile}")


def plot_auc_vs_round(bdt, outfile):
    results = bdt.evals_result()
    plt.figure(figsize=(8, 6))
    plt.plot(results['validation_0']['auc'], label='Train')
    plt.plot(results['validation_1']['auc'], label='Validation')
    plt.xlabel('Boosting Rounds')
    plt.ylabel('AUC')
    plt.legend()
    plt.title('XGBoost Overfitting Plot')
    plt.savefig(outfile, dpi=200)
    plt.close()
    print(f"Saved AUC vs round to: {outfile}")


def compare_even_odd(score_even, score_odd, y_even, y_odd, w_even, w_odd,
                     signal_classes, background_classes, outfile):
    sig_mask_even = np.isin(y_even, signal_classes)
    sig_mask_odd  = np.isin(y_odd, signal_classes)
    bkg_mask_even = np.isin(y_even, background_classes)
    bkg_mask_odd  = np.isin(y_odd, background_classes)

    bins = np.linspace(0, 1, 50)
    plt.figure(figsize=(8, 6))

    plt.hist(score_even[sig_mask_even], bins=bins, weights=w_even[sig_mask_even],
             density=True, histtype='step', linewidth=2, label="Even (signal)")
    plt.hist(score_odd[sig_mask_odd], bins=bins, weights=w_odd[sig_mask_odd],
             density=True, histtype='step', linewidth=2, linestyle="--", label="Odd (signal)")
    plt.hist(score_even[bkg_mask_even], bins=bins, weights=w_even[bkg_mask_even],
             density=True, histtype='step', linewidth=2, label="Even (bkg)")
    plt.hist(score_odd[bkg_mask_odd], bins=bins, weights=w_odd[bkg_mask_odd],
             density=True, histtype='step', linewidth=2, linestyle="--", label="Odd (bkg)")

    plt.xlabel("BDT score")
    plt.ylabel("Normalized entries")
    plt.title("Even/Odd training stability")
    plt.legend()
    plt.grid(alpha=0.3)
    plt.savefig(outfile, dpi=200)
    plt.close()
    print(f"Saved {outfile}")


def plot_score_vs_mass(score, sig_mask, bkg_mask, data, mass_col, outdir, category, suffix=""):
    mass = data[mass_col].values
    valid = ~(np.isnan(score) | np.isnan(mass))
    score, mass = score[valid], mass[valid]
    sig_mask, bkg_mask = sig_mask[valid], bkg_mask[valid]

    corr_sig = np.corrcoef(score[sig_mask], mass[sig_mask])[0, 1]
    corr_bkg = np.corrcoef(score[bkg_mask], mass[bkg_mask])[0, 1]
    print(f"corr(score,mass) signal = {corr_sig:.4f}")
    print(f"corr(score,mass) bkg    = {corr_bkg:.4f}")

    mass_range = [mass.min() - 5, mass.max() + 5]

    for name, mask in [("sig", sig_mask), ("bkg", bkg_mask)]:
        plt.figure(figsize=(8, 6))
        plt.hist2d(mass[mask], score[mask], bins=[60, 50], range=[mass_range, [0, 1]])
        plt.xlabel(f"{mass_col} [GeV]")
        plt.ylabel("BDT score")
        plt.title(f"{name}: BDT score vs mass")
        plt.colorbar(label="Events")
        plt.tight_layout()
        outfile = f"{outdir}score_vs_mass_{name}_{category}{suffix}.png"
        plt.savefig(outfile, dpi=200)
        plt.close()
        print(f"Saved {outfile}")


def plot_mass_in_score_bins(score, sig_mask, bkg_mask, data, mass_col, outdir, category, suffix=""):
    mass = data[mass_col].values
    mass_range = (mass.min() - 5, mass.max() + 5)
    score_bins = np.arange(0, 1.01, 0.1)

    for name, mask in [("sig", sig_mask), ("bkg", bkg_mask)]:
        plt.figure(figsize=(8, 6))
        for low, high in zip(score_bins[:-1], score_bins[1:]):
            sel = mask & (score >= low) & (score < high)
            if np.sum(sel) < 10:
                continue
            plt.hist(mass[sel], bins=60, range=mass_range, density=True,
                     histtype="step", linewidth=1.5, label=f"{low:.1f}-{high:.1f}")
        plt.xlabel(f"{mass_col} [GeV]")
        plt.ylabel("Normalized events")
        plt.title(f"{category} {name}: mass in score bins")
        plt.legend(ncol=2, fontsize=8)
        plt.tight_layout()
        outfile = f"{outdir}mass_in_score_bins_{name}_{category}{suffix}.png"
        plt.savefig(outfile, dpi=200)
        plt.close()
        print(f"Saved {outfile}")


def print_var_diagnostics(train_data, train_labels, variables):
    for v in variables:
        col = train_data[v]
        print(f"{v:25s} min={col.min():.4g}  max={col.max():.4g}  std={col.std():.4g}")

    for v in variables:
        sig = train_data.loc[train_labels == 1, v]
        bkg = train_data.loc[train_labels == 0, v]
        ks_stat = scipy.stats.ks_2samp(sig, bkg).statistic
        if   ks_stat < 0.05: comment = 'No separation'
        elif ks_stat < 0.15: comment = 'very weak'
        elif ks_stat < 0.30: comment = 'some discrimination'
        else:                comment = 'STRONG'
        print(f"{v:25s} <sig>={sig.mean():.4g}  <bkg>={bkg.mean():.4g}  KS={ks_stat:.4f}  => {comment}")


def build_training_weights(data,
                           doMultiClass=False,
                           flatten_background=True,
                           signal_sigma_power=1.0,
                           signal_relative_sigma=True,
                           mass_col="HiggsCandCorrMass",
                           sigma_col="HiggsCandCorrMassErr"):
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

        for k in range(nclass):
            mask = data["target"] == k
            weights[mask] *= class_weights[k]

        weights /= weights.mean()

        print("Multiclass weights per class:")
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
    
    sigma = data[sigma_col].to_numpy()
    sigma_safe = np.where((sigma <= 0) | (~np.isfinite(sigma)), 1e-3, sigma)

    if signal_relative_sigma:
        mass = data[mass_col].to_numpy()
        mass_safe = np.where((mass <= 0) | (~np.isfinite(mass)), 125., mass)
        sigma_used = sigma_safe / mass_safe
    else:
        sigma_used = sigma_safe

    weights[sig_mask] *= 1.0 / np.power(sigma_used[sig_mask], signal_sigma_power)

    if flatten_background and gaussian_filter1d is not None:
        mass = data[mass_col].to_numpy()
        bins = np.linspace(mass[bkg_mask].min() - 1, mass[bkg_mask].max() + 1, 41)
        hist, edges = np.histogram(mass[bkg_mask], bins=bins, weights=w[bkg_mask])
        hist = gaussian_filter1d(hist.astype(float), sigma=1.5)
        hist = np.maximum(hist, 1e-6)
        inv = 1. / hist
        inv *= hist.sum() / np.sum(hist * inv)
        ibin = np.clip(np.digitize(mass, edges) - 1, 0, len(inv) - 1)
        weights[bkg_mask] *= inv[ibin][bkg_mask]

    weights[sig_mask] /= weights[sig_mask].mean()
    weights[bkg_mask] /= weights[bkg_mask].mean()

    return weights


def save_xgb_model(bdt, variables, fOutName, model_name):
    ROOT.TMVA.Experimental.SaveXGBoost(bdt, model_name, fOutName, num_inputs=len(variables))
    print(f"Model exported to {fOutName} as {model_name}")

    variables_ = ROOT.TList()
    for var in variables:
        variables_.Add(ROOT.TObjString(var))
    fOut = ROOT.TFile(fOutName, "UPDATE")
    fOut.WriteObject(variables_, "variables")
    fOut.Close()
    print("Variable list saved")
