"""
sigFitDCBG.py -- signal fits with a double-sided Crystal Ball PLUS a wider
Gaussian, to test whether the oscillation in the single-CB pull goes away.

Why a second component
----------------------
`HiggsCandCorrMass` mixes muon eta topologies with genuinely different mass
resolutions, so the true lineshape is a SUM of narrow and broad components
rather than one. A single DSCB then splits the difference: too wide at the
peak, too narrow on the shoulders, which shows up as the coherent wave seen
through 121-135 GeV in the ratio panel (chi2/ndf = 4.7 for ggHcat 2024).

Model
-----
    f = fcore * DSCB(mean, sigma)  +  (1 - fcore) * Gauss(mean, rwide*sigma)

Both components share the mean, and the Gaussian's width is sigma times a
RATIO constrained to >= 1. That ordering is deliberate: with two free widths
the fit can swap the roles of the components (an exact degeneracy), and
tying the wide one to a ratio >= 1 removes it.

Because the Gaussian's width is built from the same `sigma` formula var, it
inherits the CMS_res_m nuisance automatically -- the nuisance structure is
unchanged from sigFit.py.

Eight free parameters instead of six, so thin templates will be harder to
constrain than before. Watch N_eff and the "(bound)" flags.

Everything else -- plotting, diagnostics, JSON, CLI -- is sigFit.py, reused
by overriding three module-level names. sigFit.py itself is NOT modified.
The output PDF keeps the name `crystal_ball_<tag>` so the downstream naming
contract with bwsHrare.py still holds.

Usage
-----
    python sigFitDCBG.py --selftest
    python sigFitDCBG.py -c ggHcat -y 2024 -b incl -s ggH

    # then compare against the single-CB fit
    python sigFit.py     -c ggHcat -y 2024 -b incl -s ggH
"""

import os
import sys

import ROOT

import sigFit as sf


# ---------------------------------------------------------------------------
# extra parameters
# ---------------------------------------------------------------------------
# fcore  fraction of the yield in the DSCB core
# rwide  sigma_wide / sigma_core, >= 1 by construction
EXTRA_BOUNDS = {
    "fcore": (0.85, 0.20, 1.00),
    "rwide": (2.50, 1.00, 10.0),
}

EXTRA_LABELS = {
    "fcore": ("f_{core}", ""),
    "rwide": ("#sigma_{wide}/#sigma", ""),
}

# order of the parameter block on the plot
PARAM_ORDER = ("mu", "sigma", "alphaL", "nL", "alphaR", "nR",
               "fcore", "rwide")


# ---------------------------------------------------------------------------
# the PDF
# ---------------------------------------------------------------------------

def make_signal_pdf(x, tag, scale_unc=sf.SCALE_UNC, res_unc=sf.RES_UNC,
                    correlated=True):
    """DSCB + wider Gaussian, sharing the mean and the resolution nuisance.

    Same signature and same (pdf, bundle) contract as sigFit.make_signal_pdf,
    so sigFit.fit_one can call this unchanged.
    """
    # --- shape parameters measured from MC ---
    nom = {}
    for k, (start, lo, hi) in sf.BOUNDS.items():
        nom[k] = ROOT.RooRealVar(f"cb_{k}_{tag}", f"cb_{k}", start, lo, hi)
    for k, (start, lo, hi) in EXTRA_BOUNDS.items():
        nom[k] = ROOT.RooRealVar(f"cb_{k}_{tag}", f"cb_{k}", start, lo, hi)

    # --- mass hypothesis (see sigFit) ---
    mh = ROOT.RooRealVar("MH", "Higgs mass hypothesis", sf.MH_REF, 120.0, 130.0)
    mh.setConstant(True)

    # --- nuisance parameters ---
    nsuf = "" if correlated else f"_{tag}"
    nuis = {
        "scale": ROOT.RooRealVar(f"CMS_scale_m{nsuf}", "muon scale", 0.0, -5.0, 5.0),
        "res": ROOT.RooRealVar(f"CMS_res_m{nsuf}", "muon resolution", 0.0, -5.0, 5.0),
    }

    mean = ROOT.RooFormulaVar(
        f"cb_mean_{tag}", "mean",
        f"(@0 + @1 - {sf.MH_REF})*(1 + {scale_unc}*@2)",
        ROOT.RooArgList(nom["mu"], mh, nuis["scale"]),
    )
    sigma = ROOT.RooFormulaVar(
        f"cb_sigmaEff_{tag}", "sigmaEff",
        f"@0*(1 + {res_unc}*@1)",
        ROOT.RooArgList(nom["sigma"], nuis["res"]),
    )
    # the wide component tracks sigma, so it follows CMS_res_m too
    sigma_wide = ROOT.RooFormulaVar(
        f"cb_sigmaWide_{tag}", "sigmaWide", "@0*@1",
        ROOT.RooArgList(sigma, nom["rwide"]),
    )

    core = ROOT.RooCrystalBall(
        f"dscb_core_{tag}", "dscb_core",
        x, mean, sigma,
        nom["alphaL"], nom["nL"], nom["alphaR"], nom["nR"],
    )
    wide = ROOT.RooGaussian(
        f"gauss_wide_{tag}", "gauss_wide", x, mean, sigma_wide,
    )

    # Keep the sigFit name for the SUM: bwsHrare.py looks up
    # `crystal_ball_<tag>` and `<pdfname>_norm`.
    pdf = ROOT.RooAddPdf(
        f"crystal_ball_{tag}", "dscb_plus_gauss",
        ROOT.RooArgList(core, wide),
        ROOT.RooArgList(nom["fcore"]),
    )

    bundle = {
        "nom": nom,
        "nuis": nuis,
        "mh": mh,
        "func": {"mean": mean, "sigma": sigma, "sigma_wide": sigma_wide},
        # components must stay referenced or PyROOT collects them
        "components": {"core": core, "wide": wide},
        "cfg": {"scale_unc": scale_unc, "res_unc": res_unc,
                "correlated": correlated, "mh_ref": sf.MH_REF,
                "model": "DSCB + Gaussian"},
    }
    return pdf, bundle


# ---------------------------------------------------------------------------
# annotation: sigFit._annotate iterates a hardcoded six-parameter tuple
# ---------------------------------------------------------------------------

def _annotate(nom, chi2_ndf, norm, n_eff, label):
    """sigFit._annotate, extended to the two extra parameters."""
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(sf.ANN_TEXT_SIZE)

    left = sf.PAD_LEFT_MARGIN + 0.05
    y0, dy = 0.85, sf.LINE_SPACING * sf.ANN_TEXT_SIZE
    for i, line in enumerate(sf.header_lines(label)):
        latex.DrawLatex(left, y0 - i * dy, line)

    latex.SetTextSize(sf.PARAM_TEXT_SIZE)
    latex.SetTextAlign(12)
    right = 1.0 - sf.PAD_RIGHT_MARGIN - 0.36
    y = 0.85
    lines = [f"yield = {norm:.1f}", f"N_{{eff}} = {n_eff:.0f}"]
    for k in PARAM_ORDER:
        if k not in nom:
            continue
        v = nom[k]
        lab, unit = LABELS[k]
        txt = sf.fmt_with_unc(lab, v.getVal(), v.getError(), unit)
        if sf.near_bound(v):
            txt += " #color[2]{(bound)}"
        lines.append(txt)
    lines.append(f"#chi^{{2}}/ndf = {chi2_ndf:.2f}")
    pdy = sf.LINE_SPACING * sf.PARAM_TEXT_SIZE
    for i, line in enumerate(lines):
        latex.DrawLatex(right, y - i * pdy, line)

    return latex


LABELS = dict(sf.PARAM_LABEL)
LABELS.update(EXTRA_LABELS)


# ---------------------------------------------------------------------------
# install the overrides
# ---------------------------------------------------------------------------
# sigFit.fit_one and sigFit._annotate resolve these as module globals at call
# time, so reassigning the module attributes is enough -- no edit to
# sigFit.py. This must run before any fit, which it does: import time.
sf.make_signal_pdf = make_signal_pdf
sf._annotate = _annotate
sf.PARAM_LABEL = LABELS


def _default_arg(flag, value):
    """Supply a default for one of sigFit.main's flags if absent, so the
    two models do not write over each other's workspaces and plots."""
    if flag not in sys.argv:
        sys.argv.extend([flag, value])


if __name__ == "__main__":
    _default_arg("--wsdir", "WS_DSCBG")
    _default_arg("--plotdir",
                 os.path.expanduser("~/public_html/HmumuFits/signal_fits_dscbg"))
    print("model: DSCB + wider Gaussian "
          "(8 free parameters; sigFit.py is 6)")
    sf.setup_style()
    sf.main()
