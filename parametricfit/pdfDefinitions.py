"""
pdfDefinitions.py -- the PDF definitions shared by every fit script.

Design rules, all of which the previous version broke somewhere:

  1. NO ROOT CLASSES IMPORTED AT MODULE LEVEL. `import ROOT` only, and
     classes resolved inside the functions. The old
     `from ROOT import RooDoubleCBFast` ran at import time, so the class
     had to exist before this module was imported -- in SIGfits.py the
     import is five lines BEFORE gSystem.Load, and it worked only because
     ROOT autoloads from the library's rootmap. Now load order cannot
     matter.

  2. NO COMBINE DEPENDENCY. Every background PDF here is plain ROOT. The
     signal shape picks RooCrystalBall (ROOT >= 6.26) when it exists and
     falls back to Combine's RooDoubleCBFast otherwise -- identical
     function, identical arguments. So this module imports in the plain
     conda environment; only workspace writing (RooMultiPdf) needs the
     container.

  3. EVERY CANDIDATE OWNS ITS PARAMETERS. The old create_bkg_pdfs gave
     bern2 and bern3 the same bern_c0..c2 objects, and all three
     exponentials the same exp_p1, so fitting one moved another's
     starting point and the result depended on fit order. Each candidate
     is now built by a factory with its own variables.

  4. ONE SET OF BOUNDS, ONE SET OF NAMES. Defined here and imported, not
     duplicated in each script.

The signal factory returns (pdf, bundle) with the bundle layout
sigFit.fit_one expects: nom / nuis / mh / func / cfg.
"""

import ROOT


# ---------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------

MH_REF = 125.0            # reference mass the shape is fitted at

# TODO: fractional uncertainties folded into the shape formulas.
# Replace with Run 3 muon POG numbers, or from objScaleSmear() variations.
SCALE_UNC = 0.002         # 0.2% muon momentum scale
RES_UNC = 0.05            # 5% muon momentum resolution

# Signal shape bounds. alpha/n ranges follow the ATLAS H->mumu choice; n
# must be allowed well above 10 or the fits park on the upper bound.
BOUNDS = {
    "mu":     (MH_REF, 120.0, 130.0),
    "sigma":  (2.0,      0.5,   6.0),
    "alphaL": (1.5,      0.1,   3.0),
    "nL":     (5.0,      0.1,  50.0),
    "alphaR": (1.5,      0.1,   3.0),
    "nR":     (5.0,      0.1,  50.0),
}

# With RooDoubleCBFast the tail integral's 1/(n-1) factor is NOT special-
# cased at n = 1, so a fit wandering there returns NaN with no warning.
# RooCrystalBall handles it. Raise the floor when falling back.
N_FLOOR_NO_SPECIAL_CASE = 1.05

PARAM_LABEL = {
    "mu":     ("#mu",           "GeV"),
    "sigma":  ("#sigma",        "GeV"),
    "alphaL": ("#alpha_{lo}",   ""),
    "nL":     ("n_{lo}",        ""),
    "alphaR": ("#alpha_{high}", ""),
    "nR":     ("n_{high}",      ""),
}

COMBINE_LIB = "/code/HiggsAnalysis/CombinedLimit/build/lib/" \
              "libHiggsAnalysisCombinedLimit.so"

_combine_loaded = False


# ---------------------------------------------------------------------------
# class resolution
# ---------------------------------------------------------------------------

def load_combine(path=None, quiet=False):
    """Load the Combine library, once. Only needed for RooMultiPdf (and
    RooDoubleCBFast, if this ROOT has no RooCrystalBall)."""
    global _combine_loaded
    if _combine_loaded:
        return True
    ok = ROOT.gSystem.Load(path or COMBINE_LIB) >= 0
    _combine_loaded = ok
    if not ok and not quiet:
        print(f"warning: could not load {path or COMBINE_LIB}")
    return ok


def dscb_class():
    """The double-sided Crystal Ball class available in this environment.

    Returns (class, name). Same function and same 7 arguments either way:
        (x, mean, sigma, alphaL, nL, alphaR, nR)
    """
    cls = getattr(ROOT, "RooCrystalBall", None)
    if cls is not None:
        return cls, "RooCrystalBall"
    load_combine(quiet=True)
    cls = getattr(ROOT, "RooDoubleCBFast", None)
    if cls is None:
        raise RuntimeError(
            "no double-sided Crystal Ball available: this ROOT has no "
            "RooCrystalBall (needs >= 6.26) and Combine's RooDoubleCBFast "
            f"could not be loaded from {COMBINE_LIB}")
    return cls, "RooDoubleCBFast"


def signal_bounds(cls_name=None):
    """BOUNDS, with the n floor raised if the class has no n = 1 case."""
    if cls_name is None:
        cls_name = dscb_class()[1]
    b = {k: tuple(v) for k, v in BOUNDS.items()}
    if cls_name != "RooCrystalBall":
        for k in ("nL", "nR"):
            start, lo, hi = b[k]
            b[k] = (max(start, N_FLOOR_NO_SPECIAL_CASE),
                    max(lo, N_FLOOR_NO_SPECIAL_CASE), hi)
    return b


# ---------------------------------------------------------------------------
# signal
# ---------------------------------------------------------------------------

def create_signal_pdf(x, tag, scale_unc=None, res_unc=None,
                      correlated=True, bounds=None, name=None):
    """Double-sided Crystal Ball with mass-hypothesis and nuisance hooks.

        mean  = (cb_mu + MH - MH_REF) * (1 + scale_unc * CMS_scale_m)
        sigma =  cb_sigma             * (1 + res_unc   * CMS_res_m)

    MH is the mass hypothesis: an ADDITIVE shift, so the same workspace
    can be evaluated at other masses. Named MH because that is what
    Combine's -m looks for. Not a systematic.

    CMS_scale_m / CMS_res_m are unit Gaussians. They are exactly
    degenerate with cb_mu / cb_sigma, so the caller must pin them for the
    MC fit and release them afterwards.

    Returns (pdf, bundle). Everything lives in the bundle because PyROOT
    collects any RooFit object nothing holds a reference to.
    """
    scale_unc = SCALE_UNC if scale_unc is None else scale_unc
    res_unc = RES_UNC if res_unc is None else res_unc
    cls, cls_name = dscb_class()
    bounds = signal_bounds(cls_name) if bounds is None else bounds

    nom = {k: ROOT.RooRealVar(f"cb_{k}_{tag}", f"cb_{k}", *v)
           for k, v in bounds.items()}

    mh = ROOT.RooRealVar("MH", "Higgs mass hypothesis", MH_REF, 120.0, 130.0)
    mh.setConstant(True)

    nsuf = "" if correlated else f"_{tag}"
    nuis = {
        "scale": ROOT.RooRealVar(f"CMS_scale_m{nsuf}", "muon scale",
                                 0.0, -5.0, 5.0),
        "res": ROOT.RooRealVar(f"CMS_res_m{nsuf}", "muon resolution",
                               0.0, -5.0, 5.0),
    }

    mean = ROOT.RooFormulaVar(
        f"cb_mean_{tag}", "mean",
        f"(@0 + @1 - {MH_REF})*(1 + {scale_unc}*@2)",
        ROOT.RooArgList(nom["mu"], mh, nuis["scale"]))
    sigma = ROOT.RooFormulaVar(
        f"cb_sigmaEff_{tag}", "sigmaEff",
        f"@0*(1 + {res_unc}*@1)",
        ROOT.RooArgList(nom["sigma"], nuis["res"]))

    # bwsHrare.py looks the PDF up as crystal_ball_<tag>; keep it.
    pdf = cls(name or f"crystal_ball_{tag}", "crystal_ball",
              x, mean, sigma,
              nom["alphaL"], nom["nL"], nom["alphaR"], nom["nR"])

    return pdf, {"nom": nom, "nuis": nuis, "mh": mh,
                 "func": {"mean": mean, "sigma": sigma},
                 "cfg": {"scale_unc": scale_unc, "res_unc": res_unc,
                         "correlated": correlated, "mh_ref": MH_REF,
                         "pdf_class": cls_name,
                         "model": "double-sided Crystal Ball"}}


# ---------------------------------------------------------------------------
# background
# ---------------------------------------------------------------------------
# Each factory builds ONE candidate with its own RooRealVars, so no two
# candidates share a parameter and fit order cannot matter. Coefficient
# ranges are deliberately wide: a coefficient pinned against a bound makes
# a good model look bad, and the discrete-profiling envelope needs each
# candidate fitted at its own best point, not at the edge of a box.

def _bernstein(x, tag, order):
    """Bernstein of the given order: order+1 coefficients, all in [0, 10].

    Bernstein coefficients are non-negative by construction (that is what
    keeps the PDF positive), so the lower bound of 0 is physical; the
    upper bound is loose.
    """
    pars = [ROOT.RooRealVar(f"bern{order}_c{i}_{tag}", f"c{i}",
                            1.0 if i == 0 else 0.5, 0.0, 10.0)
            for i in range(order + 1)]
    pdf = ROOT.RooBernstein(f"bern{order}_{tag}", f"bern{order}",
                            x, ROOT.RooArgList(*pars))
    return pdf, pars


def _chebychev(x, tag, order):
    """Chebychev of the given order: `order` coefficients in [-1, 1]."""
    pars = [ROOT.RooRealVar(f"cheb{order}_c{i}_{tag}", f"c{i}",
                            0.1 if i == 0 else 0.0, -1.0, 1.0)
            for i in range(order)]
    pdf = ROOT.RooChebychev(f"cheb{order}_{tag}", f"cheb{order}",
                            x, ROOT.RooArgList(*pars))
    return pdf, pars


def _sum_of(x, tag, kind, nterms):
    """Sum of `nterms` exponentials or power laws, with recursive
    fractions so the coefficients stay in [0, 1] and the sum is positive.

    RooAddPdf with recursiveFraction=True means coefficient i is the
    fraction of what is left after the first i-1 terms, which removes the
    "fractions must sum to < 1" constraint that trips up the explicit
    form. The old code emulated this with RooFormulaVars built from
    exp_c1/exp_c2, which exp2 and exp3 SHARED -- so fitting one moved the
    other's fractions, on top of the shared slope parameters.
    """
    comps, pars = [], []
    for i in range(nterms):
        if kind == "exp":
            p = ROOT.RooRealVar(f"{kind}{nterms}_p{i}_{tag}", f"p{i}",
                                -0.02 * (i + 1), -2.0, 0.0)
            c = ROOT.RooExponential(f"{kind}{nterms}_c{i}_{tag}",
                                    f"comp{i}", x, p)
        else:
            p = ROOT.RooRealVar(f"{kind}{nterms}_p{i}_{tag}", f"p{i}",
                                -2.5 * (i + 1), -20.0, 0.0)
            c = ROOT.RooGenericPdf(f"{kind}{nterms}_c{i}_{tag}", f"comp{i}",
                                   "TMath::Power(@0,@1)",
                                   ROOT.RooArgList(x, p))
        comps.append(c)
        pars.append(p)

    if nterms == 1:
        return comps[0], pars, comps

    fracs = [ROOT.RooRealVar(f"{kind}{nterms}_f{i}_{tag}", f"f{i}",
                             0.5, 0.0, 1.0) for i in range(nterms - 1)]
    pdf = ROOT.RooAddPdf(f"{kind}{nterms}_{tag}", f"{kind}{nterms}",
                         ROOT.RooArgList(*comps), ROOT.RooArgList(*fracs),
                         True)                      # recursive fractions
    return pdf, pars + fracs, comps


# name -> builder. Add a candidate here and every script sees it.
BKG_FACTORIES = {
    **{f"bern{o}": (lambda x, t, o=o: _bernstein(x, t, o)) for o in range(6)},
    **{f"cheb{o}": (lambda x, t, o=o: _chebychev(x, t, o)) for o in (1, 2, 3)},
    # no [:2] slice: the third return value holds the components, and
    # dropping it lets PyROOT delete them under the RooAddPdf (segfault)
    **{f"exp{n}": (lambda x, t, n=n: _sum_of(x, t, "exp", n))
       for n in (1, 2, 3)},
    **{f"pow{n}": (lambda x, t, n=n: _sum_of(x, t, "pow", n))
       for n in (1, 2, 3)},
}

# number of FREE parameters, for the F-test and for ndf
BKG_NPAR = {
    **{f"bern{o}": o + 1 for o in range(6)},
    **{f"cheb{o}": o for o in (1, 2, 3)},
    **{f"exp{n}": 2 * n - 1 for n in (1, 2, 3)},
    **{f"pow{n}": 2 * n - 1 for n in (1, 2, 3)},
}

# candidates that are nested within one another, for the F-test
BKG_FAMILIES = {
    "bernstein": [f"bern{o}" for o in range(6)],
    "chebychev": [f"cheb{o}" for o in (1, 2, 3)],
    "exponential": [f"exp{n}" for n in (1, 2, 3)],
    "power": [f"pow{n}" for n in (1, 2, 3)],
}


def create_bkg_pdf(x, tag, name):
    """One background candidate by name, with its own parameters.

    Returns (pdf, keep) -- `keep` holds every object that must stay
    referenced, which PyROOT otherwise collects.
    """
    if name not in BKG_FACTORIES:
        raise KeyError(f"unknown background model '{name}'; "
                       f"known: {sorted(BKG_FACTORIES)}")
    pdf, pars, *extra = BKG_FACTORIES[name](x, tag)   # extra = [comps]
    return pdf, {"params": pars, "pdf": pdf, "extra": extra}


def create_bkg_pdfs(x, tag, names=None):
    """Several candidates at once. Each is independent of the others."""
    names = sorted(BKG_FACTORIES) if names is None else names
    pdfs, keep = {}, {}
    for n in names:
        pdfs[n], keep[n] = create_bkg_pdf(x, tag, n)
    return {"pdfs": pdfs, "keep": keep}


# ---------------------------------------------------------------------------
# backwards compatibility
# ---------------------------------------------------------------------------

class PDFDefinitions:
    """The old class interface, kept so SIGfits.py still imports.

    create_signal_pdfs returned {'pdf':..., 'params':...} with no nuisance
    hooks; new code should call create_signal_pdf and use the bundle.
    """

    @staticmethod
    def create_signal_pdfs(x, tag, suffix=""):
        pdf, bundle = create_signal_pdf(
            x, tag, name=f"crystal_ball{suffix}{tag}")
        return {"pdf": pdf, "params": bundle["nom"], "bundle": bundle}

    @staticmethod
    def create_bkg_pdfs(x, tag):
        out = create_bkg_pdfs(x, tag)
        # the old code exposed 'chebychev1' etc.
        pdfs = dict(out["pdfs"])
        for o in (1, 2, 3):
            if f"cheb{o}" in pdfs:
                pdfs[f"chebychev{o}"] = pdfs[f"cheb{o}"]
        return {"pdfs": pdfs, "params": {}, "components": {},
                "formulas": {}, "keep": out["keep"]}
