# pdfDefinitions.py
import ROOT
from ROOT import RooRealVar, RooBernstein, RooChebychev
from ROOT import RooExponential, RooAddPdf, RooGenericPdf, RooArgList, RooFormulaVar

class PDFDefinitions:
    """Centralized PDF definitions for signal and background fits"""

    @staticmethod
    def create_signal_pdfs(x, tag, suffix="", scale_unc=None, res_unc=None,
                           correlated=True, bounds=None):
        """Create double Crystal Ball signal PDFs

            mean  = (cb_mu + MH - MH_REF) * (1 + scale_unc * CMS_scale_m)
            sigma =  cb_sigma             * (1 + res_unc   * CMS_res_m)

        MH is the mass hypothesis (an additive shift; Combine's -m), not a
        systematic. CMS_scale_m / CMS_res_m are unit Gaussians, exactly
        degenerate with cb_mu / cb_sigma, so the caller pins them for the MC
        fit and releases them afterwards.
        """
        scale_unc = SCALE_UNC if scale_unc is None else scale_unc
        res_unc = RES_UNC if res_unc is None else res_unc
        bounds = BOUNDS if bounds is None else bounds
        RooDSCB, cls_name = dscb_class()

        cb_mu = RooRealVar(f'cb_mu{suffix}{tag}', 'cb_mu', *bounds['mu'])
        cb_sigma = RooRealVar(f'cb_sigma{suffix}{tag}', 'cb_sigma', *bounds['sigma'])
        cb_alphaL = RooRealVar(f'cb_alphaL{suffix}{tag}', 'cb_alphaL', *bounds['alphaL'])
        cb_alphaR = RooRealVar(f'cb_alphaR{suffix}{tag}', 'cb_alphaR', *bounds['alphaR'])
        cb_nL = RooRealVar(f'cb_nL{suffix}{tag}', 'cb_nL', *bounds['nL'])
        cb_nR = RooRealVar(f'cb_nR{suffix}{tag}', 'cb_nR', *bounds['nR'])

        # mass hypothesis and the two shape nuisances
        mh = RooRealVar('MH', 'Higgs mass hypothesis', MH_REF, 120., 130.)
        mh.setConstant(True)
        nsuf = "" if correlated else f'{suffix}{tag}'
        cms_scale_m = RooRealVar(f'CMS_scale_m{nsuf}', 'muon scale', 0., -5., 5.)
        cms_res_m = RooRealVar(f'CMS_res_m{nsuf}', 'muon resolution', 0., -5., 5.)

        cb_mean = RooFormulaVar(f'cb_mean{suffix}{tag}', 'mean',
                                f'(@0 + @1 - {MH_REF})*(1 + {scale_unc}*@2)',
                                RooArgList(cb_mu, mh, cms_scale_m))
        cb_sigmaEff = RooFormulaVar(f'cb_sigmaEff{suffix}{tag}', 'sigmaEff',
                                    f'@0*(1 + {res_unc}*@1)',
                                    RooArgList(cb_sigma, cms_res_m))

        # bwsHrare.py looks the PDF up as crystal_ball_<tag>; keep it.
        pdf_cb = RooDSCB(f'crystal_ball{suffix}{tag}', 'crystal_ball', x,
                         cb_mean, cb_sigmaEff, cb_alphaL, cb_nL, cb_alphaR, cb_nR)

        return {
            'pdf': pdf_cb,
            'params': {
                'mu': cb_mu, 'sigma': cb_sigma, 'alphaL': cb_alphaL,
                'nL': cb_nL, 'alphaR': cb_alphaR, 'nR': cb_nR
            },
            'nuis': {'scale': cms_scale_m, 'res': cms_res_m},
            'mh': mh,
            'func': {'mean': cb_mean, 'sigma': cb_sigmaEff},
            'cfg': {'scale_unc': scale_unc, 'res_unc': res_unc,
                    'correlated': correlated, 'mh_ref': MH_REF,
                    'pdf_class': cls_name,
                    'model': 'double-sided Crystal Ball'},
        }

    @staticmethod
    def create_bkg_pdfs(x, tag):
        """
        Create all background PDFs - KEEP ALL OBJECTS IN MEMORY!

        Every candidate owns its parameters (keys and names start with the
        candidate name), so fitting one cannot move another's starting point.

        Returns dict with:
            - 'pdfs': final PDF objects to use
            - 'params': all RooRealVar parameters
            - 'components': intermediate PDF components
            - 'formulas': all RooFormulaVar objects
        """

        params = {}
        pdfs = {}
        components = {}  # Store intermediate PDFs!
        formulas = {}    # Store formula variables!

        # ===== BERNSTEIN PDFs =====
        # last coefficient fixed at 1 (only ratios matter), the others in [-10, 10]
        for o in range(1, 6):
            for i in range(o):
                params[f'bern{o}_c{i}'] = RooRealVar(f'bern{o}_c{i}_{tag}', f'c{i}', 1., -10., 10.)
            pdfs[f'bern{o}'] = RooBernstein(f'bern{o}_{tag}', f'bern{o}', x,
                                            RooArgList(*[params[f'bern{o}_c{i}'] for i in range(o)],
                                                       ROOT.RooFit.RooConst(1.)))

        # ===== CHEBYCHEV PDFs =====
        # `order` coefficients in [-1, 1]
        for o in (1, 2, 3):
            for i in range(o):
                params[f'cheb{o}_c{i}'] = RooRealVar(f'cheb{o}_c{i}_{tag}', f'c{i}',
                                                     0.1 if i == 0 else 0., -1., 1.)
            pdfs[f'cheb{o}'] = RooChebychev(f'cheb{o}_{tag}', f'cheb{o}', x,
                                            RooArgList(*[params[f'cheb{o}_c{i}'] for i in range(o)]))

        # ===== EXPONENTIAL PDFs =====
        # ===== POWER PDFs =====
        # sum of n terms with fractions f0, (1-f0)*f1, rest: each in [0, 1], so
        # the sum stays positive
        for kind in ('exp', 'pow'):
            for n in (1, 2, 3):
                comps = []
                for i in range(n):
                    if kind == 'exp':
                        params[f'{kind}{n}_p{i}'] = RooRealVar(f'{kind}{n}_p{i}_{tag}', f'p{i}',
                                                               -0.02 * (i + 1), -2., 0.)
                        components[f'{kind}{n}_c{i}'] = RooExponential(f'{kind}{n}_c{i}_{tag}', f'comp{i}', x,
                                                                       params[f'{kind}{n}_p{i}'])
                    else:
                        params[f'{kind}{n}_p{i}'] = RooRealVar(f'{kind}{n}_p{i}_{tag}', f'p{i}',
                                                               -2.5 * (i + 1), -20., 0.)
                        components[f'{kind}{n}_c{i}'] = RooGenericPdf(f'{kind}{n}_c{i}_{tag}', f'comp{i}',
                                                                      'TMath::Power(@0,@1)',
                                                                      RooArgList(x, params[f'{kind}{n}_p{i}']))
                    comps.append(components[f'{kind}{n}_c{i}'])

                if n == 1:
                    pdfs[f'{kind}1'] = comps[0]
                    continue

                for i in range(n - 1):
                    params[f'{kind}{n}_f{i}'] = RooRealVar(f'{kind}{n}_f{i}_{tag}', f'f{i}', 0.5, 0., 1.)
                formulas[f'{kind}{n}_frac0'] = RooFormulaVar(f'{kind}{n}_frac0_{tag}', '@0',
                                                             RooArgList(params[f'{kind}{n}_f0']))
                if n == 3:
                    formulas[f'{kind}{n}_frac1'] = RooFormulaVar(f'{kind}{n}_frac1_{tag}', '(1-@0)*@1',
                                                                 RooArgList(params[f'{kind}{n}_f0'],
                                                                            params[f'{kind}{n}_f1']))
                pdfs[f'{kind}{n}'] = RooAddPdf(f'{kind}{n}_{tag}', f'{kind}{n}',
                                               RooArgList(*comps),
                                               RooArgList(*[formulas[f'{kind}{n}_frac{i}'] for i in range(n - 1)]))

        # Return EVERYTHING to keep all references alive
        return {
            'pdfs': pdfs,
            'params': params,
            'components': components,
            'formulas': formulas
        }


# ---------------------------------------------------------------------------
# Below: configuration and helpers for the parametricfit scripts. Unlike
# analysis/FITS/pdfDefinitions.py, no Combine class is imported at module
# level: the signal picks RooCrystalBall when this ROOT has it and Combine's
# RooDoubleCBFast otherwise (same function, same arguments), so the module
# imports in plain conda and only the workspace step needs the container.
# ---------------------------------------------------------------------------

MH_REF = 125.0            # reference mass the shape is fitted at

# TODO: fractional uncertainties folded into the shape formulas.
# Replace with Run 3 muon POG numbers, or from objScaleSmear() variations.
SCALE_UNC = 0.002         # 0.2% muon momentum scale
RES_UNC = 0.05            # 5% muon momentum resolution

# Signal shape bounds; n must be allowed well above 10 or the fits park on
# the upper bound.
BOUNDS = {
    "mu":     (MH_REF, 120.0, 130.0),
    "sigma":  (2.0,      0.5,   6.0),
    "alphaL": (1.5,      0.1,   3.0),
    "nL":     (5.0,      0.1,  50.0),
    "alphaR": (1.5,      0.1,   3.0),
    "nR":     (5.0,      0.1,  50.0),
}

# number of FREE parameters, for the F-test and for ndf
BKG_NPAR = {
    **{f"bern{o}": o for o in range(1, 6)},
    **{f"cheb{o}": o for o in (1, 2, 3)},
    **{f"exp{n}": 2 * n - 1 for n in (1, 2, 3)},
    **{f"pow{n}": 2 * n - 1 for n in (1, 2, 3)},
}

# candidates that are nested within one another, for the F-test
BKG_FAMILIES = {
    "bernstein": [f"bern{o}" for o in range(1, 6)],
    "chebychev": [f"cheb{o}" for o in (1, 2, 3)],
    "exponential": [f"exp{n}" for n in (1, 2, 3)],
    "power": [f"pow{n}" for n in (1, 2, 3)],
}

COMBINE_LIB = "/code/HiggsAnalysis/CombinedLimit/build/lib/" \
              "libHiggsAnalysisCombinedLimit.so"

_combine_loaded = False


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
            "RooCrystalBall (needs >= 6.24) and Combine's RooDoubleCBFast "
            f"could not be loaded from {COMBINE_LIB}")
    return cls, "RooDoubleCBFast"


def create_signal_pdf(x, tag, scale_unc=None, res_unc=None,
                      correlated=True, bounds=None):
    """PDFDefinitions.create_signal_pdfs as (pdf, bundle), the layout
    sigFit.fit_one expects: nom / nuis / mh / func / cfg."""
    out = PDFDefinitions.create_signal_pdfs(x, f"_{tag}", scale_unc=scale_unc,
                                            res_unc=res_unc, correlated=correlated,
                                            bounds=bounds)
    return out['pdf'], {"nom": out['params'], "nuis": out['nuis'],
                        "mh": out['mh'], "func": out['func'], "cfg": out['cfg']}


def create_bkg_pdf(x, tag, name):
    """One background candidate by name, with its own parameters.

    Returns (pdf, keep) -- `keep` holds every object that must stay
    referenced, which PyROOT otherwise collects.
    """
    if name not in BKG_NPAR:
        raise KeyError(f"unknown background model '{name}'; "
                       f"known: {sorted(BKG_NPAR)}")
    out = PDFDefinitions.create_bkg_pdfs(x, tag)
    pars = [v for k, v in out['params'].items() if k.startswith(f"{name}_")]
    return out['pdfs'][name], {"params": pars, "pdf": out['pdfs'][name], "all": out}
