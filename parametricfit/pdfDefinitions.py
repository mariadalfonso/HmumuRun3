# pdfDefinitions.py
import ROOT
from ROOT import RooRealVar, RooBernstein, RooChebychev
from ROOT import RooExponential, RooGaussian, RooAddPdf, RooGenericPdf, RooArgList, RooFormulaVar

class PDFDefinitions:
    """Centralized PDF definitions for signal and background fits"""
    
    @staticmethod
    def create_signal_pdfs(x, tag, suffix=""):
        """Create double Crystal Ball signal PDFs"""
        cb_mu = RooRealVar(f'cb_mu{suffix}{tag}', 'cb_mu', 125., 120., 130.)
        cb_sigma = RooRealVar(f'cb_sigma{suffix}{tag}', 'cb_sigma', 2., 0.5, 6.)
        cb_alphaL = RooRealVar(f'cb_alphaL{suffix}{tag}', 'cb_alphaL', 1.5, 0.1, 3.)
        cb_alphaR = RooRealVar(f'cb_alphaR{suffix}{tag}', 'cb_alphaR', 1.5, 0.1, 3.)
        cb_nL = RooRealVar(f'cb_nL{suffix}{tag}', 'cb_nL', 5., 0.1, 50.)
        cb_nR = RooRealVar(f'cb_nR{suffix}{tag}', 'cb_nR', 5., 0.1, 50.)

        pdf_cb = RooDoubleCBFast(f'crystal_ball{suffix}{tag}', 'crystal_ball', x, 
                                 cb_mu, cb_sigma, cb_alphaL, cb_nL, cb_alphaR, cb_nR)

        return {
            'pdf': pdf_cb,
            'params': {
                'mu': cb_mu, 'sigma': cb_sigma, 'alphaL': cb_alphaL,
                'alphaR': cb_alphaR, 'nL': cb_nL, 'nR': cb_nR
            }
        }

    @staticmethod
    def create_bkg_pdfs(x, tag):
        """
        Create all background PDFs - KEEP ALL OBJECTS IN MEMORY!
        
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
        params['bern_c0'] = RooRealVar(f'bern_c0_{tag}', 'bern_c0', 1.)  # constant value for the normalisation; only the ratios to it are fitted
        params['bern_c1'] = RooRealVar(f'bern_c1_{tag}', 'bern_c1', 1., -10., 10.)
        params['bern_c2'] = RooRealVar(f'bern_c2_{tag}', 'bern_c2', 1., -10., 10.)
        params['bern_c3'] = RooRealVar(f'bern_c3_{tag}', 'bern_c3', 1., -10., 10.)
        params['bern_c4'] = RooRealVar(f'bern_c4_{tag}', 'bern_c4', 1., -10., 10.)
        params['bern_c5'] = RooRealVar(f'bern_c5_{tag}', 'bern_c5', 1., -10., 10.)

        pdfs['bern0'] = RooBernstein(f'bern0_{tag}', 'bern0', x, 
                                      RooArgList(params['bern_c0']))
        pdfs['bern1'] = RooBernstein(f'bern1_{tag}', 'bern1', x, 
                                      RooArgList(params['bern_c0'], params['bern_c1']))
        pdfs['bern2'] = RooBernstein(f'bern2_{tag}', 'bern2', x, 
                                      RooArgList(params['bern_c0'], params['bern_c1'], params['bern_c2']))
        pdfs['bern3'] = RooBernstein(f'bern3_{tag}', 'bern3', x, 
                                      RooArgList(params['bern_c0'], params['bern_c1'], 
                                                params['bern_c2'], params['bern_c3']))
        pdfs['bern4'] = RooBernstein(f'bern4_{tag}', 'bern4', x,
                                      RooArgList(params['bern_c0'], params['bern_c1'], 
                                                params['bern_c2'], params['bern_c3'], params['bern_c4']))
        pdfs['bern5'] = RooBernstein(f'bern5_{tag}', 'bern5', x,
                                      RooArgList(params['bern_c0'], params['bern_c1'], 
                                                params['bern_c2'], params['bern_c3'], 
                                                params['bern_c4'], params['bern_c5']))

        # ===== CHEBYCHEV PDFs =====
        params['chebychev_c0'] = RooRealVar(f'chebychev_c0_{tag}', 'chebychev_c0', 1.08, -1.1, 10.)
        params['chebychev_c1'] = RooRealVar(f'chebychev_c1_{tag}', 'chebychev_c1', 0.4, -1., 1.)
        params['chebychev_c2'] = RooRealVar(f'chebychev_c2_{tag}', 'chebychev_c2', 0.01, -0.1, 0.1)
        params['chebychev_c3'] = RooRealVar(f'chebychev_c3_{tag}', 'chebychev_c3', 0., -1., 1.)

        pdfs['chebychev1'] = RooChebychev(f'chebychev1_{tag}', 'chebychev1', x,
                                          RooArgList(params['chebychev_c0'], params['chebychev_c1']))
        pdfs['chebychev2'] = RooChebychev(f'chebychev2_{tag}', 'chebychev2', x,
                                          RooArgList(params['chebychev_c0'], params['chebychev_c1'], 
                                                    params['chebychev_c2']))
        pdfs['chebychev3'] = RooChebychev(f'chebychev3_{tag}', 'chebychev3', x,
                                          RooArgList(params['chebychev_c0'], params['chebychev_c1'], 
                                                    params['chebychev_c2'], params['chebychev_c3']))

        # ===== GAUSSIAN PDF =====
        params['gauss_mu'] = RooRealVar(f'gauss_mu{tag}', 'gauss_mu', 120, 100, 140)
        params['gauss_sigma'] = RooRealVar(f'gauss_sigma{tag}', 'gauss_sigma', 30, 20, 40)
        pdfs['gauss'] = RooGaussian(f'gauss{tag}', 'gauss', x, 
                                     params['gauss_mu'], params['gauss_sigma'])

        # ===== EXPONENTIAL PDFs =====
        params['exp_p1'] = RooRealVar(f'exp_p1{tag}', 'exp_p1', -0.0207, -0.5, 0.)
        params['exp_p2'] = RooRealVar(f'exp_p2{tag}', 'exp_p2', -1e-2, -10, 10)
        params['exp_p3'] = RooRealVar(f'exp_p3{tag}', 'exp_p3', -1e-3, -10, 0)
        params['exp_c1'] = RooRealVar(f'exp_c1{tag}', 'exp_c1', 0.5, 0., 1.)
        params['exp_c2'] = RooRealVar(f'exp_c2{tag}', 'exp_c2', 0.5, 0., 1.)

        formulas['exp_frac1'] = RooFormulaVar(f"exp_frac1{tag}", "@0",
                                              RooArgList(params['exp_c1']))
        formulas['exp_frac2'] = RooFormulaVar(f"exp_frac2{tag}", "(1-@0)*@1",
                                              RooArgList(params['exp_c1'], params['exp_c2']))

        # Store intermediate exponential PDFs
        components['pdf_exp1'] = RooExponential(f'exp1_{tag}', 'exp1', x, params['exp_p1'])
        components['pdf_single_exp2'] = RooExponential(f'single_exp2{tag}', 'single_exp2', x, params['exp_p2'])
        components['pdf_single_exp3'] = RooExponential(f'single_exp3{tag}', 'single_exp3', x, params['exp_p3'])

        pdfs['exp1'] = components['pdf_exp1']
        pdfs['exp2'] = RooAddPdf(f'exp2{tag}', 'exp2',
                                 RooArgList(components['pdf_exp1'], components['pdf_single_exp2']),
                                 RooArgList(formulas['exp_frac1']))
        pdfs['exp3'] = RooAddPdf(f'exp3{tag}', 'exp3',
                                 RooArgList(components['pdf_exp1'], components['pdf_single_exp2'], 
                                           components['pdf_single_exp3']),
                                 RooArgList(formulas['exp_frac1'], formulas['exp_frac2']))

        # ===== POWER PDFs =====
        params['frac1'] = RooRealVar(f'frac1{tag}', 'frac1', 0.3, 0., 1.)
        params['frac2'] = RooRealVar(f'frac2{tag}', 'frac2', 0.5, 0., 1.)
        formulas['pow_frac1'] = RooFormulaVar(f"pow_frac1{tag}", "@0", 
                                              RooArgList(params['frac1']))
        formulas['pow_frac2'] = RooFormulaVar(f"pow_frac2{tag}", "(1-@0)*@1", 
                                              RooArgList(params['frac1'], params['frac2']))

        params['pow_p1'] = RooRealVar(f'p1{tag}', 'p1', -2.555, -10., 0.)
        params['pow_p2'] = RooRealVar(f'p2{tag}', 'p2', -8., -10., 0.)
        params['pow_p3'] = RooRealVar(f'p3{tag}', 'p3', -10., -10., 0.)

        # Store intermediate power PDFs
        components['pdf_pow1'] = RooGenericPdf(f'pow1_{tag}', 'pow1',
                                               'TMath::Power(@0,@1)',
                                               RooArgList(x, params['pow_p1']))
        components['pdf_pow2_comp'] = RooGenericPdf(f'pow2comp_{tag}', 'pow2comp',
                                                    'TMath::Power(@0,@1)',
                                                    RooArgList(x, params['pow_p2']))
        components['pdf_pow3_comp'] = RooGenericPdf(f'pow3comp_{tag}', 'pow3comp',
                                                    'TMath::Power(@0,@1)',
                                                    RooArgList(x, params['pow_p3']))

        pdfs['pow1'] = components['pdf_pow1']
        pdfs['pow2'] = RooAddPdf(f'pow2_{tag}', 'pow2',
                                 RooArgList(components['pdf_pow1'], components['pdf_pow2_comp']),
                                 RooArgList(formulas['pow_frac1']))
        pdfs['pow3'] = RooAddPdf(f'pow3_{tag}', 'pow3',
                                 RooArgList(components['pdf_pow1'], components['pdf_pow2_comp'], 
                                           components['pdf_pow3_comp']),
                                 RooArgList(formulas['pow_frac1'], formulas['pow_frac2']))

        # Return EVERYTHING to keep all references alive
        return {
            'pdfs': pdfs,
            'params': params,
            'components': components,
            'formulas': formulas
        }

# ---------------------------------------------------------------------------
# Helpers for the parametricfit scripts (sigFit, bkgFit, makeWorkspaces).
# ---------------------------------------------------------------------------

MH_REF = 125.0            # reference mass the shape is fitted at

# TODO: fractional uncertainties folded into the shape formulas.
# Replace with Run 3 muon POG numbers, or from objScaleSmear() variations.
SCALE_UNC = 0.002         # 0.2% muon momentum scale
RES_UNC = 0.05            # 5% muon momentum resolution

# signal starting values and bounds, as in create_signal_pdfs
BOUNDS = {
    "mu":     (MH_REF, 120.0, 130.0),
    "sigma":  (2.0,      0.5,   6.0),
    "alphaL": (1.5,      0.1,   3.0),
    "nL":     (5.0,      0.1,  50.0),
    "alphaR": (1.5,      0.1,   3.0),
    "nR":     (5.0,      0.1,  50.0),
}

# background candidates and their number of FREE parameters (F-test, ndf)
BKG_NPAR = {
    **{f"bern{o}": o for o in range(1, 6)},        # bern_c0 is constant
    **{f"chebychev{o}": o + 1 for o in (1, 2, 3)},
    **{f"exp{n}": 2 * n - 1 for n in (1, 2, 3)},
    **{f"pow{n}": 2 * n - 1 for n in (1, 2, 3)},
}

# candidates that are nested within one another, for the F-test
BKG_FAMILIES = {
    "bernstein": [f"bern{o}" for o in range(1, 6)],
    "chebychev": [f"chebychev{o}" for o in (1, 2, 3)],
    "exponential": [f"exp{n}" for n in (1, 2, 3)],
    "power": [f"pow{n}" for n in (1, 2, 3)],
}

COMBINE_LIB = "/code/HiggsAnalysis/CombinedLimit/build/lib/" \
              "libHiggsAnalysisCombinedLimit.so"

_combine_loaded = False


def load_combine(path=None, quiet=False):
    """Load the Combine library, once."""
    global _combine_loaded
    if _combine_loaded:
        return True
    ok = ROOT.gSystem.Load(path or COMBINE_LIB) >= 0
    _combine_loaded = ok
    if not ok and not quiet:
        print(f"warning: could not load {path or COMBINE_LIB}")
    return ok


def dscb_class():
    """RooCrystalBall if this ROOT has it, else Combine's RooDoubleCBFast
    (same function, same arguments). Returns (class, name)."""
    cls = getattr(ROOT, "RooCrystalBall", None)
    if cls is not None:
        return cls, "RooCrystalBall"
    load_combine(quiet=True)
    cls = getattr(ROOT, "RooDoubleCBFast", None)
    if cls is None:
        raise RuntimeError("no double-sided Crystal Ball: no RooCrystalBall "
                           f"and RooDoubleCBFast not loadable from {COMBINE_LIB}")
    return cls, "RooDoubleCBFast"


def RooDoubleCBFast(*args):
    """Used by create_signal_pdfs: the class from dscb_class(), so this module
    imports without Combine."""
    return dscb_class()[0](*args)


def create_signal_pdf(x, tag, scale_unc=None, res_unc=None, correlated=True):
    """DSCB from create_signal_pdfs, with the mass hypothesis and the two
    shape nuisances in its mean and width:

        mean  = (cb_mu + MH - MH_REF) * (1 + scale_unc * CMS_scale_m)
        sigma =  cb_sigma             * (1 + res_unc   * CMS_res_m)

    CMS_scale_m / CMS_res_m are unit Gaussians, degenerate with cb_mu /
    cb_sigma: the caller pins them for the MC fit and releases them after.
    Returns (pdf, bundle) with bundle = nom / nuis / mh / func / cfg.
    """
    scale_unc = SCALE_UNC if scale_unc is None else scale_unc
    res_unc = RES_UNC if res_unc is None else res_unc
    cls, cls_name = dscb_class()

    out = PDFDefinitions.create_signal_pdfs(x, f"_{tag}")
    nom = out['params']

    mh = RooRealVar("MH", "Higgs mass hypothesis", MH_REF, 120., 130.)
    mh.setConstant(True)
    nsuf = "" if correlated else f"_{tag}"
    nuis = {"scale": RooRealVar(f"CMS_scale_m{nsuf}", "muon scale", 0., -5., 5.),
            "res": RooRealVar(f"CMS_res_m{nsuf}", "muon resolution", 0., -5., 5.)}
    mean = ROOT.RooFormulaVar(f"cb_mean_{tag}", "mean",
                              f"(@0 + @1 - {MH_REF})*(1 + {scale_unc}*@2)",
                              RooArgList(nom['mu'], mh, nuis['scale']))
    sigma = ROOT.RooFormulaVar(f"cb_sigmaEff_{tag}", "sigmaEff",
                               f"@0*(1 + {res_unc}*@1)",
                               RooArgList(nom['sigma'], nuis['res']))
    # replaces the bare DSCB of create_signal_pdfs; bwsHrare.py looks it up
    # as crystal_ball_<tag>
    pdf = cls(f"crystal_ball_{tag}", "crystal_ball", x, mean, sigma,
              nom['alphaL'], nom['nL'], nom['alphaR'], nom['nR'])

    return pdf, {"nom": nom, "nuis": nuis, "mh": mh,
                 "func": {"mean": mean, "sigma": sigma},
                 "cfg": {"scale_unc": scale_unc, "res_unc": res_unc,
                         "correlated": correlated, "mh_ref": MH_REF,
                         "pdf_class": cls_name,
                         "model": "double-sided Crystal Ball"},
                 "keep": out}


def create_bkg_pdf(x, tag, name):
    """One background candidate by name, with its own (fresh) parameters.

    Returns (pdf, keep); keep["params"] are its FLOATING parameters.
    """
    if name not in BKG_NPAR:
        raise KeyError(f"unknown background model '{name}'; "
                       f"known: {sorted(BKG_NPAR)}")
    out = PDFDefinitions.create_bkg_pdfs(x, tag)
    pdf = out['pdfs'][name]
    free = {p.GetName() for p in pdf.getParameters(ROOT.RooArgSet(x))
            if not p.isConstant()}
    pars = [v for v in out['params'].values() if v.GetName() in free]
    return pdf, {"params": pars, "pdf": pdf, "all": out}
