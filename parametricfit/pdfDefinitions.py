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
        params['chebychev_c2'] = RooRealVar(f'chebychev_c2_{tag}', 'chebychev_c2', 0.01, -1., 1.)
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
        params['exp_p1'] = RooRealVar(f'exp_p1{tag}', 'exp_p1', -0.0207, -10., 10.)
        params['exp_p2'] = RooRealVar(f'exp_p2{tag}', 'exp_p2', -1e-2, -10, 10)
        params['exp_p3'] = RooRealVar(f'exp_p3{tag}', 'exp_p3', -1e-3, -10, 10)
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

        params['pow_p1'] = RooRealVar(f'p1{tag}', 'p1', -2.555, -10., 10.)
        params['pow_p2'] = RooRealVar(f'p2{tag}', 'p2', -8., -20., 10.)
        params['pow_p3'] = RooRealVar(f'p3{tag}', 'p3', -10., -20., 10.)

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

        # ===== DY CORE x CORRECTION PDFs =====
        # core(m) is the fixed leading-order DY shape (dy_core_expr); only
        # the correction floats. All coefficients start at 0, the pure core.
        #   dyexp<n>: core * exp(a1*T1(t) + ... + an*Tn(t))
        #   dypow<n>: core * (m/125)^(a0 + a1*t + ... + a(n-1)*t^(n-1))
        # with t = the window mapped onto [-1, 1]. The core is written into
        # each formula, so members of both families can share a workspace.
        core = dy_core_expr("@0")
        lo, hi = x.getMin(), x.getMax()
        t = f"((2*@0-{lo + hi!r})/{hi - lo!r})"
        cheb = [t, f"(2*{t}*{t}-1)", f"(4*{t}*{t}*{t}-3*{t})",
                f"(8*{t}*{t}*{t}*{t}-8*{t}*{t}+1)"]
        for n in range(1, DY_MAX_ORDER + 1):
            params[f'dyexp_a{n}'] = RooRealVar(f'dyexp_a{n}_{tag}', f'dyexp_a{n}', 0., -10., 10.)
            params[f'dypow_a{n - 1}'] = RooRealVar(f'dypow_a{n - 1}_{tag}', f'dypow_a{n - 1}', 0., -10., 10.)
        for n in range(1, DY_MAX_ORDER + 1):
            coefs = RooArgList(x)
            for i in range(1, n + 1):
                coefs.add(params[f'dyexp_a{i}'])
            arg = " + ".join(f"@{i}*{cheb[i - 1]}" for i in range(1, n + 1))
            pdfs[f'dyexp{n}'] = RooGenericPdf(f'dyexp{n}_{tag}', f'dyexp{n}',
                                              f"{core}*exp({arg})", coefs)

            coefs = RooArgList(x)
            for i in range(n):
                coefs.add(params[f'dypow_a{i}'])
            expo = " + ".join(f"@{i + 1}" + f"*{t}" * i for i in range(n))
            pdfs[f'dypow{n}'] = RooGenericPdf(f'dypow{n}_{tag}', f'dypow{n}',
                                              f"{core}*TMath::Power(@0/{DY_MREF!r},{expo})",
                                              coefs)

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

# ---------------------------------------------------------------------------
# Drell-Yan core of the dyexp / dypow families
# ---------------------------------------------------------------------------
# The leading-order q qbar -> gamma*/Z -> mu mu cross section for u, d and s
# quarks, each weighted by its parton luminosity: a degree-6 Chebyshev series
# over 70-190 GeV. Ported verbatim from flashggFinalFit,
# DYCore/src/RooDYCore.cc (its RooDYCore::evaluate). That repository does not
# record the beam energy or PDF set behind the luminosity coefficients; the
# fitted correction absorbs slow differences. No detector smearing.
#
# Written as a formula string rather than a compiled class, so the same
# function works in conda, in the Combine container and inside Combine.
DY_MAX_ORDER = 4           # dyexp1-4, dypow1-4
DY_MREF = 125.0            # the core is divided by its value here
DY_LUMI_RANGE = (70.0, 190.0)
# quark: (charge Q, weak isospin T3, Chebyshev coefficients of its luminosity)
DY_QUARKS = {
    "u": (2. / 3., 0.5, [5.65683e+06, -203204, 1628.22, -7.19829, 0.0182027,
                         -2.47312e-05, 1.40276e-08]),
    "d": (-1. / 3., -0.5, [4.62365e+06, -166744, 1337.52, -5.91684, 0.0149719,
                           -2.03548e-05, 1.15522e-08]),
    "s": (-1. / 3., -0.5, [2.68735e+06, -100109, 816.659, -3.64815,
                           0.00928968, -1.26864e-05, 7.22469e-09]),
}
DY_ALPHA = 1. / 137.035999
DY_GF = 1.1663787e-5
DY_MZ = 91.1876
DY_GZ = 2.4952
DY_SW2 = 0.23126


def _dy_kappa():
    return 2 ** 0.5 * DY_GF * DY_MZ ** 2 / (4 * 3.141592653589793 * DY_ALPHA)


def _dy_couplings(q, t3):
    """(photon, interference, Z) coefficients of one quark's cross section."""
    vq, aq = t3 - 2 * DY_SW2 * q, t3
    vl, al = -0.5 + 2 * DY_SW2, -0.5
    return q * q, -2 * q * vl * vq, (al * al + vl * vl) * (aq * aq + vq * vq)


def dy_core_value(m):
    """The core at mass m, unnormalised, line by line as
    RooDYCore::evaluate(). The reference dy_core_expr is checked against."""
    import math
    lo, hi = DY_LUMI_RANGE
    u = (2 * m - (lo + hi)) / (hi - lo)
    s = m * m
    d = (s - DY_MZ ** 2) ** 2 + DY_GZ ** 2 * DY_MZ ** 2
    k = _dy_kappa()
    chi_z, chi_int = k * k * s * s / d, k * s * (s - DY_MZ ** 2) / d
    val = 0.0
    for q, t3, c in DY_QUARKS.values():
        lum = sum(ci * math.cos(i * math.acos(u)) for i, ci in enumerate(c))
        a, b, z = _dy_couplings(q, t3)
        # the common prefactor 4 pi alpha^2 / (3 Nc) is dropped: it cancels
        # in the normalisation
        val += lum * (a + b * chi_int + z * chi_z) / s
    return val if val > 0 else 1e-12


def _cheb_to_power(c):
    """Chebyshev series coefficients -> power-series coefficients."""
    t = [[1.0], [0.0, 1.0]]
    while len(t) < len(c):          # T(n+1) = 2x T(n) - T(n-1)
        nxt = [0.0] + [2 * v for v in t[-1]]
        for i, v in enumerate(t[-2]):
            nxt[i] -= v
        t.append(nxt)
    out = [0.0] * len(c)
    for ci, tn in zip(c, t):
        for k, v in enumerate(tn):
            out[k] += ci * v
    return out


def _horner(coefs, var):
    s = repr(coefs[-1])
    for c in reversed(coefs[:-1]):
        s = f"({c!r}+{var}*{s})"
    return s


def dy_core_expr(m="@0"):
    """The core as a formula string in the mass `m`, normalised to 1 at
    DY_MREF. The three quark luminosities are summed per coupling into one
    polynomial each, so the string holds three polynomials, not nine."""
    lo, hi = DY_LUMI_RANGE
    u = f"((2*{m}-{lo + hi!r})/{hi - lo!r})"
    per_coupling = [[0.0] * 7 for _ in range(3)]
    for q, t3, c in DY_QUARKS.values():
        pw = _cheb_to_power(c)
        for j, g in enumerate(_dy_couplings(q, t3)):
            for k, v in enumerate(pw):
                per_coupling[j][k] += g * v
    p_gam, p_int, p_z = (_horner(p, u) for p in per_coupling)
    mz2, k = DY_MZ ** 2, _dy_kappa()
    s = f"({m}*{m})"
    d = f"(({s}-{mz2!r})*({s}-{mz2!r})+{DY_GZ ** 2 * mz2!r})"
    chi_int = f"({k!r}*{s}*({s}-{mz2!r})/{d})"
    chi_z = f"({k * k!r}*{s}*{s}/{d})"
    norm = dy_core_value(DY_MREF)
    return f"(({p_gam}+{p_int}*{chi_int}+{p_z}*{chi_z})/{s}/{norm!r})"


# background candidates and their number of FREE parameters (F-test, ndf)
BKG_NPAR = {
    **{f"bern{o}": o for o in range(1, 6)},        # bern_c0 is constant
    **{f"exp{n}": 2 * n - 1 for n in (1, 2, 3)},
    **{f"pow{n}": 2 * n - 1 for n in (1, 2, 3)},
    **{f"dyexp{n}": n for n in range(1, DY_MAX_ORDER + 1)},
    **{f"dypow{n}": n for n in range(1, DY_MAX_ORDER + 1)},
}

# candidates that are nested within one another, for the F-test
BKG_FAMILIES = {
    "bernstein": [f"bern{o}" for o in range(1, 6)],
    "exponential": [f"exp{n}" for n in (1, 2, 3)],
    "power": [f"pow{n}" for n in (1, 2, 3)],
    "dyexp": [f"dyexp{n}" for n in range(1, DY_MAX_ORDER + 1)],
    "dypow": [f"dypow{n}" for n in range(1, DY_MAX_ORDER + 1)],
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
