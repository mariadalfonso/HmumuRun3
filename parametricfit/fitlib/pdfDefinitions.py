# pdfDefinitions.py
import ROOT
from ROOT import RooRealVar, RooBernstein
from ROOT import RooExponential, RooAddPdf, RooGenericPdf, RooArgList, RooFormulaVar

from . import agnostic, bwz, dyCore, fewz, mcCore

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
    def create_bkg_pdfs(x, tag, mc_core=None):
        """
        Create all background PDFs - KEEP ALL OBJECTS IN MEMORY!

        mc_core: the frozen MC cores ({model: {param: value}}, from
        mcCoreFit.py); without it the mcexp/mcpow/mcdy families are absent.
        
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
        # c0 = 1 fixes the normalisation (only the ratios to it are
        # fitted); c1..c5 are the SQUARES of the parameters bern_c<k>, so
        # the polynomial is never negative (bwz.BERN_RANGE)
        params['bern_c0'] = RooRealVar(f'bern_c0_{tag}', 'bern_c0', 1.)
        coefs = [params['bern_c0']]
        for k in range(1, 6):
            params[f'bern_c{k}'] = p = RooRealVar(
                f'bern_c{k}_{tag}', f'bern_c{k}', bwz.BERN_INIT,
                *bwz.BERN_RANGE)
            p.setError(bwz.BERN_STEP)       # Minuit's first step
            formulas[f'bern_sq{k}'] = RooFormulaVar(
                f'bern_sq{k}_{tag}', f'bern_c{k}^2', '@0*@0', RooArgList(p))
            coefs.append(formulas[f'bern_sq{k}'])
        for o in range(1, 6):
            args = RooArgList()
            for c in coefs[:o + 1]:
                args.add(c)
            pdfs[f'bern{o}'] = RooBernstein(f'bern{o}_{tag}', f'bern{o}', x,
                                            args)

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

        # ===== DY CORE x CORRECTION PDFs (dyCore.py) =====
        dy_pdfs, dy_params = dyCore.create_pdfs(x, tag)
        pdfs.update(dy_pdfs)
        params.update(dy_params)

        # ===== BREIT-WIGNER Z PDFs, Run 2 note eqs. 5.1-5.4 (bwz.py) =====
        bw_pdfs, bw_params = bwz.create_pdfs(x, tag)
        pdfs.update(bw_pdfs)
        params.update(bw_params)

        # ===== FEWZ SPLINE x BERNSTEIN, Run 2 note eq. 5.5 (fewz.py) =====
        fz_pdfs, fz_params = fewz.create_pdfs(x, tag)
        pdfs.update(fz_pdfs)
        params.update(fz_params)

        # ===== EXP/POW + CONSTANT, LANDAU x BERNSTEIN, Run 2 note
        # eqs. 5.9-5.11 (agnostic.py) =====
        ag_pdfs, ag_params = agnostic.create_pdfs(x, tag)
        pdfs.update(ag_pdfs)
        params.update(ag_params)

        # ===== MC-TUNED CORE x CORRECTION PDFs (mcCore.py) =====
        if mc_core:
            mc_pdfs, mc_params = mcCore.create_pdfs(x, tag, mc_core)
            pdfs.update(mc_pdfs)
            params.update(mc_params)

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
    **{f"exp{n}": 2 * n - 1 for n in (1, 2, 3)},
    **{f"pow{n}": 2 * n - 1 for n in (1, 2, 3)},
    **dyCore.NPAR,
    **bwz.NPAR,
    **fewz.NPAR,
    **agnostic.NPAR,
    **mcCore.NPAR,                      # only with an MC core file
}

# candidates that are nested within one another, for the F-test
BKG_FAMILIES = {
    "bernstein": [f"bern{o}" for o in range(1, 6)],
    "exponential": [f"exp{n}" for n in (1, 2, 3)],
    "power": [f"pow{n}" for n in (1, 2, 3)],
    **dyCore.FAMILIES,                  # dyexp, dypow
    **bwz.FAMILIES,                     # bwzbern, bwzgamma, bwzredux
    **fewz.FAMILIES,                    # fewzbern
    **agnostic.FAMILIES,                # expconst, powconst, landbern
    **mcCore.FAMILIES,                  # mcexp, mcpow, mcdy
}

# The default candidates of bkgFit.py: the families of the Run 2 H->mumu
# analysis note, Sec. 5.3, the physics-inspired 5.1-5.5 and the agnostic
# 5.6-5.11. Every parameter floats on data; FEWZ is the only fixed input.
RUN2_FAMILIES = [
    "bwzbern", "bwzgamma", "bwzredux", "fewzbern",          # 5.1-5.5
    "bernstein", "exponential", "power",                    # 5.6-5.8
    "expconst", "powconst", "landbern",                     # 5.9-5.11
]
# bkgFit.py's default: the six families of the Run 2 note's single-function
# bias study (Sec. 6.4.1, Fig. 6.20). --all-run2 adds the rest of
# RUN2_FAMILIES (BWZ x Bernstein, Bernstein, exp / power law + constant).
DEFAULT_FAMILIES = [
    "bwzredux", "bwzgamma", "fewzbern",                     # 5.3-5.5
    "exponential", "power", "landbern",                     # 5.7, 5.8, 5.11
]
assert set(DEFAULT_FAMILIES) <= set(RUN2_FAMILIES)
# Beyond the Run 2 note, added by bkgFit.py --extra-models: the LO DY core
# x a data correction, and the cores tuned to MC (mcCoreFit.py) x exp(Cheb).
EXTRA_FAMILIES = [*dyCore.FAMILIES, *mcCore.FAMILIES]
assert sorted(RUN2_FAMILIES + EXTRA_FAMILIES) == sorted(BKG_FAMILIES)


def family_members(families):
    """Every candidate of the given BKG_FAMILIES, in order."""
    return [nm for fam in families for nm in BKG_FAMILIES[fam]]


def member_order(name):
    """bern5 -> 5, fewzbern3 -> 3; a member without a number (landau,
    bwzredux) -> 0"""
    digits = name[len(name.rstrip("0123456789")):]
    return int(digits) if digits else 0


def family_member(fam, order):
    """("fewzbern", 3) -> "fewzbern3": the member of that order"""
    found = [m for m in BKG_FAMILIES[fam] if member_order(m) == order]
    if not found:
        raise ValueError(f"{fam} has no order {order}: "
                         f"{', '.join(BKG_FAMILIES[fam])}")
    return found[0]


def bern_families():
    """The families built on a Bernstein polynomial (bern, BWZ x Bern,
    FEWZ x Bern, Landau x Bern): the ones labelled with their order, and
    whose order buildWorkspace.py --members families fixes across bins."""
    return [f for f in BKG_FAMILIES if "Bern{n}" in FAMILY_LABEL.get(f, "")]


# How the Bernstein parameters enter, written to every fit JSON's
# provenance ("bernstein_coefficients"): "squared", c_k = p_k^2 (bwz.py).
# A JSON without it is from before, when the stored value WAS c_k; read
# as p_k it would give another shape, so check_bern_coefficients stops.
BERN_COEFFICIENTS = "squared"


def bern_members(names):
    """The members of `names` that have Bernstein parameters (order >= 1
    of a Bernstein family: not bwz or landau)."""
    fam_of = {m: f for f, ms in BKG_FAMILIES.items() for m in ms}
    bern = set(bern_families())
    return [n for n in names if fam_of.get(n) in bern and member_order(n)]


def check_bern_coefficients(prov, names, source):
    """Stop if `source` (a fit JSON, provenance `prov`) stores Bernstein
    parameters of `names` in another form than BERN_COEFFICIENTS."""
    used = bern_members(names)
    form = prov.get("bernstein_coefficients", "signed (before squaring)")
    if used and form != BERN_COEFFICIENTS:
        raise SystemExit(
            f"{source}: Bernstein parameters stored as {form}, but the PDFs "
            f"now use {BERN_COEFFICIENTS} coefficients ({', '.join(used)}); "
            "rerun the fit that wrote it (bkgFit.py / corePdfFit.py)")


# Plot labels, short forms of the Run 2 note's names (Fig. 6.19, Table
# 6.13), ROOT TLatex; {n} is the member's trailing number (an order, or the
# terms of a sum for S-Exp / S-Pow).
FAMILY_LABEL = {
    "bwzbern": "BWZ #times Bern{n}",
    "bwzgamma": "BWZGamma",
    "bwzredux": "BWZRedux",
    "fewzbern": "FEWZ #times Bern{n}",
    "bernstein": "Bern{n}",
    "exponential": "S-Exp{n}",
    "power": "S-Pow{n}",
    "expconst": "Exp + const",
    "powconst": "Pow + const",
    "landbern": "Landau #times Bern{n}",
    "dyexp": "DY #times ExpCheb{n}",
    "dypow": "DY #times Pow{n}",
    "mcexp": "MC 2Exp #times ExpCheb{n}",
    "mcpow": "MC 2Pow #times ExpCheb{n}",
    "mcdy": "MC DY #times ExpCheb{n}",
}
# members whose name has no number, or a different function at order 0
MEMBER_LABEL = {"bwz": "BWZ", "landau": "Landau"}


def model_label(name):
    """bern5 -> "Bern5", exp2 -> "S-Exp2", fewzbern3 -> "FEWZ #times Bern3";
    unknown names as given."""
    if name in MEMBER_LABEL:
        return MEMBER_LABEL[name]
    fam = next((f for f, ms in BKG_FAMILIES.items() if name in ms), None)
    if fam not in FAMILY_LABEL:
        return name
    return FAMILY_LABEL[fam].format(n=name[len(name.rstrip("0123456789")):])

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


def create_bkg_pdf(x, tag, name, mc_core=None):
    """One background candidate by name, with its own (fresh) parameters.

    mc_core: the frozen MC cores, needed for the mcexp/mcpow/mcdy members.

    Returns (pdf, keep); keep["params"] are its FLOATING parameters.
    """
    if name not in BKG_NPAR:
        raise KeyError(f"unknown background model '{name}'; "
                       f"known: {sorted(BKG_NPAR)}")
    out = PDFDefinitions.create_bkg_pdfs(x, tag, mc_core)
    if name not in out['pdfs']:
        raise KeyError(f"'{name}' needs its MC core (mcCoreFit.py), "
                       f"which was not given")
    pdf = out['pdfs'][name]
    free = {p.GetName() for p in pdf.getParameters(ROOT.RooArgSet(x))
            if not p.isConstant()}
    pars = [v for v in out['params'].values() if v.GetName() in free]
    return pdf, {"params": pars, "pdf": pdf, "all": out}
