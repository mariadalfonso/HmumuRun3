"""
bwz.py -- the Breit-Wigner Z families of the Run 2 H->mumu analysis note.

The "physics-inspired models" of the Run 2 note, Sec. 5.3.1, written as in
eqs. 5.1-5.4 (not as flashggFinalFit's PdfModelBuilder writes them). Every
parameter floats on the data; mZ and GammaZ are fixed (dyCore.MZ, GZ):

    BWZ(m)            = GZ e^(a m) / ((m - mZ)^2 + (GZ/2)^2)            (5.4)
    BWZ x Bern_n(m)   = BWZ(m) * Bern_n(m)                             (5.1)
    BWZGamma(m)       = f GZ e^(a m) / ((m - mZ)^2 + (GZ/2)^2)
                        + (1 - f) e^(a m) / m^2                        (5.2)
    BWZRedux(m)       = GZ e^(a t + b t^2) / ((m - mZ)^c + (GZ/2)^c)   (5.3)
                        t = (m - m0) / h, m0 and h the window's centre and
                        half-width (130 and 20 GeV): t in [-1, 1]

Bern_n is the Bernstein polynomial of order n on the fit window mapped onto
[0, 1], with its first coefficient fixed to 1 (only the shape is fitted)
and the others the squares of the fitted parameters, c_k = p_k^2, so that
it is never negative (BERN_RANGE). It is written out here rather than with
RooBernstein so the mapping never depends on the fit range. BWZ is Bern_0, so BWZ and BWZ x Bern_n are one nested
family for the F-test:

    family    members                        free parameters
    bwzbern   bwz, bwzbern1 ... bwzbern4     a, Bernstein c1..cn   (1 + n)
    bwzgamma  bwzgamma                       a, f                  (2)
    bwzredux  bwzredux                       a, b, c               (3)

The note writes the exponentials in m itself, not in m - 125; so does this
file, except BWZRedux. There e^(a m + b m^2) has a and b almost fully
correlated (m^2 ~ 2 m0 m over the window) and on scales 1e4 apart, so
Minuit forces its covariance positive-definite (status 1) on every fit,
and bkgFit.py drops it. e^(a t + b t^2) is the same family (a, b map
linearly onto the note's; the constant goes into the normalisation): the
same minimum, and toy fits converge with an accurate covariance about
70% of the time. In the window m > mZ, so (m - mZ)^c is real for any c.

Everything is a formula string in standard RooFit classes, so the same
functions work in conda, in the Combine container and in Combine, as for
dyCore.py.
"""

from ROOT import RooArgList, RooGenericPdf, RooRealVar

from .dyCore import GZ, MZ

MAX_ORDER = 4               # bwzbern1-4

# start value, range and first Minuit step of every parameter. The steps
# matter: a multiplies m ~ 130 GeV, so Minuit's default first step (a
# tenth of the range) would move the exponent by 6.
A_INIT, A_RANGE, A_STEP = 0.0, (-0.5, 0.5), 1e-3      # e^(a m), per GeV
F_INIT, F_RANGE, F_STEP = 0.1, (0.0, 1.0), 0.02       # BWZGamma: BW fraction
C_INIT, C_RANGE, C_STEP = 2.0, (0.5, 5.0), 0.05       # BWZRedux: BW exponent
# BWZRedux's exponent a t + b t^2, t in [-1, 1]. The note's form in m with
# a +/- 0.5 per GeV and b +/- 2e-3 per GeV^2 is about a +/- 10, b +/- 0.8
RA_INIT, RA_RANGE, RA_STEP = 0.0, (-20.0, 20.0), 0.02
RB_INIT, RB_RANGE, RB_STEP = 0.0, (-10.0, 10.0), 0.02
# Bernstein coefficients, every family (bern, BWZ / FEWZ / Landau x Bern):
# each parameter p_k enters SQUARED, c_k = p_k^2 (bernstein()), so the
# polynomial is >= 0 on the whole window, as flashggFinalFit's
# PdfModelBuilder::getBernstein. With signed c_k a toy fit could drive the
# background negative under the peak and trade it against r; RooFit then
# forces the PDF to 0 there and Minuit fails (bias study, Nov 2026).
# p_k in [-5, 5] is c_k in [0, 25]; p_k = 1 starts every c_k at 1 (a flat
# polynomial), away from p_k = 0, where dc_k/dp_k = 0.
BERN_INIT, BERN_RANGE, BERN_STEP = 1.0, (-5.0, 5.0), 0.05

NPAR = {"bwz": 1,
        **{f"bwzbern{n}": 1 + n for n in range(1, MAX_ORDER + 1)},
        "bwzgamma": 2,
        "bwzredux": 3}
FAMILIES = {"bwzbern": ["bwz"] + [f"bwzbern{n}"
                                  for n in range(1, MAX_ORDER + 1)],
            "bwzgamma": ["bwzgamma"],
            "bwzredux": ["bwzredux"]}

_GZ, _MZ, _HALF = repr(GZ), repr(MZ), repr(GZ / 2.0)


def _bw(m, expo, c=None):
    """GZ e^(expo) / ((m - mZ)^c + (GZ/2)^c); c = None is the square."""
    if c is None:
        den = f"(({m}-{_MZ})*({m}-{_MZ})+{_HALF}*{_HALF})"
    else:
        den = f"(TMath::Power({m}-{_MZ},{c})+TMath::Power({_HALF},{c}))"
    return f"({_GZ}*exp({expo})/{den})"


def _binom(n, k):
    out = 1
    for i in range(k):
        out = out * (n - i) // (i + 1)
    return out


def bernstein(u, coefs):
    """sum_k coefs[k]^2 C(n,k) u^k (1-u)^(n-k), n = len(coefs) - 1: every
    coefficient squared (BERN_RANGE), so the sum is >= 0 for u in [0, 1].
    coefs are formula terms ("1.0", "@2")."""
    n = len(coefs) - 1
    terms = []
    for k, ck in enumerate(coefs):
        t = ([f"{_binom(n, k)}*{ck}*{ck}"] + [u] * k
             + [f"(1-{u})"] * (n - k))
        terms.append("*".join(t))
    return "(" + "+".join(terms) + ")"


def create_pdfs(x, tag):
    """Every bwz / bwzbern / bwzgamma / bwzredux candidate in x.

    Returns (pdfs, params), keyed like PDFDefinitions.create_bkg_pdfs.
    """
    lo, hi = x.getMin(), x.getMax()
    u = f"((@0-{lo!r})/{hi - lo!r})"
    params, pdfs = {}, {}

    def var(key, init, rng, step):
        params[key] = v = RooRealVar(f"{key}_{tag}", key, init, *rng)
        v.setError(step)        # Minuit's first step
        return v

    # ---- BWZ (5.4) and BWZ x Bern_n (5.1) ------------------------------
    a = var("bwzbern_a", A_INIT, A_RANGE, A_STEP)
    bern = [var(f"bwzbern_c{k}", BERN_INIT, BERN_RANGE, BERN_STEP)
            for k in range(1, MAX_ORDER + 1)]
    pdfs["bwz"] = RooGenericPdf(f"bwz_{tag}", "bwz", _bw("@0", "@1*@0"),
                                RooArgList(x, a))
    for n in range(1, MAX_ORDER + 1):
        args = RooArgList(x, a)
        for v in bern[:n]:
            args.add(v)
        coefs = ["1.0"] + [f"@{i + 2}" for i in range(n)]
        pdfs[f"bwzbern{n}"] = RooGenericPdf(
            f"bwzbern{n}_{tag}", f"bwzbern{n}",
            f"{_bw('@0', '@1*@0')}*{bernstein(u, coefs)}", args)

    # ---- BWZGamma (5.2) ------------------------------------------------
    a = var("bwzgamma_a", A_INIT, A_RANGE, A_STEP)
    f = var("bwzgamma_f", F_INIT, F_RANGE, F_STEP)
    pdfs["bwzgamma"] = RooGenericPdf(
        f"bwzgamma_{tag}", "bwzgamma",
        f"(@2*{_bw('@0', '@1*@0')}+(1-@2)*exp(@1*@0)/(@0*@0))",
        RooArgList(x, a, f))

    # ---- BWZRedux (5.3) ------------------------------------------------
    a = var("bwzredux_a", RA_INIT, RA_RANGE, RA_STEP)
    b = var("bwzredux_b", RB_INIT, RB_RANGE, RB_STEP)
    c = var("bwzredux_c", C_INIT, C_RANGE, C_STEP)
    t = f"((@0-{0.5 * (lo + hi)!r})/{0.5 * (hi - lo)!r})"
    pdfs["bwzredux"] = RooGenericPdf(
        f"bwzredux_{tag}", "bwzredux",
        _bw("@0", f"@1*{t}+@2*{t}*{t}", "@3"), RooArgList(x, a, b, c))

    return pdfs, params
