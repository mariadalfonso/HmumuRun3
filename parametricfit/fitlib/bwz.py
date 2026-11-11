"""
bwz.py -- the Breit-Wigner Z families of the Run 2 H->mumu analysis note.

The "physics-inspired models" of the Run 2 note, Sec. 5.3.1, written as in
eqs. 5.1-5.4 (not as flashggFinalFit's PdfModelBuilder writes them). Every
parameter floats on the data; mZ and GammaZ are fixed (dyCore.MZ, GZ):

    BWZ(m)            = GZ e^(a m) / ((m - mZ)^2 + (GZ/2)^2)            (5.4)
    BWZ x Bern_n(m)   = BWZ(m) * Bern_n(m)                             (5.1)
    BWZGamma(m)       = f GZ e^(a m) / ((m - mZ)^2 + (GZ/2)^2)
                        + (1 - f) e^(a m) / m^2                        (5.2)
    BWZRedux(m)       = GZ e^(a m + b m^2) / ((m - mZ)^c + (GZ/2)^c)   (5.3)

Bern_n is the Bernstein polynomial of order n on the fit window mapped onto
[0, 1], with its first coefficient fixed to 1 (only the shape is fitted),
written out here rather than with RooBernstein so the mapping never depends
on the fit range. BWZ is Bern_0, so BWZ and BWZ x Bern_n are one nested
family for the F-test:

    family    members                        free parameters
    bwzbern   bwz, bwzbern1 ... bwzbern4     a, Bernstein c1..cn   (1 + n)
    bwzgamma  bwzgamma                       a, f                  (2)
    bwzredux  bwzredux                       a, b, c               (3)

The note writes the exponentials in m itself, not in m - 125; so does this
file. In the window m > mZ, so (m - mZ)^c is real for any c.

Everything is a formula string in standard RooFit classes, so the same
functions work in conda, in the Combine container and in Combine, as for
dyCore.py.
"""

from ROOT import RooArgList, RooGenericPdf, RooRealVar

from .dyCore import GZ, MZ

MAX_ORDER = 4               # bwzbern1-4

# start value, range and first Minuit step of every parameter. The steps
# matter: b multiplies m^2 ~ 2e4 GeV^2, so Minuit's default first step (a
# tenth of the range) would move the exponent by tens.
A_INIT, A_RANGE, A_STEP = 0.0, (-0.5, 0.5), 1e-3      # e^(a m), per GeV
F_INIT, F_RANGE, F_STEP = 0.1, (0.0, 1.0), 0.02       # BWZGamma: BW fraction
B_INIT, B_RANGE, B_STEP = 0.0, (-2e-3, 2e-3), 1e-6    # BWZRedux: per GeV^2
C_INIT, C_RANGE, C_STEP = 2.0, (0.5, 5.0), 0.05       # BWZRedux: BW exponent
BERN_INIT, BERN_RANGE, BERN_STEP = 1.0, (-10.0, 10.0), 0.05   # as bern

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
    """sum_k coefs[k] C(n,k) u^k (1-u)^(n-k), n = len(coefs) - 1."""
    n = len(coefs) - 1
    terms = []
    for k, ck in enumerate(coefs):
        t = [f"{_binom(n, k)}*{ck}"] + [u] * k + [f"(1-{u})"] * (n - k)
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
    a = var("bwzredux_a", A_INIT, A_RANGE, A_STEP)
    b = var("bwzredux_b", B_INIT, B_RANGE, B_STEP)
    c = var("bwzredux_c", C_INIT, C_RANGE, C_STEP)
    pdfs["bwzredux"] = RooGenericPdf(
        f"bwzredux_{tag}", "bwzredux",
        _bw("@0", "@1*@0+@2*@0*@0", "@3"), RooArgList(x, a, b, c))

    return pdfs, params
