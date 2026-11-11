"""
agnostic.py -- the rest of the Run 2 note's "agnostic models" (Sec. 5.3.2).

Bernstein (5.6), the sum of exponentials (5.7) and the sum of power laws
(5.8) are the bern, exp and pow families of pdfDefinitions.py. This file
adds the other three, every parameter floating on the data:

    Exponential + constant   f + (1 - f) e^(a (m - 125))              (5.9)
    Power law + constant     f + (1 - f) (m/125)^a                     (5.10)
    Landau x Bern_n          Landau(m; mZ, s) * Bern_n(m)              (5.11)

5.9 and 5.10 are written about 125 GeV rather than as the note's e^(a m) and
m^a. That is the same set of shapes: the constant only moves the f that
gives a shape, f_note = f / (f + (1 - f) e^(125 a)). Written literally, the
falling term is ~1e-11 for a power law (125^-5), so f would sit at 1e-11
against its bound at 0 and Minuit could not fit it.

Landau is TMath::Landau with its location at mZ and a floating width s, as
in the note; Bern_n is the Bernstein of bwz.py (fit window mapped onto
[0, 1], first coefficient 1). Landau alone is Bern_0, so the two make one
nested family for the F-test:

    family     members                         free parameters
    expconst   expconst                        f, a          (2)
    powconst   powconst                        f, a          (2)
    landbern   landau, landbern1 ... landbern4 s, c1..cn     (1 + n)
"""

from ROOT import RooArgList, RooGenericPdf, RooRealVar

from .bwz import BERN_INIT, BERN_RANGE, BERN_STEP, bernstein
from .dyCore import MZ

MREF = 125.0
MAX_ORDER = 4               # landbern1-4

# start value, range and first Minuit step of every parameter
F_INIT, F_RANGE, F_STEP = 0.05, (0.0, 1.0), 0.01        # the constant
EA_INIT, EA_RANGE, EA_STEP = -0.05, (-1.0, 1.0), 1e-3   # e^(a (m-125))
PA_INIT, PA_RANGE, PA_STEP = -5.0, (-50.0, 50.0), 0.05  # (m/125)^a
S_INIT, S_RANGE, S_STEP = 5.0, (0.1, 100.0), 0.1        # Landau width, GeV

NPAR = {"expconst": 2, "powconst": 2, "landau": 1,
        **{f"landbern{n}": 1 + n for n in range(1, MAX_ORDER + 1)}}
FAMILIES = {"expconst": ["expconst"],
            "powconst": ["powconst"],
            "landbern": ["landau"] + [f"landbern{n}"
                                      for n in range(1, MAX_ORDER + 1)]}


def create_pdfs(x, tag):
    """Every expconst / powconst / landau / landbern candidate in x.

    Returns (pdfs, params), keyed like PDFDefinitions.create_bkg_pdfs.
    """
    lo, hi = x.getMin(), x.getMax()
    u = f"((@0-{lo!r})/{hi - lo!r})"
    params, pdfs = {}, {}

    def var(key, init, rng, step):
        params[key] = v = RooRealVar(f"{key}_{tag}", key, init, *rng)
        v.setError(step)        # Minuit's first step
        return v

    # ---- exponential + constant (5.9) ----------------------------------
    f = var("expconst_f", F_INIT, F_RANGE, F_STEP)
    a = var("expconst_a", EA_INIT, EA_RANGE, EA_STEP)
    pdfs["expconst"] = RooGenericPdf(
        f"expconst_{tag}", "expconst",
        f"(@1+(1-@1)*exp(@2*(@0-{MREF!r})))", RooArgList(x, f, a))

    # ---- power law + constant (5.10) -----------------------------------
    f = var("powconst_f", F_INIT, F_RANGE, F_STEP)
    a = var("powconst_a", PA_INIT, PA_RANGE, PA_STEP)
    pdfs["powconst"] = RooGenericPdf(
        f"powconst_{tag}", "powconst",
        f"(@1+(1-@1)*TMath::Power(@0/{MREF!r},@2))", RooArgList(x, f, a))

    # ---- Landau, and Landau x Bern_n (5.11) ----------------------------
    s = var("landbern_s", S_INIT, S_RANGE, S_STEP)
    bern = [var(f"landbern_c{k}", BERN_INIT, BERN_RANGE, BERN_STEP)
            for k in range(1, MAX_ORDER + 1)]
    landau = f"TMath::Landau(@0,{MZ!r},@1)"
    pdfs["landau"] = RooGenericPdf(f"landau_{tag}", "landau", landau,
                                   RooArgList(x, s))
    for n in range(1, MAX_ORDER + 1):
        args = RooArgList(x, s)
        for v in bern[:n]:
            args.add(v)
        coefs = ["1.0"] + [f"@{k + 1}" for k in range(1, n + 1)]
        pdfs[f"landbern{n}"] = RooGenericPdf(
            f"landbern{n}_{tag}", f"landbern{n}",
            f"{landau}*{bernstein(u, coefs)}", args)
    return pdfs, params
