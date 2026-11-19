"""
fewz.py -- the FEWZ x Bernstein family of the Run 2 H->mumu analysis note.

    FEWZ x Bern_n(m) = FEWZ-spline(m) * Bern_n(m)                      (5.5)

The FEWZ Drell-Yan prediction (NNLO QCD, NLO EW) in bins of m_mumu,
interpolated by a cubic spline, times a Bernstein polynomial that absorbs
the non-Z/gamma* backgrounds and the detector effects (Run 2 note, Sec.
5.3.1). The spline is fixed; only the Bernstein coefficients float.

The spline is the TSpline3 `fewz_spline` in fewz_spline.root next to this
file, copied unchanged from flashggFinalFit (LPC-HH),
Background/DYshapes/fewz_spline.root, commit f7924bd: 280 knots, 0.25 GeV
apart, over 95-165 GeV. The file also holds `fewz_spline_smear`, which
flashgg uses; we take the unsmeared one, as the note leaves the detector
effects to the Bernstein. Within 110-150 GeV the two differ by < 1%.

The spline is written out as a piecewise cubic formula from its own
coefficients (exact to rounding), so the PDF is a plain RooGenericPdf: the
same code in conda, in the Combine container and in Combine, and the
workspace does not need this file.

    family    members                      free parameters
    fewzbern  fewzbern1 ... fewzbern4      Bernstein c1..cn    (n)

Bern_n is the same Bernstein as bwz.py: on the fit window mapped onto
[0, 1], first coefficient fixed to 1, the others c_k = p_k^2 (never
negative; bwz.BERN_RANGE).

The normalisation is numeric, and RooFit's default 1D integrator (Romberg,
at most 20 steps) gets the integral of this 161-piece formula wrong by
~8e-4, by an amount that changes with the coefficients. With ~1e6 events
that is enough to move the minimum: fewzbern3 fitted to its own Asimov in
ggHcat bdt0 (5 bins) gave a Bias of -1.0. Every FEWZ pdf therefore carries
its own integrator, adaptive Gauss-Kronrod (error < 1e-7); the setting is
stored with the pdf, so Combine reading the workspace uses it too.
"""

import ctypes
import functools
import os

from ROOT import RooArgList, RooGenericPdf, RooRealVar, TFile

from .bwz import BERN_INIT, BERN_RANGE, BERN_STEP, bernstein

FEWZ_FILE = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                         "fewz_spline.root")
SPLINE = "fewz_spline"
MREF = 125.0                # the spline is divided by its value here
MAX_ORDER = 4               # fewzbern1-4

# the normalisation integral: see the docstring. A tighter tolerance only
# adds "roundoff error" messages (the pieces meet to rounding), not accuracy
# that matters
INT_METHOD = "RooAdaptiveGaussKronrodIntegrator1D"
INT_EPS = 1e-5              # relative and absolute; the result is ~4e-8 off
INT_MAX_SEGMENTS = 1000     # the default 100 is too few for 161 pieces

NPAR = {f"fewzbern{n}": n for n in range(1, MAX_ORDER + 1)}
FAMILIES = {"fewzbern": list(NPAR)}


@functools.lru_cache(maxsize=None)
def spline_expr(lo, hi, m="@0"):
    """The spline over [lo, hi] as a formula string in `m`, 1 at MREF.

    One cubic per knot interval that overlaps the window,
    y + b*dx + c*dx^2 + d*dx^3 with dx = m - x_knot (TSpline3::GetCoeff);
    the first and last pieces carry on past the window edges.
    """
    f = TFile.Open(FEWZ_FILE)
    if not f or f.IsZombie():
        raise RuntimeError(f"cannot open {FEWZ_FILE}")
    sp = f.Get(SPLINE)
    if not sp:
        raise RuntimeError(f"{FEWZ_FILE} has no {SPLINE}")
    if not sp.GetXmin() <= lo < hi <= sp.GetXmax():
        raise RuntimeError(f"fit window {lo}-{hi} GeV is outside the FEWZ "
                           f"spline ({sp.GetXmin()}-{sp.GetXmax()} GeV)")
    norm = sp.Eval(MREF)
    pieces = []
    for i in range(sp.GetNp() - 1):
        x0, y, b, c, d = (ctypes.c_double() for _ in range(5))
        sp.GetCoeff(i, x0, y, b, c, d)
        x1, y1 = ctypes.c_double(), ctypes.c_double()
        sp.GetKnot(i + 1, x1, y1)
        if x1.value > lo and x0.value < hi:
            pieces.append((x0.value, x1.value,
                           *(v.value / norm for v in (y, b, c, d))))
    f.Close()

    terms = []
    for k, (x0, x1, y, b, c, d) in enumerate(pieces):
        cond = []
        if k > 0:
            cond.append(f"({m}>={x0!r})")
        if k < len(pieces) - 1:
            cond.append(f"({m}<{x1!r})")
        dx = f"({m}-{x0!r})"
        poly = f"({y!r}+{dx}*({b!r}+{dx}*({c!r}+{dx}*{d!r})))"
        terms.append("*".join(cond + [poly]))
    return "(" + "+".join(terms) + ")"


def precise_integration(pdf):
    """Give `pdf` its own numeric integrator (INT_METHOD), kept with it
    when it is written to a workspace."""
    cfg = pdf.specialIntegratorConfig(True)
    cfg.method1D().setLabel(INT_METHOD)
    cfg.setEpsRel(INT_EPS)
    cfg.setEpsAbs(INT_EPS)
    cfg.getConfigSection(INT_METHOD).setRealValue("maxSeg", INT_MAX_SEGMENTS)
    return pdf


def create_pdfs(x, tag):
    """Every fewzbern candidate in x. Returns (pdfs, params), keyed like
    PDFDefinitions.create_bkg_pdfs."""
    lo, hi = x.getMin(), x.getMax()
    core = spline_expr(lo, hi)
    u = f"((@0-{lo!r})/{hi - lo!r})"
    params, pdfs = {}, {}
    for k in range(1, MAX_ORDER + 1):
        key = f"fewzbern_c{k}"
        params[key] = v = RooRealVar(f"{key}_{tag}", key, BERN_INIT,
                                     *BERN_RANGE)
        v.setError(BERN_STEP)       # Minuit's first step
    for n in range(1, MAX_ORDER + 1):
        args = RooArgList(x)
        for k in range(1, n + 1):
            args.add(params[f"fewzbern_c{k}"])
        coefs = ["1.0"] + [f"@{k}" for k in range(1, n + 1)]
        pdfs[f"fewzbern{n}"] = precise_integration(RooGenericPdf(
            f"fewzbern{n}_{tag}", f"fewzbern{n}",
            f"{core}*{bernstein(u, coefs)}", args))
    return pdfs, params
