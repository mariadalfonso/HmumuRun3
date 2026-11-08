"""
mcCore.py -- background families on a core tuned to the background MC.

    background(m) = core_MC(m) * exp(a1*T1(t) + ... + an*Tn(t))

The flashggFinalFit ggH setup. The core is a smooth shape fitted to the
category's background MC once (mcCoreFit.py) and then FROZEN: its
parameters are numbers in the formula, never floated on data. Only the
correction is fitted on the data sideband, exactly as for dyexp. Three
cores, after flashgg DYCore/src/RooDYTexpcore, RooDYTPowcore, RooDYCoreMod:

    mcexp  f*exp(b1*t + b2*t^2 + b3*t^3) + (1-f)*exp(b4*t)          (2Exp)
    mcpow  f*x^-(b1 + b2*d + b3*d^2) + (1-f)*x^-b4                  (2Pow)
    mcdy   LO DY core (dyCore.py) * (1 + c1*T1 + c2*T2 + c3*T3)

with t = the fit window mapped onto [-1, 1], x = m/125, d = x - 1 and
f = 1/(1 + exp(-theta)). flashgg writes the 2Pow exponent in x rather than
d; a quadratic in x is a quadratic in d, so it is the same shape, but the
fit is far better conditioned around x = 1.

The correction is ours, not flashgg's 1 + 0.5*sum(a_i T_i): the exponential
keeps the PDF positive for any coefficient.

mcCoreFit.py writes the frozen numbers to WS_<year>/MCcore_<tag>.json;
bkgFit.py copies them into Bkg_<tag>_params.json, so makeWorkspaces.py
rebuilds the same PDFs without the MC.
"""

from ROOT import RooArgList, RooGenericPdf, RooRealVar

from . import dyCore

MAX_ORDER = 4               # mcexp1-4, mcpow1-4, mcdy1-4
COEF_RANGE = dyCore.COEF_RANGE

# core model -> (parameter names, start values, ranges) for the MC fit
CORES = {
    "mcexp": (["theta", "b1", "b2", "b3", "b4"],
              [0.0, -1.0, 0.0, 0.0, -0.3],
              [(-10, 10), (-20, 20), (-20, 20), (-20, 20), (-20, 20)]),
    "mcpow": (["theta", "b1", "b2", "b3", "b4"],
              [0.0, 4.0, 0.0, 0.0, 6.0],
              [(-10, 10), (-50, 50), (-200, 200), (-1000, 1000), (-50, 50)]),
    "mcdy": (["c1", "c2", "c3"],
             [0.0, 0.0, 0.0],
             [(-2, 2), (-2, 2), (-2, 2)]),
}
CORE_LABEL = {"mcexp": "2Exp", "mcpow": "2Pow", "mcdy": "DY core #times Cheb3"}

# number of FREE parameters of each candidate (the core adds none), and the
# nested families; merged into BKG_NPAR and BKG_FAMILIES in pdfDefinitions
NPAR = {f"{c}{n}": n for c in CORES for n in range(1, MAX_ORDER + 1)}
FAMILIES = {c: [f"{c}{n}" for n in range(1, MAX_ORDER + 1)] for c in CORES}


def _t(m, lo, hi):
    return f"((2*{m}-{lo + hi!r})/{hi - lo!r})"


def _cheb(t):
    """T1..T4 written out: cos(n*acos(t)) is NaN for |t| > 1."""
    return [t, f"(2*{t}*{t}-1)", f"(4*{t}*{t}*{t}-3*{t})",
            f"(8*{t}*{t}*{t}*{t}-8*{t}*{t}+1)"]


def core_expr(model, p, m, lo, hi):
    """The core as a formula string in the mass `m`.

    p: one string per parameter of CORES[model], either a reference
    ("@1") for the MC fit or a number (repr) once frozen.
    """
    t = _t(m, lo, hi)
    if model == "mcexp":
        th, b1, b2, b3, b4 = p
        f = f"(1/(1+exp(-{th})))"
        return (f"({f}*exp({b1}*{t}+{b2}*{t}*{t}+{b3}*{t}*{t}*{t})"
                f"+(1-{f})*exp({b4}*{t}))")
    if model == "mcpow":
        th, b1, b2, b3, b4 = p
        f = f"(1/(1+exp(-{th})))"
        x = f"({m}/{dyCore.MREF!r})"
        d = f"({x}-1)"
        return (f"({f}*TMath::Power({x},-({b1}+{b2}*{d}+{b3}*{d}*{d}))"
                f"+(1-{f})*TMath::Power({x},-{b4}))")
    if model == "mcdy":
        c1, c2, c3 = p
        T = _cheb(t)
        corr = f"(1+{c1}*{T[0]}+{c2}*{T[1]}+{c3}*{T[2]})"
        # guard: a linear correction can dip below zero far from the MC
        return f"({dyCore.core_expr(m)}*TMath::Max({corr},1e-9))"
    raise KeyError(f"unknown MC core '{model}' (have {sorted(CORES)})")


def make_fit_pdf(x, model, tag):
    """The core with FLOATING parameters, for the fit to the MC.

    Returns (pdf, [RooRealVar]); keep both alive.
    """
    names, starts, ranges = CORES[model]
    pars = [RooRealVar(f"{model}_core_{n}_{tag}", n, s, *r)
            for n, s, r in zip(names, starts, ranges)]
    args = RooArgList(x)
    for v in pars:
        args.add(v)
    refs = [f"@{i + 1}" for i in range(len(pars))]
    expr = core_expr(model, refs, "@0", x.getMin(), x.getMax())
    pdf = RooGenericPdf(f"{model}_core_{tag}", f"{model} core", expr, args)
    return pdf, pars


def frozen_expr(model, values, m, lo, hi):
    """The core with its MC-fit values written in as numbers.

    values: {parameter name: value} as mcCoreFit.py stores it.
    """
    names = CORES[model][0]
    missing = [n for n in names if n not in values]
    if missing:
        raise KeyError(f"MC core '{model}' lacks {missing}")
    return core_expr(model, [repr(float(values[n])) for n in names],
                     m, lo, hi)


def create_pdfs(x, tag, mc_core):
    """Every mcexp / mcpow / mcdy candidate whose core is in mc_core.

    mc_core: {model: {parameter: value}}, the "cores" block of
    MCcore_<tag>.json. A model missing there gives no candidates.
    Returns (pdfs, params), keyed like PDFDefinitions.create_bkg_pdfs.
    """
    lo, hi = x.getMin(), x.getMax()
    cheb = _cheb(_t("@0", lo, hi))
    pdfs, params = {}, {}
    for model in CORES:
        if model not in mc_core:
            continue
        core = frozen_expr(model, mc_core[model], "@0", lo, hi)
        for n in range(1, MAX_ORDER + 1):
            params[f"{model}_a{n}"] = RooRealVar(
                f"{model}_a{n}_{tag}", f"{model}_a{n}", 0., *COEF_RANGE)
        for n in range(1, MAX_ORDER + 1):
            coefs = RooArgList(x)
            for i in range(1, n + 1):
                coefs.add(params[f"{model}_a{i}"])
            arg = " + ".join(f"@{i}*{cheb[i - 1]}" for i in range(1, n + 1))
            pdfs[f"{model}{n}"] = RooGenericPdf(
                f"{model}{n}_{tag}", f"{model}{n}",
                f"{core}*exp({arg})", coefs)
    return pdfs, params


def is_member(name):
    """True for an mcexp / mcpow / mcdy candidate name."""
    return name.rstrip("0123456789") in CORES
