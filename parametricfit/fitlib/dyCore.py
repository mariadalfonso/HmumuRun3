"""
dyCore.py -- the Drell-Yan core background families, dyexp and dypow.

    background(m) = core(m) * correction(m)

The core is the leading-order q qbar -> gamma*/Z -> mu mu cross section for
u, d and s quarks, each weighted by its parton luminosity: a degree-6
Chebyshev series over 70-190 GeV. It is ported verbatim from flashggFinalFit,
DYCore/src/RooDYCore.cc (RooDYCore::evaluate). That repository does not
record the beam energy or PDF set behind the luminosity coefficients; the
fitted correction absorbs slow differences. No detector smearing.

The core has no free parameters; only the correction is fitted. With t = the
fit window mapped onto [-1, 1]:

    dyexp<n> = core * exp(a1*T1(t) + ... + an*Tn(t))         (ATLAS HixEpoly)
    dypow<n> = core * (m/125)^(a0 + a1*t + ... + a(n-1)*t^(n-1))  (HixPower)

Everything is a formula string in standard RooFit classes, not a compiled
class, so the same functions work in conda, in the Combine container and in
Combine. pdfDefinitions.py adds these families to the background candidates.
"""

import math

from ROOT import RooArgList, RooGenericPdf, RooRealVar

MAX_ORDER = 4              # dyexp1-4, dypow1-4
MREF = 125.0               # the core is divided by its value here
COEF_RANGE = (-10.0, 10.0)  # every correction coefficient; starts at 0

# number of FREE parameters of each candidate, and the nested families,
# merged into BKG_NPAR and BKG_FAMILIES in pdfDefinitions.py
NPAR = {**{f"dyexp{n}": n for n in range(1, MAX_ORDER + 1)},
        **{f"dypow{n}": n for n in range(1, MAX_ORDER + 1)}}
FAMILIES = {"dyexp": [f"dyexp{n}" for n in range(1, MAX_ORDER + 1)],
            "dypow": [f"dypow{n}" for n in range(1, MAX_ORDER + 1)]}

# ---------------------------------------------------------------------------
# the core, as RooDYCore
# ---------------------------------------------------------------------------
LUMI_RANGE = (70.0, 190.0)
# quark: (charge Q, weak isospin T3, Chebyshev coefficients of its luminosity)
QUARKS = {
    "u": (2. / 3., 0.5, [5.65683e+06, -203204, 1628.22, -7.19829, 0.0182027,
                         -2.47312e-05, 1.40276e-08]),
    "d": (-1. / 3., -0.5, [4.62365e+06, -166744, 1337.52, -5.91684, 0.0149719,
                           -2.03548e-05, 1.15522e-08]),
    "s": (-1. / 3., -0.5, [2.68735e+06, -100109, 816.659, -3.64815,
                           0.00928968, -1.26864e-05, 7.22469e-09]),
}
ALPHA = 1. / 137.035999
GF = 1.1663787e-5
MZ = 91.1876
GZ = 2.4952
SW2 = 0.23126


def _kappa():
    return 2 ** 0.5 * GF * MZ ** 2 / (4 * math.pi * ALPHA)


def _couplings(q, t3):
    """(photon, interference, Z) coefficients of one quark's cross section."""
    vq, aq = t3 - 2 * SW2 * q, t3
    vl, al = -0.5 + 2 * SW2, -0.5
    return q * q, -2 * q * vl * vq, (al * al + vl * vl) * (aq * aq + vq * vq)


def core_value(m):
    """The core at mass m, unnormalised, line by line as
    RooDYCore::evaluate(). The reference core_expr is checked against."""
    lo, hi = LUMI_RANGE
    u = (2 * m - (lo + hi)) / (hi - lo)
    s = m * m
    d = (s - MZ ** 2) ** 2 + GZ ** 2 * MZ ** 2
    k = _kappa()
    chi_z, chi_int = k * k * s * s / d, k * s * (s - MZ ** 2) / d
    val = 0.0
    for q, t3, c in QUARKS.values():
        lum = sum(ci * math.cos(i * math.acos(u)) for i, ci in enumerate(c))
        a, b, z = _couplings(q, t3)
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


def core_expr(m="@0"):
    """The core as a formula string in the mass `m`, normalised to 1 at
    MREF. The three quark luminosities are summed per coupling into one
    polynomial each, so the string holds three polynomials, not nine."""
    lo, hi = LUMI_RANGE
    u = f"((2*{m}-{lo + hi!r})/{hi - lo!r})"
    per_coupling = [[0.0] * 7 for _ in range(3)]
    for q, t3, c in QUARKS.values():
        pw = _cheb_to_power(c)
        for j, g in enumerate(_couplings(q, t3)):
            for k, v in enumerate(pw):
                per_coupling[j][k] += g * v
    p_gam, p_int, p_z = (_horner(p, u) for p in per_coupling)
    mz2, k = MZ ** 2, _kappa()
    s = f"({m}*{m})"
    d = f"(({s}-{mz2!r})*({s}-{mz2!r})+{GZ ** 2 * mz2!r})"
    chi_int = f"({k!r}*{s}*({s}-{mz2!r})/{d})"
    chi_z = f"({k * k!r}*{s}*{s}/{d})"
    norm = core_value(MREF)
    return f"(({p_gam}+{p_int}*{chi_int}+{p_z}*{chi_z})/{s}/{norm!r})"


# ---------------------------------------------------------------------------
# the families
# ---------------------------------------------------------------------------

def create_pdfs(x, tag):
    """Every dyexp / dypow candidate in the observable x.

    Returns (pdfs, params), keyed like PDFDefinitions.create_bkg_pdfs.
    The core is written into each formula rather than shared, so members of
    both families can sit in one workspace without a name clash.
    """
    core = core_expr("@0")
    lo, hi = x.getMin(), x.getMax()
    t = f"((2*@0-{lo + hi!r})/{hi - lo!r})"
    # Chebyshev polynomials written out: cos(n*acos(t)) is NaN for |t| > 1
    cheb = [t, f"(2*{t}*{t}-1)", f"(4*{t}*{t}*{t}-3*{t})",
            f"(8*{t}*{t}*{t}*{t}-8*{t}*{t}+1)"]

    params, pdfs = {}, {}
    for n in range(1, MAX_ORDER + 1):
        params[f"dyexp_a{n}"] = RooRealVar(f"dyexp_a{n}_{tag}", f"dyexp_a{n}",
                                           0., *COEF_RANGE)
        params[f"dypow_a{n - 1}"] = RooRealVar(f"dypow_a{n - 1}_{tag}",
                                               f"dypow_a{n - 1}",
                                               0., *COEF_RANGE)
    for n in range(1, MAX_ORDER + 1):
        coefs = RooArgList(x)
        for i in range(1, n + 1):
            coefs.add(params[f"dyexp_a{i}"])
        arg = " + ".join(f"@{i}*{cheb[i - 1]}" for i in range(1, n + 1))
        pdfs[f"dyexp{n}"] = RooGenericPdf(f"dyexp{n}_{tag}", f"dyexp{n}",
                                          f"{core}*exp({arg})", coefs)

        coefs = RooArgList(x)
        for i in range(n):
            coefs.add(params[f"dypow_a{i}"])
        expo = " + ".join(f"@{i + 1}" + f"*{t}" * i for i in range(n))
        pdfs[f"dypow{n}"] = RooGenericPdf(
            f"dypow{n}_{tag}", f"dypow{n}",
            f"{core}*TMath::Power(@0/{MREF!r},{expo})", coefs)
    return pdfs, params
