"""
corePdf.py -- the Run 2 note's core-PDF (Sec. 6.4.2), the ggH background
model of the Run 2 result.

One core function is shared by every BDT bin of a category, with the SAME
parameters in every bin; each bin multiplies it by its own shape modifier
(SMF), a polynomial with its own parameters:

    B_bin(m) = core_k(m; theta)  x  SMF_bin(m; c_bin)      (theta shared)

The core is itself an envelope of three functions with three parameters
each, as in the note, so the bins also share ONE discrete index:

    core      the Run 2 function           members here
    bwzredux  BWZRedux (5.3)               bwzredux  (a, b, c)
    fewzbern3 FEWZ x Bernstein (5.5)       fewzbern3 (c1, c2, c3)
    exp2      S-exponential (5.7)          exp2      (two slopes, one fraction)

With the same number of parameters in every member, the 0.5-per-parameter
penalty does not change the choice, as in the note.

The SMF is a Bernstein polynomial over the fit window, its first
coefficient fixed to 1 (the norm floats separately): order n has n free
parameters and spans the same functions as the note's order-n polynomial.
The note used order 3 in ggH-cat0 and cat1 and order 2 in the others, and
chose the orders with the bias study.

Names, so that Combine shares the core across the bins of the combined
card (text2workspace merges same-named objects of the bin workspaces):

    core parameters, core PDFs   <name>_<cat>_core_<year>   e.g. bwzredux_a_ggHcat_core_Run3
    the discrete index           pdfindex_<cat>_core_<year>
    SMF parameters, per bin      smf_<core>_c<k>_<tag>
    bin members                  core_<core>_<tag>          = core x SMF
    bin multipdf                 multipdf_<tag>_bkg          (as for the envelope)

Each (core, bin) pair has its own SMF parameters, so every member of a
bin's multipdf starts at its own fitted values.
"""

import ROOT

from . import pdfDefinitions as pdefs

CORES = ["bwzredux", "fewzbern3", "exp2"]
SMF_ORDER = 2                 # default; the note: 3 in cat0/cat1, 2 elsewhere
SMF_INIT, SMF_RANGE, SMF_STEP = 1.0, (-10.0, 10.0), 0.05


def core_tag(category, year):
    """ggHcat_core_Run3: the tag of the shared core's objects"""
    return f"{category}_core_{year}"


def index_name(category, year):
    """the discrete index the bins share"""
    return f"pdfindex_{core_tag(category, year)}"


def member_name(core, tag):
    return f"core_{core}_{tag}"


def create_core(x, category, year, core):
    """The shared core function `core`, its parameters named for the
    category. Build it ONCE per process and use it in every bin.
    Returns (pdf, keep) as pdfDefinitions.create_bkg_pdf."""
    if core not in pdefs.BKG_NPAR:
        raise KeyError(f"unknown core '{core}'")
    pdf, keep = pdefs.create_bkg_pdf(x, core_tag(category, year), core)
    if pdf.InheritsFrom("RooAddPdf"):
        # inside the product RooProdPdf evaluates the core without a
        # normalisation set, and a RooAddPdf then reads its fraction as one
        # of UNnormalised terms; pin it to the terms normalised over x, as
        # when the core is fitted alone (1.4% off for exp2 otherwise)
        pdf.fixCoefNormalization(ROOT.RooArgSet(x))
    return pdf, keep


def create_smf(x, core, tag, order):
    """Bernstein SMF of `order` for (core, bin): c0 = 1, c1..cn free.
    Returns (pdf, [free parameters], keep)."""
    c0 = ROOT.RooRealVar(f"smf_{core}_c0_{tag}", "c0", 1.0)
    c0.setConstant(True)
    coefs = ROOT.RooArgList(c0)
    free = []
    for k in range(1, order + 1):
        v = ROOT.RooRealVar(f"smf_{core}_c{k}_{tag}", f"c{k}",
                            SMF_INIT, *SMF_RANGE)
        v.setError(SMF_STEP)
        coefs.add(v)
        free.append(v)
    if order == 0:
        pdf = ROOT.RooUniform(f"smf_{core}_{tag}", "smf", ROOT.RooArgSet(x))
    else:
        pdf = ROOT.RooBernstein(f"smf_{core}_{tag}", "smf", x, coefs)
    return pdf, free, [c0, coefs]


def create_member(x, core_pdf, core, tag, order):
    """core x SMF for one bin. RooProdPdf normalises the product over x
    numerically, as the two factors share the observable.
    Returns (pdf, smf parameters, keep)."""
    smf, free, keep = create_smf(x, core, tag, order)
    prod = ROOT.RooProdPdf(member_name(core, tag), f"{core} x SMF",
                           ROOT.RooArgList(core_pdf, smf))
    return prod, free, keep + [smf]



def npar(core, order):
    """free shape parameters of one bin's member (the core counted in full)"""
    return pdefs.BKG_NPAR[core] + order
