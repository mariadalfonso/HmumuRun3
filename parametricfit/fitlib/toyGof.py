"""
toyGof.py -- the goodness of fit of a bias-study toy fit: the S+B
expectation at the fit's best point against the toy, as a binned chi2
with fitlib/gof.py's bin merging and its p-value (the same test bkgFit.py
applies to the data, over the full window: the toy fit sees every bin).

biasStudy.py --mode single runs each toy fit with Combine's
--trackParameters on every floating parameter of the fit function, so
the best-fit row of the fit file holds the point; here the workspace is
set to it, toy by toy, and the expectation compared with the toy. The
workspace has RooMultiPdfs, so this needs the Combine library.

    ndf = merged bins - (free parameters + 1 for r)

where the free parameters are the fit function's shape parameters and
the background normalisation; the constrained nuisances (the muon scale
and resolution) float in the fit too but are not counted, their
constraint holds them near nominal.
"""

import math

import ROOT

from fitlib.gof import CHI2_MIN_BIN_SIGNIFICANCE, GOF_THRESHOLD, \
    merge_chi2_bins

BIN_GEV = 1.0           # chi2 binning, as bkgFit.py's
NSUB = 5                # expectation: midpoints of NSUB slices per bin
TRACKED = "trackedParam_"
Q_BEST = -1.0           # quantileExpected of the best-fit row


def find_add(pdf, multi):
    """The extended RooAddPdf with a term that holds `multi`, under pdf
    (a RooSimultaneous of one channel, a RooProdPdf with the constraints,
    or the RooAddPdf itself)."""
    if pdf.InheritsFrom("RooAddPdf"):
        if (pdf.coefList().getSize() == pdf.pdfList().getSize()
                and any(p.dependsOn(multi) for p in pdf.pdfList())):
            return pdf
    if pdf.InheritsFrom("RooSimultaneous"):
        for label in category_labels(pdf.indexCat()):
            sub = pdf.getPdf(label)
            if sub:
                found = find_add(sub, multi)
                if found:
                    return found
    if pdf.InheritsFrom("RooProdPdf"):
        for p in pdf.pdfList():
            if p.dependsOn(multi):
                found = find_add(p, multi)
                if found:
                    return found
    return None


def category_labels(cat):
    """The state names of a RooAbsCategory."""
    try:
        return [str(s.first) for s in cat.states()]    # ROOT >= 6.22
    except AttributeError:
        it = cat.typeIterator()
        out, t = [], it.Next()
        while t:
            out.append(t.GetName())
            t = it.Next()
        return out


class Channel:
    """The S+B expectation of one channel per chi2 bin.

    x      the observable
    terms  [(pdf, yield)] of the channel's extended S+B RooAddPdf
    vars   every RooRealVar of the workspace, by name
    edges  the chi2 bin edges, BIN_GEV wide across x's range
    """

    def __init__(self, x, add, variables):
        self.x = x
        self.nset = ROOT.RooArgSet(x)
        self.terms = list(zip(add.pdfList(), add.coefList()))
        self.vars = variables
        lo, hi = x.getMin(), x.getMax()
        n = max(1, int(round((hi - lo) / BIN_GEV)))
        self.edges = [lo + k * (hi - lo) / n for k in range(n + 1)]
        self.step = (hi - lo) / float(n * NSUB)
        self.grid = [lo + (k + 0.5) * self.step for k in range(n * NSUB)]

    def set_point(self, values):
        """Set the named parameters (those the workspace has)."""
        for name, val in values.items():
            v = self.vars.get(name)
            if v is not None:
                v.setVal(val)

    def expected(self):
        """Expected events per chi2 bin at the current point."""
        dens = [0.0] * len(self.grid)
        for pdf, coef in self.terms:
            n = coef.getVal()
            for k, v in enumerate(self.grid):
                self.x.setVal(v)
                dens[k] += n * pdf.getVal(self.nset)
        return [sum(dens[b * NSUB:(b + 1) * NSUB]) * self.step
                for b in range(len(self.edges) - 1)]

    def binned(self, data, name="h_toy"):
        """The data (a RooDataSet or RooDataHist) in the chi2 bins, with
        sqrt(N) errors."""
        h = ROOT.TH1D(name, "", len(self.edges) - 1, self.edges[0],
                      self.edges[-1])
        h.SetDirectory(0)
        xname = self.x.GetName()
        for i in range(data.numEntries()):
            row = data.get(i)
            h.Fill(row.getRealValue(xname), data.weight())
        for b in range(1, h.GetNbinsX() + 1):
            h.SetBinError(b, math.sqrt(max(h.GetBinContent(b), 0.0)))
        return h


def chi2_pvalue(h_data, expected, npar, threshold=None):
    """(chi2, ndf, p) of the data against the expectation per bin, bins
    merged as bkgFit.py merges them (no blinded window: the toy fit sees
    every bin)."""
    threshold = (CHI2_MIN_BIN_SIGNIFICANCE if threshold is None
                 else threshold)
    h_exp = h_data.Clone(h_data.GetName() + "_exp")
    h_exp.SetDirectory(0)
    for b, e in enumerate(expected):
        h_exp.SetBinContent(b + 1, e)
    chi2, nbins = merge_chi2_bins(h_data, h_exp, 0.0, 0.0, threshold)
    ndf = nbins - npar
    p = ROOT.TMath.Prob(chi2, ndf) if ndf > 0 else 0.0
    return chi2, ndf, p


def read_tracked(path, names):
    """{iToy: {name: value, "r": r}} at the best fit of each toy of a
    Combine fit file written with --trackParameters `names`. Returns None
    if the file has none of those branches (a fit made without them)."""
    out = {}
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        return {}
    t = f.Get("limit")
    if not t:
        f.Close()
        return {}
    have = [n for n in names if t.GetBranch(TRACKED + n)]
    if names and not have:
        f.Close()
        return None
    for _ in t:
        if abs(float(t.quantileExpected) - Q_BEST) > 1e-3:
            continue
        row = {n: float(getattr(t, TRACKED + n)) for n in have}
        row["r"] = float(t.r)
        out[int(t.iToy)] = row
    f.Close()
    return out


def toys_gof(channel, toys_path, rows, npar):
    """{iToy: (chi2, ndf, p)} for the toys of toys_path at the points
    `rows` (read_tracked)."""
    out = {}
    f = ROOT.TFile.Open(toys_path)
    if not f or f.IsZombie():
        return out
    for k, values in sorted(rows.items()):
        ds = f.Get(f"toys/toy_{k}")
        if not ds:
            continue
        channel.set_point(values)
        h = channel.binned(ds, f"h_toy_{k}")
        out[k] = chi2_pvalue(h, channel.expected(), npar)
    f.Close()
    return out


def summarize_p(pvalues, threshold=GOF_THRESHOLD):
    """The p-values of one cell's toys: their number, median and the
    fraction above the goodness-of-fit threshold."""
    ps = sorted(pvalues)
    n = len(ps)
    if not n:
        return {"n": 0}
    med = ps[n // 2] if n % 2 else 0.5 * (ps[n // 2 - 1] + ps[n // 2])
    return {"n": n, "median_p": med, "threshold": threshold,
            "frac_pass": sum(1 for p in ps if p > threshold) / float(n)}


# ---------------------------------------------------------------------------
# the lowest r the fits may reach
# ---------------------------------------------------------------------------

FLOOR_STEP_GEV = 0.05   # grid the S+B >= 0 limit is found on


def channel_adds(pdf):
    """Every channel's extended S+B RooAddPdf under pdf (model_s: a
    RooSimultaneous of channels, or one channel with its constraints)."""
    if (pdf.InheritsFrom("RooAddPdf")
            and pdf.coefList().getSize() == pdf.pdfList().getSize()):
        return [pdf]
    out = []
    if pdf.InheritsFrom("RooSimultaneous"):
        for label in category_labels(pdf.indexCat()):
            sub = pdf.getPdf(label)
            if sub:
                out += channel_adds(sub)
    elif pdf.InheritsFrom("RooProdPdf"):
        for p in pdf.pdfList():
            out += channel_adds(p)
    return out


def r_floor(w, poi="r"):
    """The lowest r at which S + B stays >= 0 everywhere: -min of B(x) /
    S(x) over the channels and x, S the signal at r = 1 (the terms whose
    yield depends on r) and B the rest, every RooMultiPdf member in turn,
    all at the workspace's (pre-fit) values. None if there is no signal.
    Below it a fit meets a negative expectation: the likelihood is not
    defined there."""
    r = w.var(poi)
    model = w.pdf("model_s")
    obs = w.set("observables")
    if not r or not model or not obs:
        return None
    saved = r.getVal()
    r.setVal(1.0)
    worst = None
    for add in channel_adds(model):
        xs = [v for v in add.getObservables(obs)
              if v.InheritsFrom("RooRealVar")]
        if len(xs) != 1:
            continue
        x = xs[0]
        nset = ROOT.RooArgSet(x)
        sig = [(p, c) for p, c in zip(add.pdfList(), add.coefList())
               if c.dependsOn(r)]
        bkg = [(p, c) for p, c in zip(add.pdfList(), add.coefList())
               if not c.dependsOn(r)]
        if not sig or not bkg:
            continue
        n = int(round((x.getMax() - x.getMin()) / FLOOR_STEP_GEV))
        grid = [x.getMin() + (k + 0.5) * (x.getMax() - x.getMin()) / n
                for k in range(n)]
        saved_x = x.getVal()
        s_of = []
        for v in grid:
            x.setVal(v)
            s_of.append(sum(c.getVal() * p.getVal(nset) for p, c in sig))
        smax = max(s_of)
        # each background term as itself, or each member of a RooMultiPdf
        choices = [[p] if not p.InheritsFrom("RooMultiPdf") else
                   [p.getPdf(i) for i in range(p.getNumPdfs())]
                   for p, _ in bkg]
        combos = [[]]
        for ch in choices:
            combos = [c + [m] for c in combos for m in ch]
        for combo in combos:
            for v, s in zip(grid, s_of):
                if s <= 1e-6 * smax:
                    continue
                x.setVal(v)
                b = sum(c.getVal() * m.getVal(nset)
                        for m, (_, c) in zip(combo, bkg))
                worst = b / s if worst is None else min(worst, b / s)
        x.setVal(saved_x)
    r.setVal(saved)
    return None if worst is None else -worst
