"""
biasDiagnose.py -- why a cell of the biasStudy.py --mode single matrix
fails: for each (truth, fit function), one toy (or the truth's Asimov),
the S+B fit of that function to it, and the cell's pull and mu_fit
distributions from all its toy fits.

    python biasDiagnose.py -c ggHcat -b bdt0 -y Run3 --outdir BIAS_5bins
    python biasDiagnose.py -c ggHcat -b bdt0 -y Run3 --outdir BIAS_5bins \
        --truth exp2 --fit bwzredux3 --toy 3
    python biasDiagnose.py -c ggHcat -b bdt0 bdt1 -y Run3 \
        --outdir BIAS_5bins --asimov

Runs in the Combine container, from parametricfit/, after biasStudy.py
--mode single (it reads <outdir>/<tag>_single/: meta.json, the workspace,
the toys and the fit files).

Per cell, one plot <plotdir>/<cat>/bias_diag_<tag>_<r>[_asimov]/
<truth>__fit_<fit>.png with four panels:

  mass       the toy (points), the truth's background (red; the shape the
             toy was generated from), and the fit: S+B (blue) and its
             background (dashed blue)
  residual   per bin, the toy minus the S+B fit (points), the fitted
             background minus the truth's (blue line), and the injected
             signal (grey): where the fitted background misses the truth
             under the peak, mu_fit moves to make up for it
  pulls      the Bias of every toy fit of the cell (combine's), with its
             mean and median: a few huge pulls (mean >> median) are fits
             gone wrong; a shifted peak is a real bias
  mu_fit     combine's mu_fit of every toy, the r range of the fits marked:
             toys piled up at an edge had their minimum outside the range,
             and their pull is meaningless

and the numbers: mu_fit of this fit and of combine's fit of the same toy
(they should agree), the integral of (fitted - truth background) and of
the injected signal over 120-130 GeV.

The fit here is a RooFit fit of the channel's S+B pdf, the index frozen to
the fit function and the other functions' parameters frozen (as
biasStudy.py does), MH and the constrained nuisances fixed at nominal:
combine also profiles those, so mu_fit agrees with combine's only
approximately.

--asimov fits the truth's Asimov (expected S+B at mu_inj, nuisances
nominal) instead of a toy: the bias the shape difference alone gives,
without the statistical scatter. It also prints the matrix of these
expected Biases, (mu_fit - mu_inj) / sigma_mu, and writes it to
<outdir>/<tag>_single/diag_asimov_<tag>_<r>.json.
"""

import argparse
import glob
import math
import os

import ROOT

import biasStudy as bs
from fitlib import pdfDefinitions as pdefs
from fitlib import style

ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

WINDOW = (120.0, 130.0)     # the signal window the integrals are over
NSUB = 20                   # grid points per plot bin for the curves
TRUTH_RGB = (200, 30, 30)
FIT_RGB = (31, 90, 180)
SIG_RGB = (190, 190, 190)
EDGE_FRAC = 1e-3            # mu_fit this close to the range edge: at it


# ---------------------------------------------------------------------------
# the model
# ---------------------------------------------------------------------------

class Model:
    """The pieces of a text2workspace workspace the diagnosis needs.

    x        the observable
    add      the channel's extended S+B RooAddPdf
    terms    [(pdf, yield, is_background)] of `add`
    members  the RooMultiPdf's functions, in index order
    r        the signal strength
    fixed    the constrained nuisances and MH: fixed at nominal here
    """

    def __init__(self, w, x, add, terms, members, set_index, r, fixed,
                 names):
        self.w, self.x, self.add, self.terms = w, x, add, terms
        self.members, self.set_index, self.r = members, set_index, r
        self.fixed, self.names = fixed, names
        self.nset = ROOT.RooArgSet(x)
        params = bs.member_params(members, self.nset,
                                  keep={r.GetName()} | set(fixed))
        self.frozen = bs.others_frozen(params)
        self.vars = {v.GetName(): v for v in w.allVars()}
        self.start = {n: (v.getVal(), v.isConstant())
                      for n, v in self.vars.items()}

    def reset(self):
        """Every parameter back to the workspace's (pre-fit) value."""
        for n, (val, const) in self.start.items():
            v = self.vars[n]
            v.setConstant(False)
            v.setVal(val)
            v.setConstant(const)
        for n in self.fixed:
            if n in self.vars:
                self.vars[n].setConstant(True)

    def select(self, i):
        """Function i, the others' parameters frozen."""
        self.set_index(i)
        for n in self.frozen[i]:
            if n in self.vars:
                self.vars[n].setConstant(True)

    def density(self, xs, background=None):
        """Expected events per GeV at each x in xs: all terms, or only the
        background (True) or only the signal (False) ones."""
        out = [0.0] * len(xs)
        for pdf, coef, is_bkg in self.terms:
            if background is not None and is_bkg != background:
                continue
            n = coef.getVal()
            for k, v in enumerate(xs):
                self.x.setVal(v)
                out[k] += n * pdf.getVal(self.nset)
        return out


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


def load_model(t2w, tag):
    if not pdefs.load_combine(quiet=True):
        raise SystemExit("needs libHiggsAnalysisCombinedLimit: run inside "
                         "the Combine container")
    f = ROOT.TFile.Open(t2w)
    if not f or f.IsZombie():
        raise SystemExit(f"cannot open {t2w}")
    w = f.Get("w")
    ROOT.SetOwnership(f, False)         # keep the workspace alive
    multi = bs.find_multipdf(w, tag)
    members = [multi.getPdf(i) for i in range(multi.getNumPdfs())]
    obs = [v for v in (w.set("observables") or [])
           if v.InheritsFrom("RooRealVar")]
    if len(obs) != 1:
        raise SystemExit(f"{t2w}: expected one observable, found "
                         f"{[v.GetName() for v in obs]}")
    x = obs[0]
    model_s = w.pdf("model_s")
    add = find_add(model_s, multi) if model_s else None
    if add is None:
        raise SystemExit(f"{t2w}: no S+B RooAddPdf holding the multipdf "
                         "found under model_s")
    terms = [(p, c, bool(p.dependsOn(multi)))
             for p, c in zip(add.pdfList(), add.coefList())]
    index = w.cat(f"pdfindex_{tag}")
    # constant, so the minimizer never sees it: ROOT 6.22's RooMinimizerFcn
    # removes a floating category from its parameter list while iterating
    # over it, which throws for some orders of the parameters
    index.setConstant(True)
    globs = bs.model_set(w, "globalObservables")
    fixed = sorted(n for n in bs.model_set(w, "nuisances")
                   if f"{n}_In" in globs) + ["MH"]
    names = [bs.short_name(p.GetName(), tag) for p in members]
    return Model(w, x, add, terms, members, index.setIndex, w.var("r"),
                 fixed, names)


# ---------------------------------------------------------------------------
# data and fit
# ---------------------------------------------------------------------------

def read_toy(path, k, x):
    """Toy k (1-based, as combine counts) of a GenerateOnly --saveToys
    file, binned in x's binning."""
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise SystemExit(f"cannot open {path}")
    ds = f.Get(f"toys/toy_{k}")
    if not ds:
        raise SystemExit(f"{path}: no toys/toy_{k}")
    dh = ROOT.RooDataHist(f"toy_{k}_binned", "toy", ROOT.RooArgSet(x), ds)
    f.Close()
    return dh


def asimov(model):
    """The expected S+B of the current state, binned in x's binning."""
    return model.add.generateBinned(ROOT.RooArgSet(model.x),
                                    ROOT.RooFit.ExpectedData())


def fit(model, data):
    """S+B fit of the current function to data. Returns (status, mu,
    error on mu, min NLL)."""
    nll = model.add.createNLL(data, ROOT.RooFit.Extended(True),
                              ROOT.RooFit.Offset(True))
    m = ROOT.RooMinimizer(nll)
    m.setMinimizerType("Minuit2")
    m.setPrintLevel(-1)
    m.setPrintEvalErrors(-1)
    m.setStrategy(1)
    m.minimize("Minuit2", "Migrad")
    m.minimize("Minuit2", "Migrad")
    m.hesse()
    res = m.save()
    return res.status(), model.r.getVal(), model.r.getError(), res.minNll()


def hist_of(data, x, edges, name):
    """The data in the plot binning."""
    h = ROOT.TH1D(name, "", len(edges) - 1, edges[0], edges[-1])
    h.SetDirectory(0)
    for i in range(data.numEntries()):
        row = data.get(i)
        h.Fill(row.getRealValue(x.GetName()), data.weight())
    for b in range(1, h.GetNbinsX() + 1):
        h.SetBinError(b, math.sqrt(max(h.GetBinContent(b), 0.0)))
    return h


def grid(lo, hi, nbins, nsub=NSUB):
    """Midpoints of nsub sub-intervals of every plot bin."""
    step = (hi - lo) / float(nbins * nsub)
    return [lo + (k + 0.5) * step for k in range(nbins * nsub)], step


def per_bin(dens, nsub, step):
    """Events per plot bin from densities on the grid."""
    return [sum(dens[b * nsub:(b + 1) * nsub]) * step
            for b in range(len(dens) // nsub)]


def window_integral(xs, dens, step, lo=WINDOW[0], hi=WINDOW[1]):
    return sum(d for v, d in zip(xs, dens) if lo <= v < hi) * step


# ---------------------------------------------------------------------------
# what combine found for the cell
# ---------------------------------------------------------------------------

def cell_results(outdir, prefix, r_inj, rrange, requested):
    """Combine's fits of the cell: (pulls, every toy's mu_fit, summary)."""
    files = sorted(glob.glob(os.path.join(outdir, f"fit_{prefix}_s*.root")))
    pulls, _, _, drop = bs.collect_pulls(files, r_inj, rrange, None,
                                         requested)
    mus = []
    for path in files:
        mus += [t["best"] for t in bs.read_singles(path).values()
                if "best" in t]
    return pulls, mus, drop, len(files)


def combine_mu(outdir, prefix, seed, k):
    """Combine's mu_fit of toy k of job `seed`, or None."""
    path = os.path.join(outdir, f"fit_{prefix}_s{seed}.root")
    if not os.path.isfile(path):
        return None
    return bs.read_singles(path).get(k, {}).get("best")


# ---------------------------------------------------------------------------
# plot
# ---------------------------------------------------------------------------

def color(rgb):
    return ROOT.TColor.GetColor(*rgb)


def label(pad, right):
    """CMS Internal on the left and `right` on the right, above the frame."""
    t = ROOT.TLatex()
    t.SetNDC()
    y = 1 - 0.75 * pad.GetTopMargin()
    t.SetTextFont(42)
    t.SetTextSize(0.05)
    t.DrawLatex(pad.GetLeftMargin(), y, "#bf{CMS} #it{Internal}")
    t.SetTextAlign(31)
    t.SetTextSize(0.042)
    t.DrawLatex(1 - pad.GetRightMargin(), y, right)
    return t


def axis_sizes(h, size=0.05):
    for ax in (h.GetXaxis(), h.GetYaxis()):
        ax.SetTitleSize(size)
        ax.SetLabelSize(0.9 * size)
    h.GetYaxis().SetTitleOffset(1.4)


def draw_cell(d, outbase, pdf):
    """The four panels of one cell; `d` holds everything to draw."""
    style.setup_style()
    c = ROOT.TCanvas(f"c_{os.path.basename(outbase)}", "", 1300, 850)
    keep = []
    p_mass = ROOT.TPad("p_mass", "", 0.0, 0.32, 0.64, 1.0)
    p_res = ROOT.TPad("p_res", "", 0.0, 0.0, 0.64, 0.32)
    p_pull = ROOT.TPad("p_pull", "", 0.64, 0.5, 1.0, 1.0)
    p_mu = ROOT.TPad("p_mu", "", 0.64, 0.0, 1.0, 0.5)
    for p in (p_mass, p_res, p_pull, p_mu):
        p.SetLeftMargin(0.13)
        p.SetRightMargin(0.04)
        c.cd()
        p.Draw()
    p_mass.SetTopMargin(0.09)
    p_mass.SetBottomMargin(0.02)
    p_res.SetTopMargin(0.03)
    p_res.SetBottomMargin(0.30)
    for p in (p_pull, p_mu):
        p.SetTopMargin(0.09)
        p.SetBottomMargin(0.16)
        p.SetLeftMargin(0.16)

    edges, centers = d["edges"], d["centers"]
    width = edges[1] - edges[0]

    # --- mass
    p_mass.cd()
    h = d["data"]
    h.SetMarkerStyle(20)
    h.SetMarkerSize(0.6)
    h.SetLineColor(ROOT.kBlack)
    hi = max(h.GetMaximum(), max(d["truth_b"]))
    lo_y = min(h.GetBinContent(b) for b in range(1, h.GetNbinsX() + 1))
    h.SetMaximum(hi + 0.45 * (hi - lo_y))
    h.SetMinimum(max(0.0, lo_y - 0.15 * (hi - lo_y)))
    h.GetYaxis().SetTitle(f"events / {width:g} GeV")
    h.GetXaxis().SetLabelSize(0)
    h.GetYaxis().SetTitleSize(0.05)
    h.GetYaxis().SetLabelSize(0.045)
    h.GetYaxis().SetTitleOffset(1.25)
    h.Draw("E1")

    def curve(values, rgb, lstyle, name):
        g = ROOT.TGraph(len(d["xs"]))
        g.SetName(name)
        for k, (v, y) in enumerate(zip(d["xs"], values)):
            g.SetPoint(k, v, y * width)
        g.SetLineColor(color(rgb))
        g.SetLineStyle(lstyle)
        g.SetLineWidth(2)
        g.Draw("L same")
        keep.append(g)
        return g

    g_tb = curve(d["dens_truth_b"], TRUTH_RGB, 1, "g_tb")
    g_fsb = curve(d["dens_fit_sb"], FIT_RGB, 1, "g_fsb")
    g_fb = curve(d["dens_fit_b"], FIT_RGB, 2, "g_fb")
    h.Draw("E1 same")
    leg = ROOT.TLegend(0.52, 0.66, 0.95, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.042)
    leg.AddEntry(h, d["data_label"], "pe")
    leg.AddEntry(g_tb, f"truth background: {d['truth']}", "l")
    leg.AddEntry(g_fsb, f"fit S+B: {d['fit']}", "l")
    leg.AddEntry(g_fb, "fit, its background", "l")
    leg.Draw()
    keep.append(leg)
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextSize(0.036)
    x0, y0 = 0.50, 0.61
    for k, line in enumerate(d["text"]):
        lat.DrawLatex(x0, y0 - 0.048 * k, line)
    keep.append(lat)
    keep.append(label(p_mass, d["header"]))

    # --- residual
    p_res.cd()
    res = h.Clone("h_res")
    res.SetDirectory(0)
    for b in range(1, res.GetNbinsX() + 1):
        res.SetBinContent(b, h.GetBinContent(b) - d["fit_sb"][b - 1])
        res.SetBinError(b, h.GetBinError(b))
    sig = ROOT.TH1D("h_sig", "", len(edges) - 1, edges[0], edges[-1])
    db = ROOT.TH1D("h_db", "", len(edges) - 1, edges[0], edges[-1])
    for hh in (sig, db):
        hh.SetDirectory(0)
    for b in range(1, sig.GetNbinsX() + 1):
        sig.SetBinContent(b, d["truth_s"][b - 1])
        db.SetBinContent(b, d["fit_b"][b - 1] - d["truth_b"][b - 1])
    vals = ([res.GetBinContent(b) + s * res.GetBinError(b)
             for b in range(1, res.GetNbinsX() + 1) for s in (-1, 1)]
            + [db.GetBinContent(b) for b in range(1, db.GetNbinsX() + 1)]
            + [sig.GetMaximum()])
    span = max(abs(v) for v in vals) * 1.15 or 1.0
    res.SetMaximum(span)
    res.SetMinimum(-span)
    res.GetYaxis().SetTitle("events")
    res.GetXaxis().SetTitle("m_{#mu#mu} [GeV]")
    for ax, size in ((res.GetXaxis(), 0.11), (res.GetYaxis(), 0.10)):
        ax.SetTitleSize(size)
        ax.SetLabelSize(0.09)
    res.GetYaxis().SetTitleOffset(0.55)
    res.GetYaxis().SetNdivisions(505)
    res.Draw("E1")
    sig.SetFillColor(color(SIG_RGB))
    sig.SetLineColor(color(SIG_RGB))
    sig.Draw("HIST same")
    db.SetLineColor(color(FIT_RGB))
    db.SetLineWidth(2)
    db.Draw("HIST same")
    zero = ROOT.TLine(edges[0], 0, edges[-1], 0)
    zero.SetLineStyle(3)
    zero.Draw()
    res.Draw("E1 same")
    leg2 = ROOT.TLegend(0.15, 0.78, 0.95, 0.96)
    leg2.SetNColumns(3)
    leg2.SetBorderSize(0)
    leg2.SetFillStyle(0)
    leg2.SetTextFont(42)
    leg2.SetTextSize(0.075)
    leg2.AddEntry(res, "data - fit S+B", "pe")
    leg2.AddEntry(db, "fit B - truth B", "l")
    leg2.AddEntry(sig, f"signal, #mu = {d['r_inj']:g}", "f")
    leg2.Draw()
    keep += [res, sig, db, zero, leg2]

    # --- pulls of the cell
    p_pull.cd()
    pulls = d["pulls"]
    if pulls:
        lo = min(-5.0, math.floor(min(pulls)))
        hi = max(5.0, math.ceil(max(pulls)))
        hp = ROOT.TH1D("h_pulls", ";Bias = (#mu_{fit} - #mu_{inj}) / "
                       "#sigma_{#mu};toys", 50, lo, hi)
        hp.SetDirectory(0)
        for p in pulls:
            hp.Fill(p)
        hp.SetFillColor(ROOT.TColor.GetColor(*bs.PULL_FILL_RGB))
        hp.SetLineColor(ROOT.TColor.GetColor(*bs.PULL_LINE_RGB))
        hp.SetMaximum(1.45 * hp.GetMaximum())
        axis_sizes(hp)
        hp.Draw("HIST")
        keep.append(hp)
        n = len(pulls)
        mean = sum(pulls) / n
        lines = [f"{d['truth']} toys, fit {d['fit']}",
                 f"valid fits: {n} / {n + sum(d['drop'].values())}",
                 "mean %+.2f, median %+.2f" % (mean, bs.median(pulls))]
        dropped = ", ".join(f"{k} {v}" for k, v in d["drop"].items() if v)
        if dropped:
            lines.append(f"dropped: {dropped}")
        t = ROOT.TLatex()
        t.SetNDC()
        t.SetTextFont(42)
        t.SetTextSize(0.045)
        for k, line in enumerate(lines):
            t.DrawLatex(0.20, 0.85 - 0.06 * k, line)
        keep.append(t)
    else:
        t = ROOT.TLatex(0.2, 0.5, "no combine fits for this cell")
        t.SetNDC()
        t.Draw()
        keep.append(t)

    # --- mu_fit of the cell
    p_mu.cd()
    lo, hi = d["rrange"]
    hm = ROOT.TH1D("h_mu", ";#mu_{fit} (combine);toys", 60, lo, hi)
    hm.SetDirectory(0)
    for v in d["mus"]:
        hm.Fill(min(max(v, lo + 1e-9), hi - 1e-9))
    hm.SetFillColor(ROOT.TColor.GetColor(*bs.PULL_FILL_RGB))
    hm.SetLineColor(ROOT.TColor.GetColor(*bs.PULL_LINE_RGB))
    hm.SetMaximum(1.45 * max(hm.GetMaximum(), 1))
    axis_sizes(hm)
    hm.Draw("HIST")
    l_inj = ROOT.TLine(d["r_inj"], 0, d["r_inj"], hm.GetMaximum())
    l_inj.SetLineColor(color(TRUTH_RGB))
    l_inj.SetLineStyle(2)
    l_inj.SetLineWidth(2)
    l_inj.Draw()
    edge = EDGE_FRAC * (hi - lo)
    at_edge = sum(1 for v in d["mus"] if v < lo + edge or v > hi - edge)
    t = ROOT.TLatex()
    t.SetNDC()
    t.SetTextFont(42)
    t.SetTextSize(0.045)
    t.DrawLatex(0.20, 0.85, f"fit range [{lo:.3g}, {hi:.3g}], "
                            f"#mu_{{inj}} = {d['r_inj']:g} (red)")
    t.DrawLatex(0.20, 0.79, f"at a range edge: {at_edge} / {len(d['mus'])}")
    if d["mus"]:
        t.DrawLatex(0.20, 0.73, "median #mu_{fit} = %.3g"
                    % bs.median(d["mus"]))
    keep += [hm, l_inj, t]

    c.cd()
    bs.save(c, outbase, pdf)
    return at_edge


# ---------------------------------------------------------------------------
# main
# ---------------------------------------------------------------------------

def diagnose_tag(tag, args):
    outdir = os.path.abspath(os.path.join(args.outdir, tag + "_single"))
    meta_path = os.path.join(outdir, "meta.json")
    if not os.path.isfile(meta_path):
        print(f"\n=== {tag} === no {meta_path}; run biasStudy.py "
              "--mode single first")
        return None
    meta = bs.read_json(meta_path)
    key = bs.rlabel(args.inject)
    if key not in meta.get("rrange", {}):
        print(f"\n=== {tag} === no toys at r = {args.inject:g}")
        return None
    rrange = meta["rrange"][key]
    model = load_model(meta["workspace"], tag)
    names = model.names
    fits = [j for j in meta.get("fits", range(len(names)))
            if not args.fit or names[j] in args.fit]
    truths = [i for i in range(len(names))
              if not args.truth or names[i] in args.truth]
    x = model.x
    if args.fit_bins:
        x.setBins(args.fit_bins)
    lo, hi = x.getMin(), x.getMax()
    nplot = int(round((hi - lo) / args.plot_bin_width))
    edges = [lo + k * (hi - lo) / nplot for k in range(nplot + 1)]
    centers = [0.5 * (a + b) for a, b in zip(edges[:-1], edges[1:])]
    xs, step = grid(lo, hi, nplot)
    win_xs, win_step = grid(WINDOW[0], WINDOW[1], 1, nsub=200)

    category = tag.split("_")[0]
    sub = f"bias_diag_{tag}_{key}" + ("_asimov" if args.asimov else "")
    plotdir = os.path.join(args.plotdir, category, sub)
    os.makedirs(plotdir, exist_ok=True)
    print(f"\n=== {tag}, r = {args.inject:g} ===")
    print(f"  workspace : {meta['workspace']}")
    print(f"  fixed     : {', '.join(model.fixed)} (at nominal)")
    print(f"  plots     : {plotdir}/")
    table = {}
    for i in truths:
        truth = names[i]
        seed = args.seed + 1000 * i + args.job
        toys = os.path.join(outdir, f"toys_gen{i}_{key}_s{seed}.root")
        # the truth: its background and signal, at the pre-fit values
        model.reset()
        model.select(i)
        model.r.setVal(args.inject)
        truth_b = model.density(xs, True)
        truth_s = model.density(xs, False)
        s_win = window_integral(xs, truth_s, step)
        if args.asimov:
            data = asimov(model)
            data_label = f"Asimov of {truth}, #mu = {args.inject:g}"
        else:
            if not os.path.isfile(toys):
                print(f"  {truth}: no {toys}; skipped")
                continue
            data = read_toy(toys, args.toy, x)
            data_label = f"{truth} toy {args.toy} (job s{seed})"
        h = hist_of(data, x, edges, f"h_data_{i}")
        for j in fits:
            fitname = names[j]
            prefix = f"gen{i}_fit{j}_{key}"
            model.reset()
            model.select(j)
            model.r.setRange(*rrange)
            model.r.setVal(args.inject)
            status, mu, err, _ = fit(model, data)
            fit_sb = model.density(xs)
            fit_b = model.density(xs, True)
            db_win = (window_integral(xs, fit_b, step)
                      - window_integral(xs, truth_b, step))
            cmu = None if args.asimov else combine_mu(outdir, prefix, seed,
                                                       args.toy)
            pulls, mus, drop, nfiles = cell_results(
                outdir, prefix, args.inject, rrange, meta.get("ntoys", {}))
            chi2, ndf = 0.0, 0
            fsb_bins = per_bin(fit_sb, NSUB, step)
            for b, e in enumerate(fsb_bins):
                if e > 0:
                    chi2 += (h.GetBinContent(b + 1) - e) ** 2 / e
                    ndf += 1
            ndf -= 1 + sum(1 for v in model.vars.values()
                           if not v.isConstant()
                           and v.GetName() != x.GetName())
            win = "%g-%g GeV" % WINDOW
            text = [
                "#mu_{fit} = %.3f #pm %.3f (status %d)" % (mu, err, status),
                ("" if cmu is None else
                 "combine, same toy: #mu_{fit} = %.3f" % cmu),
                "B_{fit} - B_{truth} in %s: %+.0f events" % (win, db_win),
                "S(#mu_{inj}) in %s: %.0f events" % (win, s_win),
                "#chi^{2} / ndf = %.0f / %d" % (chi2, max(ndf, 1)),
            ]
            text = [t for t in text if t]
            d = {"edges": edges, "centers": centers, "xs": xs,
                 "data": h.Clone(f"h_{i}_{j}"), "data_label": data_label,
                 "truth": truth, "fit": fitname, "r_inj": args.inject,
                 "dens_truth_b": truth_b, "dens_fit_sb": fit_sb,
                 "dens_fit_b": fit_b,
                 "truth_b": per_bin(truth_b, NSUB, step),
                 "truth_s": per_bin(truth_s, NSUB, step),
                 "fit_b": per_bin(fit_b, NSUB, step), "fit_sb": fsb_bins,
                 "pulls": pulls, "mus": mus, "drop": drop,
                 "rrange": rrange, "text": text,
                 "header": f"{tag}, #mu_{{inj}} = {args.inject:g}"}
            at_edge = draw_cell(d, os.path.join(
                plotdir, f"{truth}__fit_{fitname}"), args.pdf)
            exp_bias = (mu - args.inject) / err if err > 0 else float("nan")
            table.setdefault(fitname, {})[truth] = {
                "mu": mu, "mu_err": err, "status": status,
                "bias": exp_bias, "delta_b_window": db_win,
                "s_window": s_win, "combine_mu": cmu,
                "chi2": chi2, "ndf": ndf,
                "combine_mean_pull": (sum(pulls) / len(pulls)
                                      if pulls else None),
                "combine_median_pull": bs.median(pulls) if pulls else None,
                "combine_at_r_edge": at_edge, "combine_toys": len(mus)}
            print(f"  {truth:>10} <- {fitname:<10} mu {mu:+8.3f} +- "
                  f"{err:6.3f}"
                  + ("" if cmu is None else f" (combine {cmu:+8.3f})")
                  + f"  dB(window) {db_win:+9.0f}  S {s_win:7.0f}"
                  + f"  combine pulls: mean "
                  + ("%+7.2f median %+6.2f" % (
                      sum(pulls) / len(pulls), bs.median(pulls))
                     if pulls else "   -")
                  + f"  at r edge {at_edge}/{len(mus)}")
    if args.asimov and table:
        print_expected(table, names, truths, args.threshold)
        bs.write_json(os.path.join(outdir, f"diag_asimov_{tag}_{key}.json"),
                      {"tag": tag, "r_inj": args.inject, "r_range": rrange,
                       "fixed_at_nominal": model.fixed, "cells": table})
    return table


def print_expected(table, names, truths, threshold):
    """The Asimov matrix: (mu_fit - mu_inj) / sigma_mu per fit and truth."""
    cols = [names[i] for i in truths]
    print("\n  expected Bias from the Asimov of each truth "
          "(rows: fit, columns: truth; * = above the threshold)")
    print("  " + " " * 12 + "".join(f"{c:>11}" for c in cols))
    for fitname, row in table.items():
        cells = []
        for c in cols:
            b = row.get(c, {}).get("bias")
            if b is None or b != b:
                cells.append(f"{'-':>11}")
            else:
                cells.append(f"{b:+10.2f}" + ("*" if abs(b) >= threshold
                                              else " "))
        print(f"  {fitname:<12}" + "".join(cells))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    ap.add_argument("-c", "--category", default="ggHcat")
    ap.add_argument("-b", "--bins", nargs="+", required=True)
    ap.add_argument("-y", "--year", default="Run3")
    ap.add_argument("--outdir", default="BIAS",
                    help="biasStudy.py's --outdir (default BIAS)")
    ap.add_argument("--plotdir",
                    default=os.path.expanduser("~/public_html/HmumuFits/bias"))
    ap.add_argument("--inject", type=float, default=1.0)
    ap.add_argument("--truth", nargs="+", default=None, metavar="NAME",
                    help="only these truths (default all)")
    ap.add_argument("--fit", nargs="+", default=None, metavar="NAME",
                    help="only these fit functions (default all)")
    ap.add_argument("--toy", type=int, default=1,
                    help="which toy of the job to fit (1-based, default 1)")
    ap.add_argument("--job", type=int, default=0,
                    help="which job's toys: seed = --seed + 1000*truth + job "
                         "(default 0, the first)")
    ap.add_argument("--seed", type=int, default=12345,
                    help="biasStudy.py's --seed (default 12345)")
    ap.add_argument("--asimov", action="store_true",
                    help="fit the truth's Asimov instead of a toy, and print "
                         "the matrix of expected Biases")
    ap.add_argument("--fit-bins", type=int, default=None,
                    help="bins of the observable for the fit (default: the "
                         "workspace's)")
    ap.add_argument("--plot-bin-width", type=float, default=1.0,
                    help="GeV per plot bin (default 1)")
    ap.add_argument("--threshold", type=float, default=bs.PULL_THRESHOLD)
    ap.add_argument("--pdf", action="store_true")
    args = ap.parse_args()
    for b in args.bins:
        diagnose_tag(f"{args.category}_{b}_{args.year}", args)


if __name__ == "__main__":
    main()
