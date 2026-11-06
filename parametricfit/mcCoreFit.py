"""
mcCoreFit.py -- fit the frozen background cores to the Drell-Yan MC.

For each BDT bin, the weighted Drell-Yan MC (BKG_DY_TAGS in
fitlib/prepareFits.py: QCD DY, plus EWK Z in VBF) is fitted over the full
window, unblinded since it holds no signal, with three smooth shapes
(fitlib/mcCore.py):

    mcexp   2Exp
    mcpow   2Pow
    mcdy    LO DY core x (1 + Cheb3)

Their fitted values are written as numbers into the mcexp<n> / mcpow<n> /
mcdy<n> background candidates, which bkgFit.py then fits on data with only
the correction exp(a1*T1 + ... + an*Tn) floating. This is the flashgg ggH
approach (Background/test/fitDYshapes.cpp); only the SHAPE of the MC is
used, never its normalisation.

Outputs, in WS_<year>/ unless --wsdir:

    MCcore_<tag>.json   frozen core values, fit quality, MC provenance

and a plot per bin in --plotdir. A core is frozen only if

  * its fit converged (status 1, a forced-positive covariance, is kept:
    only the central values are frozen), and
  * it passes a goodness-of-fit test against the MC: chi2 with bins merged
    to ~9 effective entries each (fitlib/gof.py), p > --gof-pmin (0.01,
    as bkgFit.py on data). Too little MC to leave a degree of freedom
    fails too, so a bin with sparse MC simply gets no MC cores.

A core left out has no candidates in bkgFit.py.

Usage (conda, ROOT >= 6.28), before bkgFit.py:

    python mcCoreFit.py -c ggHcat -b bdt0 bdt1 bdt2 -y Run3
    python mcCoreFit.py -c ggHcat -b bdt0 bdt1 bdt2 -y Run3 --all-bkg-mc

By default tt, tW and VV are left to the data-side correction.
--all-bkg-mc fits the total background MC instead (BKG_MC_TAGS), as
flashgg does; it writes the same MCcore_<tag>.json, which records which.
"""

import argparse
import json
import os
import random
from datetime import datetime, timezone

import ROOT

from fitlib import config, dyCore, mcCore, style
from fitlib.gof import (CHI2_MIN_BIN_SIGNIFICANCE, GOF_THRESHOLD,
                        merge_chi2_bins)

ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

CORE_RGB = {"mcexp": (42, 120, 214), "mcpow": (235, 104, 52),
            "mcdy": (27, 175, 122)}
PLOT_BINS_PER_GEV = 1          # coarser than the fit, for a readable plot
MC_FILL_RGB = (214, 214, 214)  # MC as a filled histogram, not data points
MC_LINE_RGB = (140, 140, 140)  # its outline, so the step shape reads
MC_BAND_STYLE = 3354           # hatched MC stat. uncertainty
LO_RGB = (0, 0, 0)             # the bare LO DY core, before any fit
LO_SUBSTEPS = 20               # midpoint rule per plot bin for the LO core
RANDOM_SPREAD = 0.05           # retry start: start +- this x the range


def rank(fit):
    """Better fit first: status 0, then status 1, then the lower NLL."""
    return (fit["status"] == 0, fit["status"] <= 1, -fit["nll"])


def fit_core(x, data, model, tag, tries, rng):
    """Fit one core to the MC; retry from random starts until status 0.

    Status 1 (Minuit forced the covariance positive definite) is kept if no
    try does better: the 2Exp and 2Pow components can be near-degenerate,
    which spoils the errors but not the minimum, and only the central
    values are frozen. Returns a dict with the values (and the live pdf,
    for the plot).
    """
    pdf, pars = mcCore.make_fit_pdf(x, model, tag)
    _, starts, ranges = mcCore.CORES[model]
    opts = [ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.Save(True),
            ROOT.RooFit.SumW2Error(True), ROOT.RooFit.PrintLevel(-1),
            ROOT.RooFit.PrintEvalErrors(-1)]
    best = None
    for itry in range(tries):
        for v, s, (lo, hi) in zip(pars, starts, ranges):
            if itry == 0:
                v.setVal(s)
            else:
                d = RANDOM_SPREAD * (hi - lo) * rng.uniform(-1, 1)
                v.setVal(min(max(s + d, lo), hi))
        strat = ROOT.RooFit.Strategy(1 if itry == 0 else 2)
        pdf.fitTo(data, strat, *opts)
        res = pdf.fitTo(data, strat, *opts)
        cand = {"status": res.status(), "cov_qual": res.covQual(),
                "nll": res.minNll(), "edm": res.edm(), "try": itry,
                "params": {v.GetTitle(): {"value": v.getVal(),
                                          "error": v.getError()}
                           for v in pars}}
        if best is None or rank(cand) > rank(best):
            best = cand
        if cand["status"] == 0:
            break
    for v in pars:                       # leave the best values in place
        v.setVal(best["params"][v.GetTitle()]["value"])
    best.update(pdf=pdf, keep=pars)
    return best


def gof_vs_mc(x, h, pdf, npar):
    """Goodness of fit of a core to the MC it was fitted to.

    Pearson chi2 with the MC's own sum-w^2 errors, after merging
    neighbouring bins until each has content/sqrt(sum w^2) >= 3, i.e.
    about 9 effective entries (fitlib/gof.py, as bkgFit.py on data). The
    expectation is normalised to the MC total, so ndf = bins - 1 - npar.
    A bin too sparse to leave any degree of freedom gets p = 0: a core
    that cannot be tested is not frozen.
    """
    nset = ROOT.RooArgSet(x)
    total = h.Integral(1, h.GetNbinsX())
    h_pdf = h.Clone(f"hpdf_{pdf.GetName()}")
    h_pdf.SetDirectory(0)
    h_pdf.Reset()
    for i in range(1, h.GetNbinsX() + 1):
        x.setVal(h.GetBinCenter(i))
        h_pdf.SetBinContent(i, total * pdf.getVal(nset) * h.GetBinWidth(i))
    chi2, nbins = merge_chi2_bins(h, h_pdf, 0.0, 0.0,
                                  CHI2_MIN_BIN_SIGNIFICANCE)
    ndf = nbins - 1 - npar
    p = ROOT.TMath.Prob(chi2, ndf) if ndf > 0 else 0.0
    return {"chi2": chi2, "ndf": ndf,
            "chi2_ndf": chi2 / ndf if ndf > 0 else -1.0,
            "pvalue": p, "n_bins": nbins}


def lo_core_hist(name, nbins, total):
    """The bare LO DY core (fitlib/dyCore.py, the base of mcdy) as a
    histogram in the plot binning, normalised to `total` in the window."""
    h = ROOT.TH1D(name, "", nbins, config.XLOW, config.XHIGH)
    h.SetDirectory(0)
    for i in range(1, nbins + 1):
        lo, w = h.GetBinLowEdge(i), h.GetBinWidth(i)
        step = w / LO_SUBSTEPS
        h.SetBinContent(i, step * sum(dyCore.core_value(lo + (j + 0.5) * step)
                                      for j in range(LO_SUBSTEPS)))
    h.Scale(total / h.Integral())
    h.SetLineColor(ROOT.TColor.GetColor(*LO_RGB))
    h.SetLineWidth(2)
    h.SetLineStyle(2)
    return h


def mc_hists(h, nbins, label):
    """The fitted MC histogram in the plot binning: filled, and its
    sum-w^2 uncertainty as a hatched band."""
    if h.GetNbinsX() % nbins:
        raise ValueError(f"{h.GetNbinsX()} fit bins do not rebin to {nbins}")
    filled = h.Rebin(h.GetNbinsX() // nbins, f"mc_fill_{label}")
    filled.SetDirectory(0)
    filled.SetFillColor(ROOT.TColor.GetColor(*MC_FILL_RGB))
    filled.SetLineColor(ROOT.TColor.GetColor(*MC_LINE_RGB))
    filled.SetLineWidth(1)
    band = filled.Clone(f"mc_band_{label}")
    band.SetDirectory(0)
    band.SetFillColor(ROOT.kGray + 2)
    band.SetFillStyle(MC_BAND_STYLE)
    band.SetLineWidth(0)
    band.SetMarkerSize(0)
    return filled, band


def legend_label(model, f):
    """'2Exp: chi2/ndf = 1.04, p = 0.27, pass' -- the GOF against the MC
    and whether the core is frozen."""
    name = mcCore.CORE_LABEL[model]
    if f["status"] > 1:
        return f"{name}: fit failed"
    verdict = "pass" if f["gof_pass"] else "#bf{fail}"
    if f["ndf"] <= 0:
        return f"{name}: low MC stats, {verdict}"
    p = "p < 0.001" if f["pvalue"] < 1e-3 else f"p = {f['pvalue']:.3f}"
    return (f"{name}: #chi^{{2}}/ndf = {f['chi2_ndf']:.2f}, {p}, "
            f"{verdict}")


def plot_cores(x, data, h, fits, label, what, year, outbase, pdf_too):
    """MC (filled histogram, hatched stat. band), the bare LO DY core
    (dashed) and the three fitted cores; pulls of the MC against each fit
    below."""
    nplot = int((config.XHIGH - config.XLOW) * PLOT_BINS_PER_GEV)
    c, outer, pad1, pad2 = style.make_canvas_pads(f"c_{label}")
    fr = x.frame(ROOT.RooFit.Bins(nplot))
    mc_fill, mc_band = mc_hists(h, nplot, label)
    fr.addTH1(mc_fill, "hist")
    fr.addTH1(mc_band, "E2")
    lo_core = lo_core_hist(f"lo_core_{label}", nplot, data.sumEntries())
    fr.addTH1(lo_core, "hist")
    # the points stay in the frame, invisible, as the reference for the pulls
    data.plotOn(fr, ROOT.RooFit.Name("mc"),
                ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2),
                ROOT.RooFit.Binning(nplot), ROOT.RooFit.Invisible())
    leg = ROOT.TLegend(1 - style.PAD_RIGHT_MARGIN - 0.52, 0.56,
                       1 - style.PAD_RIGHT_MARGIN - 0.03, 0.88)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.027)
    leg.AddEntry(mc_fill, "MC", "f")
    leg.AddEntry(mc_band, "MC stat. unc.", "f")
    leg.AddEntry(lo_core, "LO DY core (not fitted)", "l")
    pulls = []
    for model, f in fits.items():
        col = ROOT.TColor.GetColor(*CORE_RGB[model])
        f["pdf"].plotOn(fr, ROOT.RooFit.Name(model),
                        ROOT.RooFit.LineColor(col),
                        ROOT.RooFit.LineWidth(2))
        leg.AddEntry(fr.findObject(model), legend_label(model, f), "l")
        ph = fr.pullHist("mc", model)
        ph.SetLineColor(col)
        ph.SetMarkerColor(col)
        ph.SetMarkerSize(0.6)
        pulls.append(ph)

    pad1.cd()
    fr.SetTitle(f";;Weighted events / {1 / PLOT_BINS_PER_GEV:g} GeV")
    style._style_ratio_axis(fr.GetYaxis(), style.RATIO_Y_TITLE_OFFSET)
    fr.GetYaxis().SetTitleOffset(1.45)
    ROOT.TGaxis.SetMaxDigits(3)   # 10^n above 999: labels clear of the title
    ROOT.TGaxis.SetExponentOffset(-0.10, 0.01, "y")   # left of "CMS", above the title
    fr.GetXaxis().SetLabelSize(0)
    fr.GetXaxis().SetTitleSize(0)
    fr.SetMinimum(0.0)
    fr.SetMaximum(1.35 * fr.GetMaximum())
    fr.Draw()
    leg.Draw()
    style.cms_label(pad1, year, full_year=True)
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextSize(style.ANN_TEXT_SIZE)
    dy = style.LINE_SPACING * style.ANN_TEXT_SIZE
    for i, line in enumerate(["H #rightarrow #mu#mu", label, what,
                              "shape fit only"]):
        lat.DrawLatex(style.PAD_LEFT_MARGIN + 0.05, 0.85 - i * dy, line)

    pad2.cd()
    fr2 = x.frame()
    for ph in pulls:
        fr2.addPlotable(ph, "P")
    fr2.SetMinimum(-5)
    fr2.SetMaximum(5)
    fr2.SetTitle(";m_{#mu#mu} [GeV];Pull")
    fr2.GetYaxis().SetNdivisions(505)
    style._style_ratio_axis(fr2.GetXaxis(), style.RATIO_X_TITLE_OFFSET)
    style._style_ratio_axis(fr2.GetYaxis(), style.RATIO_Y_TITLE_OFFSET)
    fr2.Draw()
    zero = ROOT.TLine(config.XLOW, 0, config.XHIGH, 0)
    zero.SetLineStyle(2)
    zero.Draw()

    os.makedirs(os.path.dirname(outbase), exist_ok=True)
    c.SaveAs(outbase + ".png")
    if pdf_too:
        c.SaveAs(outbase + ".pdf")


def run_bin(args, binMVA, rng):
    from fitlib.prepareFits import bkg_mc_tags, getHistoBkgMC

    tag = f"{args.category}_{binMVA}_{args.year}"
    nbins = int((config.XHIGH - config.XLOW) * config.BINS_PER_GEV)
    what = "DY MC" if args.dy_only else "background MC"
    print(f"\n===== {tag}: {what}, {config.XLOW:g}-"
          f"{config.XHIGH:g} GeV, {nbins} bins =====")
    h, eras = getHistoBkgMC(nbins, config.XLOW, config.XHIGH,
                            args.category, args.year, binMVA or "incl",
                            dy_only=args.dy_only)
    if h is None or h.Integral() <= 0:
        raise SystemExit(f"{tag}: no background MC in the window")
    x = ROOT.RooRealVar(f"mh{args.category}", "m_{#mu#mu}",
                        config.XLOW, config.XHIGH)
    data = ROOT.RooDataHist(f"mc_{tag}", what,
                            ROOT.RooArgList(x), h)
    print(f"sum of weights {h.Integral():.1f}, "
          f"effective entries {h.GetEffectiveEntries():.0f}")

    fits, cores = {}, {}
    print(f"{'core':<7}{'npar':>5}{'NLL':>14}{'bins':>6}{'chi2/ndf':>10}"
          f"{'p(chi2)':>9}{'status':>8}{'try':>5}")
    for model in mcCore.CORES:
        f = fit_core(x, data, model, tag, args.tries, rng)
        f.update(gof_vs_mc(x, h, f["pdf"], len(f["params"])))
        converged = f["status"] <= 1
        f["gof_pass"] = f["pvalue"] > args.gof_pmin or args.gof_pmin <= 0
        ok = converged and f["gof_pass"]
        note = ("  <-- fit failed, left out" if not converged else
                f"  <-- {f['n_bins']} merged bins: low MC stats, "
                f"left out" if f["ndf"] <= 0 and not ok else
                f"  <-- fails GOF (p <= {args.gof_pmin:g}), left out"
                if not ok else
                "" if f["status"] == 0 else
                "  (covariance forced pos. def.; minimum kept)")
        print(f"{model:<7}{len(f['params']):>5}{f['nll']:>14.2f}"
              f"{f['n_bins']:>6}{f['chi2_ndf']:>10.3f}{f['pvalue']:>9.3f}"
              f"{f['status']:>8}{f['try']:>5}{note}")
        fits[model] = f
        if ok:
            cores[model] = {p: v["value"] for p, v in f["params"].items()}

    plot_cores(x, data, h, fits, f"{args.category} {binMVA or 'incl'}",
               what, args.year,
               os.path.join(args.plotdir, f"MCcore_{tag}"), args.pdf)

    rec = {
        "provenance": {
            "written_utc": datetime.now(timezone.utc)
                           .isoformat(timespec="seconds"),
            "git_commit": config.git_commit(),
            "category": args.category, "bin": binMVA, "year": args.year,
            "tag": tag, "window": [config.XLOW, config.XHIGH],
            "n_bins": nbins,
            "gof_pmin": args.gof_pmin,
            "gof_min_bin_significance": CHI2_MIN_BIN_SIGNIFICANCE,
            "dy_only": args.dy_only,
            "mc_tags": {era: bkg_mc_tags(args.category, era, args.dy_only)
                        for era in eras},
            "era_factors": {era: fac for era, (fac, _) in eras.items()},
            "n_files": {era: n for era, (_, n) in eras.items()},
        },
        "mc": {"sum_weights": h.Integral(),
               "effective_entries": h.GetEffectiveEntries()},
        "cores": cores,
        "fits": {m: {k: v for k, v in f.items() if k not in ("pdf", "keep")}
                 for m, f in fits.items()},
    }
    os.makedirs(args.wsdir, exist_ok=True)
    jf = os.path.join(args.wsdir, f"MCcore_{tag}.json")
    with open(jf, "w") as fh:
        json.dump(rec, fh, indent=2, sort_keys=True)
    print(f"wrote {jf}  (cores: {', '.join(cores) or 'none'})")


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--category", default="ggHcat")
    ap.add_argument("-b", "--bin", nargs="+", default=["incl"],
                    help="BDT bin label(s), as given to bkgFit")
    ap.add_argument("-y", "--year", default="Run3")
    ap.add_argument("--wsdir", default=None,
                    help="where MCcore_<tag>.json goes (default WS_<year>)")
    ap.add_argument("--plotdir",
                    default=os.path.expanduser(
                        "~/public_html/HmumuFits/mc_core"))
    ap.add_argument("--pdf", action="store_true")
    ap.add_argument("--all-bkg-mc", action="store_true",
                    help="fit the cores to the total background MC "
                         "(BKG_MC_TAGS) instead of the Drell-Yan MC alone "
                         "(BKG_DY_TAGS); same output file, the JSON records "
                         "which")
    ap.add_argument("--gof-pmin", type=float, default=GOF_THRESHOLD,
                    help="freeze a core only if its chi2 p-value against "
                         "the MC exceeds this (default %(default)s, as "
                         "bkgFit.py on data); 0 keeps every converged fit")
    ap.add_argument("--tries", type=int, default=5,
                    help="fits per core from random starts, until status 0")
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()
    args.wsdir = args.wsdir or config.default_wsdir(args.year)
    args.dy_only = not args.all_bkg_mc

    style.setup_style()
    rng = random.Random(args.seed)
    for b in args.bin:
        run_bin(args, b, rng)


if __name__ == "__main__":
    main()
