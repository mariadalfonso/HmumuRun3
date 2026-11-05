"""
mcCoreFit.py -- fit the frozen background cores to the background MC.

For each BDT bin, the weighted background MC (BKG_MC_TAGS in
fitlib/prepareFits.py) is fitted over the full window, unblinded since it
holds no signal, with three smooth shapes (fitlib/mcCore.py):

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

and a plot per bin in --plotdir. A core whose fit fails is left out of the
JSON (status 1, a forced-positive covariance, is kept: only the central
values are frozen), and bkgFit.py then has no candidates on it.

Usage (conda, ROOT >= 6.28), before bkgFit.py:

    python mcCoreFit.py -c ggHcat -b bdt0 bdt1 bdt2 -y Run3
"""

import argparse
import json
import os
import random
from datetime import datetime, timezone

import ROOT

from fitlib import config, mcCore, style

ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

CORE_RGB = {"mcexp": (42, 120, 214), "mcpow": (235, 104, 52),
            "mcdy": (27, 175, 122)}
PLOT_BINS_PER_GEV = 1          # coarser than the fit, for a readable plot
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


def chi2_vs_mc(x, h, pdf, npar):
    """Pearson chi2 of the MC histogram against the normalised core, with
    the MC's own (sum w^2) errors; bins without entries are skipped."""
    nset = ROOT.RooArgSet(x)
    total = h.Integral(1, h.GetNbinsX())
    chi2, nb = 0.0, 0
    for i in range(1, h.GetNbinsX() + 1):
        err = h.GetBinError(i)
        if err <= 0:
            continue
        x.setVal(h.GetBinCenter(i))
        exp = total * pdf.getVal(nset) * h.GetBinWidth(i)
        chi2 += ((h.GetBinContent(i) - exp) / err) ** 2
        nb += 1
    ndf = max(nb - npar, 1)
    return {"chi2": chi2, "ndf": ndf, "chi2_ndf": chi2 / ndf,
            "pvalue": ROOT.TMath.Prob(chi2, ndf)}


def plot_cores(x, data, fits, label, year, outbase, pdf_too):
    """MC with the three cores; pulls of the MC against each below."""
    nplot = int((config.XHIGH - config.XLOW) * PLOT_BINS_PER_GEV)
    c, outer, pad1, pad2 = style.make_canvas_pads(f"c_{label}")
    fr = x.frame(ROOT.RooFit.Bins(nplot))
    data.plotOn(fr, ROOT.RooFit.Name("mc"),
                ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2),
                ROOT.RooFit.Binning(nplot))
    leg = ROOT.TLegend(1 - style.PAD_RIGHT_MARGIN - 0.40, 0.62,
                       1 - style.PAD_RIGHT_MARGIN - 0.03, 0.88)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.030)
    leg.AddEntry(fr.findObject("mc"), "background MC", "pe")
    pulls = []
    for model, f in fits.items():
        col = ROOT.TColor.GetColor(*CORE_RGB[model])
        f["pdf"].plotOn(fr, ROOT.RooFit.Name(model),
                        ROOT.RooFit.LineColor(col),
                        ROOT.RooFit.LineWidth(2))
        leg.AddEntry(fr.findObject(model),
                     f"{mcCore.CORE_LABEL[model]}  "
                     f"#chi^{{2}}/ndf = {f['chi2_ndf']:.2f}", "l")
        ph = fr.pullHist("mc", model)
        ph.SetLineColor(col)
        ph.SetMarkerColor(col)
        ph.SetMarkerSize(0.6)
        pulls.append(ph)

    pad1.cd()
    fr.SetTitle(f";;Weighted events / {1 / PLOT_BINS_PER_GEV:g} GeV")
    style._style_ratio_axis(fr.GetYaxis(), style.RATIO_Y_TITLE_OFFSET)
    fr.GetYaxis().SetTitleOffset(1.45)
    ROOT.TGaxis.SetExponentOffset(-0.07, -0.02, "y")   # clear of "CMS"
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
    for i, line in enumerate(["H #rightarrow #mu#mu", label,
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
    from fitlib.prepareFits import BKG_MC_TAGS, getHistoBkgMC

    tag = f"{args.category}_{binMVA}_{args.year}"
    nbins = int((config.XHIGH - config.XLOW) * config.BINS_PER_GEV)
    print(f"\n===== {tag}: background MC, {config.XLOW:g}-"
          f"{config.XHIGH:g} GeV, {nbins} bins =====")
    h, eras = getHistoBkgMC(nbins, config.XLOW, config.XHIGH,
                            args.category, args.year, binMVA or "incl")
    if h is None or h.Integral() <= 0:
        raise SystemExit(f"{tag}: no background MC in the window")
    x = ROOT.RooRealVar(f"mh{args.category}", "m_{#mu#mu}",
                        config.XLOW, config.XHIGH)
    data = ROOT.RooDataHist(f"mc_{tag}", "background MC",
                            ROOT.RooArgList(x), h)
    print(f"sum of weights {h.Integral():.1f}, "
          f"effective entries {h.GetEffectiveEntries():.0f}")

    fits, cores = {}, {}
    print(f"{'core':<7}{'npar':>5}{'NLL':>14}{'chi2/ndf':>10}"
          f"{'p(chi2)':>9}{'status':>8}{'try':>5}")
    for model in mcCore.CORES:
        f = fit_core(x, data, model, tag, args.tries, rng)
        f.update(chi2_vs_mc(x, h, f["pdf"], len(f["params"])))
        ok = f["status"] <= 1
        note = ("" if f["status"] == 0 else
                "  (covariance forced pos. def.; minimum kept)" if ok else
                "  <-- failed, left out")
        print(f"{model:<7}{len(f['params']):>5}{f['nll']:>14.2f}"
              f"{f['chi2_ndf']:>10.3f}{f['pvalue']:>9.3f}{f['status']:>8}"
              f"{f['try']:>5}{note}")
        fits[model] = f
        if ok:
            cores[model] = {p: v["value"] for p, v in f["params"].items()}

    plot_cores(x, data, fits, f"{args.category} {binMVA or 'incl'} "
               f"background MC", args.year,
               os.path.join(args.plotdir, f"MCcore_{tag}"), args.pdf)

    per_era = BKG_MC_TAGS[args.category]
    rec = {
        "provenance": {
            "written_utc": datetime.now(timezone.utc)
                           .isoformat(timespec="seconds"),
            "git_commit": config.git_commit(),
            "category": args.category, "bin": binMVA, "year": args.year,
            "tag": tag, "window": [config.XLOW, config.XHIGH],
            "n_bins": nbins,
            "mc_tags": {era: per_era.get(era, per_era["*"]) for era in eras},
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
    ap.add_argument("--tries", type=int, default=5,
                    help="fits per core from random starts, until status 0")
    ap.add_argument("--seed", type=int, default=12345)
    args = ap.parse_args()
    args.wsdir = args.wsdir or config.default_wsdir(args.year)

    style.setup_style()
    rng = random.Random(args.seed)
    for b in args.bin:
        run_bin(args, b, rng)


if __name__ == "__main__":
    main()
