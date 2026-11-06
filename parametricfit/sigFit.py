"""
sigFit.py -- standalone H->mumu signal shape fits.

No Combine dependency for the fitting itself: uses ROOT's built-in
RooCrystalBall (ROOT >= 6.24), so this runs in the plain conda env.

Signal shape model
------------------
Six CB parameters are fitted to MC per (category, BDT bin, production mode)
and then frozen. The PDF does not see cb_mu / cb_sigma directly; it sees:

    mean  = (cb_mu + MH - 125) * (1 + scale_unc * CMS_scale_m)
    sigma =  cb_sigma          * (1 + res_unc   * CMS_res_m)

  MH           Higgs mass hypothesis (GeV). Additive shift, so the same
               workspace can be re-evaluated at 124, 126, ... for a mass
               scan. Combine drives this with -m. NOT a systematic.
  CMS_scale_m  unit-Gaussian nuisance for muon momentum scale
  CMS_res_m    unit-Gaussian nuisance for muon momentum resolution

The nuisances are exactly degenerate with cb_mu / cb_sigma, so they are
pinned to zero during the MC fit and only released afterwards.

Usage
-----
    python sigFit.py --selftest
    python sigFit.py -c VBFcat -b bdt0 -s qqH
    python sigFit.py -c VBFcat -b incl bdt0 bdt1 bdt2
    python sigFit.py -c ggHcat VBFcat VLcat VHcat Zinvcat TTLcat TTHcat \
                     -b bdt0 bdt1 bdt2
"""

import argparse
import json
import os
from datetime import datetime, timezone

import ROOT

# importing style (signal imports it too) sets ROOT up for batch plotting
from fitlib import config, signal, style


def run_category(category, binMVA, year, modes, plotdir, wsdir, freeze,
                 save_pdf=False):
    # imported late so --selftest runs without prepareFits
    from fitlib.prepareFits import (getHistoSignalByEra, normalize_year,
                             signal_yield_scales, expand_year,
                             DATA_YEAR_GROUPS, LUMIS)

    year = normalize_year(year)        # '_Run3' -> 'Run3'

    tag_cat = f"{category}_{binMVA}_{year}"
    x = ROOT.RooRealVar(f"mh{category}", "m_{#mu#mu}", config.XLOW, config.XHIGH)
    x.setRange("full", config.XLOW, config.XHIGH)

    w = ROOT.RooWorkspace("w", "workspace")
    nbins = int((config.XHIGH - config.XLOW) * config.BINS_PER_GEV)

    # plots go to <plotdir>/<category>/ so a run over all categories does not
    # dump ~30 files into one directory
    catdir = os.path.join(plotdir, category)

    # The shape is fitted on the MC eras; the yield is extrapolated to every
    # era of the year (2025/2026 borrow 2024's yield per fb^-1).
    scales = signal_yield_scales(year)

    dump = {
        "provenance": {
            "written_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "config.git_commit": config.git_commit(),
            "category": category, "bin": binMVA, "year": year,
            "yield_scales": scales,
            "lumi": {e: LUMIS[e]
                     for e in expand_year(year, DATA_YEAR_GROUPS)},
            "fit_range": [config.XLOW, config.XHIGH],
            "n_bins": nbins,
            "bins_per_gev": config.BINS_PER_GEV,
        },
        "fits": {},
    }

    keep = []   # hold references so PyROOT does not collect mid-loop
    for sig in modes:
        print(f"\n--- {category} {binMVA} {sig} ---")
        by_era = getHistoSignalByEra(nbins, config.XLOW, config.XHIGH, category, year,
                                     binMVA, sig)
        hist = None
        for h_era in by_era.values():       # the MC the shape is fitted on
            if hist is None:
                hist = h_era.Clone(f"h_{tag_cat}_{sig}")
                hist.SetDirectory(0)
            else:
                hist.Add(h_era)
        mc_by_era = {e: signal.range_integral(h) for e, h in by_era.items()}
        yield_info = {
            "value": sum(scales[e] * n for e, n in mc_by_era.items()),
            "by_era": {e: {"mc": n, "scale": scales[e]}
                       for e, n in mc_by_era.items()},
        }

        out = signal.fit_one(x, hist, f"{tag_cat}_{sig}",
                      f"{category} {binMVA} {sig}", catdir, freeze,
                      year=year, save_pdf=save_pdf, yield_info=yield_info)
        if out is None:
            dump["fits"][sig] = {"skipped": "empty histogram"}
            continue
        pdf, bundle, norm_var, record = out
        keep.append((pdf, bundle, norm_var))

        getattr(w, "import")(pdf)
        getattr(w, "import")(norm_var)
        dump["fits"][sig] = record

    # the datacard rates over the fit window, per mode and summed
    yields = {s: r["normalization"]["integral"]
              for s, r in dump["fits"].items() if "normalization" in r}
    dump["yields"] = {"window": [config.XLOW, config.XHIGH], "by_mode": yields,
                      "total": sum(yields.values())}
    print(f"\n  signal yield {config.XLOW:g}-{config.XHIGH:g} GeV: "
          + ", ".join(f"{s} {v:.2f}" for s, v in yields.items())
          + f"; total {dump['yields']['total']:.2f}")

    os.makedirs(wsdir, exist_ok=True)
    
    # Name by the modes being fit, so single-mode runs don't clobber each
    # other. When the full default set for the category is run, keep the
    # plain per-category name (that's the combined workspace Combine wants).
    if set(modes) == set(config.MC_LIST[category]):
        out_tag = tag_cat                       # e.g. ggHcat_incl_Run3
    else:
        out_tag = f"{tag_cat}_{'_'.join(modes)}" # e.g. ggHcat_incl_Run3_ggH

    wsfile = f"{wsdir}/Signal_{out_tag}_workspace.root"
    w.writeToFile(wsfile)

    jsonfile = f"{wsdir}/Signal_{out_tag}_params.json"
    with open(jsonfile, "w") as fh:
        json.dump(dump, fh, indent=2, sort_keys=True)

    print(f"\nwrote {wsfile}")
    print(f"wrote {jsonfile}")

    # Sanity check: the only floating variables should be the nuisances.
    floating = [v.GetName() for v in w.allVars()
                if not v.isConstant() and not v.GetName().startswith("mh")]
    print(f"floating in workspace: {floating or '(none)'}")

    n_bad = sum(1 for r in dump["fits"].values()
                if isinstance(r.get("fit_quality"), dict)
                and not r["fit_quality"]["ok"])
    if n_bad:
        print(f"WARNING: {n_bad} fit(s) flagged for review in {tag_cat}")



# ---------------------------------------------------------------------------
# selftest
# ---------------------------------------------------------------------------

TRUTH = {"mu": 125.10, "sigma": 2.00, "alphaL": 1.30,
         "nL": 5.00, "alphaR": 1.60, "nR": 8.00}


def closure(plotdir, n_gen=500000, tag="selftest", verbose=True):
    """Generate from a RooCrystalBall with known parameters and fit it back.

    Returns dict of {param: (fit, error, truth)} plus 'chi2_ndf', 'ok',
    or None if the fit failed.

    Generating from the same functional form is the point. A target that is
    not a CB (e.g. a sum of Gaussians) makes chi2/ndf meaningless and
    strongly binning-dependent, so it cannot serve as a pass criterion.
    """
    x = ROOT.RooRealVar(f"mh_toy_{tag}", "m_{#mu#mu}", config.XLOW, config.XHIGH)
    x.setRange("full", config.XLOW, config.XHIGH)
    x.setBins(int((config.XHIGH - config.XLOW) * config.BINS_PER_GEV))

    # Bounded and constant: RooRealVar(name, title, value) alone gets
    # [-inf, inf], and RooCrystalBall warns that sigma/alpha/n must be > 0.
    tv = {}
    for k, v in TRUTH.items():
        lo, hi = (v - 5.0, v + 5.0) if k == "mu" else (1e-3, 5.0 * v)
        tv[k] = ROOT.RooRealVar(f"true_{k}_{tag}", k, v, lo, hi)
        tv[k].setConstant(True)

    pdf_true = ROOT.RooCrystalBall(
        f"cb_true_{tag}", "cb_true", x,
        tv["mu"], tv["sigma"], tv["alphaL"], tv["nL"],
        tv["alphaR"], tv["nR"])

    ds = pdf_true.generateBinned(ROOT.RooArgSet(x), n_gen)
    h = ds.createHistogram(f"h_toy_{tag}", x)
    h.SetDirectory(0)

    out = signal.fit_one(x, h, tag, f"closure {tag}",
                  os.path.join(plotdir, "closure"), freeze=False)
    if out is None:
        return None

    _, bundle, _, record = out
    nom = bundle["nom"]
    res = {k: (nom[k].getVal(), nom[k].getError(), v) for k, v in TRUTH.items()}
    res["chi2_ndf"] = record["fit_quality"]["chi2_ndf"]
    res["ok"] = record["fit_quality"]["ok"]

    if verbose:
        print("\n  parameter recovery (pull = (fit - truth)/error):")
        for k in TRUTH:
            v, e, t = res[k]
            pull = (v - t) / e if e > 0 else float("inf")
            print(f"    {k:7s} fit={v:8.3f} +/- {e:6.3f}   truth={t:7.3f}   "
                  f"pull={pull:+6.2f}")
    return res


def selftest(plotdir, n_gen=500000):
    print("=== selftest: closure test against a generated CB ===")
    print("  truth: " + "  ".join(f"{k}={v:.3f}" for k, v in TRUTH.items()))

    res = closure(plotdir, n_gen)
    if res is None:
        print("SELFTEST FAILED: fit returned nothing")
        return 1

    def pull(k):
        v, e, t = res[k]
        return abs((v - t) / e) if e > 0 else float("inf")

    # mu / sigma propagate into the limit; the tail parameters are weakly
    # constrained and correlated, so they get a looser tolerance.
    core_ok = all(pull(k) < 3.0 for k in ("mu", "sigma"))
    tail_ok = all(pull(k) < 5.0 for k in ("alphaL", "nL", "alphaR", "nR"))
    chi2_ok = res["chi2_ndf"] < 2.0

    good = core_ok and tail_ok and chi2_ok
    print(f"\n  core (mu, sigma) within 3 sigma : {core_ok}")
    print(f"  tails within 5 sigma            : {tail_ok}")
    print(f"  chi2/ndf < 2                    : {chi2_ok}  ({res['chi2_ndf']:.2f})")
    print(f"\nselftest {'PASSED' if good else 'FAILED'}")
    return 0 if good else 1


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--category", nargs="+", default=["VBFcat"],
                    choices=sorted(config.MC_LIST))
    ap.add_argument("-b", "--bins", nargs="+", default=["bdt0"],
                    help="BDT bin label(s): bdt0 bdt1 bdt2 incl")
    ap.add_argument("-s", "--sig", nargs="+", default=None,
                    help="production mode(s); default = all for the category")
    ap.add_argument("-y", "--year", default="Run3",
                    help="concrete era (12022, 2024, ...) or group "
                         "(2022, 2023, Run3)")
    ap.add_argument("--wsdir", default=None,
                    help="fit JSONs and workspaces (default WS_<year>)")
    ap.add_argument("--plotdir",
                    default=os.path.expanduser(
                        "~/public_html/HmumuFits/signal_fits"),
                    help="plots land in <plotdir>/<category>/ (default "
                         "~/public_html/HmumuFits/signal_fits)")
    ap.add_argument("--pdf", action="store_true",
                    help="also write PDF alongside the PNG")
    ap.add_argument("--no-freeze", action="store_true",
                    help="leave shape parameters floating in the workspace")
    ap.add_argument("--selftest", action="store_true")
    ap.add_argument("--bins-per-gev", type=float, default=None,
                    help=f"override binning (default {config.BINS_PER_GEV})")
    ap.add_argument("--ngen", type=int, default=500000,
                    help="events generated for --selftest")
    args = ap.parse_args()
    args.wsdir = args.wsdir or config.default_wsdir(args.year)

    if args.bins_per_gev:
        config.BINS_PER_GEV = args.bins_per_gev
        print(f"binning override: {config.BINS_PER_GEV} bins/GeV "
              f"({int((config.XHIGH - config.XLOW) * config.BINS_PER_GEV)} bins)")

    if args.selftest:
        raise SystemExit(selftest(args.plotdir, args.ngen))

    for cat in args.category:
        modes = args.sig if args.sig else config.MC_LIST[cat]
        for b in args.bins:
            run_category(cat, b, args.year, modes,
                         args.plotdir, args.wsdir, not args.no_freeze,
                         save_pdf=args.pdf)


if __name__ == "__main__":
    style.setup_style()
    main()
