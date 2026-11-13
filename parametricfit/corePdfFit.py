"""
corePdfFit.py -- fit the Run 2 core-PDF (Sec. 6.4.2 of the note) to the
sidebands of all BDT bins of a category at once.

    python corePdfFit.py -c ggHcat -b bdt0 bdt1 bdt2 -y Run3
    python corePdfFit.py -c ggHcat -b bdt0 bdt1 bdt2 -y Run3 --smf-order 3 3 2

The model (fitlib/corePdf.py): in every bin, one of three cores (BWZRedux,
FEWZ x Bernstein, S-exponential; 3 parameters each) shared by all bins,
times that bin's own Bernstein shape modifier (SMF). Every fit is
simultaneous: one likelihood, the sum of the bins' masked sideband NLLs,
the core parameters common to all bins.

  1. each core alone on the summed sidebands of the bins: its start
  2. the SMF order of each bin, by F-test, for each core: every bin starts
     at order 1 (--smf-min-order, as the note); bin by bin the order goes
     up by one while 2 dNLL gives p < 0.05 for the extra parameter
     (bkgFit.py's F_THRESHOLD), and stops at the first step that does not
     or that does not converge. Each step starts from the order below,
     raised exactly by Bernstein degree elevation, so 2 dNLL >= 0.
     A bin then gets the highest order any core asked for, so every member
     of its multipdf has the same parameters and none has less freedom
     than it needs. --smf-order fixes the orders instead.
  3. every core at those orders; per bin, the sideband chi2 and p-value,
     and the fraction of the PDF in the sideband, which extrapolates the
     yield as bkgFit.py does

The F-test picks the order the DATA need. The note chose its orders (3 in
ggH-cat0 and cat1, 2 in the others) with the bias study, which may need
more: biasStudy.py --mode core checks the choice, and --smf-order raises a
bin that fails.

A core whose simultaneous fit does not converge is left out. Minuit
status 1 (a valid minimum, Hessian forced positive-definite) counts as
converged: the core and the SMFs partly mimic each other, so the Hessian
is near-singular by construction. --strict-status keeps status 0 only.
The other cores, ordered by their total NLL (no penalty: they all have the
same parameters), are the members of the shared envelope; the best one is
index 0, and gives every bin's extrapolated yield.

Runs in conda after bkgFit.py, from whose output it reads each bin's data
(<wsdir>/Bkg_<tag>_data.root) and binning (<wsdir>/Bkg_<tag>_params.json).
Writes, per bin,

    <wsdir>/Bkg_<tag>_core_params.json   the shared core, this bin's SMF,
                                         the yield, for makeWorkspaces.py
    <plotdir>/<cat>/core_<tag>_envelope.png

then buildWorkspace.py --members core makes the cards, and
biasStudy.py --mode core runs the note's bias test on them.

--smf-order: one order for every bin, or one per bin in the order of -b,
in place of the F-test. --smf-max-order (default 4) caps the F-test.
"""

import argparse
import json
import os
from datetime import datetime, timezone

import ROOT

import bkgFit
from fitlib import config, corePdf, style

ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)


def read_json(path):
    if not os.path.isfile(path):
        raise SystemExit(f"missing {path}; run bkgFit.py for this bin first")
    with open(path) as fh:
        return json.load(fh)


def load_bin(args, b, x):
    """The data of one bin, as bkgFit.py wrote it."""
    tag = f"{args.category}_{b}_{args.year}"
    rec = read_json(os.path.join(args.wsdir, f"Bkg_{tag}_params.json"))
    prov = rec["provenance"]
    if list(prov["window"]) != [config.XLOW, config.XHIGH]:
        raise SystemExit(f"{tag}: window {prov['window']} is not "
                         f"{config.XLOW}-{config.XHIGH}")
    path = os.path.join(args.wsdir, prov["data_file"])
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        raise SystemExit(f"cannot open {path}")
    h = f.Get(f"data_{tag}")
    if not h:
        raise SystemExit(f"{path} has no data_{tag}")
    h = h.Clone(f"h_{tag}")
    h.SetDirectory(0)
    f.Close()
    sb_idx = bkgFit.sideband_bins(h)
    return {"bin": b, "tag": tag, "h": h, "sb_idx": sb_idx,
            "n_sb": sum(h.GetBinContent(i) for i in sb_idx),
            "nbins": h.GetNbinsX(), "data_file": prov["data_file"],
            "data": ROOT.RooDataHist(f"datahist_{tag}", "data",
                                     ROOT.RooArgList(x), h)}


def param_record(pars):
    return {p.GetName(): {"value": p.getVal(), "error": p.getError(),
                          "min": p.getMin(), "max": p.getMax()}
            for p in pars}


# Minuit status 1: the minimum is valid, but the Hessian had to be made
# positive-definite. The core-PDF hits it by construction: the shared
# core and the bins' polynomials can partly mimic each other (a tilt of
# the core against the same tilt in every SMF), so the Hessian is close to
# singular. That leaves the minimum, which is all Combine uses, intact;
# --strict-status excludes such cores anyway, as bkgFit.py does.
ACCEPT_STATUS = (0, 1)


def minimize(nll, accept=ACCEPT_STATUS):
    """Migrad twice (as bkgFit.py's fits), then Hesse, with strategy 1: the
    FEWZ core costs ~13 ms per likelihood call (its formula sums all 161
    spline pieces at every bin), and strategy 2's extra Hessians would
    triple the time for the same minimum. If the status is not accepted,
    once more from there with strategy 2."""
    m = ROOT.RooMinimizer(nll)
    m.setMinimizerType("Minuit2")
    m.setPrintLevel(-1)
    m.setPrintEvalErrors(-1)
    for strategy in (1, 2):
        m.setStrategy(strategy)
        m.minimize("Minuit2", "Migrad")
        m.minimize("Minuit2", "Migrad")
        m.hesse()
        res = m.save()
        if res.status() in accept:
            break
    return res, strategy


def seed_core(core, args, x, sb_range, h_sum):
    """The start of a core: the core alone on the summed sidebands.
    Returns {core parameter name: value}, named as the shared core's."""
    data_sum = ROOT.RooDataHist(f"datahist_sum_{core}", "sum",
                                ROOT.RooArgList(x), h_sum)
    ctag = corePdf.core_tag(args.category, args.year)
    seed_tag = f"{ctag}_seed"
    seed = bkgFit.fit_candidate(x, data_sum, core, seed_tag, sb_range)
    return {k.replace(seed_tag, ctag): v["value"]
            for k, v in seed["params"].items()}, seed["status"]


def fit_core(core, args, x, bins, sb_range, orders, start):
    """The simultaneous fit of one core across the bins, bin b with an SMF
    of order orders[b], from the parameter values in `start` (by name;
    parameters not in it keep their defaults)."""
    core_pdf, ckeep = corePdf.create_core(x, args.category, args.year, core)
    members, nlls, keep = {}, ROOT.RooArgList(), [core_pdf, ckeep]
    for bd in bins:
        pdf, smf_pars, k = corePdf.create_member(x, core_pdf, core,
                                                 bd["tag"], orders[bd["bin"]])
        nll = pdf.createNLL(bd["data"], ROOT.RooFit.Range(sb_range))
        nlls.add(nll)
        members[bd["bin"]] = {"pdf": pdf, "smf": smf_pars}
        keep += [k, nll]
    pars = list(ckeep["params"]) + [v for m in members.values()
                                    for v in m["smf"]]
    for v in pars:
        if v.GetName() in start:
            v.setVal(min(max(start[v.GetName()], v.getMin()), v.getMax()))
    total = ROOT.RooAddition(f"nll_{core}", "total NLL", nlls)
    res, strategy = minimize(total, args.accept)

    out = {"core": core, "status": res.status(), "cov_qual": res.covQual(),
           "edm": res.edm(), "nll": res.minNll(), "strategy": strategy,
           "orders": dict(orders),
           "values": {v.GetName(): v.getVal() for v in pars},
           "core_params": param_record(ckeep["params"]), "bins": {},
           "keep": keep + [total, res, members]}
    # per bin: sideband fraction, yield, chi2
    nset = ROOT.RooArgSet(x)
    for bd in bins:
        mem = members[bd["bin"]]
        f_sb = mem["pdf"].createIntegral(
            nset, ROOT.RooFit.NormSet(nset),
            ROOT.RooFit.Range(sb_range)).getVal()
        rec = {"pdf": mem["pdf"], "sb_fraction": f_sb,
               "smf_params": param_record(mem["smf"]),
               "npar": corePdf.npar(core, orders[bd["bin"]])}
        rec.update(bkgFit.sideband_chi2(
            x, bd["h"], mem["pdf"], bd["sb_idx"],
            bkgFit.full_yield(bd["n_sb"], rec), rec["npar"]))
        out["bins"][bd["bin"]] = rec
    return out


def elevate(values, core, tag, n, m):
    """`values` with the SMF of (core, bin tag) raised from order n to m by
    Bernstein degree elevation: the same function, so the higher order
    starts exactly where the lower one ended.
        c'_k = k/(n+1) c_(k-1) + (1 - k/(n+1)) c_k,   c_0 = 1 stays 1"""
    out = dict(values)
    name = "smf_{}_c{{}}_{}".format(core, tag)
    for order in range(n, m):
        c = [1.0] + [out.get(name.format(k), corePdf.SMF_INIT)
                     for k in range(1, order + 1)] + [0.0]
        for k in range(1, order + 2):
            f = k / float(order + 1)
            out[name.format(k)] = f * c[k - 1] + (1 - f) * c[k]
    return out


def smf_ftest(core, args, x, bins, sb_range, start):
    """The SMF orders of one core by F-test, as bkgFit.py chooses a
    family's order: every bin starts at --smf-min-order; bin by bin, the
    order goes up by one while 2 dNLL of the simultaneous fit gives
    p < F_THRESHOLD for the one extra parameter, and stops at the first
    step that does not (or does not converge). The other bins keep their
    orders meanwhile; the core is refitted at every step.
    Returns (orders, steps, the last accepted fit)."""
    orders = {bd["bin"]: args.smf_min_order for bd in bins}
    cur = fit_core(core, args, x, bins, sb_range, orders, start)
    steps = []
    if cur["status"] not in args.accept:
        return orders, steps, cur
    for bd in bins:
        b = bd["bin"]
        while orders[b] < args.smf_max_order:
            trial = dict(orders)
            trial[b] += 1
            r = fit_core(core, args, x, bins, sb_range, trial,
                         elevate(cur["values"], core, bd["tag"], orders[b],
                                 trial[b]))
            step = {"bin": b, "from": orders[b], "to": trial[b]}
            if r["status"] not in args.accept:
                step.update(failed=True, status=r["status"])
                steps.append(step)
                break
            chi2 = max(2.0 * (cur["nll"] - r["nll"]), 0.0)
            p = ROOT.TMath.Prob(chi2, 1)
            step.update({"2dNLL": chi2, "pvalue": p,
                         "justified": p < bkgFit.F_THRESHOLD})
            steps.append(step)
            if not step["justified"]:
                break
            orders, cur = trial, r
    return orders, steps, cur


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--category", default="ggHcat",
                    choices=sorted(config.MC_LIST))
    ap.add_argument("-b", "--bins", nargs="+", required=True,
                    help="the BDT bins that share the core")
    ap.add_argument("-y", "--year", default="Run3")
    ap.add_argument("--wsdir", default=None,
                    help="bkgFit.py's output, and where this writes "
                         "(default WS_<year>)")
    ap.add_argument("--smf-order", nargs="+", type=int, default=None,
                    help="fix the SMF orders, one for all bins or one per "
                         "bin, instead of choosing them by F-test")
    ap.add_argument("--smf-min-order", type=int, default=1,
                    help="F-test: the order every bin starts at (default 1, "
                         "as the note)")
    ap.add_argument("--smf-max-order", type=int, default=4,
                    help="F-test: the highest order tried (default 4)")
    ap.add_argument("--cores", nargs="+", default=corePdf.CORES,
                    help=f"core functions (default {' '.join(corePdf.CORES)})")
    ap.add_argument("--strict-status", action="store_true",
                    help="keep only cores with Minuit status 0 (default: "
                         "also 1, a valid minimum whose Hessian was forced "
                         "positive-definite)")
    ap.add_argument("--plotdir", default=None,
                    help="default ~/public_html/HmumuFits/bkg_fits")
    ap.add_argument("--pdf", action="store_true")
    args = ap.parse_args()
    args.accept = (0,) if args.strict_status else ACCEPT_STATUS
    args.wsdir = args.wsdir or config.default_wsdir(args.year)
    args.plotdir = args.plotdir or os.path.expanduser(
        "~/public_html/HmumuFits/bkg_fits")
    if args.smf_order is not None:
        if len(args.smf_order) == 1:
            args.smf_order = args.smf_order * len(args.bins)
        if len(args.smf_order) != len(args.bins):
            ap.error(f"--smf-order: give one order, or one per bin "
                     f"({len(args.bins)})")
    style.setup_style()

    x = ROOT.RooRealVar(f"mh{args.category}", "m_{#mu#mu}",
                        config.XLOW, config.XHIGH)
    sb_range = bkgFit.setup_ranges(x)
    bins = [load_bin(args, b, x) for b in args.bins]
    if len({bd["nbins"] for bd in bins}) != 1:
        raise SystemExit("the bins have different binnings; rerun bkgFit.py "
                         "with the same --bins-per-gev")
    h_sum = bins[0]["h"].Clone("h_sum")
    h_sum.SetDirectory(0)
    for bd in bins[1:]:
        h_sum.Add(bd["h"])

    ctag = corePdf.core_tag(args.category, args.year)
    print(f"category : {args.category}, {args.year}\n"
          f"bins     : " + ", ".join(f"{bd['bin']} ({bd['n_sb']:.0f} "
                                     f"sideband events)" for bd in bins) +
          f"\ncores    : {', '.join(args.cores)}, shared as {ctag}\n"
          f"index    : {corePdf.index_name(args.category, args.year)}")

    # ---- the start of every core: alone on the summed sidebands ---------
    starts = {}
    for core in args.cores:
        starts[core], st = seed_core(core, args, x, sb_range, h_sum)
        if st != 0:
            print(f"  note: {core} alone on the summed bins: status {st}")

    # ---- the SMF orders -------------------------------------------------
    ftests = {}
    if args.smf_order is not None:
        orders = dict(zip(args.bins, args.smf_order))
        print("\nSMF orders: " + ", ".join(f"{b} {n}"
                                           for b, n in orders.items())
              + " (--smf-order)")
    else:
        print(f"\nSMF order by F-test (simultaneous fit; p < "
              f"{bkgFit.F_THRESHOLD} justifies one more order in a bin; "
              f"orders {args.smf_min_order}-{args.smf_max_order})")
        picks = {}
        for core in args.cores:
            o, steps, last = smf_ftest(core, args, x, bins, sb_range,
                                       starts[core])
            ftests[core] = steps
            if last["status"] not in args.accept:
                print(f"  {core:<11} the order-{args.smf_min_order} fit does "
                      f"not converge (status {last['status']}); no vote")
                continue
            picks[core] = o
            starts[core] = last["values"]
            for st in steps:
                if st.get("failed"):
                    print(f"  {core:<11} {st['bin']}: {st['from']} -> "
                          f"{st['to']}  fit failed (status {st['status']}): "
                          "stop here")
                else:
                    print(f"  {core:<11} {st['bin']}: {st['from']} -> "
                          f"{st['to']}  2dNLL={st['2dNLL']:7.2f}  "
                          f"p={st['pvalue']:.4f}  "
                          f"{'justified' if st['justified'] else 'NOT justified'}")
            print(f"  {core:<11} -> " + ", ".join(f"{b} {n}"
                                                 for b, n in o.items()))
        if not picks:
            raise SystemExit("no core converged at the lowest SMF order")
        # one order per bin for every core: each member of a bin's multipdf
        # then has the same parameters (no penalty needed, as in the note),
        # and no core is left with less freedom than it asked for
        orders = {b: max(o[b] for o in picks.values()) for b in args.bins}
        print("SMF orders: " + ", ".join(f"{b} {n}" for b, n in orders.items())
              + " (per bin, the highest any core needs)")

    # ---- every core at those orders -------------------------------------
    fits = {}
    print(f"\n{'core':<11}{'status':>7}{'edm':>10}{'total NLL':>16}   "
          + "  ".join(f"{bd['bin']} chi2/ndf (p)" for bd in bins))
    for core in args.cores:
        start = dict(starts[core])
        prev = ftests.get(core)
        if prev is not None:      # from the core's own F-test choice
            mine = {b: args.smf_min_order for b in args.bins}
            for st in prev:
                if st.get("justified"):
                    mine[st["bin"]] = st["to"]
            for bd in bins:
                start = elevate(start, core, bd["tag"], mine[bd["bin"]],
                                orders[bd["bin"]])
        r = fit_core(core, args, x, bins, sb_range, orders, start)
        fits[core] = r
        flag = ("" if r["status"] == 0 and r["cov_qual"] == 3 else
                "  (Hessian forced pos-def)" if r["status"] in args.accept
                else "  <-- not converged")
        if r["strategy"] == 2:
            flag += "  (strategy 2 retry)"
        print(f"{core:<11}{r['status']:>7}{r['edm']:>10.2g}{r['nll']:>16.2f}   "
              + "  ".join(f"{r['bins'][bd['bin']]['chi2_ndf']:6.3f} "
                          f"({r['bins'][bd['bin']]['pvalue']:.3f})"
                          for bd in bins) + flag)

    ok = sorted((c for c, r in fits.items() if r["status"] in args.accept),
                key=lambda c: fits[c]["nll"])
    dropped = [c for c in fits if c not in ok]
    if dropped:
        print("\nWARNING: not converged, left out of the envelope: "
              + ", ".join(f"{c} (status {fits[c]['status']})"
                          for c in dropped))
    if not ok:
        raise SystemExit("no core converged")
    best = ok[0]
    print(f"\nenvelope : {', '.join(ok)}   (best fit: {best}, "
          f"index 0)")
    for c in ok:
        print(f"  {c:<11} NLL {fits[c]['nll']:14.2f}  "
              f"(+{fits[c]['nll'] - fits[best]['nll']:.2f})")
    print("  core parameters (" + best + "): " + ", ".join(
        f"{k.replace('_' + ctag, '').replace(ctag, '')}="
        f"{v['value']:.4g}" for k, v in fits[best]["core_params"].items()))

    os.makedirs(args.wsdir, exist_ok=True)
    print()
    for bd in bins:
        b, tag = bd["bin"], bd["tag"]
        per = {c: fits[c]["bins"][b] for c in ok}
        norm_val = bkgFit.full_yield(bd["n_sb"], per[best])
        spread = [bkgFit.full_yield(bd["n_sb"], per[c]) for c in ok]
        rec = {
            "provenance": {
                "written_utc": datetime.now(timezone.utc)
                               .isoformat(timespec="seconds"),
                "git_commit": config.git_commit(),
                "model": "core-PDF (Run 2 note Sec. 6.4.2)",
                "category": args.category, "bin": b, "year": args.year,
                "tag": tag, "bins": args.bins,
                "window": [config.XLOW, config.XHIGH], "n_bins": bd["nbins"],
                "blind_region": [bkgFit.SR_LO, bkgFit.SR_HI],
                "observable": x.GetName(), "data_file": bd["data_file"],
                "core_tag": ctag,
                "index": corePdf.index_name(args.category, args.year),
                "smf_order": orders[b],
                "smf_orders": orders,
                "smf_order_source": ("--smf-order" if args.smf_order
                                     is not None else "F-test"),
                "f_threshold": bkgFit.F_THRESHOLD,
            },
            "smf_ftest": ftests,
            "members": ok, "best": best, "dropped": dropped,
            "accepted_status": list(args.accept),
            "fits": {c: {"status": fits[c]["status"],
                         "cov_qual": fits[c]["cov_qual"],
                         "nll_total": fits[c]["nll"], "edm": fits[c]["edm"],
                         "chi2": per[c]["chi2"], "ndf": per[c]["ndf"],
                         "pvalue": per[c]["pvalue"],
                         "sb_fraction": per[c]["sb_fraction"],
                         "npar": per[c]["npar"]} for c in ok},
            "core_params": {c: fits[c]["core_params"] for c in ok},
            "smf_params": {c: per[c]["smf_params"] for c in ok},
            "normalization": {"value": norm_val,
                              "sideband_count": bd["n_sb"],
                              "sb_fraction": per[best]["sb_fraction"],
                              "norm_range": list(bkgFit.NORM_RANGE),
                              "envelope_spread": [min(spread), max(spread)]},
        }
        path = os.path.join(args.wsdir, f"Bkg_{tag}_core_params.json")
        with open(path, "w") as fh:
            json.dump(rec, fh, indent=2, sort_keys=True)
        print(f"{tag}: norm {norm_val:.1f} (envelope spread "
              f"{min(spread):.1f}-{max(spread):.1f}), wrote {path}")

        plot = {c: {"pdf": per[c]["pdf"], "sb_fraction": per[c]["sb_fraction"]}
                for c in ok}
        bkgFit.plot_envelope(
            x, bd["data"], plot, ok, best, f"{args.category} {b} data",
            os.path.join(args.plotdir, args.category, f"core_{tag}"),
            args.year, sb_range, args.pdf, bd["nbins"], bd["n_sb"])


if __name__ == "__main__":
    main()
