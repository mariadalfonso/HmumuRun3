"""
biasStudy.py -- bias on r from the discrete-profiling background envelope.

For each tag = "<cat>_<bin>_<year>" and each member i of the envelope
(multipdf_<tag>_bkg in the card, indexed by pdfindex_<tag>):

  1. generate toys from member i (pdfindex_<tag> frozen to i), with r = r_inj
  2. fit every toy with the full model, pdfindex_<tag> FLOATING, exactly as
     the limit is computed (discrete profiling)
  3. Bias(mu_fit, mu_inj, sigma_mu) = (mu_fit - mu_inj) / sigma_mu, with mu
     Combine's r and sigma_mu the error on the side that faces mu_inj (the
     upper error if the fit is low, the lower if high)

The envelope is unbiased against member i if |median Bias| < 0.20 (--threshold):
the bias is then < 20% of the statistical error, which would inflate the total
error by < 2% (flashggFinalFit uses the tighter 0.14, < 1%). Every member of the
envelope must pass.

Toys come from the pre-fit workspace: each member at its own sideband fit, the
background norm at its rateParam value, CMS_scale_m / CMS_res_m drawn from
their constraints. --toysFrequentist is NOT used: it fits the full data first,
the blinded 120-130 GeV region included.

Step 6 of the chain, after bwsHrare.py (README). Runs in the Combine
container, from parametricfit/. --only report reads back the combine output
and needs only ROOT, so it also runs in conda.

Outputs, in <outdir>/<tag>/:
    <tag>.root                       text2workspace of the datacard
    meta.json                        envelope members, r range per r_inj
    toys_gen<i>_<r>_s<seed>.root     GenerateOnly output
    fit_gen<i>_<r>_s<seed>.root      MultiDimFit --algo singles output
    logs/                            combine stdout/stderr per job
    bias_<tag>_<r>.json              the pull summary per member
and plots in <plotdir>/<cat>/bias_<tag>_<r>_{<member>,summary}.png.

Usage (inside the Combine container, from parametricfit/)
-----
    python biasStudy.py -c ggHcat -b incl -y 2024 --ntoys 1000 -j 8
    python biasStudy.py -c ggHcat VBFcat -b bdt0 bdt1 -y Run3 --inject 0 1
    python biasStudy.py -c ggHcat -b incl -y 2024 --only report
    python biasStudy.py -c ggHcat -b incl -y 2024 --dry-run
"""

import argparse
import glob
import json
import math
import os
import shlex
import subprocess
from concurrent.futures import ThreadPoolExecutor

import ROOT

import pdfDefinitions as pdefs
import sigFit as sf

ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)

HERE = os.path.dirname(os.path.abspath(__file__))
FITS_DIR = os.path.normpath(os.path.join(HERE, "..", "analysis", "FITS"))

# ---------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------

MH = 125.0
MULTIPDF_CLASS = "RooMultiPdf"
PULL_THRESHOLD = 0.20      # |median pull| above this fails
MIN_VALID_FRAC = 0.90      # warn when fewer toys than this give a usable pull
PULL_RANGE = 5.0           # pull histogram and Gaussian fit over +/- this
PULL_BINS = 50
PULL_FILL_RGB = (166, 206, 227)    # light blue
PULL_LINE_RGB = (31, 120, 180)     # darker blue
PULL_GAUS_RGB = (153, 0, 13)       # dark red
SUMMARY_BAND_RGB = (217, 217, 217)  # light grey

# r range of the toy fits: r_inj -/+ R_WIDTH * sigma_r(Asimov). Wide enough
# that the 1 sigma interval of a toy never reaches the edge; --algo singles
# reports an interval that hits the edge as the edge itself, which is why the
# README's rMin = 0 runs lose the lower error.
R_WIDTH = 10.0
R_SEED_RANGE = 50.0        # range of the Asimov fit that sets it

# Discrete profiling, as recommended in the Combine documentation: Minuit
# strategy 0, the parameters of the inactive members frozen, and the
# constant terms hidden from the minimizer when the index changes.
COMBINE_FIT_OPTS = [
    "--cminDefaultMinimizerStrategy", "0",
    "--X-rtd", "MINIMIZER_freezeDisassociatedParams",
    "--X-rtd", "MINIMIZER_multiMin_hideConstants",
    "--X-rtd", "MINIMIZER_multiMin_maskConstraints",
    "--X-rtd", "MINIMIZER_multiMin_maskChannels=2",
    "--X-rtd", "REMOVE_CONSTANT_ZERO_POINT=1",
    "--robustFit", "1",
    "--cminFallbackAlgo", "Minuit2,Migrad,1:0.1",
    "--cminFallbackAlgo", "Minuit2,Simplex,1:0.5",
]

# quantileExpected of the three rows --algo singles writes per toy
Q_BEST, Q_LO, Q_HI = -1.0, -0.32, 0.32


def rlabel(r):
    """r_inj as a file-name token: 1 -> r1p00, -0.5 -> rm0p50."""
    return ("r%.2f" % r).replace("-", "m").replace(".", "p")


def short_name(pdf_name, tag):
    """bern3_ggHcat_incl_2024 -> bern3 (some members carry the tag without
    the underscore, e.g. exp2ggHcat_incl_2024)."""
    return pdf_name.replace(f"_{tag}", "").replace(tag, "")


def read_json(path):
    with open(path) as fh:
        return json.load(fh)


def write_json(path, obj):
    with open(path, "w") as fh:
        json.dump(obj, fh, indent=2, sort_keys=True)


# ---------------------------------------------------------------------------
# combine
# ---------------------------------------------------------------------------

def run_cmd(cmd, cwd, log, dry):
    """Run one command, its output to `log`. Returns True on exit code 0."""
    line = " ".join(shlex.quote(c) for c in cmd)
    if dry:
        print(f"[dry-run] (cd {cwd} && {line})")
        return True
    os.makedirs(os.path.dirname(log), exist_ok=True)
    with open(log, "w") as fh:
        fh.write(f"# cwd: {cwd}\n# {line}\n\n")
        fh.flush()
        rc = subprocess.call(cmd, cwd=cwd, stdout=fh, stderr=subprocess.STDOUT)
    if rc != 0:
        print(f"  !! exit code {rc}: see {log}")
    return rc == 0


def take_output(cwd, name, method, dest, dry):
    """Move higgsCombine<name>.<method>.mH*.root to dest; combine puts the
    mass and the seed in the file name."""
    if dry:
        return True
    found = glob.glob(os.path.join(cwd, f"higgsCombine{name}.{method}.mH*.root"))
    if len(found) != 1:
        print(f"  !! expected one higgsCombine{name}.{method} output, "
              f"found {len(found)}")
        return False
    os.replace(found[0], dest)
    return True


def shapes_basedir(card, candidates):
    """The directory the card's relative `shapes` paths resolve from.

    bwsHrare.py writes the path as given on its command line, relative to
    where it ran (analysis/FITS for the README commands), so text2workspace
    must run from there."""
    paths = []
    with open(card) as fh:
        for line in fh:
            f = line.split()
            if len(f) >= 4 and f[0] == "shapes":
                paths.append(f[3])
    rel = [p for p in paths if not os.path.isabs(p)]
    if not rel:
        return os.path.dirname(card)
    for c in candidates:
        if all(os.path.exists(os.path.join(c, p)) for p in rel):
            return c
    raise SystemExit(f"{card}: cannot find {rel[0]} relative to any of "
                     f"{candidates}")


def build_t2w(card, t2w, dry):
    """text2workspace the card, unless the output is newer than the card and
    the workspace file it points to."""
    base = shapes_basedir(card, [os.path.dirname(card),
                                 os.path.dirname(os.path.dirname(card)),
                                 FITS_DIR])
    if os.path.isfile(t2w) and os.path.getmtime(t2w) > os.path.getmtime(card):
        print(f"  reusing {t2w}")
        return
    cmd = ["text2workspace.py", "-m", "%g" % MH, card, "-o", t2w]
    log = os.path.join(os.path.dirname(t2w), "logs", "text2workspace.log")
    if not run_cmd(cmd, base, log, dry):
        raise SystemExit(f"text2workspace failed: see {log}")


def find_multipdf(w, tag):
    """The RooMultiPdf indexed by pdfindex_<tag>. text2workspace renames the
    card's multipdf_<tag>_bkg to shapeBkg_<process>_<bin>, so match on the
    index category, not the name."""
    index = f"pdfindex_{tag}"
    multis = [p for p in w.allPdfs() if p.InheritsFrom(MULTIPDF_CLASS)]
    found = [p for p in multis if p.findServer(index)]
    if len(found) != 1:
        raise SystemExit(f"expected one {MULTIPDF_CLASS} indexed by {index}, "
                         f"found {len(found)}; {MULTIPDF_CLASS}s in the "
                         f"workspace: {[p.GetName() for p in multis]}")
    return found[0]


def read_envelope(t2w, tag):
    """Names of the envelope members, in pdfindex order."""
    if not pdefs.load_combine(quiet=True):
        raise SystemExit("needs libHiggsAnalysisCombinedLimit: run inside "
                         "the Combine container")
    f = ROOT.TFile.Open(t2w)
    if not f or f.IsZombie():
        raise SystemExit(f"cannot open {t2w}")
    w = f.Get("w")
    if not w:
        raise SystemExit(f"{t2w}: no workspace 'w'")
    multi = find_multipdf(w, tag)
    names = [multi.getPdf(i).GetName() for i in range(multi.getNumPdfs())]
    f.Close()
    return names


def read_singles(path, index_name=None):
    """{iToy: {best, lo, hi, index}} from a MultiDimFit --algo singles file.

    Rows are grouped by iToy and identified by quantileExpected, not by their
    order, so a toy whose interval was not found simply lacks lo/hi."""
    out = {}
    f = ROOT.TFile.Open(path)
    if not f or f.IsZombie():
        print(f"  !! cannot open {path}")
        return out
    t = f.Get("limit")
    if not t:
        print(f"  !! {path} has no limit tree")
        f.Close()
        return out
    has_index = bool(index_name) and bool(t.GetBranch(index_name))
    for _ in t:
        toy = out.setdefault(int(t.iToy), {})
        q = float(t.quantileExpected)
        if abs(q - Q_BEST) < 1e-3:
            toy["best"] = float(t.r)
            if has_index:
                toy["index"] = int(getattr(t, index_name))
        elif abs(q - Q_LO) < 5e-3:
            toy["lo"] = float(t.r)
        elif abs(q - Q_HI) < 5e-3:
            toy["hi"] = float(t.r)
    f.Close()
    return out


def asimov_rrange(t2w, tag, r_inj, outdir, extra, dry):
    """r_inj -/+ R_WIDTH * sigma_r, sigma_r from an Asimov fit with the
    envelope floating; the wider of the two errors is used."""
    name = f"_asimov_{rlabel(r_inj)}"
    dest = os.path.join(outdir, f"asimov_{rlabel(r_inj)}.root")
    cmd = (["combine", "-M", "MultiDimFit", "-d", t2w, "-m", "%g" % MH,
            "-t", "-1", "--expectSignal", "%g" % r_inj,
            "-P", "r", "--algo", "singles", "--floatOtherPOIs", "1",
            "--setParameters", "r=%g" % r_inj,
            "--setParameterRanges",
            "r=%g,%g" % (r_inj - R_SEED_RANGE, r_inj + R_SEED_RANGE),
            "-n", name] + COMBINE_FIT_OPTS + extra)
    ok = run_cmd(cmd, outdir, os.path.join(outdir, "logs", f"{name[1:]}.log"),
                 dry) and take_output(outdir, name, "MultiDimFit", dest, dry)
    if dry:
        return [r_inj - R_SEED_RANGE, r_inj + R_SEED_RANGE], "fallback"
    # one fit, whatever iToy combine gives the Asimov (not 0)
    rows = list(read_singles(dest).values()) if ok else []
    toy = rows[0] if len(rows) == 1 else {}
    errs = [abs(toy[k] - toy["best"]) for k in ("lo", "hi")
            if k in toy and "best" in toy]
    errs = [e for e in errs if 0 < e < 0.9 * R_SEED_RANGE]
    if not errs:
        print(f"  !! Asimov fit gave no usable error on r (see "
              f"logs/{name[1:]}.log); using r in +/- {R_SEED_RANGE:g}. "
              f"Set --rrange to override.")
        return [r_inj - R_SEED_RANGE, r_inj + R_SEED_RANGE], "fallback"
    sig = max(errs)
    print(f"  Asimov, r = {r_inj:g}: sigma_r = {sig:.3f}")
    return [r_inj - R_WIDTH * sig, r_inj + R_WIDTH * sig], "asimov"


def toy_job(t2w, tag, igen, r_inj, rrange, ntoys, seed, outdir, extra,
            force, dry):
    """Generate `ntoys` from member igen and fit them with the envelope
    floating. Skipped if its fit output already exists (unless force)."""
    idx = f"pdfindex_{tag}"
    stem = f"gen{igen}_{rlabel(r_inj)}_s{seed}"
    toys = os.path.join(outdir, f"toys_{stem}.root")
    fit = os.path.join(outdir, f"fit_{stem}.root")
    if os.path.isfile(fit) and not force:
        return stem, "exists"
    logs = os.path.join(outdir, "logs")

    gen = ["combine", "-M", "GenerateOnly", "-d", t2w, "-m", "%g" % MH,
           "-t", str(ntoys), "-s", str(seed), "--saveToys",
           "--expectSignal", "%g" % r_inj,
           "--setParameters", f"{idx}={igen},r={r_inj:g}",
           "--freezeParameters", f"{idx},MH",
           "-n", f"_toys_{stem}"] + extra
    if not (run_cmd(gen, outdir, os.path.join(logs, f"toys_{stem}.log"), dry)
            and take_output(outdir, f"_toys_{stem}", "GenerateOnly", toys,
                            dry)):
        return stem, "generation failed"

    # the index is NOT frozen here: the fit picks the member, as in the limit
    fitc = (["combine", "-M", "MultiDimFit", "-d", t2w, "-m", "%g" % MH,
             "-t", str(ntoys), "-s", str(seed), "--toysFile", toys,
             "-P", "r", "--algo", "singles", "--floatOtherPOIs", "1",
             "--saveSpecifiedIndex", idx,
             "--setParameters", "r=%g" % r_inj,
             "--setParameterRanges", "r=%g,%g" % tuple(rrange),
             "--freezeParameters", "MH",
             "-n", f"_fit_{stem}"] + COMBINE_FIT_OPTS + extra)
    if not (run_cmd(fitc, outdir, os.path.join(logs, f"fit_{stem}.log"), dry)
            and take_output(outdir, f"_fit_{stem}", "MultiDimFit", fit, dry)):
        return stem, "fit failed"
    return stem, "done"


def run_tag(tag, args):
    outdir = os.path.abspath(os.path.join(args.outdir, tag))
    os.makedirs(os.path.join(outdir, "logs"), exist_ok=True)
    card = os.path.abspath(os.path.join(args.carddir, f"datacard_{tag}.txt"))
    if not os.path.isfile(card):
        raise SystemExit(f"missing {card}; run bwsHrare.py first (README "
                         "step 4)")
    t2w = os.path.join(outdir, f"{tag}.root")
    extra = shlex.split(args.combine_opts)

    print(f"\n=== {tag} ===\n  card     : {card}")
    build_t2w(card, t2w, args.dry_run)

    meta_path = os.path.join(outdir, "meta.json")
    meta = read_json(meta_path) if os.path.isfile(meta_path) else {}
    if args.dry_run and not os.path.isfile(t2w):
        names = meta.get("envelope") or [f"member{i}_{tag}" for i in range(3)]
    else:
        names = read_envelope(t2w, tag)
    if meta.get("envelope") not in (None, names):
        print(f"  !! envelope changed since the last run "
              f"({meta['envelope']} -> {names}); rerun with --force")
    meta.update({"tag": tag, "card": card, "workspace": t2w,
                 "envelope": names, "mh": MH,
                 "combine_opts": COMBINE_FIT_OPTS + extra})
    meta.setdefault("rrange", {})
    meta.setdefault("rrange_source", {})   # asimov, fallback or user
    meta.setdefault("ntoys", {})       # per job, to count fits that wrote nothing
    print("  envelope : " + ", ".join(f"[{i}] {short_name(n, tag)}"
                                      for i, n in enumerate(names)))

    jobs = []
    for r_inj in args.inject:
        key = rlabel(r_inj)
        if args.rrange:
            meta["rrange"][key] = list(args.rrange)
            meta["rrange_source"][key] = "user"
        elif key not in meta["rrange"] or args.force:
            meta["rrange"][key], meta["rrange_source"][key] = asimov_rrange(
                t2w, tag, r_inj, outdir, extra, args.dry_run)
        elif key not in meta["rrange_source"]:  # meta.json from before
            wide = meta["rrange"][key][1] - meta["rrange"][key][0]
            meta["rrange_source"][key] = ("fallback" if
                                          wide >= 2 * R_SEED_RANGE - 1e-6
                                          else "asimov")
        rrange = meta["rrange"][key]
        print(f"  r = {r_inj:g}: toy fits in r [{rrange[0]:.3g}, "
              f"{rrange[1]:.3g}]")
        if meta["rrange_source"][key] == "fallback":
            print(f"  !! this is the +/- {R_SEED_RANGE:g} fallback (the "
                  f"Asimov fit failed); --force retries it")
        nchunks = int(math.ceil(args.ntoys / float(args.toys_per_job)))
        for igen in range(len(names)):
            for k in range(nchunks):
                n = min(args.toys_per_job, args.ntoys - k * args.toys_per_job)
                seed = args.seed + 1000 * igen + k
                meta["ntoys"][f"gen{igen}_{key}_s{seed}"] = n
                jobs.append((t2w, tag, igen, r_inj, rrange, n, seed, outdir,
                             extra, args.force, args.dry_run))
    if not args.dry_run:
        write_json(meta_path, meta)

    print(f"  {len(jobs)} job(s), {args.jobs} at a time")
    status = {}
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        for stem, st in pool.map(lambda j: toy_job(*j), jobs):
            status.setdefault(st, []).append(stem)
            if st not in ("done", "exists"):
                print(f"  {stem}: {st}")
    for st, stems in sorted(status.items()):
        print(f"  {st:<18} {len(stems)}")


# ---------------------------------------------------------------------------
# report
# ---------------------------------------------------------------------------

def collect_pulls(files, r_inj, rrange, index_name, requested):
    """Pulls from all fit files of one (member, r_inj), with the reasons the
    other toys were dropped. `requested` maps a file's stem to the toys its
    job generated, so a toy whose fit wrote no row still counts as failed."""
    edge = 1e-3 * (rrange[1] - rrange[0])
    pulls, dr, chosen = [], [], {}
    drop = {"no_fit": 0, "no_interval": 0, "at_r_boundary": 0}
    for path in files:
        toys = read_singles(path, index_name)
        stem = os.path.basename(path)[len("fit_"):-len(".root")]
        drop["no_fit"] += max(requested.get(stem, 0) - len(toys), 0)
        for toy in toys.values():
            if "best" not in toy:
                drop["no_fit"] += 1
                continue
            if "index" in toy:
                chosen[toy["index"]] = chosen.get(toy["index"], 0) + 1
            d = toy["best"] - r_inj
            # the error on the side facing the injected value
            side = "hi" if d < 0 else "lo"
            if side not in toy:
                drop["no_interval"] += 1
                continue
            bound = toy[side]
            if ((side == "hi" and bound > rrange[1] - edge) or
                    (side == "lo" and bound < rrange[0] + edge)):
                drop["at_r_boundary"] += 1
                continue
            sig = abs(bound - toy["best"])
            if sig <= 0:
                drop["no_interval"] += 1
                continue
            pulls.append(d / sig)
            dr.append(d)
    return pulls, dr, chosen, drop


def median(v):
    s = sorted(v)
    n = len(s)
    if n == 0:
        return float("nan")
    return s[n // 2] if n % 2 else 0.5 * (s[n // 2 - 1] + s[n // 2])


def pull_hist(name, pulls):
    h = ROOT.TH1D(name, ";Bias = (#mu_{fit} - #mu_{inj}) / #sigma_{#mu};toys",
                  PULL_BINS, -PULL_RANGE, PULL_RANGE)
    h.SetDirectory(0)
    for p in pulls:
        h.Fill(p)
    return h


def summarize(pulls, dr, chosen, drop, nmembers, threshold):
    n = len(pulls)
    ntot = n + sum(drop.values())
    s = {"n_toys": ntot, "n_valid": n, "dropped": drop,
         "chosen_index": {str(i): chosen.get(i, 0) for i in range(nmembers)}}
    if n == 0:
        s.update(status="NO TOYS", **{"pass": False})
        return s
    mean = sum(pulls) / n
    rms = math.sqrt(sum((p - mean) ** 2 for p in pulls) / max(n - 1, 1))
    s.update(median_pull=median(pulls),
             # standard error of the median of a Gaussian sample
             median_pull_err=1.2533 * rms / math.sqrt(n),
             mean_pull=mean, rms_pull=rms,
             median_dr=median(dr))
    h = pull_hist("h_tmp_summarize", pulls)
    if n >= 20:
        res = h.Fit("gaus", "QSN0").Get()
        if res and res.Status() == 0:
            s.update(gaus_mean=res.Parameter(1), gaus_mean_err=res.ParError(1),
                     gaus_sigma=res.Parameter(2),
                     gaus_sigma_err=res.ParError(2))
    s["pass"] = abs(s["median_pull"]) < threshold
    s["status"] = "ok" if s["pass"] else "BIASED"
    if ntot and n < MIN_VALID_FRAC * ntot:
        s["status"] += " (few valid toys)"
    return s


def draw_label(pad, right):
    t, l, r = pad.GetTopMargin(), pad.GetLeftMargin(), pad.GetRightMargin()
    y = 1 - t + 0.2 * t
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(61)
    lat.SetTextSize(sf.CMS_TEXT_SIZE_FRAC * t)
    lat.DrawLatex(l, y, sf.CMS_TEXT)
    if sf.CMS_EXTRA_TEXT:
        lat.SetTextFont(52)
        lat.SetTextSize(sf.EXTRA_OVER_CMS_TEXT_SIZE * sf.CMS_TEXT_SIZE_FRAC * t)
        lat.DrawLatex(l + 0.14, y, sf.CMS_EXTRA_TEXT)   # full-height pad
    lat.SetTextFont(42)
    lat.SetTextAlign(31)
    lat.SetTextSize(sf.LUMI_TEXT_SIZE_FRAC * t)
    lat.DrawLatex(1 - r, y, right)
    return lat


def new_canvas(name):
    sf.setup_style()
    c = ROOT.TCanvas(name, "", sf.CANVAS_W, sf.CANVAS_H)
    c.SetTopMargin(sf.PAD_TOP_MARGIN)
    c.SetBottomMargin(sf.PAD_BOTTOM_MARGIN)
    c.SetLeftMargin(sf.PAD_LEFT_MARGIN)
    c.SetRightMargin(sf.PAD_RIGHT_MARGIN)
    return c


def save(c, outbase, pdf):
    c.SaveAs(f"{outbase}.png")
    if pdf:
        c.SaveAs(f"{outbase}.pdf")


def plot_member(pulls, s, tag, gen, r_inj, outbase, pdf):
    c = new_canvas(f"c_{os.path.basename(outbase)}")
    h = pull_hist(f"h_{os.path.basename(outbase)}", pulls)
    h.SetFillColor(ROOT.TColor.GetColor(*PULL_FILL_RGB))
    h.SetLineColor(ROOT.TColor.GetColor(*PULL_LINE_RGB))
    h.SetLineWidth(2)
    h.SetMaximum(1.6 * max(h.GetMaximum(), 1))
    h.Draw("HIST")
    keep = [h]
    if "gaus_mean" in s:
        g = ROOT.TF1(f"g_{os.path.basename(outbase)}", "gaus",
                     -PULL_RANGE, PULL_RANGE)
        g.SetParameters(h.GetMaximum() / 1.6, s["gaus_mean"], s["gaus_sigma"])
        h.Fit(g, "Q0")
        g.SetLineColor(ROOT.TColor.GetColor(*PULL_GAUS_RGB))
        g.SetLineWidth(2)
        g.Draw("same")
        keep.append(g)
    x0 = c.GetLeftMargin() + 0.03
    leg = ROOT.TLegend(x0, 0.68, x0 + 0.45, 0.89)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.035)
    leg.SetHeader(f"{tag}, #mu_{{inj}} = {r_inj:g}")
    leg.AddEntry(h, f"{gen} (n_{{toys}} = {s['n_valid']})", "f")
    if "gaus_mean" in s:
        leg.AddEntry(g, "Gaussian fit, mean = %+.2f, #sigma = %.2f"
                     % (s["gaus_mean"], s["gaus_sigma"]), "l")
    if s.get("n_valid"):
        leg.AddEntry(ROOT.nullptr, "median = %+.3f #pm %.3f"
                     % (s["median_pull"], s["median_pull_err"]), "")
    leg.Draw()
    keep += [leg, draw_label(c, "toys (%.1f TeV)" % sf.SQRT_S_TEV)]
    save(c, outbase, pdf)


def plot_summary(results, tag, r_inj, threshold, outbase, pdf):
    """Median Bias per generating member, with the +/- threshold band."""
    n = len(results)
    c = new_canvas(f"c_{os.path.basename(outbase)}")
    c.SetBottomMargin(0.18)
    frame = ROOT.TH1D(f"f_{os.path.basename(outbase)}",
                      ";background model;median Bias", n, -0.5, n - 0.5)
    frame.SetDirectory(0)
    lim = max([threshold * 2.5] +
              [abs(r.get("median_pull", 0)) + r.get("median_pull_err", 0)
               for r in results.values()]) * 1.3
    frame.SetMinimum(-lim)
    frame.SetMaximum(lim)
    for i, gen in enumerate(results):
        frame.GetXaxis().SetBinLabel(i + 1, gen)
    frame.GetXaxis().SetLabelSize(0.05)
    frame.GetXaxis().SetTitleOffset(1.6)
    frame.Draw("AXIS")
    band = ROOT.TBox(-0.5, -threshold, n - 0.5, threshold)
    band.SetFillColor(ROOT.TColor.GetColor(*SUMMARY_BAND_RGB))
    band.Draw()
    zero = ROOT.TLine(-0.5, 0, n - 0.5, 0)
    zero.SetLineStyle(2)
    zero.Draw()
    gr = ROOT.TGraphErrors()
    for i, r in enumerate(results.values()):
        if "median_pull" in r:
            k = gr.GetN()
            gr.SetPoint(k, i, r["median_pull"])
            gr.SetPointError(k, 0, r["median_pull_err"])
    gr.SetMarkerStyle(20)
    gr.Draw("P same")
    frame.Draw("AXIS same")
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(42)
    lat.SetTextSize(0.035)
    lat.DrawLatex(c.GetLeftMargin() + 0.04, 0.86,
                  f"{tag}, #mu_{{inj}} = {r_inj:g}")
    keep = [frame, band, zero, gr,  # noqa: F841 (ROOT needs them alive)
            draw_label(c, "toys (%.1f TeV)" % sf.SQRT_S_TEV)]
    save(c, outbase, pdf)


def report_tag(tag, args):
    outdir = os.path.join(args.outdir, tag)
    meta_path = os.path.join(outdir, "meta.json")
    if not os.path.isfile(meta_path):
        print(f"\n=== {tag} === no {meta_path}; run the toys first")
        return {}
    meta = read_json(meta_path)
    names = [short_name(n, tag) for n in meta["envelope"]]
    index_name = f"pdfindex_{tag}"
    category = tag.split("_")[0]
    plotdir = os.path.join(args.plotdir, category)
    os.makedirs(plotdir, exist_ok=True)

    out = {}
    for r_inj in args.inject:
        key = rlabel(r_inj)
        if key not in meta["rrange"]:
            print(f"\n=== {tag}, r = {r_inj:g} === no toys")
            continue
        rrange = meta["rrange"][key]
        results = {}
        for igen, gen in enumerate(names):
            files = sorted(glob.glob(os.path.join(
                outdir, f"fit_gen{igen}_{key}_s*.root")))
            requested = meta.get("ntoys", {})
            pulls, dr, chosen, drop = collect_pulls(files, r_inj, rrange,
                                                    index_name, requested)
            # jobs that never wrote a fit file (combine crashed): their toys
            # count as dropped, so the valid fraction is out of all requested
            have = {os.path.basename(p)[len("fit_"):-len(".root")]
                    for p in files}
            failed = sorted(st for st in requested
                            if st.startswith(f"gen{igen}_{key}_s")
                            and st not in have)
            drop["job_failed"] = sum(requested[st] for st in failed)
            s = summarize(pulls, dr, chosen, drop, len(names), args.threshold)
            s["n_files"] = len(files)
            s["failed_jobs"] = failed
            results[gen] = s
            if pulls:
                plot_member(pulls, s, tag, gen, r_inj,
                            os.path.join(plotdir, f"bias_{tag}_{key}_{gen}"),
                            args.pdf)
        if any("median_pull" in s for s in results.values()):
            plot_summary(results, tag, r_inj, args.threshold,
                         os.path.join(plotdir, f"bias_{tag}_{key}_summary"),
                         args.pdf)
        summary = {"tag": tag, "r_inj": r_inj, "r_range": rrange,
                   "threshold": args.threshold, "envelope": names,
                   "pull": "Bias = (mu_fit - mu_inj) / sigma_mu, mu = r, "
                           "sigma_mu on the side facing mu_inj",
                   "members": results,
                   "pass": all(s.get("pass", False) for s in results.values())}
        write_json(os.path.join(outdir, f"bias_{tag}_{key}.json"), summary)
        print_table(summary)
        out[key] = summary
    return out


def print_table(summary):
    names = summary["envelope"]
    print(f"\n=== {summary['tag']}, r_inj = {summary['r_inj']:g}, "
          f"r in [{summary['r_range'][0]:.3g}, {summary['r_range'][1]:.3g}]"
          f", |median| < {summary['threshold']:g} ===")
    print(f"  {'toys from':<10} {'valid':>11}  {'median':>15}  {'mean':>6}  "
          f"{'gaus mu':>7}  {'gaus sig':>8}  {'median dr':>9}  "
          f"{'fit picks ' + '/'.join(names):<24}  status")
    for gen, s in summary["members"].items():
        tot = sum(s["chosen_index"].values()) or 1
        picks = "/".join("%.2f" % (v / float(tot))
                         for v in s["chosen_index"].values())
        if "median_pull" not in s:
            print(f"  {gen:<10} {s['n_valid']:>5}/{s['n_toys']:<5}  "
                  f"{'-':>15}  {'':>6}  {'':>7}  {'':>8}  {'':>9}  "
                  f"{'':<24}  {s['status']}")
            continue
        med = "%+.3f +/- %.3f" % (s["median_pull"], s["median_pull_err"])
        gm = "%+.3f" % s["gaus_mean"] if "gaus_mean" in s else "-"
        gs = "%.3f" % s["gaus_sigma"] if "gaus_sigma" in s else "-"
        print(f"  {gen:<10} {s['n_valid']:>5}/{s['n_toys']:<5}  {med:>15}  "
              f"{s['mean_pull']:+6.3f}  {gm:>7}  {gs:>8}  "
              f"{s['median_dr']:+9.3f}  {picks:<24}  {s['status']}")
        dropped = {k: v for k, v in s["dropped"].items() if v}
        if dropped:
            print(f"  {'':<10} dropped: " +
                  ", ".join(f"{k} {v}" for k, v in dropped.items()))
    failed = [st for s in summary["members"].values()
              for st in s.get("failed_jobs", [])]
    if failed:
        print(f"  !! {len(failed)} job(s) wrote no fit: {', '.join(failed)}; "
              "rerun the same command to retry them")
    print("  envelope " + ("PASSES" if summary["pass"] else "FAILS"))


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--category", nargs="+", default=["ggHcat"])
    ap.add_argument("-b", "--bins", nargs="+", default=["incl"],
                    help="BDT bin label(s), as given to sigFit/bkgFit")
    ap.add_argument("-y", "--year", default="Run3")
    ap.add_argument("--carddir", default=os.path.join(FITS_DIR, "workspaces"),
                    help="where bwsHrare.py wrote datacard_<tag>.txt "
                         "(default analysis/FITS/workspaces)")
    ap.add_argument("--outdir", default="BIAS",
                    help="toys, fits and summaries go to <outdir>/<tag>/")
    ap.add_argument("--plotdir",
                    default=os.path.expanduser("~/public_html/HmumuFits/bias"),
                    help="plots land in <plotdir>/<category>/")
    ap.add_argument("--inject", nargs="+", type=float, default=[1.0],
                    help="injected signal strength(s) (default 1)")
    ap.add_argument("--ntoys", type=int, default=1000,
                    help="toys per envelope member and r_inj (default 1000)")
    ap.add_argument("--toys-per-job", type=int, default=100)
    ap.add_argument("-j", "--jobs", type=int, default=1,
                    help="combine processes run at once")
    ap.add_argument("--seed", type=int, default=12345,
                    help="job k of member i uses seed + 1000*i + k")
    ap.add_argument("--rrange", nargs=2, type=float, default=None,
                    metavar=("RMIN", "RMAX"),
                    help="r range of the toy fits (default: from an Asimov "
                         f"fit, r_inj -/+ {R_WIDTH:g} sigma_r)")
    ap.add_argument("--threshold", type=float, default=PULL_THRESHOLD,
                    help=f"max |median Bias| (default {PULL_THRESHOLD})")
    ap.add_argument("--combine-opts", default="",
                    help="extra options for every combine call, e.g. "
                         "\"--toysNoSystematics\"")
    ap.add_argument("--only", choices=("run", "report"), default=None)
    ap.add_argument("--force", action="store_true",
                    help="redo jobs whose output exists, and the r range")
    ap.add_argument("--dry-run", action="store_true",
                    help="print the combine commands, run nothing")
    ap.add_argument("--pdf", action="store_true",
                    help="also write PDF alongside the PNG")
    args = ap.parse_args()

    tags = [f"{c}_{b}_{args.year}" for c in args.category for b in args.bins]
    if args.only != "report":
        for tag in tags:
            run_tag(tag, args)
    if args.only != "run" and not args.dry_run:
        for tag in tags:
            report_tag(tag, args)


if __name__ == "__main__":
    main()
