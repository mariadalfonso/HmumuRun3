"""
biasStudy.py -- bias on r from the background model: the discrete-profiling
envelope (default), or each function on its own (--mode single, Run 2).

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

--mode single: the Run 2 note's single-function study (Sec. 6.4.1)
------------------------------------------------------------------
The truths are every family's F-test pick, the ones the goodness-of-fit cut
rejected included (the buildWorkspace.py --members families cards), and
each toy is fitted once per function, pdfindex FROZEN to it: a truth x fit
matrix of Bias. As in the note, the Bias is the MEAN of the pulls
(--metric), a function is usable if |Bias| < 0.20 against every truth,
and of the usable ones the one with the highest expected significance
(Asimov, background frozen to that function, one -M Significance per
function) is chosen. The toys of a truth are generated once and fitted by
every function. --fit-models restricts the fit functions (the note fitted
with six, and used Landau x Bernstein only as a truth).

Step 6 of the chain, after buildWorkspace.py (README). Runs in the Combine
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

--mode single writes to <outdir>/<tag>_single/ instead:
    toys_gen<i>_<r>_s<seed>.root          as above, truth i
    fit_gen<i>_fit<j>_<r>_s<seed>.root    those toys fitted with member j
    sig_fit<j>.root                       expected significance with j
    bias_<tag>_<r>_matrix.json            the matrix, usable fits, choice
and the plot <plotdir>/<cat>/bias_<tag>_<r>_matrix.png.

Usage (inside the Combine container, from parametricfit/)
-----
    python biasStudy.py -c ggHcat -b incl -y 2024 --ntoys 1000 -j 8
    python biasStudy.py -c ggHcat VBFcat -b bdt0 bdt1 -y Run3 --inject 0 1
    python biasStudy.py -c ggHcat -b incl -y 2024 --only report
    python biasStudy.py -c ggHcat -b incl -y 2024 --dry-run
    python biasStudy.py --mode single -c ggHcat -b bdt0 bdt1 bdt2 -y Run3 \
        --ntoys 2000 -j 10
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

from fitlib import pdfDefinitions as pdefs
from fitlib import config, style

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
# summary plots: one row per generating member, sizes in pixels
ROWS_CANVAS_W = 700
ROWS_MIN_H = 450
ROW_PX = 24                 # height of one row
ROW_HEADROOM = 2            # empty rows on top for the title and legend
ROWS_TOP_PX, ROWS_BOTTOM_PX = 48, 64
ROWS_LEFT_PX, ROWS_RIGHT_PX = 120, 150   # member names | value column
ROW_LABEL_PX, ROW_VALUE_PX = 16, 14

# r range of the toy fits: r_inj -/+ R_WIDTH * sigma_r(Asimov). Wide enough
# that the 1 sigma interval of a toy never reaches the edge; --algo singles
# reports an interval that hits the edge as the edge itself, which is why the
# README's rMin = 0 runs lose the lower error.
R_WIDTH = 10.0
R_SEED_RANGE = 50.0        # range of the Asimov fit that sets it

# Discrete profiling, as recommended in the Combine documentation: Minuit
# strategy 0, the parameters of the inactive members frozen, and the
# constant terms hidden from the minimizer when the index changes.
DISCRETE_OPTS = [
    "--cminDefaultMinimizerStrategy", "0",
    "--X-rtd", "MINIMIZER_freezeDisassociatedParams",
    "--X-rtd", "MINIMIZER_multiMin_hideConstants",
    "--X-rtd", "MINIMIZER_multiMin_maskConstraints",
    "--X-rtd", "MINIMIZER_multiMin_maskChannels=2",
    "--X-rtd", "REMOVE_CONSTANT_ZERO_POINT=1",
]
COMBINE_FIT_OPTS = DISCRETE_OPTS + [
    "--robustFit", "1",
    "--cminFallbackAlgo", "Minuit2,Migrad,1:0.1",
    "--cminFallbackAlgo", "Minuit2,Simplex,1:0.5",
]

# --mode: which cards it reads, where it writes, which Bias it quotes
#   envelope  toys from each envelope member, fitted with pdfindex floating
#   single    the Run 2 note (Sec. 6.4.1): toys from every family's F-test
#             pick, each fitted with ONE function at a time (pdfindex frozen)
CARD_SUFFIX = {"envelope": "", "single": "_families"}   # buildWorkspace.py
OUT_SUFFIX = {"envelope": "", "single": "_single"}
DEFAULT_METRIC = {"envelope": "median", "single": "mean"}
MATRIX_PASS_RGB = (232, 232, 232)   # light grey
MATRIX_FAIL_RGB = (246, 196, 190)   # light red

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
    where it ran (parametricfit/ for buildWorkspace.py, analysis/FITS/ for
    older cards), so text2workspace must run from there."""
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
    base = shapes_basedir(card, [HERE, os.path.dirname(card),
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


def gen_toys(t2w, tag, igen, r_inj, ntoys, seed, outdir, extra, force, dry):
    """`ntoys` toys from member igen (pdfindex frozen to it), r = r_inj.
    Returns (stem, toys file or None). An existing file is reused, unless
    force (the envelope may have changed since)."""
    idx = f"pdfindex_{tag}"
    stem = f"gen{igen}_{rlabel(r_inj)}_s{seed}"
    toys = os.path.join(outdir, f"toys_{stem}.root")
    if os.path.isfile(toys) and not force:
        return stem, toys
    gen = ["combine", "-M", "GenerateOnly", "-d", t2w, "-m", "%g" % MH,
           "-t", str(ntoys), "-s", str(seed), "--saveToys",
           "--expectSignal", "%g" % r_inj,
           "--setParameters", f"{idx}={igen},r={r_inj:g}",
           "--freezeParameters", f"{idx},MH",
           "-n", f"_toys_{stem}"] + extra
    log = os.path.join(outdir, "logs", f"toys_{stem}.log")
    if not (run_cmd(gen, outdir, log, dry)
            and take_output(outdir, f"_toys_{stem}", "GenerateOnly", toys,
                            dry)):
        return stem, None
    return stem, toys


def fit_toys(t2w, tag, toys, stem, r_inj, rrange, ntoys, seed, outdir,
             extra, dry, ifit=None):
    """Fit the toys of `stem`: pdfindex floating (ifit None), as in the
    limit, or frozen to member ifit. Returns the fit output's stem."""
    idx = f"pdfindex_{tag}"
    fstem = stem if ifit is None else stem.replace("_", f"_fit{ifit}_", 1)
    if ifit is None:
        index = ["--saveSpecifiedIndex", idx,
                 "--setParameters", "r=%g" % r_inj,
                 "--freezeParameters", "MH"]
    else:
        index = ["--setParameters", f"r={r_inj:g},{idx}={ifit}",
                 "--freezeParameters", f"MH,{idx}"]
    fitc = (["combine", "-M", "MultiDimFit", "-d", t2w, "-m", "%g" % MH,
             "-t", str(ntoys), "-s", str(seed), "--toysFile", toys,
             "-P", "r", "--algo", "singles", "--floatOtherPOIs", "1",
             "--setParameterRanges", "r=%g,%g" % tuple(rrange)]
            + index + ["-n", f"_fit_{fstem}"] + COMBINE_FIT_OPTS + extra)
    fit = os.path.join(outdir, f"fit_{fstem}.root")
    log = os.path.join(outdir, "logs", f"fit_{fstem}.log")
    ok = (run_cmd(fitc, outdir, log, dry)
          and take_output(outdir, f"_fit_{fstem}", "MultiDimFit", fit, dry))
    return fstem, ok


def toy_job(t2w, tag, igen, r_inj, rrange, ntoys, seed, outdir, extra,
            force, dry, fits=None):
    """One chunk of toys from member igen, then its fits: with the envelope
    floating (fits None), or once per member in `fits` with the index
    frozen. Fits whose output exists are skipped (unless force). Returns
    [(fit stem, status)]."""
    stem = f"gen{igen}_{rlabel(r_inj)}_s{seed}"
    todo = [None] if fits is None else list(fits)

    def out(ifit):
        f = stem if ifit is None else stem.replace("_", f"_fit{ifit}_", 1)
        return f, os.path.join(outdir, f"fit_{f}.root")
    if not force:
        done = [(f, "exists") for f, path in map(out, todo)
                if os.path.isfile(path)]
        todo = [i for i in todo if not os.path.isfile(out(i)[1])]
        if not todo:
            return done
    else:
        done = []
    stem, toys = gen_toys(t2w, tag, igen, r_inj, ntoys, seed, outdir, extra,
                          force, dry)
    if toys is None:
        return done + [(out(i)[0], "generation failed") for i in todo]
    for ifit in todo:
        fstem, ok = fit_toys(t2w, tag, toys, stem, r_inj, rrange, ntoys,
                             seed, outdir, extra, dry, ifit)
        done.append((fstem, "done" if ok else "fit failed"))
    return done


def sig_job(t2w, tag, ifit, outdir, extra, force, dry):
    """Expected significance with the background frozen to member ifit:
    Asimov (-t -1) from the pre-fit model, r = 1. The figure of merit the
    Run 2 note chooses the function by (its Table 6.14)."""
    idx = f"pdfindex_{tag}"
    name = f"_sig_fit{ifit}"
    dest = os.path.join(outdir, f"sig_fit{ifit}.root")
    if os.path.isfile(dest) and not force:
        return ifit, "exists"
    cmd = (["combine", "-M", "Significance", "-d", t2w, "-m", "%g" % MH,
            "-t", "-1", "--expectSignal", "1",
            "--setParameters", f"{idx}={ifit}",
            "--freezeParameters", f"MH,{idx}", "-n", name]
           + DISCRETE_OPTS + extra)
    ok = (run_cmd(cmd, outdir, os.path.join(outdir, "logs", f"{name[1:]}.log"),
                  dry)
          and take_output(outdir, name, "Significance", dest, dry))
    return ifit, "done" if ok else "failed"


def read_significance(path):
    """The expected significance in a -M Significance output, or None."""
    if not os.path.isfile(path):
        return None
    f = ROOT.TFile.Open(path)
    t = f.Get("limit") if f and not f.IsZombie() else None
    z = None
    if t and t.GetEntries():
        t.GetEntry(0)
        z = float(t.limit)
    if f:
        f.Close()
    return z


def fit_indices(names, tag, wanted):
    """pdfindex of every member to fit with in --mode single: all, or those
    named in --fit-models (members this bin does not have are skipped)."""
    short = [short_name(n, tag) for n in names]
    if not wanted:
        return list(range(len(names)))
    missing = [w for w in wanted if w not in short]
    if missing:
        print(f"  note: no {', '.join(missing)} in this bin; not fitted")
    return [short.index(w) for w in wanted if w in short]


def tag_outdir(tag, args):
    """<outdir>/<tag>, or <outdir>/<tag>_single for --mode single"""
    return os.path.abspath(os.path.join(args.outdir,
                                        tag + OUT_SUFFIX[args.mode]))


def run_tag(tag, args):
    outdir = tag_outdir(tag, args)
    os.makedirs(os.path.join(outdir, "logs"), exist_ok=True)
    card = os.path.abspath(os.path.join(
        args.carddir, f"datacard_{tag}{CARD_SUFFIX[args.mode]}.txt"))
    if not os.path.isfile(card):
        hint = (" --members families" if args.mode == "single" else "")
        raise SystemExit(f"missing {card}; run buildWorkspace.py{hint} "
                         "first, or point --carddir at the cards")
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
    fits = None
    if args.mode == "single":
        fits = fit_indices(names, tag, args.fit_models)
        meta["fits"] = fits
    meta.update({"tag": tag, "card": card, "workspace": t2w, "mode": args.mode,
                 "envelope": names, "mh": MH,
                 "combine_opts": COMBINE_FIT_OPTS + extra})
    meta.setdefault("rrange", {})
    meta.setdefault("rrange_source", {})   # asimov, fallback or user
    meta.setdefault("ntoys", {})       # per job, to count fits that wrote nothing
    print(("  envelope : " if args.mode == "envelope" else "  truths   : ")
          + ", ".join(f"[{i}] {short_name(n, tag)}"
                      for i, n in enumerate(names)))
    if fits is not None:
        print("  fitted   : " + ", ".join(short_name(names[i], tag)
                                         for i in fits)
              + "   (one at a time, pdfindex frozen)")

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
                for ifit in ([None] if fits is None else fits):
                    fs = (f"gen{igen}_{key}_s{seed}" if ifit is None else
                          f"gen{igen}_fit{ifit}_{key}_s{seed}")
                    meta["ntoys"][fs] = n
                jobs.append((t2w, tag, igen, r_inj, rrange, n, seed, outdir,
                             extra, args.force, args.dry_run, fits))
    if not args.dry_run:
        write_json(meta_path, meta)

    nfit = 1 if fits is None else len(fits)
    print(f"  {len(jobs)} job(s) of {nfit} fit(s) each, {args.jobs} at a time")
    status = {}
    with ThreadPoolExecutor(max_workers=max(1, args.jobs)) as pool:
        if fits is not None:     # the significances first: one fit each
            for ifit, st in pool.map(
                    lambda i: sig_job(t2w, tag, i, outdir, extra, args.force,
                                      args.dry_run), fits):
                if st not in ("done", "exists"):
                    print(f"  significance with "
                          f"{short_name(names[ifit], tag)}: {st}")
        for res in pool.map(lambda j: toy_job(*j), jobs):
            for stem, st in res:
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


def summarize(pulls, dr, chosen, drop, nmembers, threshold, metric="median"):
    """Pull statistics of one truth (and fit function). The quoted Bias,
    "bias" +/- "bias_err", is the median (flashggFinalFit) or the mean (the
    Run 2 note) of the pulls, as `metric` says; pass is |bias| < threshold."""
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
             mean_pull=mean, mean_pull_err=rms / math.sqrt(n), rms_pull=rms,
             median_dr=median(dr), metric=metric)
    s["bias"] = s[f"{metric}_pull"]
    s["bias_err"] = s[f"{metric}_pull_err"]
    h = pull_hist("h_tmp_summarize", pulls)
    if n >= 20:
        res = h.Fit("gaus", "QSN0").Get()
        if res and res.Status() == 0:
            s.update(gaus_mean=res.Parameter(1), gaus_mean_err=res.ParError(1),
                     gaus_sigma=res.Parameter(2),
                     gaus_sigma_err=res.ParError(2))
    s["pass"] = abs(s["bias"]) < threshold
    s["status"] = "ok" if s["pass"] else "BIASED"
    if ntot and n < MIN_VALID_FRAC * ntot:
        s["status"] += " (few valid toys)"
    return s


def draw_label(pad, right, t=None):
    """CMS label and `right` text in the top band of height t (NDC; default
    the top margin)."""
    l, r = pad.GetLeftMargin(), pad.GetRightMargin()
    t = pad.GetTopMargin() if t is None else t
    y = 1 - t + 0.2 * t
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(61)
    lat.SetTextSize(style.CMS_TEXT_SIZE_FRAC * t)
    lat.DrawLatex(l, y, style.CMS_TEXT)
    if style.CMS_EXTRA_TEXT:
        lat.SetTextFont(52)
        lat.SetTextSize(style.EXTRA_OVER_CMS_TEXT_SIZE * style.CMS_TEXT_SIZE_FRAC * t)
        # 0.14 of the 600 px wide canvas, kept at 84 px on wider ones
        lat.DrawLatex(l + 0.14 * style.CANVAS_W / pad.GetWw(), y,
                      style.CMS_EXTRA_TEXT)
    lat.SetTextFont(42)
    lat.SetTextAlign(31)
    lat.SetTextSize(style.LUMI_TEXT_SIZE_FRAC * t)
    lat.DrawLatex(1 - r, y, right)
    return lat


def new_canvas(name):
    style.setup_style()
    c = ROOT.TCanvas(name, "", style.CANVAS_W, style.CANVAS_H)
    c.SetTopMargin(style.PAD_TOP_MARGIN)
    c.SetBottomMargin(style.PAD_BOTTOM_MARGIN)
    c.SetLeftMargin(style.PAD_LEFT_MARGIN)
    c.SetRightMargin(style.PAD_RIGHT_MARGIN)
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
    keep += [leg, draw_label(c, "toys (%.1f TeV)" % style.SQRT_S_TEV)]
    save(c, outbase, pdf)


def plot_rows(groups, title, threshold, outbase, pdf):
    """Median Bias as one row per generating member, grouped under a header
    row per BDT bin. groups is a list of (header or None, {member: summary}).

    Rows rather than columns, so a long envelope never crowds the labels:
    the member is the axis label, its value is printed in a column on the
    right, and a failing member is a red square with "fails" next to its
    value. All text is in pixel fonts (precision 3), so the canvas grows
    with the number of rows without the text changing size.
    """
    rows = []                                   # (label, summary or None)
    for header, members in groups:
        if header is not None:
            rows.append((f"#bf{{{header}}}", None))
        rows += [(gen, s) for gen, s in members.items()]
    n, head = len(rows), ROW_HEADROOM
    ny = n + head
    w = ROWS_CANVAS_W
    h = max(ROWS_MIN_H, ROWS_TOP_PX + ROWS_BOTTOM_PX + ROW_PX * ny)
    style.setup_style()
    c = ROOT.TCanvas(f"c_{os.path.basename(outbase)}", "", w, h)
    c.SetTopMargin(ROWS_TOP_PX / h)
    c.SetBottomMargin(ROWS_BOTTOM_PX / h)
    c.SetLeftMargin(ROWS_LEFT_PX / w)
    c.SetRightMargin(ROWS_RIGHT_PX / w)
    c.SetTicky(0)

    vals = [abs(s["bias"]) + s["bias_err"] for _, s in rows
            if s and "bias" in s]
    metric = next((s["metric"] for _, s in rows if s and "metric" in s),
                  "median")
    lim = max([threshold * 2.5] + vals) * 1.12
    frame = ROOT.TH2F(f"f_{os.path.basename(outbase)}", f";{metric} Bias;",
                      100, -lim, lim, ny, 0, ny)
    frame.SetDirectory(0)
    frame.SetStats(0)
    ya = frame.GetYaxis()
    for k, (label, _) in enumerate(rows):
        ya.SetBinLabel(n - k, label)
    for b in range(n + 1, ny + 1):
        ya.SetBinLabel(b, "")
    ya.SetLabelFont(43)
    ya.SetLabelSize(ROW_LABEL_PX)
    ya.SetTickLength(0)
    xa = frame.GetXaxis()
    xa.SetLabelFont(43)
    xa.SetLabelSize(ROW_LABEL_PX)
    xa.SetTitleFont(43)
    xa.SetTitleSize(ROW_LABEL_PX + 2)
    xa.SetTitleOffset(1.3)
    xa.SetNdivisions(505)
    frame.Draw("AXIS")

    band = ROOT.TBox(-threshold, 0, threshold, n)
    band.SetFillColor(ROOT.TColor.GetColor(*SUMMARY_BAND_RGB))
    band.Draw()
    zero = ROOT.TLine(0, 0, 0, n)
    zero.SetLineStyle(2)
    zero.Draw()
    keep = [c, frame, band, zero]
    # a thin line above every header row after the first
    for k, (label, s) in enumerate(rows):
        if s is None and k > 0:
            ln = ROOT.TLine(-lim, n - k, lim, n - k)
            ln.SetLineColor(ROOT.kGray)
            ln.Draw()
            keep.append(ln)

    def ndc_y(y):
        b, t = c.GetBottomMargin(), c.GetTopMargin()
        return b + y / float(ny) * (1 - b - t)

    ok_col = ROOT.TColor.GetColor(*PULL_LINE_RGB)
    bad_col = ROOT.TColor.GetColor(*PULL_GAUS_RGB)
    ok, bad = ROOT.TGraphErrors(), ROOT.TGraphErrors()
    val = ROOT.TLatex()
    val.SetNDC()
    val.SetTextFont(43)
    val.SetTextSize(ROW_VALUE_PX)
    val.SetTextAlign(12)
    xv = 1 - c.GetRightMargin() + 0.015
    for k, (label, s) in enumerate(rows):
        if s is None:
            continue
        y = n - k - 0.5
        if "bias" not in s:
            val.SetTextColor(ROOT.kGray + 2)
            val.DrawLatex(xv, ndc_y(y), s.get("status", "no toys"))
            continue
        gr = ok if s.get("pass") else bad
        i = gr.GetN()
        gr.SetPoint(i, s["bias"], y)
        gr.SetPointError(i, s["bias_err"], 0)
        val.SetTextColor(ROOT.kBlack if s.get("pass") else bad_col)
        val.DrawLatex(xv, ndc_y(y), "%+.2f #pm %.2f%s" % (
            s["bias"], s["bias_err"],
            "" if s.get("pass") else "  fails"))
    for gr, col, mk in ((ok, ok_col, 20), (bad, bad_col, 21)):
        gr.SetMarkerStyle(mk)
        gr.SetMarkerSize(1.1)
        gr.SetMarkerColor(col)
        gr.SetLineColor(col)
        gr.SetLineWidth(2)
        if gr.GetN():
            gr.Draw("P same")
    frame.Draw("AXIS same")

    # headroom: title on the left, legend on the right, value column header
    yh = ndc_y(n + head / 2.0)
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(43)
    lat.SetTextSize(ROW_LABEL_PX)
    lat.SetTextAlign(12)
    lat.DrawLatex(c.GetLeftMargin() + 0.02, yh, title)
    leg = ROOT.TLegend(0.47, yh - 0.6 * ROW_PX / h,
                       1 - c.GetRightMargin() - 0.01, yh + 0.6 * ROW_PX / h)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(43)
    leg.SetTextSize(ROW_VALUE_PX)
    leg.SetNColumns(3)
    leg.SetMargin(0.3)
    leg.AddEntry(ok, "passes", "pe")
    leg.AddEntry(bad, "fails", "pe")
    leg.AddEntry(band, "|%s| < %g" % (metric, threshold), "f")
    leg.Draw()
    val.SetTextColor(ROOT.kGray + 2)
    val.DrawLatex(xv, yh, f"{metric} #pm err")
    keep += [ok, bad, val, lat, leg,
             draw_label(c, "toys (%.1f TeV)" % style.SQRT_S_TEV)]
    save(c, outbase, pdf)


def plot_summary(results, tag, r_inj, threshold, outbase, pdf):
    """Median Bias per generating member of one tag's envelope."""
    plot_rows([(None, results)], f"{tag}, #mu_{{inj}} = {r_inj:g}",
              threshold, outbase, pdf)


def plot_summary_bins(by_bin, title, r_inj, threshold, outbase, pdf):
    """Median Bias per generating member, one block of rows per BDT bin.
    by_bin maps a bin label to that bin's report summary; each bin has its
    own envelope, so the members differ from block to block."""
    groups = [(style.BIN_LABEL.get(b, b), summary["members"])
              for b, summary in by_bin.items()]
    plot_rows(groups, f"{title}, #mu_{{inj}} = {r_inj:g}",
              threshold, outbase, pdf)


def collect_member(outdir, prefix, r_inj, rrange, index_name, requested,
                   nmembers, threshold, metric):
    """Summary of the fit files fit_<prefix>_s*.root (one truth, and one fit
    function in --mode single), with the toys of jobs that wrote nothing
    counted as dropped. Returns (summary, pulls)."""
    files = sorted(glob.glob(os.path.join(outdir, f"fit_{prefix}_s*.root")))
    pulls, dr, chosen, drop = collect_pulls(files, r_inj, rrange,
                                            index_name, requested)
    # jobs that never wrote a fit file (combine crashed): their toys count
    # as dropped, so the valid fraction is out of all requested
    have = {os.path.basename(p)[len("fit_"):-len(".root")] for p in files}
    failed = sorted(st for st in requested
                    if st.startswith(f"{prefix}_s") and st not in have)
    drop["job_failed"] = sum(requested[st] for st in failed)
    s = summarize(pulls, dr, chosen, drop, nmembers, threshold, metric)
    s["n_files"] = len(files)
    s["failed_jobs"] = failed
    return s, pulls


def read_meta(tag, args):
    outdir = tag_outdir(tag, args)
    meta_path = os.path.join(outdir, "meta.json")
    if not os.path.isfile(meta_path):
        print(f"\n=== {tag} === no {meta_path}; run the toys first")
        return outdir, None
    return outdir, read_json(meta_path)


def report_tag(tag, args):
    if args.mode == "single":
        return report_single(tag, args)
    outdir, meta = read_meta(tag, args)
    if meta is None:
        return {}
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
            s, pulls = collect_member(outdir, f"gen{igen}_{key}", r_inj,
                                      rrange, index_name,
                                      meta.get("ntoys", {}), len(names),
                                      args.threshold, args.metric)
            results[gen] = s
            if pulls:
                plot_member(pulls, s, tag, gen, r_inj,
                            os.path.join(plotdir, f"bias_{tag}_{key}_{gen}"),
                            args.pdf)
        if any("bias" in s for s in results.values()):
            plot_summary(results, tag, r_inj, args.threshold,
                         os.path.join(plotdir, f"bias_{tag}_{key}_summary"),
                         args.pdf)
        summary = {"tag": tag, "mode": "envelope", "r_inj": r_inj,
                   "r_range": rrange, "threshold": args.threshold,
                   "metric": args.metric, "envelope": names,
                   "pull": "Bias = (mu_fit - mu_inj) / sigma_mu, mu = r, "
                           "sigma_mu on the side facing mu_inj",
                   "members": results,
                   "pass": all(s.get("pass", False) for s in results.values())}
        write_json(os.path.join(outdir, f"bias_{tag}_{key}.json"), summary)
        print_table(summary)
        out[key] = summary
    return out


def report_single(tag, args):
    """The Run 2 note's matrix: Bias of every fit function (pdfindex frozen)
    on the toys of every truth. A function is usable if |Bias| < threshold
    against every truth; of those, the one with the best expected
    significance is chosen."""
    outdir, meta = read_meta(tag, args)
    if meta is None:
        return {}
    names = [short_name(n, tag) for n in meta["envelope"]]
    fits = meta.get("fits", list(range(len(names))))
    category = tag.split("_")[0]
    plotdir = os.path.join(args.plotdir, category)
    os.makedirs(plotdir, exist_ok=True)
    signif = {names[j]: read_significance(
        os.path.join(outdir, f"sig_fit{j}.root")) for j in fits}

    out = {}
    for r_inj in args.inject:
        key = rlabel(r_inj)
        if key not in meta["rrange"]:
            print(f"\n=== {tag}, r = {r_inj:g} === no toys")
            continue
        rrange = meta["rrange"][key]
        matrix = {}               # fit -> truth -> summary
        for j in fits:
            row = {}
            for igen, gen in enumerate(names):
                s, _ = collect_member(outdir, f"gen{igen}_fit{j}_{key}",
                                      r_inj, rrange, None,
                                      meta.get("ntoys", {}), len(names),
                                      args.threshold, args.metric)
                s.pop("chosen_index", None)     # the index is frozen
                row[gen] = s
            matrix[names[j]] = row
        functions = {}
        for fit, row in matrix.items():
            vals = [abs(s["bias"]) for s in row.values() if "bias" in s]
            complete = all("bias" in s for s in row.values())
            functions[fit] = {
                "max_abs_bias": max(vals) if vals else None,
                "worst_truth": (max((g for g in row if "bias" in row[g]),
                                    key=lambda g: abs(row[g]["bias"]))
                                if vals else None),
                "pass": complete and all(s["pass"] for s in row.values()),
                "significance": signif.get(fit)}
        usable = [f for f, v in functions.items() if v["pass"]]
        ranked = sorted((f for f in usable
                         if functions[f]["significance"] is not None),
                        key=lambda f: -functions[f]["significance"])
        chosen = ranked[0] if ranked else None
        summary = {"tag": tag, "mode": "single", "r_inj": r_inj,
                   "r_range": rrange, "threshold": args.threshold,
                   "metric": args.metric, "truths": names,
                   "fits": [names[j] for j in fits],
                   "pull": "Bias = (mu_fit - mu_inj) / sigma_mu, mu = r, "
                           "sigma_mu on the side facing mu_inj; each fit "
                           "with pdfindex frozen to one function",
                   "matrix": matrix, "functions": functions,
                   "usable": usable, "chosen": chosen}
        write_json(os.path.join(outdir, f"bias_{tag}_{key}_matrix.json"),
                   summary)
        print_matrix(summary)
        plot_matrix(summary, os.path.join(
            plotdir, f"bias_{tag}_{key}_matrix"), args.pdf)
        out[key] = summary
    return out


def plot_matrix(summary, outbase, pdf):
    """Truth x fit matrix of the Bias, one cell per pair: the value in the
    cell, light red where |Bias| >= threshold. On the right, per fit
    function: its worst |Bias|, its expected significance, and which one is
    chosen. Pixel fonts, so the canvas grows with the matrix."""
    truths, fits = summary["truths"], summary["fits"]
    thr, metric = summary["threshold"], summary["metric"]
    nx, ny = len(truths), len(fits)
    cw, ch = 64, 30                                  # cell, pixels
    band = 48                    # CMS label, as on the 600 px canvases
    left, right, top, bottom = 120, 210, band + 40, 120
    w, h = left + right + cw * nx, top + bottom + ch * ny
    style.setup_style()
    c = ROOT.TCanvas(f"c_{os.path.basename(outbase)}", "", w, h)
    c.SetLeftMargin(left / w)
    c.SetRightMargin(right / w)
    c.SetTopMargin(top / h)
    c.SetBottomMargin(bottom / h)
    frame = ROOT.TH2F(f"f_{os.path.basename(outbase)}", "",
                      nx, 0, nx, ny, 0, ny)
    frame.SetDirectory(0)
    frame.SetStats(0)
    for i, t in enumerate(truths):
        frame.GetXaxis().SetBinLabel(i + 1, t)
    for k, f in enumerate(fits):
        frame.GetYaxis().SetBinLabel(ny - k, f)       # first fit on top
    for ax in (frame.GetXaxis(), frame.GetYaxis()):
        ax.SetLabelFont(43)
        ax.SetLabelSize(14)
        ax.SetTickLength(0)
    frame.GetXaxis().LabelsOption("v")
    frame.Draw("AXIS")
    keep = [c, frame]
    ok_col = ROOT.TColor.GetColor(*MATRIX_PASS_RGB)
    bad_col = ROOT.TColor.GetColor(*MATRIX_FAIL_RGB)
    txt = ROOT.TLatex()
    txt.SetTextFont(43)
    txt.SetTextSize(13)
    txt.SetTextAlign(22)
    for k, f in enumerate(fits):
        y = ny - k - 1
        for i, t in enumerate(truths):
            s = summary["matrix"][f][t]
            box = ROOT.TBox(i + 0.03, y + 0.06, i + 0.97, y + 0.94)
            box.SetFillColor(ok_col if s.get("pass") else
                             bad_col if "bias" in s else ROOT.kWhite)
            box.Draw()
            keep.append(box)
            txt.DrawLatex(i + 0.5, y + 0.5,
                          "%+.2f" % s["bias"] if "bias" in s else "-")
    frame.Draw("AXIS same")

    # per fit function, right of the matrix
    def ndc_y(y):
        return c.GetBottomMargin() + y / float(ny) * (
            1 - c.GetBottomMargin() - c.GetTopMargin())
    col = ROOT.TLatex()
    col.SetNDC()
    col.SetTextFont(43)
    col.SetTextSize(13)
    col.SetTextAlign(12)
    x0 = 1 - c.GetRightMargin() + 0.012
    xs = [x0, x0 + 62.0 / w, x0 + 120.0 / w]
    col.SetTextColor(ROOT.kGray + 2)
    for x, head in zip(xs, ("max |Bias|", "Z_{exp}", "")):
        col.DrawLatex(x, ndc_y(ny + 0.6), head)
    for k, f in enumerate(fits):
        v = summary["functions"][f]
        y = ndc_y(ny - k - 0.5)
        col.SetTextColor(ROOT.kBlack)
        col.DrawLatex(xs[0], y, "%.2f" % v["max_abs_bias"]
                      if v["max_abs_bias"] is not None else "-")
        col.DrawLatex(xs[1], y, "%.2f" % v["significance"]
                      if v["significance"] is not None else "-")
        label = ("#bf{chosen}" if f == summary["chosen"] else
                 "usable" if v["pass"] else "fails")
        col.SetTextColor(ROOT.kBlack if v["pass"] else
                         ROOT.TColor.GetColor(*PULL_GAUS_RGB))
        col.DrawLatex(xs[2], y, label)
    lat = ROOT.TLatex()
    lat.SetNDC()
    lat.SetTextFont(43)
    lat.SetTextSize(15)
    lat.SetTextAlign(12)
    lat.DrawLatex(c.GetLeftMargin(), 1 - (band + 20.0) / h,
                  f"{summary['tag']}, #mu_{{inj}} = {summary['r_inj']:g}: "
                  f"{metric} Bias, red where |Bias| #geq {thr:g}")
    # axis titles by hand: pixel-font title offsets do not scale
    lat.SetTextSize(14)
    lat.SetTextColor(ROOT.kGray + 2)
    lat.SetTextAlign(21)
    lat.DrawLatex(c.GetLeftMargin() + 0.5 * cw * nx / w, 8.0 / h,
                  "toys generated from (truth)")
    lat.SetTextAlign(32)
    lat.DrawLatex(c.GetLeftMargin() - 8.0 / w, ndc_y(ny + 0.6), "fitted with")
    keep += [txt, col, lat, draw_label(c, "toys (%.1f TeV)" % style.SQRT_S_TEV,
                                       t=band / float(h))]
    save(c, outbase, pdf)


def print_matrix(summary):
    truths, fits = summary["truths"], summary["fits"]
    print(f"\n=== {summary['tag']}, r_inj = {summary['r_inj']:g}, "
          f"single functions, {summary['metric']} Bias, "
          f"|Bias| < {summary['threshold']:g} ===")
    wid = max(9, max(len(t) for t in truths) + 1)
    print("  " + "fit / truth".ljust(11)
          + "".join(t.rjust(wid) for t in truths)
          + "max|B|".rjust(8) + "Z_exp".rjust(7) + "  status")
    for f in fits:
        row = summary["matrix"][f]
        v = summary["functions"][f]
        cells = "".join(
            ("%+.2f%s" % (row[t]["bias"], "" if row[t]["pass"] else "*")
             if "bias" in row[t] else "-").rjust(wid) for t in truths)
        mb = "%.2f" % v["max_abs_bias"] if v["max_abs_bias"] is not None else "-"
        z = "%.3f" % v["significance"] if v["significance"] is not None else "-"
        st = ("CHOSEN" if f == summary["chosen"] else
              "usable" if v["pass"] else f"fails (worst: {v['worst_truth']})")
        print(f"  {f:<11}{cells}{mb:>8}{z:>7}  {st}")
    low = [(f, t, s["n_valid"], s["n_toys"]) for f in fits
           for t, s in summary["matrix"][f].items()
           if s["n_toys"] and s["n_valid"] < MIN_VALID_FRAC * s["n_toys"]]
    for f, t, nv, nt in low:
        print(f"  !! {f} on {t} toys: only {nv}/{nt} valid fits")
    failed = sum(len(s.get("failed_jobs", [])) for f in fits
                 for s in summary["matrix"][f].values())
    if failed:
        print(f"  !! {failed} fit job(s) wrote nothing; rerun the same "
              "command to retry them")
    print("  * = |Bias| above the threshold")
    print("  chosen: " + (summary["chosen"] or
                          "none (no function passes against every truth)"))


def print_table(summary):
    names = summary["envelope"]
    print(f"\n=== {summary['tag']}, r_inj = {summary['r_inj']:g}, "
          f"r in [{summary['r_range'][0]:.3g}, {summary['r_range'][1]:.3g}]"
          f", |{summary.get('metric', 'median')}| < "
          f"{summary['threshold']:g} ===")
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
    ap.add_argument("--mode", choices=("envelope", "single"),
                    default="envelope",
                    help="envelope (default): toys from each envelope member, "
                         "fitted with pdfindex floating. single: the Run 2 "
                         "note, toys from every family's F-test pick, each "
                         "fitted with one function at a time; reads the "
                         "buildWorkspace.py --members families cards")
    ap.add_argument("--fit-models", nargs="+", default=None, metavar="NAME",
                    help="--mode single: fit with these members only "
                         "(default all); every member stays a truth")
    ap.add_argument("--metric", choices=("median", "mean"), default=None,
                    help="the quoted Bias: median (default for envelope, as "
                         "flashggFinalFit) or mean (default for single, as "
                         "the Run 2 note)")
    ap.add_argument("--carddir", default=None,
                    help="where buildWorkspace.py wrote datacard_<tag>.txt "
                         "(datacard_<tag>_families.txt for --mode single; "
                         "default WS_<year>)")
    ap.add_argument("--outdir", default="BIAS",
                    help="toys, fits and summaries go to <outdir>/<tag>/ "
                         "(<outdir>/<tag>_single/ for --mode single)")
    ap.add_argument("--plotdir",
                    default=os.path.expanduser("~/public_html/HmumuFits/bias"),
                    help="plots land in <plotdir>/<category>/")
    ap.add_argument("--inject", nargs="+", type=float, default=[1.0],
                    help="injected signal strength(s) (default 1)")
    ap.add_argument("--ntoys", type=int, default=1000,
                    help="toys per truth and r_inj (default 1000; the "
                         "Run 2 note used 2000)")
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
                    help=f"max |Bias| (default {PULL_THRESHOLD})")
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
    if args.metric is None:
        args.metric = DEFAULT_METRIC[args.mode]
    if args.fit_models and args.mode != "single":
        ap.error("--fit-models needs --mode single")
    if args.carddir is None:
        args.carddir = os.path.join(HERE, config.default_wsdir(args.year))

    tags = [f"{c}_{b}_{args.year}" for c in args.category for b in args.bins]
    if args.only != "report":
        for tag in tags:
            run_tag(tag, args)
    if args.only != "run" and not args.dry_run:
        reports = {tag: report_tag(tag, args) for tag in tags}
        if len(args.bins) > 1 and args.mode == "envelope":
            for cat in args.category:
                for r_inj in args.inject:
                    key = rlabel(r_inj)
                    by_bin = {b: reports[f"{cat}_{b}_{args.year}"][key]
                              for b in args.bins
                              if key in reports[f"{cat}_{b}_{args.year}"]}
                    if len(by_bin) > 1:
                        plot_summary_bins(
                            by_bin, f"{cat}, {args.year}", r_inj,
                            args.threshold,
                            os.path.join(args.plotdir, cat,
                                         f"bias_{cat}_{args.year}_{key}"
                                         "_summary_bins"), args.pdf)


if __name__ == "__main__":
    main()
