"""
submitBias.py -- biasStudy.py on Slurm (SubMIT): one job per test.

    python submitBias.py -c ggHcat -b bdt0 bdt1 bdt2 bdt3 bdt4 -y Run3 \
        --carddir WS_Run3_5bins --outdir BIAS_5bins_v2 --ntoys 1000

A test is one truth and one fit function in one bin (--mode single, the
default), or one truth with the envelope floating (--mode envelope). Three
stages, chained with Slurm dependencies:

  1. prepare  one job per bin: biasStudy.py --only prepare. The workspace,
              the r range (Asimov), the expected significances and every
              toy (GenerateOnly, in --toys-per-job chunks); no fits.
  2. fit      one array task per test: biasStudy.py --only fit --truth I
              --fit J. All --ntoys toys of truth I, fitted with J. Starts
              when every prepare job has succeeded.
  3. report   one job: biasStudy.py --only report, the matrices and the
              plots. Starts when the array has ended, however its tasks did.

The tests come from each bin's Bkg_<tag>_params.json in --carddir: with
--members families the multipdf has one member per family_choice entry,
so a bin with n families has n x n tests (n for --mode envelope). A task
whose member does not exist exits at once.

Everything is resumable as biasStudy.py is: a fit whose output exists is
skipped, so submitting again redoes only what is missing (prepare reuses
its toys). The job scripts, the task list and the Slurm logs go to
<outdir>/slurm/. --dry-run writes them and prints the sbatch commands
without submitting.

Runs on the submit node in any environment (no ROOT needed); the jobs run
biasStudy.py inside the Combine container (--container).
"""

import argparse
import json
import os
import shlex
import subprocess
import sys

from fitlib import config

HERE = os.path.dirname(os.path.abspath(__file__))
CONTAINER = ("/work/submit/jaeyserm/software/docker/"
             "combine-standalone_v9.2.1.sif")
BINDS = ["/work", "/ceph"]


def members(carddir, tag, mode):
    """How many members the bin's multipdf has, from bkgFit.py's JSON"""
    path = os.path.join(carddir, f"Bkg_{tag}_params.json")
    if not os.path.isfile(path):
        raise SystemExit(f"missing {path}; run bkgFit.py (and "
                         "buildWorkspace.py) for this bin first")
    with open(path) as fh:
        rec = json.load(fh)
    return len(rec["family_choice"] if mode == "single" else rec["envelope"])


def bias_args(args, cats, bins):
    """The biasStudy.py options every stage shares, shell-quoted; `cats`
    and `bins` go in as given, so a job's shell variables ("$CAT") stay
    unquoted and expand there."""
    out = ["--mode", args.mode, "-y", args.year, "--carddir", args.carddir, "--outdir", args.outdir,
           "--inject", *("%g" % r for r in args.inject),
           "--ntoys", str(args.ntoys), "--toys-per-job",
           str(args.toys_per_job), "--seed", str(args.seed)]
    if args.plotdir:
        out += ["--plotdir", args.plotdir]
    if args.combine_opts:
        out += ["--combine-opts", args.combine_opts]
    if args.r_width is not None:
        out += ["--r-width", "%g" % args.r_width]
    if args.rrange:
        out += ["--rrange", *("%g" % r for r in args.rrange)]
    return " ".join(["-c", *cats, "-b", *bins] + [shlex.quote(o) for o in out])


def in_container(opts, args, tail=""):
    """The shell lines that run biasStudy.py `opts` (a string) from
    parametricfit/ inside the Combine container"""
    binds = " ".join(f"-B {b}" for b in BINDS)
    return [f"cd {shlex.quote(HERE)}",
            f"singularity exec --cleanenv {binds} "
            f"{shlex.quote(args.container)} python -u biasStudy.py {opts}"
            + (f" {tail}" if tail else "")]


def script(path, name, logs, args, body, time, cpus=1, array=None):
    """Write an sbatch script; returns its path"""
    lines = ["#!/bin/bash",
             f"#SBATCH --job-name={name}",
             f"#SBATCH --partition={args.partition}",
             f"#SBATCH --time={time}",
             f"#SBATCH --mem-per-cpu={args.mem}",
             f"#SBATCH --cpus-per-task={cpus}",
             f"#SBATCH --output={logs}/{name}_%A" + ("_%a" if array else "")
             + ".log"]
    if array:
        lines.append(f"#SBATCH --array={array}")
    lines += ["set -e", 'echo "host $(hostname), $(date)"', *body,
              'echo "done $(date)"', ""]
    with open(path, "w") as fh:
        fh.write("\n".join(lines))
    os.chmod(path, 0o755)
    return path


def sbatch(path, args, after=None):
    """Submit; returns the job id (a placeholder with --dry-run)"""
    cmd = ["sbatch", "--parsable"]
    if after:
        cmd.append(f"--dependency={after}")
    cmd.append(path)
    print("$ " + " ".join(cmd))
    if args.dry_run:
        return f"<{os.path.basename(path)}>"
    out = subprocess.run(cmd, check=True, capture_output=True, text=True)
    return out.stdout.strip().split(";")[0]


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--category", nargs="+", default=["ggHcat"])
    ap.add_argument("-b", "--bins", nargs="+", required=True)
    ap.add_argument("-y", "--year", default="Run3")
    ap.add_argument("--mode", choices=("single", "envelope"),
                    default="single")
    ap.add_argument("--carddir", default=None,
                    help="the cards and fit JSONs (default WS_<year>)")
    ap.add_argument("--outdir", default="BIAS")
    ap.add_argument("--plotdir", default=None)
    ap.add_argument("--inject", nargs="+", type=float, default=[1.0])
    ap.add_argument("--ntoys", type=int, default=1000,
                    help="toys per test (default 1000)")
    ap.add_argument("--toys-per-job", type=int, default=100,
                    help="toys per combine call and file (default 100)")
    ap.add_argument("--seed", type=int, default=12345)
    ap.add_argument("--combine-opts", default="")
    ap.add_argument("--r-width", type=float, default=None, metavar="N",
                    help="biasStudy.py --r-width: the r range's half-width "
                         "in sigma_r (its default if not given)")
    ap.add_argument("--rrange", nargs=2, type=float, default=None,
                    metavar=("RMIN", "RMAX"),
                    help="biasStudy.py --rrange: a fixed r range instead")
    ap.add_argument("--partition", default="submit")
    ap.add_argument("--time-prepare", default="04:00:00")
    ap.add_argument("--time-fit", default="12:00:00",
                    help="per test (default 12 h: the FEWZ fits are slow)")
    ap.add_argument("--time-report", default="01:00:00")
    ap.add_argument("--mem", default="2G", help="per CPU (default 2G)")
    ap.add_argument("--prepare-cpus", type=int, default=8,
                    help="combine processes at once in a prepare job")
    ap.add_argument("--max-running", type=int, default=None,
                    help="at most this many fit tasks at once (array %%N)")
    ap.add_argument("--container", default=CONTAINER)
    ap.add_argument("--no-report", action="store_true")
    ap.add_argument("--dry-run", action="store_true",
                    help="write the scripts, print the sbatch commands, "
                         "submit nothing")
    args = ap.parse_args()
    # absolute, from parametricfit/ (where the jobs run), as biasStudy.py
    args.carddir = os.path.join(HERE, args.carddir
                                or config.default_wsdir(args.year))
    args.outdir = os.path.join(HERE, args.outdir)
    if args.plotdir:
        args.plotdir = os.path.abspath(os.path.expanduser(args.plotdir))

    slurm = os.path.join(args.outdir, "slurm")
    logs = os.path.join(slurm, "logs")
    os.makedirs(logs, exist_ok=True)
    # the tests
    tasks = []
    for cat in args.category:
        for b in args.bins:
            n = members(args.carddir, f"{cat}_{b}_{args.year}", args.mode)
            fits = range(n) if args.mode == "single" else [None]
            tasks += [(cat, b, i, j) for i in range(n) for j in fits]
    task_file = os.path.join(slurm, "tasks.txt")
    with open(task_file, "w") as fh:
        for cat, b, i, j in tasks:
            fh.write(f"{cat} {b} {i} {'-' if j is None else j}\n")
    print(f"{len(tasks)} test(s) over {len(args.category) * len(args.bins)} "
          f"bin(s), {args.ntoys} toys each; task list {task_file}")

    # 1. prepare, one job per bin
    prep_ids = []
    for cat in args.category:
        for b in args.bins:
            opts = bias_args(args, [cat], [b])
            path = script(os.path.join(slurm, f"prepare_{cat}_{b}.sh"),
                          f"bias_prepare_{cat}_{b}", logs, args,
                          in_container(opts, args, "--only prepare -j "
                                       f"{args.prepare_cpus}"),
                          args.time_prepare, cpus=args.prepare_cpus)
            prep_ids.append(sbatch(path, args))

    # 2. fit, one array task per test
    body = [f'read CAT BIN I J <<< "$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" '
            f'{shlex.quote(task_file)})"',
            'FITOPT=""',
            'if [ "$J" != "-" ]; then FITOPT="--fit $J"; fi',
            'echo "test: $CAT $BIN truth $I fit $J"',
            *in_container(bias_args(args, ['"$CAT"'], ['"$BIN"']), args,
                          '--only fit --truth "$I" $FITOPT')]
    array = f"0-{len(tasks) - 1}" + (f"%{args.max_running}"
                                     if args.max_running else "")
    path = script(os.path.join(slurm, "fit.sh"), "bias_fit", logs, args,
                  body, args.time_fit, array=array)
    fit_id = sbatch(path, args, after="afterok:" + ":".join(prep_ids))

    # 3. report
    if not args.no_report:
        path = script(os.path.join(slurm, "report.sh"), "bias_report", logs,
                      args, in_container(bias_args(args, args.category,
                                                   args.bins), args,
                                         "--only report"), args.time_report)
        sbatch(path, args, after=f"afterany:{fit_id}")

    print(f"\nlogs: {logs}/\nwatch: squeue -u $USER\n"
          "resubmit the same command to redo what is missing")


if __name__ == "__main__":
    sys.exit(main())
