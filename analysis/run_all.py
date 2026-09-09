#!/usr/bin/env python3
"""
run_all.py — run Hmm.py one sample per process instead of all in sequence.

Samples are independent and each writes its own snapshot, so this
parallelises both the event loop and the per-RDataFrame JIT compilation
(which a single process cannot, since Cling is serialised).

    # the sample list comes from datasets/VBF_2024.txt unless you say otherwise
    python run_all.py 2024 isVBF --dry-run

    # split the list into 8 balanced jobs instead of one job per sample
    python run_all.py 2024 isVBF -n 8 --slurm slurm/vbf_2024.sh

    # Slurm job array: the right backend on SubMIT.
    # Bare --slurm writes slurm/<category>_<year>.sh
    python run_all.py 2024 isVBF --slurm --ncores 2
    sbatch slurm/VBF_2024.sh

    # ...or generate and submit in one go
    python run_all.py 2024 isVBF --slurm slurm/vbf_2024.sh --ncores 2 --submit

    # smoke test on the login node (NOT for real work)
    python run_all.py 2024 isVBF -s 10,11 -j 2 --ncores 2

    # rerun only what is missing
    python run_all.py 2024 isVBF --slurm slurm/retry.sh --skip-existing

    # an alternative list, or an ad-hoc selection
    python run_all.py 2024 isVBF datasets/VBF_2024_test.txt --dry-run
    python run_all.py 2024 isVBF -s signal_hmm --dry-run

Sample lists live in datasets/ and are committed to git. Regenerate them
only when the datasets change:  python tools/datasets.py --write-all

HTCondor is deliberately not supported: its worker nodes cannot see /ceph,
/scratch, /work or /home, so the inputs and the conda environment would have
to be shipped in via CVMFS or file transfer. Slurm workers see all of them.
"""

import argparse
import getpass
import os
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import tools.datasets as datasets
from tools.datasets import (resolve_ids, sample_label, read_list, list_path,
                            CATEGORY, MODE_MAP, VALID_YEARS, VALID_MODES)
from tools.utilsAna import SwitchSample

HERE = Path(__file__).resolve().parent
HMM = HERE / "Hmm.py"

def parse_args():
    p = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("year", choices=VALID_YEARS)
    p.add_argument("mode", choices=VALID_MODES)
    p.add_argument("samplelist", nargs="?", default=None,
                   help="sample list to run (default: the committed "
                        "datasets/<category>_<year>.txt)")
    p.add_argument("-o", "--outdir", default=None, help="passed through to Hmm.py -o")
    p.add_argument("-s", "--samples", default=None,
                   help="restrict to this selection: ids (11, -41), group names from "
                        "samples.yaml (signal_hmm, vv), 'mc', 'data', or globs on the "
                        "dataset name ('DYto2Mu*'). Comma-separated.")
    p.add_argument("-l", "--list", action="store_true",
                   help="print the id / dataset / group table and exit")
    p.add_argument("--mc-only", action="store_true")
    p.add_argument("--data-only", action="store_true")
    p.add_argument("--stream", default="prompt", choices=["prompt", "parking"])
    p.add_argument("-n", "--split", type=int, default=0, metavar="K",
                   help="split the sample list into K jobs, balanced by file count "
                        "(default 0 = one job per sample)")
    p.add_argument("--ncores", type=int, default=2,
                   help="IMT threads per Hmm.py process (default: %(default)s)")
    p.add_argument("--skip-existing", action="store_true",
                   help="drop samples whose snapshot already exists")
    p.add_argument("--logdir", default="logs")
    p.add_argument("--dry-run", action="store_true", help="print the plan and exit")

    g = p.add_argument_group("local mode (login node, smoke tests only)")
    g.add_argument("-j", "--jobs", type=int, default=2,
                   help="concurrent processes (default: %(default)s). "
                        "Keep jobs x ncores well below the core count.")

    g = p.add_argument_group("Slurm mode")
    g.add_argument("--slurm", nargs="?", const="AUTO", metavar="FILE",
                   help="write a job-array script instead of running. With no value, "
                        "defaults to slurm/<category>_<year>.sh")
    g.add_argument("--submit", action="store_true",
                   help="run sbatch on the generated script straight away")
    g.add_argument("--partition", default="submit")
    g.add_argument("--time", default="04:00:00",
                   help="walltime per task, max 6 days (default: %(default)s)")
    g.add_argument("--mem-per-cpu", type=int, default=4000, metavar="MB",
                   help="memory per cpu (default: %(default)s). Measure a real "
                        "sample with /usr/bin/time -v before trusting this.")
    g.add_argument("--max-concurrent", type=int, default=0, metavar="N",
                   help="cap simultaneously running tasks (0 = uncapped)")
    g.add_argument("--conda-env", default=os.environ.get("CONDA_DEFAULT_ENV", "pyenv"),
                   help="env to activate on the worker (default: %(default)s). "
                        "Slurm does NOT inherit the submitting environment.")
    g.add_argument("--conda-root", default=None,
                   help="conda prefix; auto-detected from CONDA_EXE")
    return p.parse_args()


def outdir_for(args):
    base = args.outdir or f"/work/submit/{getpass.getuser()}/HmumuRun3/ROOTFILES/"
    return base if base.endswith("/") else base + "/"


def snapshot_path(args, sid):
    cat = MODE_MAP[args.mode]
    return f"{outdir_for(args)}{cat}/snapshot_mc_{sid}_{args.year}_{cat}.root"


def build_units(args):
    """Enumerate work units for the array, largest first.

    A unit is a list of sample ids that one Hmm.py process will handle. With
    --split K the samples are packed into K balanced bundles; otherwise each
    sample is its own unit.

    Resolving the file lists here is the one place the driver touches the
    filesystem; the counts are used to order and balance the work.
    """
    thisdict = datasets.BuildDict(args.year)

    if args.samples:
        ids = resolve_ids(args.year, args.mode, args.samples, args.stream)
        source = f"-s {args.samples}"
    else:
        source = args.samplelist or list_path(args.year, args.mode)
        ids = read_list(source, args.year, args.mode, args.stream)

    if args.mc_only:
        ids = [i for i in ids if i > 0]
    if args.data_only:
        ids = [i for i in ids if i < 0]

    print(f"sample list: {source}  ->  {len(ids)} samples")

    sized = []
    for sid in ids:
        if args.skip_existing and os.path.exists(snapshot_path(args, sid)):
            print(f"  skip {sid}: output exists")
            continue
        n = len(SwitchSample(thisdict, sid)[0])
        if n == 0:
            print(f"  skip {sid}: no files")
            continue
        sized.append((sid, n))

    sized.sort(key=lambda u: -u[1])          # heaviest first

    if not args.split or args.split >= len(sized):
        return [([sid], n) for sid, n in sized]

    # longest-processing-time-first packing into K bundles
    bundles = [([], 0) for _ in range(args.split)]
    for sid, n in sized:
        k = min(range(len(bundles)), key=lambda j: bundles[j][1])
        bundles[k] = (bundles[k][0] + [sid], bundles[k][1] + n)
    units = [(b, w) for b, w in bundles if b]
    units.sort(key=lambda u: -u[1])
    return units


def hmm_args(args, ids):
    out = [args.year, args.mode, "-s", ",".join(str(i) for i in ids),
           "--ncores", str(args.ncores)]
    if args.outdir:
        out += ["-o", args.outdir]
    return out


def unit_tag(args, ids, k):
    return f"{args.year}_{args.mode}_" + (str(ids[0]) if len(ids) == 1 else f"part{k:02d}")


def run_one(args, ids, k, logdir):
    log = Path(logdir) / f"{unit_tag(args, ids, k)}.log"
    t0 = time.time()
    with open(log, "w") as fh:
        rc = subprocess.call([sys.executable, str(HMM)] + hmm_args(args, ids),
                             stdout=fh, stderr=subprocess.STDOUT, cwd=str(HERE))
    return ids, rc, time.time() - t0, log


def conda_root(args):
    if args.conda_root:
        return args.conda_root
    exe = os.environ.get("CONDA_EXE")   # /work/submit/USER/miniforge3/bin/conda
    if exe:
        return str(Path(exe).parent.parent)
    return str(Path.home() / "miniforge3")


def write_slurm(args, units, path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    joblist = path.with_suffix(".jobs")

    with open(joblist, "w") as fh:
        for ids, weight in units:
            fh.write(" ".join(hmm_args(args, ids))
                     + f"   # {len(ids)} sample(s), ~{weight} files\n")

    cap = f"%{args.max_concurrent}" if args.max_concurrent else ""
    logdir = Path(args.logdir).resolve()
    logdir.mkdir(parents=True, exist_ok=True)
    tag = f"{args.year}_{args.mode}"

    with open(path, "w") as fh:
        fh.write(f"""#!/bin/bash
#
# Generated by run_all.py -- one array task per work unit.
#   sbatch {path.name}
#
#SBATCH --job-name=hmm_{tag}
#SBATCH --partition={args.partition}
#SBATCH --array=1-{len(units)}{cap}
#SBATCH --cpus-per-task={args.ncores}
#SBATCH --mem-per-cpu={args.mem_per_cpu}
#SBATCH --time={args.time}
#SBATCH --output={logdir}/{tag}_%A_%a.out
#SBATCH --error={logdir}/{tag}_%A_%a.err

set -eo pipefail

# Slurm does not inherit the submitting conda environment, so activate it here.
# -u must be off for this: conda's ROOT deactivate hook reads
# CONDA_BACKUP_ROOTSYS, which is unset, and that is fatal under set -u.
set +u
source {conda_root(args)}/etc/profile.d/conda.sh
conda activate {args.conda_env}
set -u

python -c "import ROOT, yaml" || {{ echo "environment {args.conda_env} is not usable" >&2; exit 1; }}

cd {HERE}

ARGS=$(sed -n "${{SLURM_ARRAY_TASK_ID}}p" {joblist.resolve()} | sed 's/#.*//')
if [ -z "$ARGS" ]; then
    echo "no work unit on line ${{SLURM_ARRAY_TASK_ID}} of {joblist.name}" >&2
    exit 1
fi

echo "host    : $(hostname)"
echo "task    : ${{SLURM_ARRAY_TASK_ID}} of {len(units)}"
echo "args    : $ARGS"
echo "started : $(date)"

/usr/bin/time -v python Hmm.py $ARGS

echo "finished: $(date)"
""")
    path.chmod(0o755)
    print(f"wrote {path} and {joblist} ({len(units)} tasks)")
    print(f"  cpus-per-task {args.ncores}, mem-per-cpu {args.mem_per_cpu} MB, "
          f"time {args.time}, env {args.conda_env}")
    if not args.submit:
        print(f"submit with:  sbatch {path}")
        print(f"monitor with: squeue -u $USER     then     seff <jobid>")
        return

    try:
        out = subprocess.run(["sbatch", str(path)], capture_output=True, text=True)
    except FileNotFoundError:
        sys.exit(f"sbatch not found; submit by hand:  sbatch {path}")
    if out.returncode != 0:
        sys.exit(f"sbatch failed ({out.returncode}):\n{out.stderr.strip()}\n"
                 f"the script is still there: {path}")
    print(out.stdout.strip())          # "Submitted batch job 12345"
    print("monitor with: squeue -u $USER     then     seff <jobid>")


def main():
    args = parse_args()
    Path(args.logdir).mkdir(parents=True, exist_ok=True)

    if args.list:
        return print_table(args)

    units = build_units(args)
    if not units:
        sys.exit("nothing to do")

    total = sum(len(ids) for ids, _ in units)
    print(f"{len(units)} job(s) covering {total} samples")

    if args.slurm:
        path = (f"slurm/{CATEGORY[args.mode]}_{args.year}.sh"
                if args.slurm == "AUTO" else args.slurm)
        return write_slurm(args, units, path)

    if args.dry_run:
        for k, (ids, n) in enumerate(units, 1):
            print(f"  [{n:>6} files] python Hmm.py " + " ".join(hmm_args(args, ids)))
        return

    print(f"local mode: {args.jobs} concurrent x {args.ncores} threads "
          f"(login node -- smoke tests only)")
    t0, failed = time.time(), []
    with ThreadPoolExecutor(max_workers=args.jobs) as pool:
        futures = [pool.submit(run_one, args, ids, j, args.logdir)
                   for j, (ids, _) in enumerate(units, 1)]
        for k, fut in enumerate(as_completed(futures), 1):
            ids, rc, dt, log = fut.result()
            status = "ok  " if rc == 0 else f"FAIL({rc})"
            what = str(ids[0]) if len(ids) == 1 else f"{len(ids)} samples"
            print(f"[{k}/{len(units)}] {status} {what}  {dt/60:.1f} min  -> {log}")
            if rc != 0:
                failed.append((ids, log))

    print(f"\ntotal wall time {(time.time()-t0)/60:.1f} min")
    if failed:
        print(f"{len(failed)} failed:")
        for ids, log in failed:
            print(f"  samples {ids}  see {log}")
        sys.exit(1)
    print("all samples completed")


if __name__ == "__main__":
    main()
