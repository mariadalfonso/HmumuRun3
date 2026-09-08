
## create the evironment
First activate your conda environment
```
conda activate pyenv
```

If is it the first time creating the environment, run
```
conda env create --name pyenv --file=environment.yml
```

## Layout

```
analysis/
  Hmm.py                  event loop, one category at a time
  run_all.py              parallel submission (Slurm on SubMIT)
  make_datasets.py        regenerates the datasets/*.txt sample lists
  config/
    samples.yaml          dataset definitions: paths, xsecs, BRs, groups
    selection.json        object selections
    ...                   functions.h, POG JSONs, cert JSONs
  datasets/
    VBF_2024.txt          committed sample lists, one per category+year
    ggH_2024.txt
    ...
  tools/
    datasets.py           reads samples.yaml: paths, xsecs, selections
    utilsAna.py           corrections, weights, data quality
    helper_tmva.py        MVA inference helpers
    validate_datasets.py  checks samples.yaml against the legacy datasets.py
```

The three scripts at the top level are the entry points you run;
everything in `tools/` is imported by them.

Run everything from the `analysis/` directory.

## Running `Hmm.py`

```
    python Hmm.py <year> <mode> [samplelist] [-o OUTDIR] [-s SELECTION] [-j NCORES]
```
- **year**: `12022`, `22022`, `12023`, `22023`, `2024`, `2025`, `2026`
- **mode**: `isVBF`, `isGGH`, `isZinv`, `isVlep`, `isVhad`, `isTTlep`, `isTThad`
- **samplelist**: which samples to process. Defaults to the committed
  `datasets/<category>_<year>.txt` for this mode, so you normally omit it.
- **-o/--outdir**: output dir (default `/work/submit/$USER/HmumuRun3/ROOTFILES/`)
- **-s/--samples**: override the list with an explicit selection (see below)
- **-j/--ncores**: implicit-MT threads; 0 = all cores

### Choosing what to run

There are three ways to say which samples to process. They take precedence
in this order: **`-s` beats a positional list file, which beats the
default.** `Hmm.py` prints which one it used at startup, e.g.
`sample list: datasets/VBF_2024.txt  ->  50 MC + 14 data`.

**1. Nothing — use the committed list (the normal case)**
```
python Hmm.py 2024 isVBF          # reads datasets/VBF_2024.txt
python Hmm.py 12022 isGGH         # reads datasets/ggH_12022.txt
```
The filename comes from the mode and year, so there is nothing to remember.
Missing file raises and tells you how to generate it.

**2. A different list file, as a positional argument**
```
python Hmm.py 2024 isVBF datasets/VBF_2024_test.txt
python Hmm.py 2024 isVBF /tmp/rerun_these.txt
```
Use this for a whole alternative list: a cut-down test set, or a rerun list
of the samples that failed. Note the list is *not* checked against the mode
— running `isVBF` with `datasets/ggH_2024.txt` will happily push those
samples through the VBF selection, so watch the startup line.

**3. `-s` for an explicit selection**
```
python Hmm.py 2024 isVBF -s 11                 # one sample id
python Hmm.py 2024 isVBF -s 10,11,-41          # several ids (negative = data)
python Hmm.py 2024 isVBF -s signal_hmm         # a group from samples.yaml
python Hmm.py 2024 isVBF -s vv,vvv             # several groups
python Hmm.py 2024 isVBF -s 'DYto2Mu*'         # glob on the dataset name
python Hmm.py 2024 isVBF -s mc                 # all MC for this mode
python Hmm.py 2024 isVBF -s data               # all data for this mode
python Hmm.py 2024 isVBF -s @datasets/VBF_2024.txt,141   # a list plus extras
```
Quote globs so the shell does not expand them first. `-s` **replaces** the
list rather than filtering it, so `-s 141` runs sample 141 even though it is
not in `datasets/VBF_2024.txt` — handy for the training samples that are not
in any category list. An unrecognised token raises rather than silently
matching nothing.

`run_all.py` accepts the same three forms, so anything above can be
submitted in parallel by swapping `Hmm.py` for `run_all.py` and adding
`--slurm FILE`.

Writes snapshots to `<outdir>/<category>/snapshot_mc_<mc>_<year>_<category>.root`.
One file per sample, so parallel jobs never collide.

### Selection syntax

Anywhere a selection is accepted (`-s`, or a line in a list file):

| token | meaning |
|---|---|
| `11`, `-41` | a sample id |
| `signal_hmm`, `vv` | a group from `groups:` in `samples.yaml` |
| `mc`, `data` | everything of that kind for this year and mode |
| `'DYto2Mu*'` | glob on the dataset name |
| `@path/to/list.txt` | another list file |

## Sample lists: `datasets/*.txt`

One file per category and year, generated from `config/samples.yaml` and
**committed to git**. Nobody needs to regenerate them unless the datasets
change.

```
# VBF 2024 -- sample list for isVBF
# generated 2026-09-08 03:41 from config/samples.yaml by tools/datasets.py

# --- MC (50) ---
10     # VBF*Hto2Mu_*M-125
11     # GluGlu*Hto2Mu_*M-125
...
# --- data (14) ---
-41    # Run2024C/Muon0
```

Comment a line out to skip that sample without touching `samples.yaml`.

### Regenerating

```
python make_datasets.py --write-all              # every category and year
python make_datasets.py --write 2024 isVBF       # just one
python make_datasets.py --show 2024 isVBF        # print, write nothing
```

Do this after adding a sample to `samples.yaml` or changing a group, then
commit the updated `datasets/` files.


## Running in parallel: `run_all.py`

Samples are independent and each writes its own snapshot, so one process per
sample parallelises both the event loop and the per-RDataFrame JIT compilation.

`run_all.py` is a warpper for sending multiple `Hmm.py` jobs. 
It creates the lists to pass to `Hmm.py` and runs them locally or writes a Slurm job array for you to submit.

```
python run_all.py <year> <mode> [samplelist] [options]
```

### How to run it

From `analysis/`, with `pyenv` active:

```
python run_all.py 2024 isVBF --dry-run            # what would run, largest first

python run_all.py 2024 isVBF --slurm slurm/VBF_2024.sh \
    --ncores 2 --mem-per-cpu 1500 --time 04:00:00 --max-concurrent 20
sbatch slurm/VBF_2024.sh

squeue -u $USER                                   # pending / running
seff <jobid>                                      # efficiency, once finished
```

Then rerun whatever is missing — `--skip-existing` queues only the gaps:

```
python run_all.py 2024 isVBF --slurm slurm/VBF_2024_retry.sh \
    --skip-existing --ncores 2 --mem-per-cpu 1500
sbatch slurm/VBF_2024_retry.sh
```

Those flag values are a reasonable starting point for `isVBF` 2024
(~1.1 GB per job measured, longest sample ~1 h). Read on to size them for
other categories, and see the notes below on what the generated script does.

### More detail

`--slurm FILE` writes two files and submits nothing:

- `slurm/VBF_2024.sh` — the `sbatch` script, one array task per work unit
- `slurm/VBF_2024.jobs` — the exact `Hmm.py` arguments for each task, which
  doubles as the record of what ran

Worth reading the `.sh` once: it is short, and shows the resource request,
the conda activation and the working directory.

**Before the first submission of a category**, measure one sample rather
than guessing at memory and walltime:

```
srun --partition=submit --cpus-per-task=2 --mem-per-cpu=4000 \
     --time=01:00:00 --pty bash
conda activate pyenv                  # Slurm does not inherit your env
cd ~/HmumuRun3/analysis
/usr/bin/time -v python Hmm.py 2024 isVBF -s 11 --ncores 2
exit
```

"Maximum resident set size" plus a cushion is your `--mem-per-cpu`; scale the
wall time by the largest sample's file count from `--dry-run`.

**Monitoring and cleanup:**

```
squeue -u $USER -t running            # just the running tasks
scancel <jobid>                       # kill the whole array
scancel <jobid>_5                     # kill one task
tail -20 logs/2024_isVBF_<jobid>_1.out
ls -l /work/submit/$USER/HmumuRun3/ROOTFILES/VBFcat/
```

Per-task stdout and stderr land in `logs/<year>_<mode>_<arrayid>_<task>.out`
and `.err`. Every task is wrapped in `/usr/bin/time -v`, so peak memory is
at the end of each log whether the task succeeded or not.

### Sizing the request

Keep `--ncores` small. Slurm on SubMIT hands out cores, not whole nodes, and
the fair-share system gives lower priority to users requesting more, so 50
jobs at 2 cores schedule sooner than 8 jobs at 16 cores. Small jobs also
parallelise the JIT compilation, which threads within one job cannot.

| option | meaning | note |
|---|---|---|
| `--ncores N` | IMT threads per job | maps to `--cpus-per-task`; 2 is a good start |
| `--mem-per-cpu MB` | memory per core | measure it first (see above); default 4000 is a guess |
| `--time HH:MM:SS` | walltime per task | max 6 days on `submit`; task is killed at the limit |
| `--partition` | Slurm partition | `submit` (default), `submit-gpu` for GPUs |
| `--max-concurrent N` | cap running tasks | becomes `--array=1-M%N` |
| `-n/--split K` | pack the list into K jobs | balanced by file count, longest first |

`-n/--split K` is for when one job per sample is too many tasks. With 64
samples, `-n 8` gives 8 jobs of roughly equal total file count instead of 64
tasks. Without it, each sample is its own task.

### Running locally

Omitting `--slurm` runs the jobs here, as subprocesses:
```
python run_all.py 2024 isVBF -s signal_hmm -j 4 --ncores 2
```
`-j` sets how many run at once. This is fine for a smoke test, but the
login nodes are shared — for real work either submit an array, or take an
interactive allocation and use local mode inside it:
```
srun --partition=submit --cpus-per-task=16 --mem-per-cpu=4000 \
     --time=04:00:00 --pty bash
conda activate pyenv && cd ~/HmumuRun3/analysis
python run_all.py 2024 isVBF -j 8 --ncores 2
```