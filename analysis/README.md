# Higgs to dimuon Run 3 analysis code

## Environment

```
conda activate pyenv
conda env create --name pyenv --file=environment.yml    # first time only
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
Run from the `analysis/` directory.

---

## Run one category: `Hmm.py`

```
python Hmm.py <year> <mode> [samplelist] [options]
```

| argument | meaning |
|---|---|
| `year` | `12022`, `22022`, `12023`, `22023`, `2024`, `2025`, `2026` |
| `mode` | `isVBF`, `isGGH`, `isZinv`, `isVlep`, `isVhad`, `isTTlep`, `isTThad` |
| `samplelist` | optional list file; defaults to `datasets/<category>_<year>.txt` |
| `-o, --outdir` | output dir (default `/work/submit/$USER/HmumuRun3/ROOTFILES/`) |
| `-s, --samples` | explicit selection, overrides the list file |
| `-j, --ncores` | implicit-MT threads; 0 = all cores |
| `--maxfiles N` | process only the first N files of each sample |

Writes `<outdir>/<category>/snapshot_mc_<mc>_<year>_<category>.root`, one
file per sample. Prints which selection it used at startup:
`sample list: datasets/VBF_2024.txt -> 50 MC + 14 data`.

### Choosing samples

Precedence: `-s` > positional list file > default list.

```
python Hmm.py 2024 isGGH                                 # datasets/ggH_2024.txt
python Hmm.py 2024 isGGH datasets/ggH_2024_test.txt      # another list
python Hmm.py 2024 isGGH -s 11                           # explicit selection
```

Selection tokens, comma-separated, usable in `-s` or as a line in a list file:

| token | meaning |
|---|---|
| `11`, `-41` | a sample id (negative = data) |
| `signal_hmm`, `vv` | a group from `groups:` in `samples.yaml` |
| `mc`, `data` | everything of that kind for this year and mode |
| `'DYto2Mu*'` | glob on the dataset name (quote it) |
| `@path/to/list.txt` | another list file |

`-s` replaces the list rather than filtering it, so `-s 141` runs a sample
that is in no category list. Unrecognised tokens raise. A positional list is
not checked against the mode.

### `--maxfiles` for quick tests

```
python Hmm.py 2024 isGGH -s 103 --maxfiles 12 -o /tmp/test/    # ~2 min
python Hmm.py 2024 isGGH -s 103                                # ~35 min
```

`sumW` is computed from the same subset, so the output is correctly weighted with reduced statistics. 

Pass `-o` somewhere temporary: the filename does not record the truncation, so it will overwrite a full snapshot.

---

## `datasets/*.txt` — the sample lists

One file per category and year, generated from `config/samples.yaml` and
committed to git. Comment a line out to skip that sample.

```
# VBF 2024 -- sample list for isVBF
# generated 2026-09-08 05:26 from config/samples.yaml by make_datasets.py

# --- MC (50) ---
10     # VBF*Hto2Mu_*M-125
...
# --- data (14) ---
-41    # Run2024C/Muon0
```

Regenerate after editing `samples.yaml`, then commit `datasets/`:

```
python make_datasets.py --write-all              # every category and year
python make_datasets.py --write 2024 isVBF       # just one
python make_datasets.py --show 2024 isVBF        # print, write nothing
```

---

## `run_all.py` — run many samples in parallel

A wrapper that splits the sample list into independent `Hmm.py` jobs, then
either runs them locally or writes a Slurm array to submit.

```
python run_all.py <year> <mode> [samplelist] [options]
```

| option | meaning |
|---|---|
| `--dry-run` | print the plan, run nothing |
| `-l, --list` | print the id / dataset / group table |
| `-s, --samples` | same selection syntax as `Hmm.py` |
| `-o, --outdir` | passed through to `Hmm.py` |
| `--mc-only`, `--data-only` | restrict by kind |
| `--stream prompt\|parking` | which data stream |
| `-n, --split K` | pack the list into K jobs (default: one per sample) |
| `--ncores N` | threads per job → `--cpus-per-task` |
| `--skip-existing` | drop samples whose snapshot already exists |
| `--logdir` | per-job logs (default `logs/`) |
| `-j, --jobs` | local mode: how many processes at once |
| `--slurm FILE` | write a Slurm array script instead of running |
| `--submit` | run `sbatch` on the generated script straight away |
| `--partition` | Slurm partition (default `submit`) |
| `--time HH:MM:SS` | walltime per task (max 6 days) |
| `--mem-per-cpu MB` | memory per core |
| `--max-concurrent N` | cap running array tasks → `--array=1-M%N` |
| `--conda-env`, `--conda-root` | env to activate on the worker |

Note `-j` differs between the scripts: `--ncores` in `Hmm.py`, `--jobs` in
`run_all.py`.

### Slurm

```
python run_all.py 2024 isGGH --dry-run

python run_all.py 2024 isGGH --slurm slurm/ggH_2024.sh \
    --ncores 8 --mem-per-cpu 400 --time 04:00:00 --max-concurrent 12

sbatch slurm/ggH_2024.sh

squeue -u $USER                       # pending / running
squeue -u $USER -t running
scancel <jobid>                       # whole array
scancel <jobid>_5                     # one task
seff <jobid>                          # efficiency, once finished
```

If using `--slurm`, it writes two files and does not submit, unless also using `--submit`.

The two files are the `sbatch` script, and a `.jobs` file with the exact `Hmm.py` arguments per task, which doubles as the record of what ran. 
Logs go to `logs/<year>_<mode>_<arrayid>_<task>.out`.
Each task is wrapped in `/usr/bin/time -v` so peak memory is in every log.

Rerun the gaps:

```
python run_all.py 2024 isGGH --slurm slurm/retry.sh --skip-existing --ncores 4
sbatch slurm/retry.sh
```

### Local

```
python run_all.py 2024 isGGH -s signal_hmm -j 4 --ncores 2
```

Or use local mode inside an interactive allocation:
```
srun --partition=submit --cpus-per-task=16 --mem-per-cpu=4000 \
     --time=04:00:00 --pty bash
conda activate pyenv && cd ~/HmumuRun3/analysis
python run_all.py 2024 isGGH -j 8 --ncores 2
```

