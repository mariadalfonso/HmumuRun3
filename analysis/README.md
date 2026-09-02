
## create the evironment
First activate your conda environment
```
conda activate myenv
```

If is it the first time creating the environment, run
```
conda env create --name myenv --file=environment.yml
```

## Running `Hmm.py`

Run from the `analysis/` directory:
```
    python Hmm.py <year> <mode> [-o OUTDIR]
```
- **year**: `12022`, `22022`, `12023`, `22023`, `2024`, `2025`, `2026`
- **mode**: `isVBF`, `isGGH`, `isZinv`, `isVlep`, `isVhad`, `isTTlep`, `isTThad`
- **-o/--outdir**: output dir (default `/work/submit/$USER/HmumuRun3/ROOTFILES/`)

Example: `python Hmm.py 12022 isVBF`

Writes snapshots to `<outdir>/<category>/snapshot_mc_<mc>_<year>_<category>.root`.

