
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

## Sample definitions: `config/samples.yaml`

Datasets, cross-sections and branching ratios live in `config/samples.yaml`.
`datasets.py` reads it; nothing is globbed until a sample is actually
selected, so only the samples a run needs touch the filesystem.

Blocks in the file:
- **`periods`** — campaign glob and base directories per era (v12 vs v15)
- **`branching_ratios`** — BRs and decay fractions, in one place
- **`xsecs`** — named cross-sections **in fb**, split into `run3` (13.6 TeV)
  and `run2` (13 TeV reference values)
- **`groups`** — sample id lists used by `getMCList`
- **`samples`** / **`data`** — id -> path pattern + cross-section
- **`data_selection`** — which data ids to run per period

### Adding a sample

```yaml
  151:
    path: "{ceph}/{year}/MyNewSample_TuneCP5_13p6TeV_powheg-pythia8/NANOAODSIM/{campaign}"
    xsec: {ref: TTH, br: H_to_mumu}
```

Path placeholders `{ceph}` `{scratch}` `{year}` `{campaign}` are filled from
`periods`. An `xsec` is written one of these ways:

```yaml
    xsec: 2219000                             # a plain value in fb
    xsec: {ref: W}                            # named process from xsecs.run3
    xsec: {ref: VBFH, br: H_to_mumu}          # times a branching ratio
    xsec: {ref: Wm, br: [H_to_WW, W_to_qq]}   # times several
    xsec: {value: 21650, factor: 0.6}         # times an empirical scaling
```

`br` takes physics constants from `branching_ratios`; `factor` is for everything else, so empirical rescalings stay visible as such.

To include the sample in a mode's selection, add its id to a list under `groups`. The per-period and per-mode branching lives in `datasets.getMCList`.

Helpers: `getXsec(name, run)`, `getBR(name)`, `resolve_xsec(spec)`.