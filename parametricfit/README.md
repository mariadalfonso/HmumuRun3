# Parametric fit framework

Parametric fits of the dimuon mass for H→μμ: a signal shape from MC, a
background envelope from the data sidebands, and the Combine workspaces and
datacard built from them.

## What runs where

The fits need ROOT ≥ 6.28. The Combine container on SubMIT
(`analysis/FITS/runSing`, combine v9.2.1) ships ROOT 6.22, which normalises
each sideband of a blinded `Range("sbLow,sbHigh")` fit separately and has no
`RooCrystalBall`. So the chain is split in two:

| step | script | where | reads | writes (in `--wsdir`) |
|---|---|---|---|---|
| 1 | `sigFit.py` | conda | signal snapshots | `Signal_<tag>_params.json` (+ a conda-only workspace) |
| 2 | `bkgFit.py` | conda | data snapshots | `Bkg_<tag>_params.json`, `Bkg_<tag>_data.root` |
| 3 | `makeWorkspaces.py` | container | the two JSONs + TH1D | `Signal_<tag>_workspace.root`, `Bkg_<tag>_workspace.root` |
| 4 | `analysis/FITS/bwsHrare.py` | container | the two workspaces | datacard + combined workspace |
| 5 | `text2workspace.py`, `combine` | container | the datacard | limits, fits |
| 6 | `biasStudy.py` | container | the datacard | `BIAS/<tag>/`: toys, fits, bias summary |

`<tag>` is `<category>_<bin>_<year>`, e.g. `ggHcat_incl_2024`. Only numbers
(JSON) and a plain TH1D cross from conda to the container, and nothing is
refitted there.

- **`sigFit.py`** fits a double-sided Crystal Ball to the signal MC of each
  production mode, freezes the shape, and leaves two shape nuisances for
  Combine: `CMS_scale_m` (muon momentum scale) and `CMS_res_m` (resolution),
  plus `MH` for mass scans. The shape is fitted on the MC eras (2022–2024),
  but the yield, which becomes the datacard rate, covers every era of the
  requested year. 2025 and 2026 have no signal MC, so they take 2024's yield
  per fb⁻¹: for `-y Run3` the 2024 yield is scaled by
  (L₂₀₂₄ + L₂₀₂₅ + L₂₀₂₆)/L₂₀₂₄ (`SIGNAL_PROXY` and `signal_yield_scales` in
  `prepareFits.py`, luminosities from `LUMIS`). The JSON keeps the unscaled MC
  yield and the per-era factors under `normalization.extrapolation`.
- **`bkgFit.py`** fits the background candidates (Bernstein, exponential
  and power-law sums) on the sidebands 110–120 and 130–150 GeV, with
  120–130 blinded. It picks one order per family with an F-test, keeps those
  passing a goodness-of-fit cut, and extrapolates the yield into the blinded
  region. The survivors form the envelope Combine profiles over (discrete
  profiling).
- **`makeWorkspaces.py`** rebuilds both models with the container's classes
  (`RooDoubleCBFast` for the signal) and writes `crystal_ball_<tag>_<mode>`,
  `multipdf_<tag>_bkg`, `pdfindex_<tag>`, `multipdf_<tag>_bkg_norm` and
  `datahist_<tag>`, the names `bwsHrare.py` expects.
- **`bwsHrare.py`** writes the datacard: signal rates from the MC, the
  background yield as a free `rateParam`, `CMS_scale_m` / `CMS_res_m` as
  `param 0 1`, and `pdfindex_<tag> discrete`.

## How to run: ggHcat, inclusive, 2024

**1–2. Fits, in conda**, from `parametricfit/`:
```bash
conda activate pyenv
python sigFit.py -c ggHcat -b incl -y 2024 --wsdir WS_2024
python bkgFit.py -c ggHcat -b incl -y 2024 --wsdir WS_2024
```
Run `sigFit.py` without `-s`, so the signal JSON covers every production mode
of the category. Check the printouts for `CHECK THIS FIT` / `<-- check`, and the
plots in `~/public_html/HmumuFits/{signal_fits,bkg_fits}/ggHcat/`.

**3. Workspaces, in the container.** Leave conda first, and start the container
with `--cleanenv`. Otherwise the conda ROOT headers leak into it and the
Combine classes fail to load ("cannot instantiate incomplete class
`RooDoubleCBFast`").
```bash
conda deactivate
singularity shell --cleanenv -B /work -B /ceph --nv /work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif
cd <path to>/parametricfit
python makeWorkspaces.py -c ggHcat -b incl -y 2024 --wsdir WS_2024
```

**4. Datacard**, still in the container, from `analysis/FITS/`:
```bash
python bwsHrare.py --whichCat=ggHcat --whichBIN=incl --whichYear=2024 \
    --inputFileSig=../../parametricfit/WS_2024/Signal_ggHcat_incl_2024_workspace.root \
    --inputFileBKG=../../parametricfit/WS_2024/Bkg_ggHcat_incl_2024_workspace.root \
    --output=workspaces/workspace_ggHcat_incl_2024_workspace.root \
    --datCardName=workspaces/datacard_ggHcat_incl_2024.txt
```
The card should end with `CMS_res_m param 0 1`, `CMS_scale_m param 0 1` and
`pdfindex_ggHcat_incl_2024 discrete`.

**5. Combine**, from `analysis/FITS/`:
```bash
text2workspace.py -m 125 workspaces/datacard_ggHcat_incl_2024.txt -o workspaces/ggHcat_incl_2024.root

# expected limit
combine -M AsymptoticLimits -m 125 -t -1 --run expected --rMin 0 \
    -d workspaces/ggHcat_incl_2024.root --X-rtd MINIMIZER_freezeDisassociatedParams -n _ggH2024

# Asimov fit with r = 1: r ~ 1, CMS_scale_m and CMS_res_m ~ 0 +/- 1
combine -M MultiDimFit -m 125 -t -1 --expectSignal=1 --algo singles \
    -P r -P CMS_scale_m -P CMS_res_m --floatOtherPOIs 1 \
    -d workspaces/ggHcat_incl_2024.root --X-rtd MINIMIZER_freezeDisassociatedParams -n _ggH2024_check
```
**Keep `-t -1` on every combine call.** The workspace holds the full data,
the blinded region included.

**6. Bias study**, in the container, from `parametricfit/`:
```bash
python biasStudy.py -c ggHcat -b incl -y 2024 --ntoys 1000 -j 8
python biasStudy.py -c ggHcat -b incl -y 2024 --only report   # re-read only
```
For each envelope member it generates toys from that member alone
(`pdfindex_<tag>` frozen), with r = 1, and fits them with the full model, the
index floating, as in the limit. The bias is
Bias(μ_fit, μ_inj, σ_μ) = (μ_fit − μ_inj)/σ_μ, with μ Combine's r and σ_μ the
error on the side facing μ_inj. A member passes if |median Bias| < 0.20, and
every member must pass. The table also gives how often the fit picked each
member, and why toys were dropped (fit failed, no interval, interval at the
r range edge). The r range of the toy fits is r_inj ± 10 σ_r from an Asimov
fit, so the 1σ interval does not reach the edge (unlike `--rMin 0` above).

Toys come from the pre-fit workspace, never from a fit to data:
`--toysFrequentist` would fit the blinded region. Output goes to
`BIAS/<tag>/` (`bias_<tag>_r1p00.json`), plots to
`~/public_html/HmumuFits/bias/<category>/`. Jobs whose output exists are
skipped, so a crashed run resumes; add `--seed` for more toys,
`--force` to redo. `--dry-run` prints the combine commands.

For another category, bin or year, change `-c` / `-b` / `-y` and the
`--which*` options; the same values must be used in every step.
`analysis/FITS/combineCommand.sh` loops steps 4–5 over categories and bins
once its `wsDIR`, `years`, `cases` and `cats` are edited.

## Status (2026-10-20)

- The full chain runs for ggHcat incl 2024: expected limit r < 1.76 (median).
- **Open: background envelope.** Failed fits are not yet excluded from the
  envelope, and `bkgFit.py` prints the full-window event count. Bernstein
  orders stop at 5 by choice.
- In `--algo singles`, the low-side errors of `r` and the nuisances are not
  found (reported to the range edge); use `--algo grid` scans for those.
- Not wired into this chain: the studies `sigFitCond.py`, `sigFitDCBG.py`,
  `fitEtaBins.py`, `fitMVABins.py`, `massResVsMass.py`.

## Signal fit options

- **`-c`, `--category`** (default `VBFcat`): one or more categories
- **`-b`, `--bins`** (default `bdt0`): BDT bin labels, e.g. `incl bdt0 bdt1 bdt2`
- **`-s`, `--sig`** (default: all modes of the category): `ggH`, `qqH`, `VH`, `ttH`
- **`-y`, `--year`** (default `Run3`): era or group (`2022`, `2023`, `2024`, `Run3`, ...)
- **`--wsdir`** (default `WS_LOCAL`): JSON + workspace output dir
- **`--plotdir`** (default `~/public_html/HmumuFits/signal_fits`): plots land in `<plotdir>/<cat>/`
- **`--pdf`**: also write PDFs
- **`--no-freeze`**: leave the shape parameters floating in the workspace
- **`--selftest`**: closure test (500k generated events, fitted back), then exit

## Bias study options

- **`-c`, `-b`, `-y`**: as for the fits; one study per `<cat>_<bin>_<year>`
- **`--carddir`** (default `analysis/FITS/workspaces`): where `datacard_<tag>.txt` is
- **`--outdir`** (default `BIAS`), **`--plotdir`** (default `~/public_html/HmumuFits/bias`)
- **`--inject`** (default `1`): injected r, one or more, e.g. `--inject 0 1`
- **`--ntoys`** (default 1000) per member and r, in jobs of **`--toys-per-job`** (100),
  **`-j`** at a time
- **`--rrange RMIN RMAX`**: fixed r range instead of the Asimov one
- **`--threshold`** (default 0.20): max |median Bias|
- **`--combine-opts`**: passed to every combine call, e.g. `"--toysNoSystematics"`
- **`--only run|report`**, **`--force`**, **`--dry-run`**, **`--pdf`**
