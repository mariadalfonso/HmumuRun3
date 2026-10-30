# Parametric fit framework

H→μμ dimuon-mass fits: signal shapes from MC, a background envelope from the
data sidebands, and the Combine inputs built from them.

## Environments

| | fits (`sigFit.py`, `bkgFit.py`) | everything else |
|---|---|---|
| where | conda (`pyenv`) | Combine container (v9.2.1) |
| why | needs ROOT ≥ 6.28 | needs the Combine libraries (ROOT 6.22) |

Enter the container from a clean shell. Otherwise conda's ROOT leaks in and
`RooDoubleCBFast` fails to load.
```bash
conda deactivate
singularity shell --cleanenv -B /work -B /ceph --nv /work/submit/jaeyserm/software/docker/combine-standalone_v9.2.1.sif
```
Run everything from `parametricfit/`.

## Workflow

Example: ggHcat, BDT bins, full Run 3. For one bin use e.g. `-b incl`.

**1. Fits** (conda)
```bash
python sigFit.py -c ggHcat -b bdt0 bdt1 bdt2 -y Run3 --wsdir WS_Run3
for b in bdt0 bdt1 bdt2; do python bkgFit.py -c ggHcat -b $b -y Run3 --wsdir WS_Run3; done
```
Look for `CHECK THIS FIT` / `<-- check` in the output.

**2. Cards and workspace** (container)
```bash
python buildWorkspace.py -c ggHcat -b bdt0 bdt1 bdt2 -y Run3 --wsdir WS_Run3
```
Runs `makeWorkspaces.py`, `bwsHrare.py`, `combineCards.py` and
`text2workspace.py`, prints the yield table, and prints the commands below.
`--dry-run` only prints the commands.

**3. Limit and significance** (container). Always use `-t -1`: the workspace
holds the unblinded data.
```bash
OPTS="--cminDefaultMinimizerStrategy 0 --X-rtd MINIMIZER_freezeDisassociatedParams \
  --X-rtd MINIMIZER_multiMin_hideConstants --X-rtd MINIMIZER_multiMin_maskConstraints \
  --X-rtd MINIMIZER_multiMin_maskChannels=2"
combine -M AsymptoticLimits -m 125 -t -1 --run expected --rMin 0 -d WS_Run3/ggHcat_bdt_Run3.root -n _ggHcat_bdt_Run3 $OPTS
combine -M Significance -m 125 -t -1 --expectSignal 1 -d WS_Run3/ggHcat_bdt_Run3.root -n _ggHcat_bdt_Run3_sig $OPTS
```

**4. Bias study** (container, hours: use tmux)
```bash
python -u biasStudy.py -c ggHcat -b bdt0 bdt1 bdt2 -y Run3 --ntoys 1000 -j 10
python biasStudy.py -c ggHcat -b bdt0 bdt1 bdt2 -y Run3 --only report   # replot only; conda is fine
```

## Outputs

In `WS_<year>/`, with `<tag>` = `<category>_<bin>_<year>`:

| file | from | content |
|---|---|---|
| `Signal_<tag>_params.json` | `sigFit.py` | shape parameters, yield per mode (`yields`) |
| `Bkg_<tag>_params.json` | `bkgFit.py` | envelope, fits, background yield (`normalization.value`) |
| `Signal_/Bkg_<tag>_workspace.root` | `makeWorkspaces.py` | models built with Combine's classes |
| `datacard_<tag>.txt`, `workspace_<tag>_workspace.root` | `bwsHrare.py` | per-bin card |
| `datacard_<name>.txt`, `<name>.root` | `buildWorkspace.py` | combined card and workspace, e.g. `ggHcat_bdt_Run3` |
| `Yields_<name>.json` | `buildWorkspace.py` | signal per mode and total, background, per bin and summed |

- Combine writes `higgsCombine*.root` to `parametricfit/`.
- The bias study writes to `BIAS/<tag>/`.
- Plots go to `~/public_html/HmumuFits/{signal_fits,bkg_fits,bias}/<category>/`.

## What the scripts do

- **`sigFit.py`**: fits a double-sided Crystal Ball to each production mode's
  MC and freezes the shape.
  - `CMS_scale_m` and `CMS_res_m` stay floating as nuisances.
  - The shape is fitted on the 2022–2024 MC. The yield covers every data era:
    2025–2026 use 2024's yield per fb⁻¹ (`signal_yield_scales` in
    `prepareFits.py`).
- **`bkgFit.py`**: fits the sidebands (110–120 and 130–150 GeV) with
  Bernstein, exponential and power-law sums.
  - The F-test picks one order per family; a goodness-of-fit cut builds the
    envelope.
  - The yield is extrapolated over 110–150 GeV.
- **`makeWorkspaces.py`**: rebuilds both models from the JSONs with Combine's
  classes. Nothing is refitted.
- **`bwsHrare.py`**: a copy of `analysis/FITS/bwsHrare.py`.
  - Signal rates come from the workspace norms; the background is a free
    `rateParam`.
  - It adds the lumi and BR uncertainties, `CMS_scale_m` / `CMS_res_m` as
    `param 0 1`, and `pdfindex_<tag> discrete`.
  - Shapes paths are relative to `parametricfit/`, so run combine from there.
- **`biasStudy.py`**: for each envelope member, generates toys from that
  member and fits them with the full envelope.
  - Bias = (μ_fit − μ_inj)/σ_μ. Every member needs |median Bias| < 0.20.
  - Toys come from the pre-fit model, so the blinded region is never fitted.
  - Finished jobs are skipped, so a crashed run resumes.

## Status (2026-10-29)

Expected results for ggHcat, all from Asimov:

| | incl 2024 | BDT bins 2024 | BDT bins Run 3 |
|---|---|---|---|
| median expected limit | r < 1.76 | r < 1.39 | r < 0.85 |
| expected significance | | 1.57σ | 2.33σ |

The significance scales less than √L between 2024 and Run 3; still to check.

- Failed fits are not yet excluded from the envelope.
- Bernstein orders stop at 5.
- `--algo singles` misses the low-side errors; use `--algo grid`.
- Not part of the chain: `sigFitCond.py`, `sigFitDCBG.py`, `fitEtaBins.py`,
  `fitMVABins.py`, `massResVsMass.py`.
- `analysis/FITS/combineCommand.sh` still uses the old `analysis/FITS/` layout.

## Options

**`sigFit.py`**
- `-c`, `-b`, `-s`: categories, bins and modes. All modes by default; keep it
  that way for the chain.
- `-y`: era or group, e.g. `2024`, `Run3`.
- `--wsdir`, `--plotdir`, `--pdf`.
- `--no-freeze`: leave the shape floating.
- `--selftest`: closure test.

**`buildWorkspace.py`**: `-c`, `-b`, `-y`, `--wsdir` (default `WS_<year>`),
`--name`, `--dry-run`.

**`biasStudy.py`**
- `-c`, `-b`, `-y`.
- `--carddir`: default `WS_<year>`. For old cards use
  `../analysis/FITS/workspaces`.
- `--outdir`, `--plotdir`.
- `--inject`: injected r, default 1.
- `--ntoys` (default 1000), `--toys-per-job` (100), `-j`.
- `--rrange`: fixed r range.
- `--threshold`: default 0.20.
- `--combine-opts`, `--seed`, `--force`, `--dry-run`, `--only run|report`.
