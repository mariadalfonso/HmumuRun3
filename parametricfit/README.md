# Parametric fit framework

H→μμ dimuon-mass fits: signal shapes from MC, a background envelope from the
data sidebands, and the Combine inputs built from them.

## Environments

| | fits (`sigFit.py`, `mcCoreFit.py`, `bkgFit.py`) | everything else |
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

## Layout

| where | what |
|---|---|
| top level | the chain: `sigFit.py`, `mcCoreFit.py`, `bkgFit.py`, `buildWorkspace.py`, `biasStudy.py` |
| `fitlib/` | shared code: `config` (fit window, modes, lumi), `style` (CMS plot style), `signal` (the signal fit), `pdfDefinitions`, `dyCore`, `bwz`, `fewz` (with `fewz_spline.root`), `agnostic`, `mcCore`, `gof` (the χ² goodness of fit), `prepareFits`, and the steps `buildWorkspace.py` runs (`makeWorkspaces`, `bwsHrare`) |
| `studies/` | side studies of the signal shape, run as modules: `python -m studies.fitEtaBins 2024` |

## Workflow

Example: ggHcat, BDT bins, full Run 3. For one bin use e.g. `-b incl`.

**1. Fits** (conda)
```bash
python sigFit.py -c ggHcat -b bdt0 bdt1 bdt2 -y Run3
for b in bdt0 bdt1 bdt2; do python bkgFit.py -c ggHcat -b $b -y Run3; done
```
By default `bkgFit.py` fits the background families of the Run 2 analysis
note (Sec. 5.3, eqs. 5.1–5.11; `RUN2_FAMILIES` in `fitlib/pdfDefinitions.py`),
every parameter on data. `--extra-models` adds the families beyond it
(`EXTRA_FAMILIES`): `dyexp`/`dypow`, and `mcexp`/`mcpow`/`mcdy` once
`mcCoreFit.py` has run:
```bash
python mcCoreFit.py -c ggHcat -b bdt0 bdt1 bdt2 -y Run3
for b in bdt0 bdt1 bdt2; do python bkgFit.py -c ggHcat -b $b -y Run3 --extra-models; done
```
Both write `WS_<year>/Bkg_<tag>_params.json`; to keep the two side by side,
give one of them its own `--wsdir` (and the signal there: `sigFit.py
--wsdir`), and `biasStudy.py` its own `--carddir` and `--outdir`, since it
skips toy jobs whose output already exists. Look for `CHECK THIS FIT` /
`<-- check` in the output.

**2. Cards and workspace** (container)
```bash
python buildWorkspace.py -c ggHcat -b bdt0 bdt1 bdt2 -y Run3
```
Runs `fitlib/makeWorkspaces.py`, `fitlib/bwsHrare.py`, `combineCards.py` and
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
The Run 2 note's single-function study (Sec. 6.4.1) instead: every family's
F-test pick as a truth, each fitted with one function at a time.
```bash
python buildWorkspace.py -c ggHcat -b bdt0 bdt1 bdt2 -y Run3 --members families
python -u biasStudy.py --mode single -c ggHcat -b bdt0 bdt1 bdt2 -y Run3 --ntoys 2000 -j 10
```

## Outputs

In `WS_<year>/`, with `<tag>` = `<category>_<bin>_<year>`. Every script
reads and writes there unless `--wsdir` says otherwise.

| file | from | content |
|---|---|---|
| `Signal_<tag>_params.json` | `sigFit.py` | shape parameters, yield per mode (`yields`) |
| `MCcore_<tag>.json` | `mcCoreFit.py` | frozen MC cores, fit quality, MC tags per era |
| `Bkg_<tag>_params.json` | `bkgFit.py` | envelope, fits, background yield (`normalization.value`) |
| `Signal_/Bkg_<tag>_workspace.root` | `fitlib/makeWorkspaces.py` | models built with Combine's classes |
| `datacard_<tag>.txt`, `workspace_<tag>_workspace.root` | `fitlib/bwsHrare.py` | per-bin card |
| `datacard_<name>.txt`, `<name>.root` | `buildWorkspace.py` | combined card and workspace, e.g. `ggHcat_bdt_Run3` |
| `Yields_<name>.json` | `buildWorkspace.py` | signal per mode and total, background, per bin and summed |

- Combine writes `higgsCombine*.root` to `parametricfit/`.
- The bias study writes to `BIAS/<tag>/`.
- Plots go to `~/public_html/HmumuFits/{signal_fits,bkg_fits,bias}/<category>/`,
  and `~/public_html/HmumuFits/mc_core/` for `mcCoreFit.py`.

## What the scripts do

- **`sigFit.py`**: fits a double-sided Crystal Ball to each production mode's
  MC and freezes the shape.
  - `CMS_scale_m` and `CMS_res_m` stay floating as nuisances.
  - The shape is fitted on the 2022–2024 MC. The yield covers every data era:
    2025–2026 use 2024's yield per fb⁻¹ (`signal_yield_scales` in
    `fitlib/prepareFits.py`).
- **`mcCoreFit.py`**: fits the weighted Drell–Yan MC of each bin
  (`BKG_DY_TAGS` in `fitlib/prepareFits.py`: QCD DY, plus EWK Z in VBF,
  from samples generated with 105 < m_μμ < 160 GeV) over 110–150 GeV with three smooth cores, as flashgg does for ggH: 2Exp,
  2Pow, and the LO DY core × (1 + Cheb3). The values are frozen; only the
  shape is used, and only by `bkgFit.py --extra-models`.
  - `--all-bkg-mc` fits the total background MC (`BKG_MC_TAGS`) instead.
  - A core is frozen only if it passes a goodness-of-fit test against the
    MC: χ² with bins merged to ~9 effective entries each, p > 0.01 (the
    same test and cut as `bkgFit.py` on data, `fitlib/gof.py`). A bin with
    too little MC to test gets no MC cores.
- **`bkgFit.py`**: fits the sidebands (110–120 and 130–150 GeV). By
  default with the families of the Run 2 note, Sec. 5.3, every parameter
  floating on data:
  - eqs. 5.1–5.4 (`fitlib/bwz.py`): BWZ and BWZ × Bernstein as one nested
    family (`bwz`, `bwzbern1`–`4`), `bwzgamma`, `bwzredux`;
  - eq. 5.5 (`fitlib/fewz.py`): FEWZ × Bernstein (`fewzbern1`–`4`). The
    FEWZ spline (`fitlib/fewz_spline.root`, copied from flashggFinalFit)
    is fixed; the Bernstein coefficients float;
  - eqs. 5.6–5.8: Bernstein, exponential and power-law sums (`bern`,
    `exp`, `pow`);
  - eqs. 5.9–5.11 (`fitlib/agnostic.py`): exponential + constant,
    power law + constant, Landau × Bernstein (`expconst`, `powconst`,
    `landau`, `landbern1`–`4`).

  `--extra-models` adds the LO Drell–Yan core times exp(Chebyshev) or a
  power law (`dyexp`, `dypow`), and each frozen MC core times
  exp(Chebyshev) (`mcexp`, `mcpow`, `mcdy`).
  - The order of each family, the goodness-of-fit cut and the envelope
    are ours; the Run 2 note excerpt only asks for χ²/dof ≈ 1.
  - The F-test picks one order per family; a goodness-of-fit cut builds the
    envelope. Failed fits (status ≠ 0) are dropped, and a failed order ends
    its family: the order below it is kept, never one above the gap.
  - Each order is fitted twice, from its default start and from the fit of
    the order below (a Bernstein raised one degree reproduces it exactly),
    and the better fit is kept, so a higher order does not end in a worse
    local minimum. `(from the order below)` in the output marks where that
    start won.
  - The member with the lowest NLL + 0.5 × npar, as Combine compares them,
    is the default (`pdfindex` 0): the Asimov and the yield come from it.
  - The yield is extrapolated over 110–150 GeV.
- **`fitlib/makeWorkspaces.py`**: rebuilds both models from the JSONs with Combine's
  classes. Nothing is refitted.
- **`fitlib/bwsHrare.py`**: a copy of `analysis/FITS/bwsHrare.py`.
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
  - `--mode single` is the Run 2 note's version. The truths are every
    family's F-test pick, the ones the goodness-of-fit cut rejected included
    (`buildWorkspace.py --members families`, `*_families` files). Each truth's
    toys are fitted once per function with `pdfindex` frozen: a truth × fit
    matrix of the mean Bias. A function is usable if |Bias| < 0.20 against
    every truth; of those, the one with the highest expected significance
    (Asimov, background frozen to it) is chosen. Output in
    `<outdir>/<tag>_single/`, plot `bias_<tag>_<r>_matrix.png`.

## Status (2026-10-29)

Expected results for ggHcat, all from Asimov:

| | incl 2024 | BDT bins 2024 | BDT bins Run 3 |
|---|---|---|---|
| median expected limit | r < 1.76 | r < 1.39 | r < 0.85 |
| expected significance | | 1.57σ | 2.33σ |

The significance scales less than √L between 2024 and Run 3; still to check.

- Bernstein orders stop at 5.
- `--algo singles` misses the low-side errors; use `--algo grid`.
- `analysis/FITS/combineCommand.sh` still uses the old `analysis/FITS/` layout.

## Options

**`sigFit.py`**
- `-c`, `-b`, `-s`: categories, bins and modes. All modes by default; keep it
  that way for the chain.
- `-y`: era or group, e.g. `2024`, `Run3`.
- `--wsdir` (default `WS_<year>`), `--plotdir`, `--pdf`.
- `--no-freeze`: leave the shape floating.
- `--selftest`: closure test.

**`bkgFit.py`**
- `-c`, `-b`, `-y`, `--wsdir` (default `WS_<year>`), `--plotdir`, `--pdf`.
- `--extra-models`: add `dyexp`/`dypow` and, with `mcCoreFit.py`'s
  output, `mcexp`/`mcpow`/`mcdy` to the Run 2 families.
- `--candidates`: fit exactly these members instead.
- `--mc-core`, `--no-mc-core` (with `--extra-models`), `--bins-per-gev`.

**`buildWorkspace.py`**: `-c`, `-b`, `-y`, `--wsdir` (default `WS_<year>`),
`--name`, `--dry-run`, `--members envelope|families` (`families`: every
family's F-test pick in the multipdf, as separate `*_families` files for
`biasStudy.py --mode single`; the limit cards are not touched).

**`biasStudy.py`**
- `-c`, `-b`, `-y`.
- `--carddir`: default `WS_<year>`. For old cards use
  `../analysis/FITS/workspaces`.
- `--outdir`, `--plotdir`.
- `--inject`: injected r, default 1.
- `--ntoys` (default 1000), `--toys-per-job` (100), `-j`.
- `--rrange`: fixed r range.
- `--threshold`: default 0.20.
- `--mode envelope|single`: the envelope floating (default), or the Run 2
  single-function matrix.
- `--fit-models`: with `--mode single`, fit with these members only (all
  stay truths). The Run 2 note used Landau × Bernstein only as a truth.
- `--metric median|mean`: the quoted Bias. Median by default for the
  envelope (flashggFinalFit), mean for `single` (the Run 2 note).
- `--combine-opts`, `--seed`, `--force`, `--dry-run`, `--only run|report`.
