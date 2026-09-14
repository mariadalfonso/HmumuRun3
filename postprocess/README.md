# postprocess

Turns the snapshots from `analysis/` into plots, in two steps.

```
NanoAOD  --Hmm.py-->  snapshots  --makeHistos.py-->  histogram file  --SummaryPlots*.py-->  PNGs/PDFs
```

`makeHistos.py` reads events and fills every histogram in one event loop.

The histogram files are saved in `/work/submit/$USER/HmumuRun3/HISTOS/`
(override with `-o`), one file per category+year: `histos_<category>_<year>.root`.

The plotting scripts read the histogram files to produce plots.

## Running

```bash
# once per category+year
python makeHistos.py ggHcat 2024

# then as often as you like
python SummaryPlots.py ggHcat 2024          # stacked data/MC + ratio
python SummaryPlotsShapes.py ggHcat 2024    # normalized shape overlays
```

Categories: `VBFcat`, `ggHcat`, `VLcat`, `TTLcat`, `TTHcat`, `VHcat`, `Zinvcat`

Add `-i` / `-o` to override input and output directories, `--unblind` on
`makeHistos.py` to show blinded data. `-h` lists the rest.

## Files

| File | What it does |
|---|---|
| `makeHistos.py` | fills all histograms in one event loop |
| `SummaryPlots.py` | stacked data/MC plots |
| `SummaryPlotsShapes.py` | shape-comparison plots |
| `utils/prepareHisto.py` | RDataFrame booking |
| `utils/histo_config.py` | regions, process groups, variable lists, plot groups |
| `utils/plot_vars.py` | per-variable expression, binning, axis label |
| `utils/plot_style.py` | CMS style, colours, canvas geometry |
| `utils/LoadTree.py` | builds the snapshot TChain |

Run the three entry points from this directory: `utils/` is imported as a
package, so `python makeHistos.py ...` works but calling it by a longer path
from elsewhere does not.

## Adding a variable

One entry in `plot_vars.py` (expression, binning, label), one line in
`histo_config.get_active_vars()`. Then re-run `makeHistos.py`.

New variables land in the `objects` group; override in
`histo_config.VAR_GROUPS` only for `mass` or `mva`. The variable's branch must
be written to the snapshot for that category — see `config/branches.yaml` in
`analysis/`, whose `base` + `per_mode[mode]` is what the gating in
`get_active_vars()` follows.

## Notes

**Regions.** Every variable is histogrammed in `Inclusive`, `SR_sideband`, and
`Unrestricted` (the last supplies the full signal prediction for sideband
plots).

**Groups.** `SummaryPlots.py` writes to `<cat>_<year>/<group>/<region>/`, with
groups split by how data is treated:

| Group | Data treatment | Region split |
|---|---|---|
| `mass` | blind window cut out of the spectrum | no |
| `mva` | blinded above the per-category score cut | yes |
| `objects` | not blinded — muon, jet and MET kinematics | yes |

Use `--groups` to draw a subset.

**Currently skipped.** `mva` and `dimuon_eta` are booked but their branches
(`discrMVA0`, `HiggsCandCorrEta`) are not in the snapshots, so `makeHistos.py`
warns and skips them in every category. See `TODO.md`.