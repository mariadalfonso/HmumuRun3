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
| `prepareHisto.py` | RDataFrame booking |
| `histo_config.py` | regions, process groups, variable lists |
| `plot_vars.py` | per-variable expression, binning, axis label |
| `plot_style.py` | CMS style, colours, canvas geometry |
| `LoadTree.py` | builds the snapshot TChain |

## Adding a variable

One entry in `plot_vars.py` (expression, binning, label), one line in
`histo_config.get_active_vars()`. Then re-run `makeHistos.py`.

## Notes

Regions: every variable is histogrammed in `Inclusive`, `SR_sideband`, and
`Unrestricted` (the last supplies the full signal prediction for sideband
plots).