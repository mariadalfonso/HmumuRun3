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
| `bkgMassInBDT.py` | background m_mumu in BDT slices, reads the snapshots directly (see below) |
| `utils/prepareHisto.py` | RDataFrame booking |
| `utils/histo_config.py` | regions, process groups, variable lists, plot groups |
| `utils/plot_vars.py` | per-variable expression, binning, axis label |
| `utils/plot_style.py` | CMS style, colours, canvas geometry |
| `utils/LoadTree.py` | builds the snapshot TChain |

Run the entry points from this directory: `utils/` is imported as a
package, so `python makeHistos.py ...` works but calling it by a longer path
from elsewhere does not.

## Background mass in BDT slices: `bkgMassInBDT.py`

A standalone study, not part of the two-step chain: it reads the snapshots
directly and writes PNGs, with no histogram file in between. It checks whether
the BDT sculpts the background m_mumu, and why.

```bash
python bkgMassInBDT.py 2024 --cat ggHcat                    # default samples
python bkgMassInBDT.py 2024 --cat ggHcat -s 103,104,101     # chosen samples
python bkgMassInBDT.py 2024 --cat VLcat --mass-range 70 200
```

Input `/work/submit/$USER/HmumuRun3/bkg/<cat>/` (`-d`), output
`~/public_html/HmumuRun3/bkgshape/<cat>/` (`-o`). The snapshots have to be
written there, e.g. `run_all.py ... -o /work/submit/$USER/HmumuRun3/bkg/`.
MC only, weighted with `w_allSF` (times `boson_ptWeight` for DY).

Two controls, both printed as tables and drawn:

| Control | Categories | Plots |
|---|---|---|
| process composition: is the mix changing, or one process's shape? | all | `bkgMass_BDTslices_all`, `bkgMass_BDTbins`, `bkgMass_BDTslices_<group>`, `bkgMass_composition`, `bkgMass_noBDT_byProcess` |
| good-muon multiplicity: does the muon pairing sculpt the mass? | `VLcat`, `TTLcat` (`--mult` elsewhere) | `bkgMass_pairingControl`, `bkgMass_noBDT_pairing`, `bkgMass_noBDT_bySample_*`, `bkgMass_noBDT_dominant`, `bkgMass_BDTslices_n2`/`_n3ormore` |

The multiplicity control is off where the category vetoes extra muons, since
every event there has exactly two. The per-group plots and the composition
are skipped when only one process group is present. The script's docstring
(`-h`) has the reasoning behind both controls.

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
plots), and in `SR` (the snapshot's `isH` flag, m in [115,135]), which is
MC-only: data is never booked there.

**Groups.** `SummaryPlots.py` writes to `<cat>_<year>/<group>/<region>/`, with
groups split by how data is treated:

| Group | Data treatment | Region split |
|---|---|---|
| `mass` | blind window cut out of the spectrum | no |
| `mva` | blinded above the per-category score cut | yes |
| `objects` | not blinded — muon, jet and MET kinematics | yes |

Use `--groups` to draw a subset.
