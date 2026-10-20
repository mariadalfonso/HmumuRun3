# Fitting framwork

Performing parametric fits in dimuon mass.

## Signal fit

`sigFit.py` fits a double-sided Crystal Ball to signal MC in each `(category, BDT bin, production mode)`, freezes the shape, and writes a RooFit workspace plus a JSON record. `makeWorkspaces.py` rebuilds the workspace from that JSON for `bwsHrare.py` (see below).

### To run

Selftest to validate code setup:
```bash
python sigFit.py --selftest
```
It does a closure test: generates 500k events from a `RooCrystalBall` 
with known parameters, fits them back, and checks the parameters are recovered.

Fit a single signal:
```bash
python sigFit.py -c ggHcat -y Run3 -b incl -s ggH
```


### Options

- **`-c`, `--category`** (default `VBFcat`) — one or more categories
- **`-b`, `--bins`** (default `bdt0`) — BDT bin labels: `bdt0`, `bdt1`, `bdt2`, `incl`
- **`-s`, `--sig`** (default: all for the category) — production modes: `ggH`, `qqH`, `VH`, `ttH`
- **`-y`, `--year`** (default `Run3`) — year tag passed to `getHisto`
- **`--wsdir`** (default `WS_LOCAL`) — workspace + JSON output dir
- **`--plotdir`** (default `~/public_html/HmumuFits/signal_fits`) — plots land in `<plotdir>/<cat>/`
- **`--pdf`** (default off) — also write PDF alongside the PNG
- **`--no-freeze`** (default off) — leave shape parameters floating in the workspace
- **`--selftest`** (default off) — CB closure test, then exit

Output directories are created automatically.

## Full chain: fits in conda, workspaces in the Combine container

The fits need ROOT >= 6.28, and the Combine container has ROOT 6.22: it
normalises each sideband of a blinded `Range("sbLow,sbHigh")` fit separately,
and has no `RooCrystalBall`. So the fits run in conda and write their results
as JSON (+ the data as a TH1D), and only the workspace building runs in the
container, rebuilding the PDFs from those numbers. Nothing is refitted there.

In conda, from `parametricfit/`:
```bash
python sigFit.py -c ggHcat -b incl -y 2024 --wsdir WS_2024
python bkgFit.py -c ggHcat -b incl -y 2024 --wsdir WS_2024
```

Inside the Combine container (`analysis/FITS/runSing`), from `parametricfit/`:
```bash
python makeWorkspaces.py -c ggHcat -b incl -y 2024 --wsdir WS_2024
```
This writes `WS_2024/Signal_ggHcat_incl_2024_workspace.root` (replacing the
conda one) and `WS_2024/Bkg_ggHcat_incl_2024_workspace.root`. Then, still in
the container, `analysis/FITS/bwsHrare.py` makes the datacard from them, and
`text2workspace.py` / `combine` take it from there.

Run `sigFit.py` without `-s`, so the signal JSON covers every mode of the
category under the name `makeWorkspaces.py` expects.

