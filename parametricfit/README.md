# Fitting framwork

Performing parametric fits in dimuon mass.

## Signal fit

`sigFit.py` fits a double-sided Crystal Ball to signal MC in each `(category, BDT bin, production mode)`, freezes the shape, and writes a RooFit workspace plus a JSON record for `bwsHrare.py` to consume.

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
