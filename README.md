# HmumuRun3

Code base for H -> μμ analysis CMS Run 3.

```
NanoAOD  --analysis/-->  snapshots  --postprocess/-->  plots
```

## Setup

```bash
conda env create --name myenv --file=analysis/environment.yml   # first time
conda activate myenv
```

## Running

```bash
cd analysis
python Hmm.py 2024 isGGH          # event selection -> snapshot ntuples

cd ../postprocess
python makeHistos.py ggHcat 2024  # snapshots -> histograms
python SummaryPlots.py ggHcat 2024
```

## Layout

- **`analysis/`** — event selection. Applies the analysis channel's cuts and
  writes one slim tree per sample.
- **`postprocess/`** — histogramming and plotting.

Each directory has its own README with the details.
