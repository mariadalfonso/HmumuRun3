"""
optimizeBins.py -- the BDT bin boundaries of a category that maximise the
expected significance, as the Run 2 note chose its ggH categories.

    python optimizeBins.py -c ggHcat -y Run3
    python optimizeBins.py -c ggHcat -y Run3 --bkg data --min-bkg 2000

The Run 2 note (AN-2019/124, Sec. 6.2, "Category optimization"):
  * the figure of merit is the expected significance of a combined fit of
    the categories, signal (all modes) and background modelled over
    110-150 GeV, the background from simulation
  * the boundary is scanned in signal-efficiency quantiles, 1% steps
  * categories are added one at a time until the next gains < 1%
  * odd event numbers of the signal MC for the optimisation, even ones for
    the final signal model
It found 5 ggH categories, at 70, 35, 20 and 5% signal efficiency.

Here, per candidate bin (a range of discrMVA):
  * s_k: the signal MC in 0.25 GeV mass bins over 110-150 GeV, every mode
    (mc 10-15), with sigFit.py's per-era yield scaling
  * b_k: the background over 110-150 GeV, a smooth fit evaluated in the
    same bins. --bkg mc (default, as the note): the background MC of
    mcCoreFit.py (prepareFits.BKG_MC_TAGS, the same era scaling), a
    Poisson fit of exp(cubic in m) to its 1 GeV histogram (the note used
    BWZRedux; any smooth fit serves here). --bkg data: the data SIDEBANDS
    only (110-120, 130-150 GeV), exp(quadratic in m), extrapolated under
    the peak; the 120-130 GeV data are never read.
  * Z^2 = sum_k 2((s_k+b_k) ln(1+s_k/b_k) - s_k), the Asimov significance
    of the binned shape fit with the shapes fixed: what the note's
    combined fit gives without shape uncertainties. The total is
    sqrt(sum over bins of Z^2).

The search is exact: dynamic programming over the grid finds the best split
for each number of bins (the note's greedy search adds one boundary at a
time and can only do as well). --grid quantile (default) puts the grid at
1% steps of signal efficiency, as the note; --grid score at --step in
discrMVA. The recommended number of bins is the last one before a bin that
gains less than --min-gain (default 1%, the note's rule).

Every bin needs --min-bkg background events over 110-150 GeV (the
background MC's weighted sum, or the sideband data count with --bkg data;
default 1000) and, with --bkg mc, --min-mc-raw unweighted MC events
(default 200), so its background fit is stable. The purest bin tends to end
on that limit, so compare a few values.

--signal-parity odd uses odd event numbers only (weights x2), as the note;
the final signal model should then use the even ones, which sigFit.py does
not do yet, so the default is all events.

A guide to the boundaries, not the result: no background-shape uncertainty
(the envelope or core-PDF profiling) and no systematics. Check the best
candidates with the full chain (sigFit, bkgFit, the cards, combine -M
Significance).

Adapted from analysis/MACRO/optimize_categories.py, which counted signal and
data in |m - 125| < 25 GeV (the blinded region included), with s^2/b, at
least 20 background events per bin, and kept only the overall best.
Runs in conda, from parametricfit/. Writes <wsdir>/BinScan_<cat>_<year>.json.
"""

import argparse
import json
import math
import os
import re

import numpy as np
import ROOT

from fitlib import config, prepareFits as pf

ROOT.gROOT.SetBatch(True)

MASS = "HiggsCandCorrMass"
SCORE = "discrMVA"
SB_LO = (config.XLOW, 120.0)          # the data sidebands, --bkg data
SB_HI = (130.0, config.XHIGH)
FIT_BIN = 1.0                          # GeV, the background fit's histogram
Z_BIN = 0.25                           # GeV, the significance's mass bins
N_FINE = 1000                          # score histogram bins over [0, 1]
POLY = {"mc": 3, "data": 2}            # degree of ln(b) in m


# ---------------------------------------------------------------------------
# input: one pass per file group, (score x mass) histograms
# ---------------------------------------------------------------------------

def hist2d_array(h):
    """TH2 contents as an array [score bin, mass bin], the score's under-
    and overflow folded into the first and last bin."""
    nx, ny = h.GetNbinsX(), h.GetNbinsY()
    a = np.zeros((nx, ny))
    for i in range(nx + 2):
        ii = min(max(i - 1, 0), nx - 1)
        for j in range(1, ny + 1):
            a[ii, j - 1] += h.GetBinContent(i, j)
    return a


def fill(files, filters, nmass, name, weight=True):
    df = ROOT.RDataFrame("events", files)
    for f in filters:
        df = df.Filter(f)
    df = df.Define("wgt", "w_allSF" if weight else "1.0")
    h = df.Histo2D((name, "", N_FINE, 0.0, 1.0, nmass, config.XLOW,
                    config.XHIGH), SCORE, MASS, "wgt")
    return hist2d_array(h.GetValue())


def parity_filter(parity):
    return {"all": [], "odd": ["event % 2 == 1"],
            "even": ["event % 2 == 0"]}[parity]


def load_signal(args):
    """[fine score bin, Z_BIN mass bin], summed over modes and MC eras."""
    nz = int(round((config.XHIGH - config.XLOW) / Z_BIN))
    tags = sorted({t for s in args.modes for t in pf.SIGNAL_TAGS[s]})
    scale_parity = 1.0 if args.signal_parity == "all" else 2.0
    total = np.zeros((N_FINE, nz))
    for era, scale in pf.signal_yield_scales(args.year).items():
        files = pf._collect_files([(era, tags)], args.category,
                                  args.rootfiles_dir)
        total += scale * scale_parity * fill(
            files, [f"mc >= 10 && mc <= 15 && !isnan({MASS})"]
            + parity_filter(args.signal_parity), nz, f"s_{era}")
        print(f"  signal {era}: {len(files)} file(s), x {scale:.3f}"
              + ("" if scale_parity == 1 else
                 f", {args.signal_parity} events x 2"))
    return total


def load_background(args):
    """[fine score bin, FIT_BIN mass bin]: background MC over 110-150 GeV
    (weighted, and unweighted for the raw count), or the data sidebands."""
    nm = int(round((config.XHIGH - config.XLOW) / FIT_BIN))
    if args.bkg == "data":
        files = pf.get_data_files(args.year, args.category,
                                  args.rootfiles_dir)
        sb = (f"({MASS} >= {SB_LO[0]} && {MASS} < {SB_LO[1]}) || "
              f"({MASS} >= {SB_HI[0]} && {MASS} < {SB_HI[1]})")
        a = fill(files, [f"mc < 0 && mc > -100 && !isnan({MASS})", sb], nm,
                 "b_data", weight=False)
        print(f"  data: {len(files)} file(s), {a.sum():.0f} sideband events "
              "(120-130 GeV not read)")
        return a, a
    w = np.zeros((N_FINE, nm))
    raw = np.zeros((N_FINE, nm))
    for era, scale in pf.signal_yield_scales(args.year).items():
        files = pf.get_bkg_mc_files(args.category, era, args.rootfiles_dir)
        cut = [f"mc > 15 && !isnan({MASS})"]
        w += scale * fill(files, cut, nm, f"b_{era}")
        raw += fill(files, cut, nm, f"braw_{era}", weight=False)
        print(f"  background MC {era}: {len(files)} file(s), x {scale:.3f}")
    print(f"  background MC: {w.sum():.0f} expected, {raw.sum():.0f} raw "
          "events over 110-150 GeV")
    return w, raw


# ---------------------------------------------------------------------------
# the significance of one candidate bin
# ---------------------------------------------------------------------------

def background_density(counts, mode):
    """Poisson fit of exp(poly(u)), u = m - 130, to the FIT_BIN histogram:
    over 110-150 (mc) or the sidebands only (data). Returns the fitted
    density per GeV as a function, or None."""
    centres = config.XLOW + FIT_BIN * (np.arange(len(counts)) + 0.5)
    use = np.ones(len(counts), bool) if mode == "mc" else \
        (centres < SB_LO[1]) | (centres >= SB_HI[0])
    u, n = centres[use] - 130.0, counts[use]
    deg = POLY[mode]
    if (n > 0).sum() < deg + 2:
        return None
    X = np.vstack([u ** p for p in range(deg + 1)]).T
    pos = n > 0
    w0 = np.sqrt(n[pos])
    beta, *_ = np.linalg.lstsq(X[pos] * w0[:, None], np.log(n[pos]) * w0,
                               rcond=None)
    # Poisson regression (IRLS); quasi-Poisson for the weighted MC
    for _ in range(30):
        mu = np.exp(np.clip(X @ beta, -50, 50))
        z = X @ beta + (n - mu) / mu
        sw = np.sqrt(mu)
        new, *_ = np.linalg.lstsq(X * sw[:, None], z * sw, rcond=None)
        done = np.max(np.abs(new - beta)) < 1e-10
        beta = new
        if done:
            break

    def dens(m):
        uu = np.asarray(m) - 130.0
        return np.exp(sum(beta[p] * uu ** p for p in range(deg + 1))) / FIT_BIN
    return dens


def bin_stats(sig, bkg, raw, args):
    """(s, b over 110-150, Z^2, s in 120-130, b in 120-130) of one bin, or
    None if it fails the minimum statistics."""
    nb = bkg.sum()
    if nb < args.min_bkg or (args.bkg == "mc" and raw.sum() < args.min_mc_raw):
        return None
    dens = background_density(bkg, args.bkg)
    if dens is None:
        return None
    edges = config.XLOW + Z_BIN * np.arange(len(sig) + 1)
    centres = 0.5 * (edges[:-1] + edges[1:])
    b = dens(centres) * Z_BIN
    s = sig
    good = (b > 0) & (s > 0)
    z2 = float(np.sum(2.0 * ((s[good] + b[good]) * np.log1p(s[good] / b[good])
                             - s[good])))
    win = (centres >= 120.0) & (centres < 130.0)
    return {"s": float(s.sum()), "b_full": float(b.sum()), "Z2": z2,
            "s_window": float(s[win].sum()), "b_window": float(b[win].sum()),
            "n_bkg": float(nb), "n_raw": float(raw.sum())}


# ---------------------------------------------------------------------------
# the grid and the search
# ---------------------------------------------------------------------------

def make_grid(sig, args):
    """Cell edges as fine score-bin indices, from the top down: 1% steps in
    signal efficiency (quantile) or --step in score."""
    if args.grid == "score":
        step = int(round(args.step * N_FINE))
        return list(range(0, N_FINE + 1, step))
    per_bin = sig.sum(axis=1)
    above = np.cumsum(per_bin[::-1])[::-1]          # signal at or above bin i
    total = above[0]
    edges = {0, N_FINE}
    for q in np.arange(args.quantile_step, 1.0 - 1e-9, args.quantile_step):
        # the lowest fine edge with at most a fraction q of the signal above
        i = int(np.searchsorted(-above, -q * total, side="left"))
        if 0 < i < N_FINE:
            edges.add(i)
    return sorted(edges)


def best_partitions(grid, stats, max_bins):
    """Exact best split of the grid into k = 1..max_bins ranges maximising
    the sum of Z^2. stats(i, j) gives a range's statistics or None.
    Returns {k: (Z, [grid edge indices])}."""
    n = len(grid) - 1
    table = {}
    for i in range(n):
        for j in range(i + 1, n + 1):
            r = stats(grid[i], grid[j])
            if r is not None:
                table[i, j] = r
    best = [[-1.0] * (n + 1) for _ in range(max_bins + 1)]
    prev = [[None] * (n + 1) for _ in range(max_bins + 1)]
    best[0][0] = 0.0
    for k in range(1, max_bins + 1):
        for j in range(1, n + 1):
            for i in range(j):
                if best[k - 1][i] < 0 or (i, j) not in table:
                    continue
                v = best[k - 1][i] + table[i, j]["Z2"]
                if v > best[k][j]:
                    best[k][j], prev[k][j] = v, i
    out = {}
    for k in range(1, max_bins + 1):
        if best[k][n] < 0:
            continue
        idx, j = [n], n
        for kk in range(k, 0, -1):
            j = prev[kk][j]
            idx.append(j)
        out[k] = (math.sqrt(best[k][n]), [grid[i] for i in idx[::-1]])
    return out


def current_edges(category):
    """The cut values of fitlib/prepareFits.SELMVA for the category."""
    vals = set()
    for b, cut in pf.SELMVA.get(category, {}).items():
        if re.fullmatch(r"bdt\d+", b):
            vals.update(float(v) for v in re.findall(r"[0-9]*\.[0-9]+", cut))
    return sorted(vals)


def describe(fine_edges, stats, sig_total_above):
    rows = []
    for lo, hi in zip(fine_edges[:-1], fine_edges[1:]):
        r = dict(stats(lo, hi))
        r.update(lo=lo / N_FINE, hi=hi / N_FINE,
                 eff_lo=float(sig_total_above(hi)),
                 eff_hi=float(sig_total_above(lo)), Z=math.sqrt(r["Z2"]))
        rows.append(r)
    return rows


def print_rows(rows):
    print(f"    {'score range':<16}{'signal eff.':>14}{'s':>9}{'b (110-150)':>13}"
          f"{'s/sqrt(b) 120-130':>19}{'Z':>8}")
    for r in rows:
        swb = r["s_window"] / math.sqrt(r["b_window"]) if r["b_window"] else 0
        eff = f"{100 * r['eff_hi']:.0f}-{100 * r['eff_lo']:.0f}%"
        print(f"    [{r['lo']:.3f}, {r['hi']:.3f}){'':<2}{eff:>14}{r['s']:>9.1f}"
              f"{r['b_full']:>13.0f}{swb:>19.3f}{r['Z']:>8.3f}")


def selmva(edges):
    """Score edges (fractions of 1) as a prepareFits.SELMVA entry."""
    inner = edges[1:-1]
    cuts = {}
    for k in range(len(inner) + 1):
        parts = ([f"discrMVA>={inner[k - 1]:g}"] if k > 0 else []) + \
                ([f"discrMVA<{inner[k]:g}"] if k < len(inner) else [])
        cuts[f"bdt{k}"] = " && ".join(parts)
    return cuts


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--category", default="ggHcat",
                    choices=sorted(pf.SELMVA))
    ap.add_argument("-y", "--year", default="Run3")
    ap.add_argument("--bkg", choices=("mc", "data"), default="mc",
                    help="background from the background MC (default, as "
                         "the note) or the data sidebands")
    ap.add_argument("--grid", choices=("quantile", "score"),
                    default="quantile",
                    help="boundaries at signal-efficiency steps (default) "
                         "or at --step in discrMVA")
    ap.add_argument("--quantile-step", type=float, default=0.01,
                    help="signal-efficiency step of --grid quantile "
                         "(default 0.01, as the note)")
    ap.add_argument("--step", type=float, default=0.01,
                    help="discrMVA step of --grid score (default 0.01)")
    ap.add_argument("--max-bins", type=int, default=6)
    ap.add_argument("--min-gain", type=float, default=0.01,
                    help="recommend the last number of bins before one that "
                         "gains less than this (default 0.01, the note)")
    ap.add_argument("--min-bkg", type=float, default=1000,
                    help="minimum background events per bin over 110-150 "
                         "GeV (MC weighted sum, or sideband data count; "
                         "default 1000)")
    ap.add_argument("--min-mc-raw", type=float, default=200,
                    help="--bkg mc: minimum unweighted MC events per bin "
                         "(default 200)")
    ap.add_argument("--signal-parity", choices=("all", "odd", "even"),
                    default="all",
                    help="signal MC events used (odd: the note's choice, "
                         "weights x2; default all)")
    ap.add_argument("--modes", nargs="+", default=sorted(pf.SIGNAL_TAGS),
                    choices=sorted(pf.SIGNAL_TAGS),
                    help="signal modes (default all)")
    ap.add_argument("--rootfiles-dir", default=pf.ROOTFILES_DIR)
    ap.add_argument("--wsdir", default=None,
                    help="where the JSON goes (default WS_<year>)")
    args = ap.parse_args()
    args.wsdir = args.wsdir or config.default_wsdir(args.year)
    if args.grid == "score":
        step = args.step * N_FINE
        if abs(step - round(step)) > 1e-9 or N_FINE % int(round(step)):
            ap.error(f"--step must be a multiple of {1 / N_FINE:g} dividing 1")

    print(f"category : {args.category}, {args.year}\n"
          f"background: {'background MC' if args.bkg == 'mc' else 'data sidebands'}"
          f", fitted over 110-150 GeV\n"
          f"grid     : "
          + (f"{100 * args.quantile_step:g}% signal efficiency"
             if args.grid == "quantile" else f"{args.step:g} in discrMVA")
          + f", at least {args.min_bkg:g} background events per bin, up to "
          f"{args.max_bins} bins")
    sig = load_signal(args)
    bkg, raw = load_background(args)
    print(f"  signal over 110-150 GeV: {sig.sum():.1f}")

    sig_c = np.vstack([np.zeros(sig.shape[1]), np.cumsum(sig, axis=0)])
    bkg_c = np.vstack([np.zeros(bkg.shape[1]), np.cumsum(bkg, axis=0)])
    raw_c = np.vstack([np.zeros(raw.shape[1]), np.cumsum(raw, axis=0)])
    cache = {}

    def stats(lo, hi):
        if (lo, hi) not in cache:
            cache[lo, hi] = bin_stats(sig_c[hi] - sig_c[lo],
                                      bkg_c[hi] - bkg_c[lo],
                                      raw_c[hi] - raw_c[lo], args)
        return cache[lo, hi]

    per_bin = sig.sum(axis=1)
    above = np.concatenate([np.cumsum(per_bin[::-1])[::-1], [0.0]])

    def eff_above(i):
        return above[i] / above[0]

    grid = make_grid(sig, args)
    results = best_partitions(grid, stats, args.max_bins)

    out = {"category": args.category, "year": args.year, "bkg": args.bkg,
           "grid": args.grid, "min_bkg": args.min_bkg,
           "signal_parity": args.signal_parity, "modes": args.modes,
           "best": {}}
    cur = current_edges(args.category)
    if cur:
        fe = [0] + [int(round(v * N_FINE)) for v in cur] + [N_FINE]
        if all(stats(a, b) is not None for a, b in zip(fe[:-1], fe[1:])):
            rows = describe(fe, stats, eff_above)
            z = math.sqrt(sum(r["Z2"] for r in rows))
            print(f"\ncurrent  : edges {', '.join(f'{v:g}' for v in cur)} "
                  f"(prepareFits.SELMVA)   Z = {z:.3f}")
            print_rows(rows)
            out["current"] = {"edges": [0.0] + cur + [1.0], "Z": z,
                              "bins": rows}
        else:
            print(f"\ncurrent  : edges {cur}: a bin fails the minimum "
                  "statistics")

    prevz, recommended = None, None
    for k, (z, fe) in sorted(results.items()):
        edges = [i / N_FINE for i in fe]
        rows = describe(fe, stats, eff_above)
        gain = None if prevz is None else z / prevz - 1
        if gain is not None and gain < args.min_gain and recommended is None:
            recommended = k - 1
        print(f"\n{k} bin(s): Z = {z:.3f}"
              + ("" if gain is None else f"   (+{100 * gain:.1f}%)")
              + "   edges " + ", ".join(f"{e:g}" for e in edges[1:-1])
              + "   signal eff. "
              + ", ".join(f"{100 * eff_above(i):.0f}%" for i in fe[1:-1]))
        print_rows(rows)
        out["best"][k] = {"edges": edges, "Z": z, "bins": rows,
                          "signal_eff_edges": [eff_above(i) for i in fe],
                          "selmva": selmva(edges)}
        prevz = z
    if results:
        if recommended is None:
            recommended = max(results)
            print(f"\nevery added bin gains >= {100 * args.min_gain:g}% up to "
                  f"{recommended}; raise --max-bins to look further")
        out["recommended"] = recommended
        print(f"\nrecommended: {recommended} bin(s) (the next gains < "
              f"{100 * args.min_gain:g}%), as a prepareFits.SELMVA entry:\n  "
              + json.dumps(out["best"][recommended]["selmva"]))
    os.makedirs(args.wsdir, exist_ok=True)
    path = os.path.join(args.wsdir,
                        f"BinScan_{args.category}_{args.year}.json")
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2)
    print(f"\nwrote {path}")


if __name__ == "__main__":
    main()
