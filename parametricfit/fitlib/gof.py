"""
gof.py -- the binned chi2 goodness-of-fit shared by bkgFit.py (sidebands
vs data) and mcCoreFit.py (full window vs the background MC).
"""

# Goodness of fit: keep a candidate whose chi2 p-value exceeds this.
GOF_THRESHOLD = 0.01

# Minimum significance content/sqrt(sum err^2) a chi2 bin must reach; bins
# below it are merged with their neighbours first. Without this, the
# sparse 130-150 tail contributes bins of a handful of events whose
# residuals are not chi2-distributed, so the p-value, and therefore the
# GOF cut, is least reliable exactly where the candidates differ most.
# 0 disables merging.
CHI2_MIN_BIN_SIGNIFICANCE = 3.0
# For weighted MC this significance is sqrt(N_eff), so 3.0 asks for
# about 9 effective entries per merged bin.


def merge_chi2_bins(h_data, h_pdf, blind_lo, blind_hi, threshold):
    """Merge low-significance bins, then chi2 over what is left.

      * Bins are accumulated until data/sqrt(sum err^2) reaches
        `threshold`, then emitted as one merged bin, so no chi2 bin is
        built from a handful of events.
      * The blinded window SPLITS the histogram into segments and merging
        never spans it -- otherwise a bin from the low sideband could be
        merged with one from the high sideband, 20 GeV away.
      * A leftover at the end of a segment that does not reach the
        threshold is merged BACKWARDS into the previous bin of that same
        segment, rather than emitted under-populated or dropped.

    Returns (chi2, n_merged_bins).
    """
    bins = []                      # (data, pdf, err2)
    cur = [0.0, 0.0, 0.0]
    segment_start = 0

    def has_current():
        return cur[2] > 0.0 or abs(cur[0]) > 1e-6 or abs(cur[1]) > 1e-6

    def passes():
        if cur[2] <= 0.0:
            return False
        if threshold <= 0.0:
            return True
        if abs(cur[0]) < 1e-6:
            return False
        return cur[0] / cur[2] ** 0.5 >= threshold

    def emit():
        if has_current() and cur[2] > 0.0:
            bins.append(tuple(cur))
        cur[0] = cur[1] = cur[2] = 0.0

    def flush_segment():
        if not has_current():
            return
        if not passes() and len(bins) > segment_start:
            d, p, e = bins[-1]
            bins[-1] = (d + cur[0], p + cur[1], e + cur[2])
            cur[0] = cur[1] = cur[2] = 0.0
            return
        emit()

    ax = h_data.GetXaxis()
    go_blind = (blind_lo < blind_hi
                and (blind_lo > ax.GetXmin() or blind_hi < ax.GetXmax()))

    for i in range(1, h_data.GetNbinsX() + 1):
        c = h_data.GetBinCenter(i)
        if go_blind and blind_lo < c < blind_hi:
            flush_segment()
            segment_start = len(bins)
            continue
        cur[0] += h_data.GetBinContent(i)
        cur[1] += h_pdf.GetBinContent(i)
        cur[2] += h_data.GetBinError(i) ** 2
        if passes():
            emit()

    flush_segment()

    chi2 = sum(((d - p) / e ** 0.5) ** 2 for d, p, e in bins if e > 0)
    return chi2, len(bins)
