"""
bkgFit.py -- blinded background fit, with an F-test-selected envelope.

Blinding
--------
    Higgs SR   115 - 135 GeV        BLINDED: never fitted, never plotted
    sideband   110-115 + 135-150    what the fit sees

Half the window is masked, so the shape is constrained by 20 GeV of data
and extrapolated across 20 GeV. Two consequences worth keeping in mind:
high-order candidates are much less well determined than they would be
with a narrow blind window, and the normalisation has to be extrapolated
rather than counted.

The full-range data DOES go into the workspace, as `datahist_<tag>`,
because the datacard's `observation -1` needs it and an observed result
would otherwise be impossible. The blinding is procedural: the fit uses
only the sideband and the plots omit the SR points, but anyone opening
the file can look.

Model selection
---------------
Not a hardcoded pair. Two stages, both on the MASKED likelihood:

  1. F-test WITHIN each nested family (bernstein, chebychev,
     exponential, power): take the lowest order where adding the next
     parameter is not justified at p < F_THRESHOLD.
  2. Goodness-of-fit ACROSS families: keep every family's choice whose
     chi2 p-value exceeds GOF_THRESHOLD.

Everything that survives goes into the RooMultiPdf envelope. That is the
point of discrete profiling -- the envelope covers the ignorance about
the shape, so picking two models by hand defeats it.

Combine contract
----------------
`bwsHrare.py` looks these up by name, with tag = "<cat>_<bin>_<year>":

    multipdf_<tag>_bkg        RooMultiPdf
    pdfindex_<tag>            RooCategory
    multipdf_<tag>_bkg_norm   RooRealVar, floating
    datahist_<tag>            RooDataHist, FULL range

Note the tag is built exactly as bwsHrare does, so an empty --bin gives
the double underscore ("VBFcat__Run3") that the existing datacards have.

RooMultiPdf is the only Combine class needed, so everything except
--multipdf runs in the plain conda environment.

Usage
-----
    # fits, plots, JSON -- conda
    python bkgFit.py -c ggHcat -y Run3 -b incl

    # ...and the workspace -- inside the Combine container
    python bkgFit.py -c ggHcat -y Run3 -b incl --multipdf

    python bkgFit.py -c ggHcat -y Run3 -b "" --multipdf   # VBFcat__Run3 style
"""

import argparse
import json
import os
from datetime import datetime, timezone

import ROOT

import sigFit as sf
import pdfDefinitions as pdefs


ROOT.gROOT.SetBatch(True)
ROOT.RooMsgService.instance().setGlobalKillBelow(ROOT.RooFit.WARNING)
for _i in range(ROOT.RooMsgService.instance().numStreams()):
    ROOT.RooMsgService.instance().getStream(_i).removeTopic(ROOT.RooFit.Plotting)

# ---------------------------------------------------------------------------
# blinding
# ---------------------------------------------------------------------------
SR_LO, SR_HI = 115.0, 135.0        # blinded Higgs signal region
R_SB_LOW = "sbLow"                 # 110 - 115
R_SB_HIGH = "sbHigh"               # 135 - 150
R_SB = "sideband"                  # the union, what the fit uses
R_FULL = "full"
R_SR = "sr"

# ---------------------------------------------------------------------------
# model selection
# ---------------------------------------------------------------------------
# F-test: add a parameter only if 2*dNLL against chi2_(dnpar) gives
# p < F_THRESHOLD. 0.05 is the usual choice.
F_THRESHOLD = 0.05
# Goodness of fit: keep a candidate whose chi2 p-value exceeds this.
GOF_THRESHOLD = 0.01
# Candidates to consider at all. bern0 is a flat line and bern1 a straight
# one -- kept because the F-test should be allowed to choose them.
CANDIDATES = None                  # None -> everything in pdfDefinitions

# normalisation is extrapolated from the sideband, so give it room
NORM_RANGE = (0.2, 5.0)            # x the extrapolated value

# ---------------------------------------------------------------------------
# plotting
# ---------------------------------------------------------------------------
# y-title offset for the TALL upper pad; the ratio-pad offset is too small
# there and puts the title on top of the tick labels
MAIN_Y_TITLE_OFFSET = 1.45


# ---------------------------------------------------------------------------
# input
# ---------------------------------------------------------------------------

def load_data(category, binMVA, year, nbins):
    """Full-range data histogram, from prepareFits."""
    from prepareFits import getHisto
    h = getHisto(nbins, sf.XLOW, sf.XHIGH, False, category, year,
                 False, binMVA or "incl")
    h.SetDirectory(0)
    return h


def setup_ranges(x):
    """Name the sideband pieces, their union, the SR and the full range."""
    x.setRange(R_FULL, sf.XLOW, sf.XHIGH)
    x.setRange(R_SB_LOW, sf.XLOW, SR_LO)
    x.setRange(R_SB_HIGH, SR_HI, sf.XHIGH)
    x.setRange(R_SR, SR_LO, SR_HI)
    # RooFit takes a comma-separated list of named ranges, so the union is
    # expressed at use time rather than as its own range
    return f"{R_SB_LOW},{R_SB_HIGH}"


def sideband_bins(h):
    """Bin indices inside the sideband, for counting and for chi2."""
    out = []
    for i in range(1, h.GetNbinsX() + 1):
        c = h.GetBinCenter(i)
        if sf.XLOW <= c < SR_LO or SR_HI <= c < sf.XHIGH:
            out.append(i)
    return out


def full_yield(n_sb, fit):
    """Full-window yield extrapolated from the sideband count.

    Every PDF is evaluated normalised over the FULL window (the fit's
    Range() does not persist on x), so this one number is what the chi2,
    the plot and the workspace norm all use.
    """
    f_sb = fit["sb_fraction"]
    return n_sb / f_sb if f_sb > 0 else n_sb


# ---------------------------------------------------------------------------
# fitting
# ---------------------------------------------------------------------------

def fit_candidate(x, data, name, tag, sb_range):
    """Fit one candidate on the sideband only.

    Returns a dict with the masked NLL (for the F-test), the sideband
    chi2 and its p-value (for the GOF cut), and the fraction of the PDF
    that lies in the sideband (for the normalisation extrapolation).
    """
    pdf, keep = pdefs.create_bkg_pdf(x, tag, name)
    npar = pdefs.BKG_NPAR[name]

    opts = [ROOT.RooFit.Minimizer("Minuit2"), ROOT.RooFit.Strategy(2),
            ROOT.RooFit.Save(True),
            # the masked likelihood: only the sideband contributes, and the
            # PDF is normalised over the sideband
            ROOT.RooFit.Range(sb_range),
            ROOT.RooFit.SumW2Error(True),
            ROOT.RooFit.PrintLevel(-1)]
    pdf.fitTo(data, *opts)
    res = pdf.fitTo(data, *opts)

    # sideband fraction: what part of the PDF lies in the sideband once
    # normalised over the full range. Needed to turn the observed sideband
    # count into a full-range yield.
    nset = ROOT.RooArgSet(x)
    itg_sb = pdf.createIntegral(nset, ROOT.RooFit.NormSet(nset),
                                ROOT.RooFit.Range(sb_range)).getVal()

    return {"name": name, "pdf": pdf, "keep": keep, "npar": npar,
            "nll": res.minNll(), "status": res.status(),
            "cov_qual": res.covQual(), "edm": res.edm(),
            "sb_fraction": itg_sb, "result": res}


def sideband_chi2(x, h, pdf, sb_idx, norm_full, npar):
    """chi2 over SIDEBAND BINS ONLY.

    A goodness-of-fit that included bins the fit never saw would not be a
    goodness-of-fit. pdf.getVal(nset) is normalised over the FULL window,
    so the expectation in a bin is norm_full * density * binwidth, with
    norm_full = full_yield(n_sb, fit) = n_sb / f_sb.
    """
    nset = ROOT.RooArgSet(x)
    saved = x.getVal()
    chi2 = 0.0
    used = 0
    for i in sb_idx:
        e = h.GetBinError(i)
        o = h.GetBinContent(i)
        if e <= 0:
            continue
        x.setVal(h.GetBinCenter(i))
        dens = pdf.getVal(nset)          # normalised over the full window
        exp = norm_full * dens * h.GetBinWidth(i)
        if exp <= 0:
            continue
        chi2 += (o - exp) ** 2 / e ** 2
        used += 1
    x.setVal(saved)
    ndf = used - npar
    p = ROOT.TMath.Prob(chi2, ndf) if ndf > 0 else 0.0
    return {"chi2": chi2, "ndf": ndf, "chi2_ndf": chi2 / ndf if ndf > 0 else -1,
            "pvalue": p, "n_bins": used}


def ftest(fits, family):
    """Lowest order in `family` whose successor is not justified.

    2*(NLL_low - NLL_high) is compared with chi2 at the difference in
    parameter count. Both NLLs come from the same masked likelihood, so
    the comparison is valid; mixing a blinded and an unblinded fit here
    would not be.
    """
    names = [n for n in family if n in fits]
    if not names:
        return None, []
    steps = []
    chosen = names[0]
    for lo, hi in zip(names, names[1:]):
        dnll = fits[lo]["nll"] - fits[hi]["nll"]
        dnpar = fits[hi]["npar"] - fits[lo]["npar"]
        if dnpar <= 0:
            continue
        # a higher order cannot fit worse; guard against a failed fit
        chi2 = max(2.0 * dnll, 0.0)
        p = ROOT.TMath.Prob(chi2, dnpar)
        steps.append({"from": lo, "to": hi, "2dNLL": chi2,
                      "dnpar": dnpar, "pvalue": p,
                      "justified": p < F_THRESHOLD})
        if p < F_THRESHOLD:
            chosen = hi
        else:
            break
    return chosen, steps


# ---------------------------------------------------------------------------
# plots
# ---------------------------------------------------------------------------

def plot_envelope(x, data, fits, envelope, best, label, outbase, year,
                  sb_range, save_pdf, nbins, n_sb):
    """Data with the SR suppressed, every envelope member drawn across the
    full range, and a ratio panel against the best model.

    Two RooFit details matter here:
      * data.plotOn(..., CutRange(sideband)) ZEROES the SR bins but still
        draws them, so the SR points are removed from the RooHist
        afterwards. That removal is the blinding on the plot.
      * Each curve gets its yield explicitly, Normalization(n_sb / f_sb,
        NumEvent), AND NormRange(full). Without an explicit NormRange,
        RooFit sees the "fitrange" attribute fitTo(Range=sideband) left on
        x and normalises the curve to the sideband instead, putting it
        1/f_sb too high. The explicit yield is the same extrapolation the
        chi2 and the workspace norm use, so all three agree.
    """
    frame = x.frame(ROOT.RooFit.Title(""))
    data.plotOn(frame, ROOT.RooFit.Name("dat"),
                ROOT.RooFit.CutRange(sb_range),
                ROOT.RooFit.MarkerStyle(sf.DATA_MARKER_STYLE),
                ROOT.RooFit.MarkerSize(sf.DATA_MARKER_SIZE),
                ROOT.RooFit.MarkerColor(ROOT.kBlack),
                ROOT.RooFit.LineColor(ROOT.kBlack),
                ROOT.RooFit.DataError(ROOT.RooAbsData.SumW2))

    # blinding: drop the zeroed SR points CutRange leaves behind
    hdat = frame.getHist("dat")
    for i in reversed(range(hdat.GetN())):
        if SR_LO <= hdat.GetPointX(i) < SR_HI:
            hdat.RemovePoint(i)

    cols = [sf.fit_color(), ROOT.kAzure + 2, ROOT.kGreen + 2,
            ROOT.kOrange + 7, ROOT.kViolet + 1, ROOT.kTeal + 3]
    for i, nm in enumerate(envelope):
        fits[nm]["pdf"].plotOn(
            frame, ROOT.RooFit.Name(nm),
            ROOT.RooFit.Range(R_FULL),
            # explicit, or RooFit silently normalises to the last fit
            # range (the sideband) and the curve comes out 1/f_sb too high
            ROOT.RooFit.NormRange(R_FULL),
            ROOT.RooFit.Normalization(full_yield(n_sb, fits[nm]),
                                      ROOT.RooAbsReal.NumEvent),
            ROOT.RooFit.LineColor(cols[i % len(cols)]),
            ROOT.RooFit.LineStyle(ROOT.kSolid if nm == best
                                  else ROOT.kDashed),
            ROOT.RooFit.LineWidth(2))

    bw = (sf.XHIGH - sf.XLOW) / nbins
    frame.SetTitle(f";;Events / {bw:.2f} GeV")
    sf._style_ratio_axis(frame.GetYaxis(), sf.RATIO_Y_TITLE_OFFSET)
    frame.GetYaxis().SetTitleOffset(MAIN_Y_TITLE_OFFSET)
    frame.GetXaxis().SetLabelSize(0)
    frame.GetXaxis().SetTitleSize(0)
    frame.SetMinimum(0.0)
    frame.SetMaximum(1.35 * frame.GetMaximum())

    # ratio of data and of each member to the best, over the full range
    hdat = frame.getHist("dat")
    cbest = frame.getCurve(best)
    rgraphs = []
    g = ROOT.TGraphErrors()
    k = 0
    for i in range(hdat.GetN()):
        xv, yv = hdat.GetPointX(i), hdat.GetPointY(i)
        ey = hdat.GetErrorY(i)
        f = cbest.average(xv - 0.5 * bw, xv + 0.5 * bw)
        if yv <= 0 or f <= 0:
            continue
        g.SetPoint(k, xv, yv / f)
        g.SetPointError(k, 0.0, ey / f)
        k += 1
    g.SetMarkerStyle(sf.DATA_MARKER_STYLE)
    g.SetMarkerSize(sf.DATA_MARKER_SIZE)
    g.SetMarkerColor(ROOT.kBlack)
    g.SetLineColor(ROOT.kBlack)
    rgraphs.append(g)

    for i, nm in enumerate(envelope):
        if nm == best:
            continue
        c = frame.getCurve(nm)
        gr = ROOT.TGraph()
        npt = 0
        xv = sf.XLOW + 0.5 * bw
        while xv < sf.XHIGH:
            fb = cbest.average(xv - 0.5 * bw, xv + 0.5 * bw)
            fn = c.average(xv - 0.5 * bw, xv + 0.5 * bw)
            if fb > 0:
                gr.SetPoint(npt, xv, fn / fb)
                npt += 1
            xv += bw
        gr.SetLineColor(cols[(i) % len(cols)])
        gr.SetLineStyle(ROOT.kDashed)
        gr.SetLineWidth(2)
        rgraphs.append(gr)

    rframe = ROOT.TH1F("rframe", "", 1, sf.XLOW, sf.XHIGH)
    rframe.SetDirectory(0)
    rframe.SetStats(0)
    rframe.SetMinimum(0.80)
    rframe.SetMaximum(1.20)
    rframe.SetTitle(f";m_{{#mu#mu}} [GeV];Data / {best}")
    sf._style_ratio_axis(rframe.GetXaxis(), sf.RATIO_X_TITLE_OFFSET)
    sf._style_ratio_axis(rframe.GetYaxis(), sf.RATIO_Y_TITLE_OFFSET)

    canvas, outer, pad1, pad2 = sf.make_canvas_pads("c_bkg")
    pad1.cd()
    frame.Draw()

    # shade the blinded region so it cannot be mistaken for empty data
    box = ROOT.TBox(SR_LO, 0.0, SR_HI, frame.GetMaximum())
    box.SetFillColorAlpha(ROOT.kGray, 0.35)
    box.SetLineWidth(0)
    box.Draw("same")
    frame.Draw("same")

    keep = [box]
    latex = ROOT.TLatex()
    latex.SetNDC()
    latex.SetTextFont(42)
    latex.SetTextSize(sf.ANN_TEXT_SIZE)
    dy = sf.LINE_SPACING * sf.ANN_TEXT_SIZE
    for i, ln in enumerate(sf.header_lines(label) +
                           [f"blinded {SR_LO:g}-{SR_HI:g} GeV"]):
        latex.DrawLatex(sf.PAD_LEFT_MARGIN + 0.05, 0.85 - i * dy, ln)
    keep.append(latex)

    leg = ROOT.TLegend(1 - sf.PAD_RIGHT_MARGIN - 0.32,
                       1 - sf.PAD_TOP_MARGIN - 0.04
                       - 0.052 * (len(envelope) + 1),
                       1 - sf.PAD_RIGHT_MARGIN - 0.02,
                       1 - sf.PAD_TOP_MARGIN - 0.04)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.SetTextFont(42)
    leg.SetTextSize(0.030)
    leg.AddEntry(frame.getHist("dat"), "Data (sideband)", "pe")
    for nm in envelope:
        leg.AddEntry(frame.getCurve(nm),
                     nm + ("  (best)" if nm == best else ""), "l")
    leg.Draw()
    keep.append(leg)
    keep.append(sf.cms_label(pad1, year))

    pad2.cd()
    rframe.Draw()
    for gr in rgraphs[1:]:
        gr.Draw("L SAME")
    rgraphs[0].Draw("P SAME")
    one = ROOT.TLine(sf.XLOW, 1.0, sf.XHIGH, 1.0)
    one.SetLineColor(sf.fit_color())
    one.SetLineWidth(2)
    one.Draw("same")
    rbox = ROOT.TBox(SR_LO, 0.80, SR_HI, 1.20)
    rbox.SetFillColorAlpha(ROOT.kGray, 0.35)
    rbox.SetLineWidth(0)
    rbox.Draw("same")
    rframe.Draw("axis same")
    keep += [one, rbox] + rgraphs

    os.makedirs(os.path.dirname(outbase), exist_ok=True)
    canvas.Update()
    canvas.SaveAs(f"{outbase}_envelope.png")
    if save_pdf:
        canvas.SaveAs(f"{outbase}_envelope.pdf")
    del keep, pad1, pad2, outer, canvas


# ---------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("-c", "--category", default="ggHcat",
                    choices=sorted(sf.MC_LIST))
    ap.add_argument("-b", "--bin", default="incl",
                    help="BDT bin label; pass \"\" to reproduce the "
                         "double-underscore tag the existing datacards use")
    ap.add_argument("-y", "--year", default="Run3")
    ap.add_argument("--candidates", nargs="+", default=None,
                    help="restrict the candidate list (default: all in "
                         "pdfDefinitions)")
    ap.add_argument("--multipdf", action="store_true",
                    help="build the RooMultiPdf and write the workspace. "
                         "Needs the Combine container; everything else "
                         "runs in conda.")
    ap.add_argument("--wsdir", default="WS_BKG")
    ap.add_argument("--plotdir",
                    default=os.path.expanduser(
                        "~/public_html/HmumuFits/bkg_fits"))
    ap.add_argument("--pdf", action="store_true")
    ap.add_argument("--bins-per-gev", type=float, default=None)
    args = ap.parse_args()

    sf.setup_style()
    if args.bins_per_gev:
        sf.BINS_PER_GEV = args.bins_per_gev
    nbins = int((sf.XHIGH - sf.XLOW) * sf.BINS_PER_GEV)

    # exactly bwsHrare's tag, double underscore and all
    tag = f"{args.category}_{args.bin}_{args.year}"
    label = f"{args.category} {args.bin or 'incl'} data"

    print(f"category : {args.category}"
          f"\nbin      : {args.bin or '(empty)'}"
          f"\nyear     : {args.year}"
          f"\ntag      : {tag}"
          f"\nwindow   : {sf.XLOW:g}-{sf.XHIGH:g} GeV, {nbins} bins"
          f"\nblinded  : {SR_LO:g}-{SR_HI:g} GeV"
          f"\nsideband : {sf.XLOW:g}-{SR_LO:g} + {SR_HI:g}-{sf.XHIGH:g} GeV"
          f"\nmultipdf : {args.multipdf}")

    # ---- data ----------------------------------------------------------
    h = load_data(args.category, args.bin, args.year, nbins)
    x = ROOT.RooRealVar(f"mh{args.category}", "m_{#mu#mu}",
                        sf.XLOW, sf.XHIGH)
    sb_range = setup_ranges(x)
    data = ROOT.RooDataHist(f"datahist_{tag}", "data",
                            ROOT.RooArgList(x), h)

    sb_idx = sideband_bins(h)
    n_sb = sum(h.GetBinContent(i) for i in sb_idx)
    n_full = h.Integral(1, h.GetNbinsX())
    print(f"\nevents   : {n_full:.0f} in the window, {n_sb:.0f} in the "
          f"sideband ({100 * n_sb / n_full:.1f} %)")
    print(f"sideband bins used for chi2: {len(sb_idx)} of {nbins}")

    # ---- fit every candidate -------------------------------------------
    names = args.candidates or sorted(pdefs.BKG_FACTORIES)
    print(f"\nfitting {len(names)} candidates on the sideband")
    print(f"{'model':<8}{'npar':>5}{'NLL':>14}{'chi2/ndf':>10}"
          f"{'p(chi2)':>9}{'f_sb':>8}{'status':>8}")
    fits = {}
    for nm in names:
        try:
            f = fit_candidate(x, data, nm, tag, sb_range)
        except Exception as exc:
            print(f"{nm:<8}  failed: {exc}")
            continue
        f.update(sideband_chi2(x, h, f["pdf"], sb_idx,
                               full_yield(n_sb, f), f["npar"]))
        fits[nm] = f
        flag = "" if f["status"] == 0 and f["cov_qual"] == 3 else "  <-- check"
        print(f"{nm:<8}{f['npar']:>5}{f['nll']:>14.2f}"
              f"{f['chi2_ndf']:>10.3f}{f['pvalue']:>9.3f}"
              f"{f['sb_fraction']:>8.3f}{f['status']:>8}{flag}")

    if not fits:
        raise SystemExit("no candidate fitted")

    # ---- F-test within families ----------------------------------------
    print(f"\nF-test within families (p < {F_THRESHOLD} justifies the "
          f"extra parameter)")
    chosen, ftests = {}, {}
    for fam, members in pdefs.BKG_FAMILIES.items():
        pick, steps = ftest(fits, members)
        if pick is None:
            continue
        chosen[fam] = pick
        ftests[fam] = steps
        for st in steps:
            print(f"  {fam:<12} {st['from']:>6} -> {st['to']:<6} "
                  f"2dNLL={st['2dNLL']:7.2f} dnpar={st['dnpar']} "
                  f"p={st['pvalue']:.4f}  "
                  f"{'justified' if st['justified'] else 'NOT justified'}")
        print(f"  {fam:<12} -> {pick}")

    # ---- GOF across families -------------------------------------------
    envelope = [nm for nm in chosen.values()
                if fits[nm]["pvalue"] > GOF_THRESHOLD]
    rejected = [nm for nm in chosen.values() if nm not in envelope]
    envelope = sorted(set(envelope), key=lambda n: -fits[n]["pvalue"])
    print(f"\nGOF cut (p > {GOF_THRESHOLD})")
    for nm in sorted(set(chosen.values())):
        print(f"  {nm:<8} p={fits[nm]['pvalue']:.4f}  "
              f"{'IN' if nm in envelope else 'rejected'}")
    if not envelope:
        raise SystemExit("no candidate passed the goodness-of-fit cut; "
                         "inspect the fits before lowering GOF_THRESHOLD")
    best = envelope[0]
    print(f"\nenvelope : {', '.join(envelope)}   (best: {best})")

    # ---- normalisation, extrapolated -----------------------------------
    # The fit only saw the sideband, so the full-range yield is the
    # observed sideband count divided by the fraction of the PDF that lies
    # in the sideband. NOT the observed total, which would defeat the
    # blinding, and not the sideband count, which would be too low.
    f_sb = fits[best]["sb_fraction"]
    norm_val = full_yield(n_sb, fits[best])
    spread = [full_yield(n_sb, fits[nm]) for nm in envelope
              if fits[nm]["sb_fraction"] > 0]
    print(f"\nnorm     : {norm_val:.1f} from {best} "
          f"(sideband {n_sb:.0f} / f_sb {f_sb:.4f})")
    if len(spread) > 1:
        print(f"           envelope spread {min(spread):.1f} - "
              f"{max(spread):.1f}  "
              f"({100 * (max(spread) - min(spread)) / norm_val:.1f} %)")
        print("           this spread is the extrapolation uncertainty; "
              "Combine\n           profiles it through the floating norm "
              "and pdfindex")

    # ---- plots ---------------------------------------------------------
    plotdir = os.path.join(args.plotdir, args.category)
    outbase = os.path.join(plotdir, f"bkg_{tag}")
    plot_envelope(x, data, fits, envelope, best, label, outbase,
                  args.year, sb_range, args.pdf, nbins, n_sb)

    # ---- workspace ------------------------------------------------------
    os.makedirs(args.wsdir, exist_ok=True)
    wsfile = os.path.join(args.wsdir, f"Bkg_{tag}_workspace.root")
    if args.multipdf:
        if not pdefs.load_combine():
            raise SystemExit("--multipdf needs libHiggsAnalysisCombinedLimit; "
                             "run inside the Combine container")
        if not hasattr(ROOT, "RooMultiPdf"):
            raise SystemExit("RooMultiPdf not available after loading the "
                             "Combine library")

        cat = ROOT.RooCategory(f"pdfindex_{tag}", f"pdfindex_{tag}")
        store = ROOT.RooArgList()
        for nm in envelope:
            store.add(fits[nm]["pdf"])
        multi = ROOT.RooMultiPdf(f"multipdf_{tag}_bkg", "multipdf",
                                 cat, store)
        nrm = ROOT.RooRealVar(f"multipdf_{tag}_bkg_norm",
                              f"multipdf_{tag}_bkg_norm", norm_val,
                              NORM_RANGE[0] * norm_val,
                              NORM_RANGE[1] * norm_val)
        nrm.setConstant(False)          # Combine profiles it

        ws = ROOT.RooWorkspace("w", "w")
        getattr(ws, "import")(cat)
        getattr(ws, "import")(multi)
        getattr(ws, "import")(nrm)
        getattr(ws, "import")(data)     # FULL range, for data_obs
        ws.writeToFile(wsfile)
        print(f"\nwrote {wsfile}")
        print(f"  multipdf_{tag}_bkg        ({len(envelope)} pdfs)")
        print(f"  pdfindex_{tag}")
        print(f"  multipdf_{tag}_bkg_norm   = {norm_val:.1f}, floating")
        print(f"  datahist_{tag}            full range, "
              f"{n_full:.0f} events")
    else:
        print(f"\nno workspace written (--multipdf not given): RooMultiPdf "
              f"needs the Combine container")

    # ---- json ----------------------------------------------------------
    rec = {
        "provenance": {
            "written_utc": datetime.now(timezone.utc)
                           .isoformat(timespec="seconds"),
            "git_commit": sf.git_commit(),
            "category": args.category, "bin": args.bin,
            "year": args.year, "tag": tag,
            "window": [sf.XLOW, sf.XHIGH], "n_bins": nbins,
            "blind_region": [SR_LO, SR_HI],
            "f_threshold": F_THRESHOLD, "gof_threshold": GOF_THRESHOLD,
            "multipdf_written": bool(args.multipdf),
        },
        "data": {"n_full": n_full, "n_sideband": n_sb,
                 "n_sideband_bins": len(sb_idx)},
        "candidates": {nm: {k: v for k, v in f.items()
                            if k not in ("pdf", "keep", "result")}
                       for nm, f in fits.items()},
        "ftest": ftests,
        "family_choice": chosen,
        "envelope": envelope,
        "rejected_by_gof": rejected,
        "best": best,
        "normalization": {"value": norm_val, "sideband_count": n_sb,
                          "sb_fraction": f_sb,
                          "envelope_spread": [min(spread), max(spread)]
                          if len(spread) > 1 else None},
    }
    jf = os.path.join(args.wsdir, f"Bkg_{tag}_params.json")
    with open(jf, "w") as fh:
        json.dump(rec, fh, indent=2, sort_keys=True)
    print(f"wrote {jf}")
    print(f"plots {outbase}_envelope.png")


if __name__ == "__main__":
    main()
