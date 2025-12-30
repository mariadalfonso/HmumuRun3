import ROOT

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

outDir = "~/public_html/HMUMU_OCT/vbfStudies"
#filename = "/work/submit/mariadlf/HmumuRun3/ROOTFILES/snapshot_mc_10_2024_VBFcat.root"  # EWK-H 
#filename = "/work/submit/mariadlf/HmumuRun3/ROOTFILES/snapshot_mc_101_2024_VBFcat.root"  # EWK-Z
filename = "/work/submit/mariadlf/HmumuRun3/ROOTFILES/snapshot_mc_103_2024_VBFcat.root"  # EWK-Z

#Siglabel="H"
#Siglabel="Z"
Siglabel="Zqcd"


# ---------------------------------------
# Variables for 1D plots
# ---------------------------------------
VARS = {
    100: "Mjj",
    101: "dEtaJJ",
    102: "RPt",
    103: "ZepVar",
    104: "jetVBF1_Pt",
    105: "jetVBF2_Pt",
    106: "jetVBF1_Eta",
    107: "jetVBF2_Eta",
    108: "jetVBF1_Phi",
    109: "jetVBF2_Phi",
#    110: "dPhiJJ",
#    111: "minDR_jetVBF1_Mu",
#    112: "minDR_jetVBF2_Mu",        
}

# ---------------------------------------
# Valid 2D combinations
# ---------------------------------------
TWOD_PAIRS = [
    ("Mjj", "dEtaJJ"),
    ("Mjj", "jetVBF1_Pt"),
    ("Mjj", "jetVBF2_Pt"),        
]

# ---------------------------------------
# Cuts
# ---------------------------------------
CUT_LHE     = "jetVBF1_LHE != -1 && jetVBF2_LHE != -1"
CUT_NONLHE  = "jetVBF1_LHE == -1 || jetVBF2_LHE == -1"

# ---------------------------------------
# Binning map per variable
# ---------------------------------------
BINNING = {
    "Mjj":       (60, 0, 1500),
    "dEtaJJ":    (60, 0, 8),
    "dPhiJJ":    (60, -3.2, 3.2),
    "RPt":       (60, 0, 1),
    "ZepVar":    (60, -5, 5),
    "jetVBF1_Pt": (60, 0, 500),
    "jetVBF2_Pt": (60, 0, 500),
    "jetVBF1_Eta": (60, -5, 5),
    "jetVBF2_Eta": (60, -5, 5),
    "jetVBF1_Phi": (60, -3.2, 3.2),
    "jetVBF2_Phi": (60, -3.2, 3.2),
    "minDR_jetVBF1_Mu": (111, 0., 6.3),
    "minDR_jetVBF2_Mu": (112, 0., 6.3),
}

# ===================================================
# Draw 1D overlays
# ===================================================

def draw_single_var(tree, var):

    print(f"→ Drawing 1D overlay for {var}")

    nb, xmin, xmax = BINNING[var]
    h_all = ROOT.TH1F(f"h_all_{var}", f"{var};{var};Events", nb, xmin, xmax)
    h_lhe = ROOT.TH1F(f"h_lhe_{var}", f"{var};{var};Events", nb, xmin, xmax)

    tree.Draw(f"{var} >> h_all_{var}", "", "goff")
    tree.Draw(f"{var} >> h_lhe_{var}", CUT_LHE, "goff")

    # Style histograms
    h_all.SetLineColor(ROOT.kGray + 2)
    h_all.SetLineWidth(2)
    h_lhe.SetLineColor(ROOT.kRed)
    h_lhe.SetLineWidth(2)

    # Canvas with two pads: top for histograms, bottom for ratio
    c = ROOT.TCanvas(f"c_{var}", "", 800, 700)

    pad1 = ROOT.TPad("pad1", "pad1", 0, 0.3, 1, 1.0)
    pad2 = ROOT.TPad("pad2", "pad2", 0, 0.0, 1, 0.3)

    pad1.SetBottomMargin(0.02)  # reduce bottom margin for ratio
    pad2.SetTopMargin(0.05)
    pad2.SetBottomMargin(0.3)
    pad1.Draw()
    pad2.Draw()

    # -------------------------
    # Top pad: overlay histograms
    # -------------------------
    pad1.cd()
    h_all.SetTitle(f"{var};{var};Events")
    h_all.Draw("HIST")
    h_lhe.Draw("HIST SAME")

    leg = ROOT.TLegend(0.6, 0.75, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.AddEntry(h_all, "All jets", "l")
    leg.AddEntry(h_lhe, "LHE matched", "l")
    leg.Draw()

    # -------------------------
    # Bottom pad: ratio
    # -------------------------
    pad2.cd()
    h_ratio = h_lhe.Clone(f"h_ratio_{var}")
    h_ratio.Divide(h_all)
    h_ratio.SetMarkerStyle(20)
    h_ratio.SetTitle("")
    h_ratio.GetYaxis().SetRangeUser(0.5,1.5)
    h_ratio.GetYaxis().SetTitle("LHE / All")
    h_ratio.GetYaxis().SetNdivisions(505)
    h_ratio.GetYaxis().SetTitleSize(0.15)
    h_ratio.GetYaxis().SetTitleOffset(0.35)
    h_ratio.GetXaxis().SetTitle(var)
    h_ratio.GetXaxis().SetTitleSize(0.15)
    h_ratio.GetXaxis().SetLabelSize(0.12)
    h_ratio.GetYaxis().SetLabelSize(0.12)
    h_ratio.Draw("P")

    lineZero = ROOT.TLine(h_ratio.GetXaxis().GetXmin(), 1.,  h_ratio.GetXaxis().GetXmax(), 1.)
    lineZero.SetLineColor(11)
    lineZero.Draw("same")

    
    c.SaveAs(f"{outDir}/1D_{var}_{Siglabel}.png")
    print(f"✅ Saved {var} overlay with ratio.")

def draw_single_varOLD(tree, var):

    print(f"→ Drawing 1D overlay for {var}")

    nb, xmin, xmax = BINNING[var]
    h_all = ROOT.TH1F(f"h_all_{var}", f"{var};{var};Events", nb, xmin, xmax)
    h_lhe = ROOT.TH1F(f"h_lhe_{var}", f"{var};{var};Events", nb, xmin, xmax)

    tree.Draw(f"{var} >> h_all_{var}", "", "goff")
    tree.Draw(f"{var} >> h_lhe_{var}", CUT_LHE, "goff")

    h_all.SetLineColor(ROOT.kGray + 2)
    h_all.SetLineWidth(2)
    h_lhe.SetLineColor(ROOT.kRed)
    h_lhe.SetLineWidth(2)

    c = ROOT.TCanvas(f"c_{var}", "", 800, 600)
    h_all.Draw("HIST")
    h_lhe.Draw("HIST SAME")

    leg = ROOT.TLegend(0.6, 0.75, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.AddEntry(h_all, "All jets", "l")
    leg.AddEntry(h_lhe, "LHE matched", "l")
    leg.Draw()

    c.SaveAs(f"{outDir}/1D_{var}_{Siglabel}.png")


# ===================================================
# Draw 2D
# ===================================================
def draw_2D(tree, xvar, yvar, cut, tag):

    nbx, xmin, xmax = BINNING[xvar]
    nby, ymin, ymax = BINNING[yvar]

    hname = f"h2_{xvar}_{yvar}_{tag}"
    htitle = f"{yvar} vs {xvar} ({tag});{xvar};{yvar}"

    h2 = ROOT.TH2F(hname, htitle, nbx, xmin, xmax, nby, ymin, ymax)

    tree.Draw(f"{yvar}:{xvar} >> {hname}", cut, "goff")

    return h2


# ===================================================
# Main
# ===================================================
def main():

    f = ROOT.TFile.Open(filename)
    tree = f.Get("events")

    print("\n=== Producing 1D plots ===\n")
    for var in VARS.values():
        draw_single_var(tree, var)

    print("\n=== Producing 2D plots ===\n")
    for (xvar, yvar) in TWOD_PAIRS:
        print(f"2D: {yvar} vs {xvar}")

        h_LHE = draw_2D(tree, xvar, yvar, CUT_LHE, "LHE")
        h_NON = draw_2D(tree, xvar, yvar, CUT_NONLHE, "nonLHE")

        # LHE plot
        c1 = ROOT.TCanvas("", "", 900, 800)
        c1.SetRightMargin(0.15)
        h_LHE.Draw("COLZ")
        c1.SaveAs(f"{outDir}/2D_{yvar}_vs_{xvar}_LHE_{Siglabel}.png")

        # non-LHE plot
        h_NON.Draw("COLZ")
        c1.SaveAs(f"{outDir}/2D_{yvar}_vs_{xvar}_nonLHE_{Siglabel}.png")

        # Comparison
        comp = ROOT.TCanvas("", "", 1800, 800)
        comp.Divide(2,1)
        comp.cd(1); ROOT.gPad.SetRightMargin(0.15); h_LHE.Draw("COLZ")
        comp.cd(2); ROOT.gPad.SetRightMargin(0.15); h_NON.Draw("COLZ")
        comp.SaveAs(f"{outDir}/2D_{yvar}_vs_{xvar}_compare_{Siglabel}.png")

    print("\n✅ All plots created.\n")

main()
