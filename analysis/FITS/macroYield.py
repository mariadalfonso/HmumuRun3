import ROOT

ROOT.gROOT.SetBatch()

redDark  = (191, 34, 41)
redMed   = (237, 41, 57)
redLight = (255, 82, 82)

# =============================
# CREATE GLOBAL HISTOGRAMS
# =============================
h_ggH = ROOT.TH1F("h_ggH", "ggH", 7, 0.5, 7.5)
h_qqH = ROOT.TH1F("h_qqH", "qqH", 7, 0.5, 7.5)
h_VH  = ROOT.TH1F("h_VH",  "VH",  7, 0.5, 7.5)
h_ttH = ROOT.TH1F("h_ttH", "ttH", 7, 0.5, 7.5)

# Set colors once
for h, color in zip(
        [h_qqH, h_ggH, h_VH, h_ttH],
        [redMed, redDark, redLight, redDark]):
    h.SetLineWidth(2)
    h.SetLineColor(ROOT.TColor.GetColor(*color))
    h.SetFillColor(ROOT.TColor.GetColor(*color))


# =============================
# STORAGE OF FRACTIONS (for text)
# =============================
fractions = {1: {}, 2: {}, 3: {}, 4: {}, 5: {}, 6: {}, 7: {}}
    
# =============================
# FUNCTION TO FILL ONE BIN
# =============================
def plotByCat(cat, binN, fileToOpen, myVar):

#    f = ROOT.TFile.Open(f"WS/Signal_{cat}_2024_workspace.root")
    f = ROOT.TFile.Open(fileToOpen)    
    w = f.Get("w")
    if not w:
        raise RuntimeError("Workspace 'w' not found!")

    w.Print()
    # Retrieve normalization variables
    if cat=='TTLcat': n_ggH = 0
    else: n_ggH = w.var(f"{myVar}_ggH_N_SM")
    n_qqH = w.var(f"{myVar}_qqH_N_SM")
    n_VH  = w.var(f"{myVar}_VH_N_SM")
    n_ttH = w.var(f"{myVar}_ttH_N_SM")
    
    if cat=='TTLcat': n_ggHval = 0;
    else: n_ggHval = n_ggH.getVal()
    print('n_qqH = ', n_qqH.getVal())
    print('n_VH = ', n_VH.getVal())
    print('n_ttH = ', n_ttH.getVal())

    # Compute fractions
    total = n_ggHval + n_qqH.getVal() + n_VH.getVal() + n_ttH.getVal()
    f_ggH = n_ggHval / total
    f_qqH = n_qqH.getVal() / total
    f_VH  = n_VH.getVal()  / total
    f_ttH = n_ttH.getVal() / total

    # Fill global histograms (one bin each)
    h_ggH.SetBinContent(binN, f_ggH)
    h_qqH.SetBinContent(binN, f_qqH)
    h_VH.SetBinContent(binN,  f_VH)
    h_ttH.SetBinContent(binN, f_ttH)

    print(f"Filled bin {binN} for {cat} :",
          f"ggH={f_ggH:.3f}, qqH={f_qqH:.3f}, VH={f_VH:.3f}, ttH={f_ttH:.3f}")

    # Save fractions for text later
    fractions[binN]["ggH"] = f_ggH
    fractions[binN]["qqH"] = f_qqH
    fractions[binN]["VH"]  = f_VH
    fractions[binN]["ttH"] = f_ttH

    
# =============================
# MAIN: FILL BOTH BINS
# =============================
if __name__ == "__main__":

    plotByCat("ggHcat", 1, "workspaces_GGH_DEC16/workspace_ggHcat__Run3_workspace.root", "crystal_ball_ggHcat__Run3")
    plotByCat("VBFcat", 2, "workspaces_VBF_DEC16/workspace_VBFcat__Run3_workspace.root", "crystal_ball_VBFcat__Run3")
    plotByCat("VLcat",  3, "workspaces_MINOR_DEC16/workspace_VLcat__Run3_workspace.root", "crystal_ball_VLcat__Run3")
    plotByCat("VHcat",  4, "workspaces_MINOR_DEC16/workspace_VHcat__Run3_workspace.root", "crystal_ball_VHcat__Run3")
    plotByCat("Zinvcat",5, "workspaces_MINOR_DEC16/workspace_Zinvcat__Run3_workspace.root", "crystal_ball_Zinvcat__Run3")
    plotByCat("TTLcat", 6, "workspaces_MINOR_DEC16/workspace_TTLcat__Run3_workspace.root", "crystal_ball_TTLcat__Run3")
    plotByCat("TTHcat", 7, "workspaces_MINOR_DEC16/workspace_TTHcat__Run3_workspace.root", "crystal_ball_TTHcat__Run3")

    # =============================
    # MAKE ONE STACKED PLOT
    # =============================
    stack = ROOT.THStack("stack", "Signal Composition;Category;Fraction")

    stack.Add(h_ggH)
    stack.Add(h_qqH)
    stack.Add(h_VH)
    stack.Add(h_ttH)

    c = ROOT.TCanvas("c", "Stacked Fractions", 900, 700)
    stack.Draw("bar")

    # Fix axes
    stack.GetXaxis().SetBinLabel(1, "ggHcat")
    stack.GetXaxis().SetBinLabel(2, "VBFcat")
    stack.GetXaxis().SetBinLabel(3, "VLcat")
    stack.GetXaxis().SetBinLabel(4, "VHcat")
    stack.GetXaxis().SetBinLabel(5, "Zinvcat")
    stack.GetXaxis().SetBinLabel(6, "TTLcat")
    stack.GetXaxis().SetBinLabel(7, "TTHcat")
    stack.GetYaxis().SetRangeUser(0, 1)
    stack.GetYaxis().SetTitle("Fraction (sum = 100%)")

    # =============================
    # DRAW FRACTION TEXT
    # =============================
    t = ROOT.TLatex()
    t.SetTextFont(42)
    t.SetTextSize(0.035)

    x_positions = {1: 1.0, 2: 2.0, 3: 3.0, 4: 4.0, 5: 5.0, 6: 6.0, 7: 7.0}  # bin centers

    for binN in [1, 2, 3, 4, 5, 6, 7 ]:
        x = x_positions[binN]

        f_ggH = fractions[binN]["ggH"]
        f_qqH = fractions[binN]["qqH"]
        f_VH  = fractions[binN]["VH"]
        f_ttH = fractions[binN]["ttH"]

        # vertical stacking
        y0 = 0
        y1 = y0 + f_ggH
        y2 = y1 + f_qqH
        y3 = y2 + f_VH
        y4 = y3 + f_ttH

        # Draw text at the middle of each segment
        if f_ggH > 0.01:
            t.DrawLatex(x - 0.15, y0 + f_ggH/2, f"ggH {100*f_ggH:.1f}%")
        if f_qqH > 0.01:
            t.DrawLatex(x - 0.15, y1 + f_qqH/2, f"qqH {100*f_qqH:.1f}%")
        if f_VH > 0.01:
            t.DrawLatex(x - 0.15, y2 + f_VH/2,  f"VH {100*f_VH:.1f}%")
        if f_ttH > 0.01:
            t.DrawLatex(x - 0.15, y3 + f_ttH/2, f"ttH {100*f_ttH:.1f}%")
    
    c.SaveAs("~/public_html/HMUMU_FITS/relative_norms_ALLCATS.png")
