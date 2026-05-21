import ROOT

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

# -------------------------
# Input
# -------------------------

xvals = [0, 1, 2]
xlabels = ["bdt0", "bdt1", "bdt2"]

data = {
    "ggH" : [1.621, 1.415, 1.087],
    "VBF" : [1.600, 1.554, 1.517],
    #
    "TTH" : [1.517, 1.561, 1.739],
    "TTL" : [1.599, 1.643, 1.457],
    #
    "VH"  : [1.838, 1.962, 1.904],
    "VL"  : [1.676, 1.603, 1.723],
    "Zinv": [1.647, 1.705, 1.912],
    
}

colors = [
    ROOT.kGray+1,
    ROOT.kGreen+2,
    #
    ROOT.kBlue+1,
    ROOT.kCyan+2,
    #
    ROOT.kRed+1,
    ROOT.kMagenta+1,
    ROOT.kOrange+7,
    #    
]

# -------------------------
# Canvas
# -------------------------

c = ROOT.TCanvas("c", "Mass width", 900, 700)

frame = ROOT.TH1F("frame", ";BDT category;Mass width [GeV]", 3, -0.5, 2.5)
frame.SetMinimum(0.9)
frame.SetMaximum(2.1)

for i, lab in enumerate(xlabels):
    frame.GetXaxis().SetBinLabel(i+1, lab)

frame.Draw()

# -------------------------
# Graphs
# -------------------------

leg = ROOT.TLegend(0.8, 0.65, 0.9, 0.88)
leg.SetBorderSize(0)

graphs = []

for ic, (name, vals) in enumerate(data.items()):

    g = ROOT.TGraph(3)

    for i in range(3):
        g.SetPoint(i, xvals[i], vals[i])

    g.SetLineColor(colors[ic])
    g.SetMarkerColor(colors[ic])
    g.SetLineWidth(3)
    g.SetMarkerStyle(20)
    g.SetMarkerSize(1.3)

    g.Draw("LP SAME")

    leg.AddEntry(g, name, "lp")
    graphs.append(g)

leg.Draw()

# CMS-style text
txt = ROOT.TLatex()
txt.SetNDC()
txt.SetTextSize(0.04)
txt.DrawLatex(0.15, 0.92, "Run3")

c.SaveAs("~/public_html/HMUMU_FITS/MAY/massWidth_vs_BDT.png")
c.SaveAs("~/public_html/HMUMU_FITS/MAY/massWidth_vs_BDT.pdf")

#input("enter")
