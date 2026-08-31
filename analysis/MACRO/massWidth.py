import ROOT

ROOT.gStyle.SetOptStat(0)
ROOT.gROOT.SetBatch(True)

# -------------------------
# Input
# -------------------------

xvals = [0, 1, 2]
xlabels = ["bdt0", "bdt1", "bdt2"]

data = {
    "ggH" : [1.77, 1.369, 1.31],
    "VBF" : [1.7, 1.56, 1.46],
    #
    "TTH" : [1.6, 1.6, 1.8],
    "TTL" : [1.7, 1.6],
    #
    "VH"  : [1.9, 1.6, 1.9],
    "VL"  : [2.0, 1.6, 1.5],
    "Zinv": [1.8, 1.8, 1.9],
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
frame.SetMinimum(1.)
frame.SetMaximum(2.5)

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

    n = len(vals)
    g = ROOT.TGraph(n)

    for i in range(n):
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

c.SaveAs("~/public_html/HMUMU_FITS/JUL14/massWidth_vs_BDT_jul14.png")
c.SaveAs("~/public_html/HMUMU_FITS/JUL14/massWidth_vs_BDT_jul14.pdf")

#input("enter")
