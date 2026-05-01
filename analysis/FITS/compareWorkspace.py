import ROOT

ROOT.gROOT.SetBatch()

signal_suffix = {
    "VBFcat":  "qqH",
    "VHcat":   "VH",
    "VLcat":   "VH",
    "Zinvcat": "VH",
#    "ggHcat":  "ggHcat",
    "TTHcat":  "ttH",
    "TTLcat":  "ttH",
}


signal_binning = {
    "VBFcat":  ["bdt0","bdt1"],
    "TTLcat":  ["bdt0","bdt1"],
    "VLcat":   ["bdt0","bdt1"],    
    "VHcat":   ["bdt0","bdt1","bdt2"],
    "TTHcat":  ["bdt0","bdt1","bdt2"],
    "Zinvcat": ["bdt0","bdt1","bdt2"],        
#    "ggHcat":  "ggHcat",
}

#2 bins VBFcat , VLcat, TTLcat
#3 bins Zinv, VHcat, TTHcat

myWSdir="WS_APR27"

def compare(category):

    files = {}
    workspaces = {}

    for b in signal_binning[category]:
        fname = f"{myWSdir}/Signal_{category}_{b}_Run3_workspace.root"
        f = ROOT.TFile.Open(fname)
        files[b] = f
        workspaces[b] = f.Get("w")
    
        #w.Print()
    
    # Observable
    x = workspaces[signal_binning[category][0]].var("mh"+category)
    
    # Define range
    x.setRange("plotRange", 110, 140)
    frame = x.frame(ROOT.RooFit.Range("plotRange"),ROOT.RooFit.Title(f"m_#mu#mu {category}"))

    colors = [ROOT.kRed, ROOT.kBlue, ROOT.kGreen+2]

    plot_objects = []
    labels = []
    
    for i, b in enumerate(signal_binning[category]):
        w = workspaces[b]

        pdf = w.pdf(f"crystal_ball_{category}_{b}_Run3_{signal_suffix[category]}")
        
        obj = pdf.plotOn(frame,
                         ROOT.RooFit.LineColor(colors[i]),
                         ROOT.RooFit.Range("plotRange"),
                         ROOT.RooFit.NormRange("plotRange"))
        
        plot_objects.append(obj)

        # --- get parameters ---
        mu = w.var(f"cb_mu_{category}_{b}_Run3_{signal_suffix[category]}")
        sigma = w.var(f"cb_sigma_{category}_{b}_Run3_{signal_suffix[category]}")
        
        if not mu or not sigma:
            labels.append(f"{b} (missing)")
            continue

        val_mu = mu.getVal()
        val_sigma = sigma.getVal()
        
        print(b,
              "mu =", val_mu,
              "sigma =", val_sigma,
              "sigma/mu =", val_sigma/val_mu)
        
        labels.append(f"{b}: #sigma = {val_sigma:.3f}")
    
        
        # Canvas
        c = ROOT.TCanvas("c","CB comparison",800,600)
        frame.Draw()
        
        line = ROOT.TLine(125, 0, 125, frame.GetMaximum())
        line.SetLineColor(ROOT.kBlack)
        line.SetLineStyle(2)  # dashed
        line.SetLineWidth(2)
        line.Draw("same")
        
        # Legend
        leg = ROOT.TLegend(0.65,0.7,0.88,0.88)
        for i in range(len(labels)):
            leg.AddEntry(frame.getObject(i), labels[i], "l")
        leg.Draw()
        
        c.SaveAs("~/public_html/HMUMU_FITS/APRgood/compareCB_"+category+"_Run3.png")

if __name__ == "__main__":

    compare("VBFcat")
    compare("Zinvcat")
    compare("VLcat")
    compare("VHcat")

    compare("TTHcat")
    compare("TTLcat")    
