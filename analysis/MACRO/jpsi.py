import ROOT
from array import array


ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetOptStat(0)

# Input file and tree
treename = "events"

lumis = {
    "2022C": [5.0104],
    "2022D": [2.9700],
    "2022E": [5.8070],
    "2022F": [17.7819],
    "2022G": [3.0828],
    "2023C": [17.794],
    "2023D": [9.451],
    "2024C": [7.24],
    "2024D": [7.96],
    "2024E": [11.32],
    "2024F": [27.76],
    "2024G": [37.77],
    "2024H": [5.44],
    "2024I": [11.47],
#    "2025B": [0.26],
    "2025C": [20.78],
    "2025D": [25.29],
    "2025E": [14.00],
    "2025F": [3.04], ### recorde 16.33 whil3 golden 3.04
    "2025G": [1],    ## to update
}

# Mapping between label (year-period) and sample code
labels = {
    #    "2022C": [-11],
    #    "2022D": [2.9700],
    #    "2022E": [-115_22022],
    "2022F": "-116_22022",
    "2022G": "-117_22022",
    "2023C": "-123_12023",
    "2023D": "-131_22023",
    "2024C": "-141_2024",
    "2024D": "-142_2024",
    "2024E": "-143_2024",
    "2024F": "-144_2024",
    "2024G": "-145_2024",
    "2024H": "-146_2024",
    "2024I": "-147_2024",
#    "2025B": "-161_2025",
    "2025C": "-162_2025",
    "2025D": "-163_2025",
    "2025E": "-164_2025",
    "2025F": "-165_2025", ### recorde 16.33 whil3 golden 3.04
#    "2025G": [1],    ## to update
}

def plot(mymc, region=""):

    filename = "/work/submit/mariadlf/HmumuRun3/ROOTFILES/snapshot_mc_"+mymc+"_ggHcat.root"

    # Open the ROOT file
    file = ROOT.TFile.Open(filename)
    if not file or file.IsZombie():
        raise RuntimeError(f"❌ Could not open {filename}")
    
    tree = file.Get(treename)
    if not tree:
        raise RuntimeError(f"❌ Could not find tree '{treename}' in file")

  # Define cuts
    if region == "central":
        selection = "abs(Muon1_eta) < 0.9 && abs(Muon2_eta) < 0.9"
    elif region == "forward":
        selection = "abs(Muon1_eta) > 1.4 && abs(Muon2_eta) > 1.4"
    elif region == "middle":
        selection = "(abs(Muon1_eta) > 0.9 && abs(Muon1_eta) < 1.4) && (abs(Muon2_eta) > 0.9 && abs(Muon2_eta) < 1.4)"
    else:
        selection = ""  # no cut
    
    # Create histogram
    hist = ROOT.TH1F("hHiggsMass", "Higgs Candidate Corrected Mass;H_{corr} [GeV];Events", 400, 2.8, 3.4)

    # Fill histogram
    tree.Draw("HiggsCandCorrMass >> hHiggsMass", selection, "goff")

    # Define Gaussian fit function in the desired range (2–4 GeV)
    gaus = ROOT.TF1("gaus", "gaus", 2.7, 3.2)

    ## Parameters are: [0]=amplitude, [1]=mean, [2]=sigma
    #gaus.SetParameters(hist.GetMaximum(), 3.1, 0.1)

    # Define combined signal + background model:
    # f(x) = Gaussian (signal) + 2nd-degree polynomial (background)
    fit_func = ROOT.TF1("fit_func", "gaus(0) + pol2(3)", 2.8, 3.4)
    
    # Set initial Gaussian parameters (amplitude, mean, sigma)
    fit_func.SetParameter(0, hist.GetMaximum())  # amplitude
    fit_func.SetParameter(1, 3.0969)                # mean
    fit_func.SetParameter(2, 0.1)                # sigma
    
    # Initialize polynomial coefficients (background)
    fit_func.SetParameter(3, 10)   # p0
    fit_func.SetParameter(4, 0)    # p1
    fit_func.SetParameter(5, 0)    # p2
    
    # Optional: give names for clarity in output
    fit_func.SetParName(0, "Amp")
    fit_func.SetParName(1, "Mean")
    fit_func.SetParName(2, "Sigma")
    fit_func.SetParName(3, "p0")
    fit_func.SetParName(4, "p1")
    fit_func.SetParName(5, "p2")
    
    # Fit histogram
    fit_result = hist.Fit(fit_func, "SR")  # S=return result, R=respect range
    
    # Extract parameters
    amp = fit_result.Parameter(0)
    mean = fit_result.Parameter(1)
    sigma = fit_result.Parameter(2)
    mean_err = fit_result.ParError(1)
    sigma_err = fit_result.ParError(2)
    
    # Print results
    print(f"✅ Gaussian Fit Results:")
    print(f"   Mean  = {mean:.3f} ± {mean_err:.3f}")
    print(f"   Sigma = {sigma:.3f} ± {sigma_err:.3f}")
    
    # Draw result
    canvas = ROOT.TCanvas("canvas", "Signal+Background Fit", 800, 600)
    hist.SetLineColor(ROOT.kBlack)
    hist.Draw("")
    
    fit_func.SetLineColor(ROOT.kRed)
    fit_func.SetLineWidth(2)
    fit_func.Draw("SAME")

    # --- Draw the polynomial component ---
    poly = ROOT.TF1("poly_part", "[0] + [1]*x + [2]*x*x", 2.5, 4.0)
    for i in range(3):
        poly.SetParameter(i, fit_func.GetParameter(3 + i))
    poly.SetLineColor(ROOT.kGreen + 2)
    poly.SetLineStyle(3)
    poly.SetLineWidth(2)
    poly.Draw("SAME")
    
    # Add legend
    legend = ROOT.TLegend(0.35, 0.85, 0.85, 0.88)
    legend.AddEntry(hist, "Data", "l")
    legend.AddEntry(fit_func, f"Gauss+Poly Fit (μ={mean:.3f}, σ={sigma:.3f})", "l")
    legend.Draw()
    
    # Save output
    canvas.SaveAs("~/public_html/HMUMU_OCT/muonStudies/fitQual/HiggsCandCorrMass_fit"+mymc+"_"+region+".png")

    print("📊 Plot saved to 'HiggsCandCorrMass_fit.png'")

    return mean,sigma,amp

if __name__ == "__main__":

    means, sigmas, xvals, xlabels = [], [], [], []
    
    for i, (label, code) in enumerate(labels.items()):
        print(f"Processing {label} ({code}) ...")
        #        mean, sigma, amp = plot(code)
        mean, sigma, amp = plot(code,"middle")        
        #        mean, sigma, amp = plot(code,"central")
        #        mean, sigma, amp = plot(code,"forward")        
        means.append(mean)
        sigmas.append(sigma)
        xvals.append(i)
        xlabels.append(label)
        
    # ----------------------------------------------------------
    # Plot mean vs. sample label
    # ----------------------------------------------------------
    graph_mean = ROOT.TGraph(len(xvals), array('d', xvals), array('d', means))
    graph_mean.SetTitle("Gaussian Mean vs Data Period;Data period;Mean [GeV]")
    graph_mean.SetMarkerStyle(21)
    graph_mean.SetMarkerColor(ROOT.kBlue)
    
    canvas_mean = ROOT.TCanvas("cmean", "", 900, 600)
    graph_mean.Draw("APL")
    
    # Add custom labels on X axis
    xaxis = graph_mean.GetXaxis()
    for i, lbl in enumerate(xlabels):
        xaxis.SetBinLabel(xaxis.FindBin(i), lbl)
    xaxis.SetLabelSize(0.04)

    # Fix Y-axis range
    graph_mean.GetYaxis().SetRangeUser(3.1-0.015, 3.1+0.01)

    canvas_mean.SetBottomMargin(0.15)

    # --- Add horizontal reference line at y = 3.0969 ---
    y_ref = 3.0969
    line = ROOT.TLine(xvals[0] - 0.5, y_ref, xvals[-1] + 0.5, y_ref)
    line.SetLineColor(ROOT.kRed)
    line.SetLineStyle(2)
    line.SetLineWidth(2)
    line.Draw("same")

    # legend for clarity
    legend = ROOT.TLegend(0.45, 0.75, 0.70, 0.85)
    legend.AddEntry(graph_mean, "Fit mean", "p")
    legend.AddEntry(line, "PDG J/#psi mass = 3.0969 GeV", "l")
    legend.SetBorderSize(0)
    legend.SetFillStyle(0)
    legend.SetTextSize(0.045)  # <-- Larger font
    legend.Draw()

#    canvas_mean.SaveAs("~/public_html/HMUMU_OCT/muonStudies/mean_vs_sample.png")    
    canvas_mean.SaveAs("~/public_html/HMUMU_OCT/muonStudies/mean_vs_sample_middle.png")
#    canvas_mean.SaveAs("~/public_html/HMUMU_OCT/muonStudies/mean_vs_sample_central.png")
#    canvas_mean.SaveAs("~/public_html/HMUMU_OCT/muonStudies/mean_vs_sample_forward.png")

    # ----------------------------------------------------------
    # Plot sigma vs. sample label
    # ----------------------------------------------------------
    graph_sigma = ROOT.TGraph(len(xvals), array('d', xvals), array('d', sigmas))
    graph_sigma.SetTitle("Gaussian sigmas vs Data Period;Data period; sigmas [GeV]")
    graph_sigma.SetMarkerStyle(21)
    graph_sigma.SetMarkerColor(ROOT.kBlue)

    canvas_sigma = ROOT.TCanvas("csigma", "", 900, 600)
    graph_sigma.Draw("APL")
    
    # Add custom labels on X axis
    xaxis = graph_sigma.GetXaxis()
    for i, lbl in enumerate(xlabels):
        xaxis.SetBinLabel(xaxis.FindBin(i), lbl)
    xaxis.SetLabelSize(0.04)

    # Fix Y-axis range
    graph_sigma.GetYaxis().SetRangeUser(0.0929-0.1, 0.0929+0.1)    
    
    canvas_sigma.SetBottomMargin(0.15)

#    canvas_sigma.SaveAs("~/public_html/HMUMU_OCT/muonStudies/sigma_vs_period.png")    
    canvas_sigma.SaveAs("~/public_html/HMUMU_OCT/muonStudies/sigma_vs_period_middle.png")
#    canvas_sigma.SaveAs("~/public_html/HMUMU_OCT/muonStudies/sigma_vs_period_central.png")
#    canvas_sigma.SaveAs("~/public_html/HMUMU_OCT/muonStudies/sigma_vs_period_forward.png")
