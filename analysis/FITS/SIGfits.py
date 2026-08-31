1;95;0c binX in ["bdt0", "bdt1", "bdt2","bdt3"]:        import ROOT

from prepareFits import getHisto
from pdfDefinitions import *
from fitUtils import *

ROOT.gROOT.SetBatch()
combine_base = "/code/HiggsAnalysis/CombinedLimit"
ROOT.gSystem.Load("/code/HiggsAnalysis/CombinedLimit/build/lib/libHiggsAnalysisCombinedLimit.so")
# NB before loading the singularity, do 'conda deactivate'

blinded=False
doMultiPdf=True
do_bias_study = True

workspaceName = 'WS_JUL23'

## for GF and VBF
xlowRange = 110
xhighRange = 150

category_suffix = {
    "ggHcat":  "ggHcat",
    "VBFcat":  "VBFcat",
    "VHcat":   "VHcat",
    "VLcat":   "VLcat",
    "Zinvcat": "Zinvcat",
    "TTHcat":  "TTHcat",
    "TTLcat":  "TTLcat",
}

# TO DO BETTER
mc_list = {
    "VLcat":  ["VH","ttH"], #qqH and ggH are no relevant stat
    "TTLcat":  ["VH","ttH"],
    "Zinvcat":  ["qqH","VH","ttH"],
    "VHcat":  ["VH","ttH"],
    "TTHcat":  ["VH","ttH"],
    #
    "ggHcat":  ["ggH","qqH","VH","ttH"],
    "VBFcat":  ["ggH","qqH","VH","ttH"],
}

bkg_model_selection = {
    "ggHcat": {"model1": "pow3", "model2": "bern3"},
    "VBFcat": {"model1": "exp3", "model2": "bern2"},
    "VHcat": {"model1": "bern2", "model2": "exp1"},
    "VLcat": {"model1": "bern2", "model2": "exp1"},
    "Zinvcat": {"model1": "bern2", "model2": "exp1"},
    "TTHcat": {"model1": "bern2", "model2": "exp1"},
    "TTLcat": {"model1": "bern2", "model2": "exp1"},
}

def finalWorkspace(w, data, model, norm, year, doMultiPdf=False, tag='',storedPdfs=''):

    if doMultiPdf:
        pdf_cat = ROOT.RooCategory("pdfindex_"+tag,"pdfindex_"+tag)
        pdf_bkg = ROOT.RooMultiPdf("multipdf_"+tag+"_bkg","multipdf",pdf_cat,storedPdfs)
        getattr(w,'import')(pdf_bkg)
    else:
        # Import model and all its components into the workspace
        getattr(w,'import')(model)

    getattr(w,'import')(norm)
    print('integral signal/BKG = ',norm.Print()) # these are temporaty workspace so can do this

    # Import data into the workspace
    getattr(w,'import')(data)

    # Print workspace contents
    w.Print()

    return w


def setVar(tag, lowBlind='-1', highBlind='-1'):

    suffix = category_suffix[tag]
    x = ROOT.RooRealVar(f"mh{suffix}", "m_{#mu,#mu}", xlowRange, xhighRange)

    x.setRange("full", xlowRange, xhighRange)
    if lowBlind!='-1':
        x.setRange("left", xlowRange, lowBlind)
        x.setRange("right", highBlind, xhighRange)

    print('RooRealVar DONE')

    return x

def  fitSig(tag_ , year, binMVA):

    tag = tag_+"_"+binMVA+"_"+year

    x = setVar(tag_)

    doSignal = True
    doLog = False

    # Create a empty workspace (one for all signal)
    w = ROOT.RooWorkspace("w", "workspace")

    pdfDef = PDFDefinitions()

    for sig in mc_list[tag_]:
    
        data_full = getHisto(10*int(xhighRange - xlowRange), xlowRange, xhighRange, doLog, tag_, year, doSignal, binMVA, sig)
        print('getHisto  DONE')

        data = ROOT.RooDataHist('datahist'+tag, 'data', ROOT.RooArgList(x), data_full)

        # -----------------------------------------------------------------------------

        # Create signal PDF
        sig_pdfs = pdfDef.create_signal_pdfs(x, f'_{tag}_{sig}')
        model = sig_pdfs['pdf']

        # -----------------------------------------------------------------------------

        model.fitTo(data,ROOT.RooFit.Minimizer("Minuit2"),ROOT.RooFit.Strategy(2),ROOT.RooFit.Range("full"))

        # Here we will plot the results
        if True:
            # Create canvas with two pads
            canvas = ROOT.TCanvas("canvas", "canvas", 800, 800)
            pad1 = ROOT.TPad("pad1", "Top Pad", 0, 0.35, 1, 1)
            pad2 = ROOT.TPad("pad2", "Bottom Pad", 0, 0, 1, 0.35)

            pad1.SetBottomMargin(0.03)
            pad2.SetTopMargin(0.05)
            pad2.SetBottomMargin(0.3)

            pad1.Draw()
            pad2.Draw()

            # ========== TOP PAD ==========
            pad1.cd()

            titleSTR = "mH_"+tag+'_'+sig+"_"+str(year)
            plotFrameWithNormRange = x.frame(ROOT.RooFit.Title(titleSTR))

            # Plot only the blinded data, and then plot the PDF over the full range as well as both sidebands
            data.plotOn(plotFrameWithNormRange, ROOT.RooFit.Name("fit_data"))
            model.plotOn(plotFrameWithNormRange, ROOT.RooFit.LineColor(2), ROOT.RooFit.Range("full"), ROOT.RooFit.NormRange("full"), ROOT.RooFit.LineStyle(10), ROOT.RooFit.Name("fit_curve"))
            model.paramOn(plotFrameWithNormRange, ROOT.RooFit.Layout(0.75,0.99,0.85))
            plotFrameWithNormRange.getAttText().SetTextSize(0.02);

            plotFrameWithNormRange.Draw()

            # ========== BOTTOM PAD (ratio plot) ==========
            pad2.cd()

            # Compute residuals = (data - fit) / fit
            # Get RooHist of data and RooCurve of model
            plotFrameWithNormRange.Print("v")
            dataHist = plotFrameWithNormRange.getHist("fit_data")
            fitCurve = plotFrameWithNormRange.getCurve("fit_curve")

            # Create a new histogram for the ratio
            ratioHist = ROOT.TH1D("ratio", "", dataHist.GetN(), xlowRange, xhighRange)

            for i in range(dataHist.GetN()):
                x_val = dataHist.GetX()[i]
                y_data = dataHist.GetY()[i]
                y_fit = fitCurve.Eval(x_val)

                if y_fit != 0:
                    ratio = (y_data - y_fit) / y_fit
                else:
                    ratio = 0

                ratioHist.SetBinContent(i + 1, ratio)
                ratioHist.SetBinError(i + 1, dataHist.GetErrorY(i) / y_fit if y_fit != 0 else 0)

            ratioHist.GetYaxis().SetTitle("Ratio")
            ratioHist.GetYaxis().SetTitleSize(0.08)
            ratioHist.GetYaxis().SetLabelSize(0.08)
            ratioHist.GetYaxis().SetTitleOffset(0.4)
            ratioHist.GetXaxis().SetTitle("m_{#mu,#mu}")
            ratioHist.GetXaxis().SetTitleSize(0.1)
            ratioHist.GetXaxis().SetLabelSize(0.08)
            ratioHist.GetYaxis().SetRangeUser(-1.0,1.0)
            ratioHist.SetLineColor(ROOT.kBlack)
            ratioHist.SetMarkerStyle(20)
            ratioHist.Draw("EP")

            # Draw horizontal line at 0
            line = ROOT.TLine(xlowRange, 0.0, xhighRange, 0.0)
            line.SetLineColor(ROOT.kRed)
            line.SetLineWidth(2)
            line.SetLineStyle(2)
            line.Draw()

            # ========== Save the canvas ==========
            canvas.Draw()
            htmldir = "~/public_html/HMUMU_FITS/JUL23"
            canvas.SaveAs(htmldir+"/signal_"+tag+'_'+sig+"_"+str(year)+".png")
            chi2_ndf = plotFrameWithNormRange.chiSquare()
            print("Chi² / ndf =", chi2_ndf)

        # -----------------------------------------------------------------------------
        # -----------------------------------------------------------------------------

        binLow = data_full.GetBin(1) #contains the first bin with low-edge
        binUp = data_full.GetBin(int(xhighRange-xlowRange)*10)  # second to last bin contains the upper-edge

        norm_SR = data_full.Integral(binLow, binUp)
        print("--------------------------")
        print("SIG norm",norm_SR)
        print(' binX1 = ',data_full.GetXaxis().GetBinLowEdge(binLow)," - ",data_full.GetXaxis().GetBinUpEdge(binLow))
        print(' binX2 = ',data_full.GetXaxis().GetBinLowEdge(binUp)," - ",data_full.GetXaxis().GetBinUpEdge(binUp))
        print("--------------------------")

        # Sig_norm constant (the r is added by the text to workspace)
        Sig_norm = ROOT.RooRealVar(model.GetName()+ "_norm", model.GetName()+ "_norm", norm_SR) # no range means contants

        # -----------------------------------------------------------------------------
        # -----------------------------------------------------------------------------
        # Create workspace, import data and model

        # Set parameters constant
        for param in sig_pdfs['params'].values():
            param.setConstant()
#        cb_mu.setConstant()
#        cb_sigma.setConstant()
#        cb_alphaL.setConstant()
#        cb_alphaR.setConstant()
#        cb_nL.setConstant()
#        cb_nR.setConstant()
        Sig_norm.setConstant()

        w = finalWorkspace(w, data, model, Sig_norm, year)

    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # Save workspace in file

    w.writeToFile(workspaceName+"/Signal_"+tag+"_workspace.root")

def  fitBkg(tag_, year, binMVA):

    tag = tag_+"_"+binMVA+"_"+year

    print("tag = ", tag)
    # Create a empty workspace (one for all signal)
    w = ROOT.RooWorkspace("w", "workspace")

    lowBlind = 120
    highBlind = 130
    x = setVar(tag_, lowBlind, highBlind)
    nBins = 10*int(xhighRange-xlowRange)

    # Uniform binning for the full range
    x.setBins(nBins)

    doSignal = False
    doLog = False
    data_full = getHisto(nBins, xlowRange, xhighRange, doLog, tag_, year, doSignal, binMVA)

    data = ROOT.RooDataHist('datahist_'+tag, 'data', ROOT.RooArgList(x), data_full)
    blindedData = data.reduce(ROOT.RooFit.CutRange("left,right"))

    data_reduced_manual = data_full.Clone()

    # get bin indices for the blinded region
    bin_low = data_reduced_manual.FindBin(lowBlind)
    bin_high = data_reduced_manual.FindBin(highBlind)

    # zero out all bins in [lowBlind, highBlind]
    for i in range(bin_low, bin_high + 1):
        data_reduced_manual.SetBinContent(i, 0.0)
        data_reduced_manual.SetBinError(i, 0.0)

    data_reduced = ROOT.RooDataHist('datahistReduce_'+tag, 'dataReduced', ROOT.RooArgList(x), data_reduced_manual)

    # -----------------------------------------------------------------------------

    # Create all background PDFs
    pdfDef = PDFDefinitions()
    bkg_pdfs_result = pdfDef.create_bkg_pdfs(x, tag)
    all_pdfs = bkg_pdfs_result['pdfs']
    all_params = bkg_pdfs_result['params']
    all_components = bkg_pdfs_result['components']
    all_formulas = bkg_pdfs_result['formulas']

    # Select models for this category
    model1_name = bkg_model_selection[tag_]["model1"]
    model2_name = bkg_model_selection[tag_]["model2"]

    model1 = all_pdfs[model1_name]
    model2 = all_pdfs[model2_name]

    # -----------------------------------------------------------------------------

    storedPdfs = ROOT.RooArgList("store_"+tag)

    cmdList = ROOT.RooLinkedList()
    cmdList.Add(ROOT.RooFit.Minimizer("Minuit2"))
    cmdList.Add(ROOT.RooFit.Strategy(2))
    cmdList.Add(ROOT.RooFit.Range("full"))
    cmdList.Add(ROOT.RooFit.Save(True))

    if blinded: fitresults = model1.fitTo(blindedData,cmdList)
    else: fitresults = model1.fitTo(data,cmdList)

    if doMultiPdf:
        storedPdfs.add(model1)
        if blinded: fitresults2 = model2.fitTo(blindedData,cmdList)
        else: fitresults2 = model2.fitTo(data,cmdList)
        storedPdfs.add(model2)  # extra PDF

    # -----------------------------------------------------------------------------

    binLow = data_full.GetBin(1) #contains the first bin with low-edge
    #    binUp = data_full.GetBin(int(xhighRange-xlowRange))  # second to last bin contains the upper-edge
    binUp = data_full.GetNbinsX()

    norm_range = data_full.Integral( binLow, binUp )
    print("--------------------------")
    print("NORM BKG",norm_range)
    print(' binX1 = ',data_full.GetXaxis().GetBinLowEdge(binLow)," - ",data_full.GetXaxis().GetBinUpEdge(binLow))
    print(' binX2 = ',data_full.GetXaxis().GetBinLowEdge(binUp)," - ",data_full.GetXaxis().GetBinUpEdge(binUp))
    print("--------------------------")

    if doMultiPdf:
        BKG_norm = ROOT.RooRealVar("multipdf_"+tag+"_bkg"+"_norm", model1.GetName()+"_norm", norm_range, 0.5*norm_range, 2*norm_range)
    else:
        BKG_norm = ROOT.RooRealVar(model1.GetName()+ "_norm", model1.GetName()+ "_norm", norm_range, 0.5*norm_range, 2*norm_range)

    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # Make plot out of the frame

    # Here we will plot the results
    canvas = ROOT.TCanvas("canvas", "canvas", 800, 800)
    #canvas.Divide(2, 1)

    titleSTR = "mH_"+tag+"_"+str(year)+" -- "
    plotFrameWithNormRange = x.frame(ROOT.RooFit.Title(titleSTR))

    # Plot only the blinded data, and then plot the PDF over the full range as well as both sidebands

    if blinded:
#        data.plotOn(plotFrameWithNormRange)
        data.plotOn(plotFrameWithNormRange, ROOT.RooFit.MarkerColor(ROOT.kWhite), ROOT.RooFit.LineColor(ROOT.kWhite))

        model1.plotOn(plotFrameWithNormRange, ROOT.RooFit.Components(model1.GetName()), ROOT.RooFit.Range("full"), ROOT.RooFit.NormRange("full"), ROOT.RooFit.LineColor(ROOT.kRed)) ;
        model2.plotOn(plotFrameWithNormRange, ROOT.RooFit.Components(model2.GetName()), ROOT.RooFit.Range("full"), ROOT.RooFit.NormRange("full"), ROOT.RooFit.LineColor(ROOT.kBlue)) ;
        data_reduced.plotOn(plotFrameWithNormRange)
    else:
        data.plotOn(plotFrameWithNormRange)
        model1.plotOn(plotFrameWithNormRange, ROOT.RooFit.Components(model1.GetName()), ROOT.RooFit.Range("full"), ROOT.RooFit.NormRange("full"), ROOT.RooFit.LineColor(ROOT.kRed)) ;
        model2.plotOn(plotFrameWithNormRange, ROOT.RooFit.Components(model2.GetName()), ROOT.RooFit.Range("full"), ROOT.RooFit.NormRange("full"), ROOT.RooFit.LineColor(ROOT.kBlue)) ;
        varName = x.GetName()
        name1 = model1.GetName()+"_Norm["+varName+"]_Comp["+model1.GetName()+"]_Range[full]_NormRange[full]"
        name2 = model2.GetName()+"_Norm["+varName+"]_Comp["+model2.GetName()+"]_Range[full]_NormRange[full]"
        print(name1)
        chi2_1 = plotFrameWithNormRange.chiSquare(name1,"h_"+data.GetName(),fitresults.floatParsFinal().getSize())
        if doMultiPdf: chi2_2 = plotFrameWithNormRange.chiSquare(name2,"h_"+data.GetName(),fitresults2.floatParsFinal().getSize())
#        plotFrameWithNormRange.Print("v")
        print('--------------------')
        if doMultiPdf: print(model2.GetName(),"    chi2/ndof=",round(chi2_2,2)," ndof",fitresults2.floatParsFinal().getSize())
        print(model1.GetName(),"    chi2/ndof=",round(chi2_1,2)," ndof",fitresults.floatParsFinal().getSize())
        print('--------------------')

        fileToWrite="preselection_"+tag+"_"+"_"+str(year)+".txt"
        with open(fileToWrite, "a") as f:
            str1 = model1.GetName()+"    chi2/ndof="+str(round(chi2_1,2))+" ndof"+str(fitresults.floatParsFinal().getSize())+"\n"
            f.write(str1)
            if doMultiPdf: str2 = model2.GetName()+"    chi2/ndof="+str(round(chi2_2,2))+" ndof"+str(fitresults2.floatParsFinal().getSize())+"\n"
            if doMultiPdf: f.write(str2)

#    model.paramOn(plotFrameWithNormRange, RooFit.Layout(0.6,0.99,0.95))
#    plotFrameWithNormRange.getAttText().SetTextSize(0.02);

    plotFrameWithNormRange.Draw()
#    hresid = plotFrameWithNormRange.residHist()

    offsetY = 0.75*data_full.GetMaximum()
    latex = ROOT.TLatex()
    latex.SetTextColor(ROOT.kRed)
    latex.SetTextSize(0.04)
    latex.DrawLatex(130 ,offsetY + 0.10*data_full.GetMaximum(), model1.GetName())
    latex.SetTextColor(ROOT.kBlue)
    latex.DrawLatex(130 ,offsetY + 0.20*data_full.GetMaximum(), model2.GetName())

    canvas.Draw()
    htmldir = "~/public_html/HMUMU_FITS/JUL23"
    if blinded: canvas.SaveAs(htmldir+"/bkg_"+tag+"_"+str(year)+"blinded.png")
    else: canvas.SaveAs(htmldir+"/bkg_"+tag+"_"+str(year)+".png")

    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # Create workspace, import data and model
    # Save workspace in file

    w = finalWorkspace(w, data, model1, BKG_norm, year, doMultiPdf, tag, storedPdfs)

    if not blinded:
        w.writeToFile(workspaceName+"/Bkg_"+tag+"_workspace.root")

if __name__ == "__main__":

    for year in ['Run3']:

        for binX in ["bdt0", "bdt1", "bdt2","bdt3"]:
            fitSig('VBFcat',year,binX)
            fitBkg('VBFcat',year,binX)

        for binX in ["bdt0", "bdt1", "bdt2"]:
            fitSig('VHcat',year,binX)
            fitBkg('VHcat',year,binX)

        for binX in ["bdt0", "bdt1", "bdt2"]:
            fitSig('ggHcat',year,binX)
            fitBkg('ggHcat',year,binX)

            fitSig('Zinvcat',year,binX)
            fitBkg('Zinvcat',year,binX)

            fitSig('TTHcat',year,binX)
            fitBkg('TTHcat',year,binX)

            fitSig('VLcat',year,binX)
            fitBkg('VLcat',year,binX)

#        for binX in ["bdt0", "bdt1","incl"]:
        for binX in ["bdt0", "bdt1"]:

            fitSig('TTLcat',year,binX)
            fitBkg('TTLcat',year,binX)
