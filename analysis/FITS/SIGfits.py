import ROOT

from prepareFits import getHisto

ROOT.gROOT.SetBatch()
combine_base = "/code/HiggsAnalysis/CombinedLimit"
ROOT.gSystem.Load("/code/HiggsAnalysis/CombinedLimit/build/lib/libHiggsAnalysisCombinedLimit.so")
# NB before loading the singularity, do 'conda deactivate'

blinded=False
doMultiPdf=True

workspaceName = 'WS'

## for GF and VBF
xlowRange = 110
xhighRange = 150

category_suffix = {
    "_ggHcat":  "ggHcat",
    "_VBFcat":  "VBFcat",
    "_VHcat":   "VHcat",
    "_VLcat":   "VLcat",
    "_Zinvcat": "Zinvcat",
    "_TTHcat":  "TTHcat",
    "_TTLcat":  "TTLcat",
}

# TO DO BETTER
mc_list = {
    "_ggHcat":  ["ggH","qqH"],
    "_VBFcat":  ["qqH"],
    "_VHcat":   ["VH"],
    "_VLcat":   ["VH"],
    "_Zinvcat": ["Zinvcat"],
    "_TTHcat":  ["TTHcat"],
    "_TTLcat":  ["TTLcat"],
}

def finalWorkspace(w, data, model, norm, doMultiPdf=False, tag='',storedPdfs=''):


    if doMultiPdf:
        pdf_cat = ROOT.RooCategory("pdfindex"+tag,"pdfindex"+tag)
        pdf_bkg = ROOT.RooMultiPdf("multipdf"+tag+"_bkg","multipdf",pdf_cat,storedPdfs)
        getattr(w,'import')(pdf_bkg)
    else:
        # Import model and all its components into the workspace
        getattr(w,'import')(model)

    getattr(w,'import')(norm)
    print('integral signal/BKG = ',norm.Print())

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

def  fitSig(tag , year):

    x = setVar(tag)

    doSignal = True
    doLog = False

    # Create a empty workspace (one for all signal)
    w = ROOT.RooWorkspace("w", "workspace")

    for sig in mc_list[tag]:
    
        data_full = getHisto(10*int(xhighRange - xlowRange), xlowRange, xhighRange, doLog, tag, year, doSignal, sig)
        print('getHisto  DONE')

        data = ROOT.RooDataHist('datahist', 'data', ROOT.RooArgList(x), data_full)

        # -----------------------------------------------------------------------------

        cb_mu = ROOT.RooRealVar('cb_mu'+tag+'_'+sig, 'cb_mu', 125., 125-10. , 125+10.)
        cb_sigma = ROOT.RooRealVar('cb_sigma'+tag+'_'+sig, 'cb_sigma', 3, 0.5, 6.)
        cb_alphaL = ROOT.RooRealVar('cb_alphaL'+tag+'_'+sig, 'cb_alphaL', 2., 0., 5.)
        cb_alphaR = ROOT.RooRealVar('cb_alphaR'+tag+'_'+sig, 'cb_alphaR', 2., 0., 5.)
        cb_nL = ROOT.RooRealVar('cb_nL'+tag+'_'+sig, 'cb_nL', 0., 5.)
        cb_nR = ROOT.RooRealVar('cb_nR'+tag+'_'+sig, 'cb_nR', 0., 10.)

        pdf_crystalball = ROOT.RooDoubleCBFast('crystal_ball'+tag+'_'+sig, 'crystal_ball', x, cb_mu, cb_sigma, cb_alphaL, cb_nL, cb_alphaR, cb_nR)
        model = pdf_crystalball

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

            titleSTR = "mH"+tag+'_'+sig+"_"+str(year)
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
            htmldir = "~/public_html/HMUMU_FITS/OCT"
            canvas.SaveAs(htmldir+"/signal_"+tag+'_'+sig+"_"+str(year)+".png")
            chi2_ndf = plotFrameWithNormRange.chiSquare()
            print("Chi² / ndf =", chi2_ndf)

        # -----------------------------------------------------------------------------
        # -----------------------------------------------------------------------------

        binLow = data_full.GetBin(1) #contains the first bin with low-edge
        binUp = data_full.GetBin(int(xhighRange-xlowRange)*10)  # second to last bin contains the upper-edge

        norm_SR = data_full.Integral(binLow, binUp)

        Sig_norm = ROOT.RooRealVar(model.GetName()+ "_norm", model.GetName()+ "_norm", norm_SR) # no range means contants

        # -----------------------------------------------------------------------------
        # -----------------------------------------------------------------------------
        # Create workspace, import data and model

        cb_mu.setConstant()
        cb_sigma.setConstant()
        cb_alphaL.setConstant()
        cb_alphaR.setConstant()
        cb_nL.setConstant()
        cb_nR.setConstant()
        Sig_norm.setConstant()

        w = finalWorkspace(w, data, model, Sig_norm)

    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # Save workspace in file

    w.writeToFile(workspaceName+"/Signal"+tag+"_"+str(year)+"_workspace.root")

def  fitBkg(tag, year):

    # Create a empty workspace (one for all signal)
    w = ROOT.RooWorkspace("w", "workspace")

    lowBlind = 120
    highBlind = 130
    x = setVar(tag, lowBlind, highBlind)
    nBins = 10*int(xhighRange-xlowRange)

    # Uniform binning for the full range
    x.setBins(nBins)

    doSignal = False
    doLog = False
    data_full = getHisto(nBins, xlowRange, xhighRange, doLog, tag, year, doSignal)
    
    data = ROOT.RooDataHist('datahist'+tag, 'data', ROOT.RooArgList(x), data_full)
    blindedData = data.reduce(ROOT.RooFit.CutRange("left,right"))

    data_reduced_manual = data_full.Clone()

    # get bin indices for the blinded region
    bin_low = data_reduced_manual.FindBin(lowBlind)
    bin_high = data_reduced_manual.FindBin(highBlind)
    # zero out all bins in [lowBlind, highBlind]
    for i in range(bin_low, bin_high + 1):
        data_reduced_manual.SetBinContent(i, 0.0)
        data_reduced_manual.SetBinError(i, 0.0)

    data_reduced = ROOT.RooDataHist('datahistReduce'+tag, 'dataReduced', ROOT.RooArgList(x), data_reduced_manual)

    # -----------------------------------------------------------------------------
    # BERN law
    bern_c0 = ROOT.RooRealVar('bern_c0'+tag, 'bern_c0', 0.5, 0., 1.)
    bern_c1 = ROOT.RooRealVar('bern_c1'+tag, 'bern_c1', 0.1, 0., 1.)
    bern_c2 = ROOT.RooRealVar('bern_c2'+tag, 'bern_c2', 0.1, 0., 1.)
    bern_c3 = ROOT.RooRealVar('bern_c3'+tag, 'bern_c3', 0.1, 0., 1.)
    bern_c4 = ROOT.RooRealVar('bern_c4'+tag, 'bern_c4', 0.5, 0., 5.)
    bern_c5 = ROOT.RooRealVar('bern_c5'+tag, 'bern_c5', 1e-2, 0., 0.1)

    pdf_bern0 = ROOT.RooBernstein('bern0'+tag, 'bern0', x,
                                  ROOT.RooArgList(bern_c0))
    pdf_bern1 = ROOT.RooBernstein('bern1'+tag, 'bern1', x,
                                  ROOT.RooArgList(bern_c0, bern_c1))
    pdf_bern2 = ROOT.RooBernstein('bern2'+tag, 'bern2', x,
                                  ROOT.RooArgList(bern_c0, bern_c1, bern_c2))
    pdf_bern3 = ROOT.RooBernstein('bern3'+tag, 'bern3', x,
                                  ROOT.RooArgList(bern_c0, bern_c1, bern_c2, bern_c3))
    pdf_bern4 = ROOT.RooBernstein('bern4'+tag, 'bern4', x,
                                  ROOT.RooArgList(bern_c0, bern_c1, bern_c2, bern_c3, bern_c4))
    pdf_bern5 = ROOT.RooBernstein('bern5'+tag, 'bern5', x,
                                  ROOT.RooArgList(bern_c0, bern_c1, bern_c2, bern_c3, bern_c4, bern_c5))

#    if blinded: pdf_bern4.selectNormalizationRange(ROOT.RooFit.CutRange("left,right"))

    # -----------------------------------------------------------------------------
    # chebychev law
    chebychev_c0 = ROOT.RooRealVar('chebychev_c0'+tag, 'chebychev_c0', 1.08, -1.1, 10.)
    chebychev_c1 = ROOT.RooRealVar('chebychev_c1'+tag, 'chebychev_c1', 0.4, -1., 1.)
    chebychev_c2 = ROOT.RooRealVar('chebychev_c2'+tag, 'chebychev_c2', 0.01, -0.1, 0.1)
    chebychev_c3 = ROOT.RooRealVar('chebychev_c3'+tag, 'chebychev_c3', 0., -1., 1.) # limit this for the GF

    pdf_chebychev1 = ROOT.RooChebychev("chebychev1"+tag, "chebychev1",x,
                                       ROOT.RooArgList(chebychev_c0,chebychev_c1))

    pdf_chebychev2 = ROOT.RooChebychev("chebychev2"+tag, "chebychev2",x,
                                       ROOT.RooArgList(chebychev_c0,chebychev_c1,chebychev_c2))

    pdf_chebychev3 = ROOT.RooChebychev("chebychev3"+tag,"chebychev3",x,
                                       ROOT.RooArgList(chebychev_c0,chebychev_c1,chebychev_c2,chebychev_c3))

    # -----------------------------------------------------------------------------
    # GAUS law
    gauss_mu = ROOT.RooRealVar('gauss_mu', 'gauss_mu', 120, 100, 140)
    gauss_sigma = ROOT.RooRealVar('gauss_sigma', 'gauss_sigma', 30, 20, 40) #starting value, min max
    pdf_gauss = ROOT.RooGaussian('gauss', 'gauss', x , gauss_mu, gauss_sigma)

    # -----------------------------------------------------------------------------
    # POW law
    # str expressions of formulae
    formula_pow1 = 'TMath::Power(@0, @1)'
    formula_pow2 = '(1.-@1)*TMath::Power(@0,@2) + @1*TMath::Power(@0,@3)'
    formula_pow3 = '(1.-@1-@2)*TMath::Power(@0,@3) + @1*TMath::Power(@0,@4) + @2*TMath::Power(@0,@5)'

    # Variables
    pow_frac1 = ROOT.RooRealVar('frac1', 'frac1', 0.01, 0., 1.)
    pow_frac2 = ROOT.RooRealVar('frac2', 'frac2', 0.01, 0., 1.)
    pow_p1 = ROOT.RooRealVar('p1', 'p1', -2.555, -10., 0.)
    pow_p2 = ROOT.RooRealVar('p2', 'p2', -8., -10., 0.)
    pow_p3 = ROOT.RooRealVar('p3', 'p3', -10., -10., 0.)

    # Power Law PDFs
    pdf_pow1 = ROOT.RooGenericPdf('pow1', 'pow1', formula_pow1,
                                  ROOT.RooArgList(x, pow_p1))
    pdf_pow2 = ROOT.RooGenericPdf('pow2', 'pow2', formula_pow2,
                                  ROOT.RooArgList(x, pow_frac1, pow_p1, pow_p2))
    pdf_pow3 = ROOT.RooGenericPdf('pow3', 'pow3', formula_pow3,
                                  ROOT.RooArgList(x, pow_frac1, pow_frac2, pow_p1, pow_p2, pow_p3))
    
    # -----------------------------------------------------------------------------
    # EXP law
    # these two from Davit
    exp_p1 = ROOT.RooRealVar('exp_p1', 'exp_p1', -0.0207, -0.022, -0.018)
    exp_p2 = ROOT.RooRealVar('exp_p2', 'exp_p2', -1e-2, -10, 10)
    # old param
#    exp_p1 = ROOT.RooRealVar('exp_p1'+tag, 'exp_p1', -0.1, -10, 0)
#    exp_p2 = ROOT.RooRealVar('exp_p2'+tag, 'exp_p2', -1e-2, -10, 0)
    exp_p3 = ROOT.RooRealVar('exp_p3', 'exp_p3', -1e-3, -10, 0)
    exp_c1 = ROOT.RooRealVar('exp_c1', 'exp_c1', 0., 1.)
    exp_c2 = ROOT.RooRealVar('exp_c2', 'exp_c2', 0., 1.)
    exp_c3 = ROOT.RooRealVar('exp_c3', 'exp_c3', 0., 1.)

    pdf_exp1 = ROOT.RooExponential('exp1'+tag, 'exp1', x, exp_p1)
    pdf_single_exp2 = ROOT.RooExponential('single_exp2', 'single_exp2', x, exp_p2)
    pdf_single_exp3 = ROOT.RooExponential('single_exp3', 'single_exp3', x, exp_p3)

    pdf_exp2 = ROOT.RooAddPdf('exp2', 'exp2',
                         ROOT.RooArgList(pdf_exp1, pdf_single_exp2),
                         ROOT.RooArgList(exp_c1, exp_c2))

    pdf_exp3 = ROOT.RooAddPdf('exp3', 'exp3',
                         ROOT.RooArgList(pdf_exp1, pdf_single_exp2, pdf_single_exp3),
                         ROOT.RooArgList(exp_c1, exp_c2, exp_c3))

    pdf_exp1_conv_gauss = ROOT.RooFFTConvPdf('exp1_conv_gauss', 'exp1 (X) gauss', x, pdf_exp1, pdf_gauss)

#   Set #bins to be used for FFT sampling to 10000
#    x.setBins(10000, "cache");
    # Construct landau (x) gauss
#    model = RooFFTConvPdf ("bxg", "bernstein (X) gauss", x, pdf_bern5, pdf_gauss);
#    model = RooFFTConvPdf ("bxg", "bernstein (X) gauss", x, pdf_bern1, pdf_gauss);

    storedPdfs = ROOT.RooArgList("store_"+tag)

    models = {
        "_ggHcat":  {"model1": pdf_chebychev2,   "model2": pdf_bern3},
        "_VBFcat":  {"model1": pdf_bern2,        "model2": pdf_pow3},
        "_VHcat":   {"model1": pdf_bern2,        "model2": pdf_exp1},
        "_VLcat":   {"model1": pdf_bern2,        "model2": pdf_exp1,},
        "_Zinvcat": {"model1": pdf_bern2,        "model2": pdf_exp1},
        "_TTHcat":  {"model1": pdf_chebychev2,   "model2": pdf_exp1},
        "_TTLcat":  {"model1": pdf_chebychev3,   "model2": pdf_exp1},
    }

    if blinded: models[tag]["model1"].fitTo(blindedData,ROOT.RooFit.Minimizer("Minuit2"),ROOT.RooFit.Strategy(2),ROOT.RooFit.Range("full"))
    else: fitresults = models[tag]["model1"].fitTo(data,ROOT.RooFit.Minimizer("Minuit2"),ROOT.RooFit.Strategy(2),ROOT.RooFit.Range("full"),ROOT.RooFit.Save(ROOT.kTRUE))

    if doMultiPdf:
        storedPdfs.add(models[tag]["model1"])
        if blinded: models[tag]["model2"].fitTo(blindedData,ROOT.RooFit.Minimizer("Minuit2"),ROOT.RooFit.Strategy(2),ROOT.RooFit.Range("full"))
        else: fitresults2 = models[tag]["model2"].fitTo(data,ROOT.RooFit.Minimizer("Minuit2"),ROOT.RooFit.Strategy(2),ROOT.RooFit.Range("full"),ROOT.RooFit.Save(ROOT.kTRUE))
        storedPdfs.add(models[tag]["model2"])  # extra PDF

    # -----------------------------------------------------------------------------


    binLow = data_full.GetBin(1) #contains the first bin with low-edge
    binUp = data_full.GetBin(int(xhighRange-xlowRange))  # second to last bin contains the upper-edge

    norm_range = data_full.Integral( binLow, binUp )
    print("--------------------------")
    print("NORM BKG",norm_range)
    print(' binX1 = ',data_full.GetXaxis().GetBinLowEdge(binLow)," - ",data_full.GetXaxis().GetBinUpEdge(binLow))
    print(' binX2 = ',data_full.GetXaxis().GetBinLowEdge(binUp)," - ",data_full.GetXaxis().GetBinUpEdge(binUp))
    print("--------------------------")

    if doMultiPdf:
        BKG_norm = ROOT.RooRealVar("multipdf"+tag+"_bkg"+"_norm", models[tag]["model1"].GetName()+"_norm", norm_range, 0.5*norm_range, 2*norm_range)
    else:
        BKG_norm = ROOT.RooRealVar(models[tag]["model1"].GetName()+ "_norm", models[tag]["model1"].GetName()+ "_norm", norm_range, 0.5*norm_range, 2*norm_range)

    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # Make plot out of the frame

    # Here we will plot the results
    canvas = ROOT.TCanvas("canvas", "canvas", 800, 800)
    #canvas.Divide(2, 1)

    titleSTR = "mH"+tag+"_"+str(year)+" -- "
    plotFrameWithNormRange = x.frame(ROOT.RooFit.Title(titleSTR))

    # Plot only the blinded data, and then plot the PDF over the full range as well as both sidebands

    if blinded:
        data_reduced.plotOn(plotFrameWithNormRange, ROOT.RooFit.MarkerColor(ROOT.kWhite), ROOT.RooFit.LineColor(ROOT.kWhite))

        models[tag]["model1"].plotOn(plotFrameWithNormRange, ROOT.RooFit.Components(models[tag]["model1"].GetName()), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Range("left"), ROOT.RooFit.NormRange("left,right"))
        models[tag]["model1"].plotOn(plotFrameWithNormRange, ROOT.RooFit.Components(models[tag]["model1"].GetName()), ROOT.RooFit.LineColor(ROOT.kRed), ROOT.RooFit.Range("right"), ROOT.RooFit.NormRange("left,right"))
        models[tag]["model2"].plotOn(plotFrameWithNormRange, ROOT.RooFit.Components(models[tag]["model2"].GetName()), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Range("left"), ROOT.RooFit.NormRange("left,right"))
        models[tag]["model2"].plotOn(plotFrameWithNormRange, ROOT.RooFit.Components(models[tag]["model2"].GetName()), ROOT.RooFit.LineColor(ROOT.kBlue), ROOT.RooFit.Range("right"), ROOT.RooFit.NormRange("left,right"))

        data_reduced.plotOn(plotFrameWithNormRange, ROOT.RooFit.Binning("left"))
        data_reduced.plotOn(plotFrameWithNormRange, ROOT.RooFit.Binning("right"))

    else:
        data.plotOn(plotFrameWithNormRange)
        models[tag]["model1"].plotOn(plotFrameWithNormRange, ROOT.RooFit.Components(models[tag]["model1"].GetName()), ROOT.RooFit.Range("full"), ROOT.RooFit.NormRange("full"), ROOT.RooFit.LineColor(ROOT.kRed)) ;
        models[tag]["model2"].plotOn(plotFrameWithNormRange, ROOT.RooFit.Components(models[tag]["model2"].GetName()), ROOT.RooFit.Range("full"), ROOT.RooFit.NormRange("full"), ROOT.RooFit.LineColor(ROOT.kBlue)) ;
        name1 = models[tag]["model1"].GetName()+"_Norm[mh]_Comp["+models[tag]["model1"].GetName()+"]_Range[full]_NormRange[full]"
        name2 = models[tag]["model2"].GetName()+"_Norm[mh]_Comp["+models[tag]["model2"].GetName()+"]_Range[full]_NormRange[full]"
        chi2_1 = plotFrameWithNormRange.chiSquare(name1,"h_"+data.GetName(),fitresults.floatParsFinal().getSize())
        if doMultiPdf: chi2_2 = plotFrameWithNormRange.chiSquare(name2,"h_"+data.GetName(),fitresults2.floatParsFinal().getSize())
#        plotFrameWithNormRange.Print("v")
        print('--------------------')
        if doMultiPdf: print(models[tag]["model2"].GetName(),"    chi2/ndof=",round(chi2_2,2)," ndof",fitresults2.floatParsFinal().getSize())
        print(models[tag]["model1"].GetName(),"    chi2/ndof=",round(chi2_1,2)," ndof",fitresults.floatParsFinal().getSize())
        print('--------------------')

        fileToWrite="preselection_"+tag+"_"+"_"+str(year)+".txt"
        with open(fileToWrite, "a") as f:
            str1 = models[tag]["model1"].GetName()+"    chi2/ndof="+str(round(chi2_1,2))+" ndof"+str(fitresults.floatParsFinal().getSize())+"\n"
            f.write(str1)
            if doMultiPdf: str2 = models[tag]["model2"].GetName()+"    chi2/ndof="+str(round(chi2_2,2))+" ndof"+str(fitresults2.floatParsFinal().getSize())+"\n"
            if doMultiPdf: f.write(str2)

#    model.paramOn(plotFrameWithNormRange, RooFit.Layout(0.6,0.99,0.95))
#    plotFrameWithNormRange.getAttText().SetTextSize(0.02);

    plotFrameWithNormRange.Draw()
#    hresid = plotFrameWithNormRange.residHist()


    offsetY = 0.75*data_full.GetMaximum()
    latex = ROOT.TLatex()
    latex.SetTextColor(ROOT.kRed)
    latex.SetTextSize(0.04)
    latex.DrawLatex(130 ,offsetY + 0.10*data_full.GetMaximum(), models[tag]["model1"].GetName())
    latex.SetTextColor(ROOT.kBlue)
    latex.DrawLatex(130 ,offsetY + 0.20*data_full.GetMaximum(), models[tag]["model2"].GetName())

    canvas.Draw()
    htmldir = "~/public_html/HMUMU_FITS/OCT"
    canvas.SaveAs(htmldir+"/bkg_"+tag+"_"+str(year)+".png")

    # -----------------------------------------------------------------------------
    # -----------------------------------------------------------------------------
    # Create workspace, import data and model
    # Save workspace in file

    w = finalWorkspace(w, data, models[tag]["model1"], BKG_norm, doMultiPdf, tag, storedPdfs)

    if not blinded:
        w.writeToFile(workspaceName+"/Bkg"+tag+"_"+str(year)+"_workspace.root")

if __name__ == "__main__":

    for year in ['2024']:
        fitBkg('_ggHcat',year)
        fitBkg('_VBFcat',year)
        fitBkg('_TTLcat',year)
        fitBkg('_TTHcat',year)
        fitBkg('_Zinvcat',year)
        fitBkg('_VLcat',year)
        fitBkg('_VHcat',year)

    for year in ['12022','22022','12023','22023','2024']:
        fitSig('_ggHcat',year)
        fitSig('_VBFcat',year)
        fitSig('_TTLcat',year)
        fitSig('_TTHcat',year)
        fitSig('_Zinvcat',year)
        fitSig('_VLcat',year)
        fitSig('_VHcat',year)
