import ROOT

from prepareFits import getHisto
from pdfDefinitions import *
from fitUtils import *
from SIGfits import bkg_model_selection, setVar
import time

ROOT.gROOT.SetBatch()
combine_base = "/code/HiggsAnalysis/CombinedLimit"
ROOT.gSystem.Load("/code/HiggsAnalysis/CombinedLimit/build/lib/libHiggsAnalysisCombinedLimit.so")
# NB before loading the singularity, do 'conda deactivate'

do_bias_study = True

## for GF and VBF
xlowRange = 110
xhighRange = 150

scaleSig = 10

htmldir = "/home/submit/mariadlf/public_html/HMUMU_FITS/MAY/BIAS"

cmdList = ROOT.RooLinkedList()
cmdList.Add(ROOT.RooFit.Minimizer("Minuit"))   # to speed up
cmdList.Add(ROOT.RooFit.Strategy(1))           # to speed up
#cmdList.Add(ROOT.RooFit.Minimizer("Minuit2")) 
#cmdList.Add(ROOT.RooFit.Strategy(2))
cmdList.Add(ROOT.RooFit.Range("full"))
cmdList.Add(ROOT.RooFit.Save(True))
cmdList.Add(ROOT.RooFit.PrintLevel(-1))
cmdList.Add(ROOT.RooFit.Optimize(2))           # add to speed up
cmdList.Add(ROOT.RooFit.SumW2Error(False))     # add to speed up
cmdList.Add(ROOT.RooFit.Hesse(False))          # add to speed up
cmdList.Add(ROOT.RooFit.Timer(True))

def visualize_single_toy(tag_, binMVA, year, ntoys=1):
    """
    Generate a single toy dataset and visualize the fit
    This helps debug the bias study setup
    """
    if not do_bias_study:
        print("Bias study disabled")
        return

    print(f"\n{'='*70}")
    print(f"VISUALIZING SINGLE TOY: {tag_}  {binMVA}  {year}")
    print(f"{'='*70}\n")

    if tag_ not in bkg_model_selection:
        print(f"ERROR: {tag_} not in bkg_model_selection")
        return

    model1_name = bkg_model_selection[tag_]["model1"]
    model2_name = bkg_model_selection[tag_]["model2"]

    print(f"Generation model: Signal + {model1_name}")
    print(f"Fit model: Signal + {model2_name}\n")

    tag = tag_+"_"+binMVA+"_"+year
    x = setVar(tag_, xlowRange, xhighRange)

    pdfDef = PDFDefinitions()

    # ===== LOAD SIGNAL MC =====
    print("Loading signal MC data...")
    # Get signal process, default to "ggH" if category not in map
    sig_process_map = {
        "ggHcat": "ggH",
        "VBFcat": "qqH",
        "VHcat": "VH",
        "VLcat": "VH",
        "Zinvcat": "qqH",
        "TTHcat": "ttH",
        "TTLcat": "ttH",
    }

    sig_process = sig_process_map.get(tag_, "ggH")
    data_full_sig = getHisto(10*int(xhighRange - xlowRange), xlowRange, xhighRange, 
                             False, tag_, year, True, binMVA, sig_process)

    binLow = data_full_sig.GetBin(1)
    binUp = data_full_sig.GetNbinsX()
    N_sig = scaleSig*data_full_sig.Integral(binLow, binUp)

    print(f"Signal ({sig_process}): {N_sig:.0f} events")

    if N_sig <= 0:
        print(f"ERROR: Signal MC is empty!")
        return

    data_sig = ROOT.RooDataHist(f'datahist_sig_{tag}', 'signal_data', 
                                ROOT.RooArgList(x), data_full_sig)

    # Fit signal
    print(f"Fitting signal PDF...")
    sig_pdfs = pdfDef.create_signal_pdfs(x, f'_{tag}')
    sig_pdf = sig_pdfs['pdf']

    try:
        sig_pdf.fitTo(data_sig, cmdList)
        print(f"✓ Signal fitted\n")
    except Exception as e:
        print(f"ERROR: {e}")
        return

    # ===== LOAD BACKGROUND =====
    print("Loading background data...")
    nBins = 10*int(xhighRange-xlowRange)
    data_full_bkg = getHisto(nBins, xlowRange, xhighRange, False, tag_, year, False, binMVA)    
    data_bkg = ROOT.RooDataHist(f'datahist_bkg_{tag}', 'bkg_data', 
                                ROOT.RooArgList(x), data_full_bkg)

    N_bkg = data_full_bkg.Integral()
    print(f"Background: {N_bkg:.0f} events\n")

    if N_bkg <= 0:
        print(f"ERROR: Background data is empty!")
        return

    # Create background PDFs
    bkg_pdfs_result = pdfDef.create_bkg_pdfs(x, tag)
    all_bkg_pdfs = bkg_pdfs_result['pdfs']

    gen_pdf_bkg = all_bkg_pdfs[model1_name]
    fit_pdf_bkg = all_bkg_pdfs[model2_name]

    # Fit background
    print(f"Fitting {gen_pdf_bkg.GetName()}...")
    try:
        gen_pdf_bkg.fitTo(data_bkg, cmdList)
        print(f"✓ {gen_pdf_bkg.GetName()} fitted")
    except Exception as e:
        print(f"ERROR: {e}")
        return

    print(f"Fitting {fit_pdf_bkg.GetName()}...")
    try:
        fit_pdf_bkg.fitTo(data_bkg, cmdList)
        print(f"✓ {fit_pdf_bkg.GetName()} fitted\n")
    except Exception as e:
        print(f"ERROR: {e}")
        return

    # ===== SETUP MODELS =====
    print("Setting up generation and fit models...")

    # Set signal parameters constant
    for param in sig_pdfs['params'].values():
        param.setConstant()

    ntot = int(N_sig + N_bkg)
    print(f"Total events: {ntot}\n")

    # Normalization variables
    Nsig = ROOT.RooRealVar("Nsig", "signal_yield", N_sig, 0, N_sig*3)
    Nbkg = ROOT.RooRealVar("Nbkg", "bkg_yield", N_bkg)

    # Generation PDF (FIXED)
    Nsig.setConstant(True)
    Nbkg.setConstant(True)

    genPDF = ROOT.RooAddPdf("totPDF_gen", "gen",
                           ROOT.RooArgList(sig_pdf, gen_pdf_bkg),
                           ROOT.RooArgList(Nsig, Nbkg))

    # Fit PDF (signal FLOATS)
    Nsig.setConstant(False)
    Nbkg.setConstant(True)

    fitPDF = ROOT.RooAddPdf("totPDF_fit", "fit",
                           ROOT.RooArgList(sig_pdf, fit_pdf_bkg),
                           ROOT.RooArgList(Nsig, Nbkg))

    # ===== GENERATE SINGLE TOY =====
    print("Generating single toy dataset...")
    print(f"  Drawing {ntot} events from generation model")

    try:
        # Generate toy data
        start_time = time.time()
        toy_data = genPDF.generate(ROOT.RooArgSet(x), ntot)
        elapsed = time.time() - start_time
        print(f"✓ Toy generated with {toy_data.sumEntries()} events")
        print(f"  Time to generate toy: {elapsed:.2f}s\n")
    except Exception as e:
        print(f"ERROR generating toy: {e}")
        return

    # ===== FIT TOY WITH FIT PDF =====
    print("Fitting toy with fit model...")
    try:
        start_time = time.time()
        fitresult = fitPDF.fitTo(toy_data, cmdList)
        elapsed = time.time() - start_time
        print(f"✓ Fit completed")
        print(f"  Status: {fitresult.status()}")
        print(f"  Nsig measured: {Nsig.getVal():.2f} ± {Nsig.getError():.2f} (true: {N_sig:.2f})")
        print(f"  Nbkg measured: {Nbkg.getVal():.2f} ± {Nbkg.getError():.2f} (true: {N_bkg:.2f})")
        print(f"  Time to fit: {elapsed:.2f}s\n")
    except Exception as e:
        print(f"ERROR fitting toy: {e}")
        return

    # ===== CREATE PLOT =====
    print("Creating visualization...")

    canvas = ROOT.TCanvas("toy_viz", "Single Toy Visualization", 1200, 600)
    canvas.Divide(2, 1)

    # Pad 1: Full plot with data and fits
    canvas.cd(1)
    frame = x.frame(ROOT.RooFit.Title(f"{tag_} - Single Toy Visualization"))

    # Plot data
    toy_data.plotOn(frame, ROOT.RooFit.Name("data"), ROOT.RooFit.MarkerStyle(20))

    # Plot generation model (dashed red)
    genPDF.plotOn(frame, 
                 ROOT.RooFit.Name("gen_model"),
                 ROOT.RooFit.LineColor(ROOT.kRed), 
                 ROOT.RooFit.LineStyle(2),
                 ROOT.RooFit.LineWidth(2),
                 ROOT.RooFit.Range("full"),
                 ROOT.RooFit.NormRange("full"))

    # Plot fit model (solid blue)
    fitPDF.plotOn(frame, 
                 ROOT.RooFit.Name("fit_model"),
                 ROOT.RooFit.LineColor(ROOT.kBlue), 
                 ROOT.RooFit.LineWidth(2),
                 ROOT.RooFit.Range("full"),
                 ROOT.RooFit.NormRange("full"))

    # Plot signal component in fit (green dashed)
    fitPDF.plotOn(frame,
                 ROOT.RooFit.Components(sig_pdf.GetName()),
                 ROOT.RooFit.Name("fit_signal"),
                 ROOT.RooFit.LineColor(ROOT.kGreen),
                 ROOT.RooFit.LineStyle(3),
                 ROOT.RooFit.LineWidth(2),
                 ROOT.RooFit.Range("full"),
                 ROOT.RooFit.NormRange("full"))

    # Plot background component in fit (purple dashed)
    fitPDF.plotOn(frame,
                 ROOT.RooFit.Components(fit_pdf_bkg.GetName()),
                 ROOT.RooFit.Name("fit_bkg"),
                 ROOT.RooFit.LineColor(ROOT.kMagenta),
                 ROOT.RooFit.LineStyle(3),
                 ROOT.RooFit.LineWidth(2),
                 ROOT.RooFit.Range("full"),
                 ROOT.RooFit.NormRange("full"))

    frame.Draw()

    # Add legend
    legend = ROOT.TLegend(0.6, 0.6, 0.95, 0.95)
    legend.AddEntry("data", "Generated Data", "EP")
    legend.AddEntry("gen_model", f"Gen: {model1_name} (RED, dashed)", "L")
    legend.AddEntry("fit_model", f"Fit: {model2_name} (BLUE, solid)", "L")
    legend.AddEntry("fit_signal", "Fit Signal Component (GREEN)", "L")
    legend.AddEntry("fit_bkg", "Fit Bkg Component (MAGENTA)", "L")
    legend.Draw()

    # Pad 2: Residuals
    canvas.cd(2)

    # Create residuals
    residuals_frame = x.frame(ROOT.RooFit.Title("Residuals (Data - Fit) / Error"))
    toy_data.plotOn(residuals_frame)
    fitPDF.plotOn(residuals_frame, 
                 ROOT.RooFit.LineColor(ROOT.kBlue),
                 ROOT.RooFit.Range("full"),
                 ROOT.RooFit.NormRange("full"))

    # Extract and compute residuals
    hresid = residuals_frame.residHist()
    hresid.SetName("residuals")

    residuals_frame.addPlotable(hresid, "P")
    residuals_frame.Draw()

    # Add zero line
    line = ROOT.TLine(xlowRange, 0, xhighRange, 0)
    line.SetLineColor(ROOT.kRed)
    line.SetLineStyle(2)
    line.Draw()

    canvas.Draw()

    # Save
    canvas.SaveAs(f"{htmldir}/toy_viz_{tag}_{model1_name}_vs_{model2_name}_scaleSig{scaleSig}.png")

    print(f"✓ Visualization saved\n")

    # ===== PRINT SUMMARY =====
    print(f"{'='*70}")
    print("SINGLE TOY SUMMARY")
    print(f"{'='*70}")
    print(f"Generation Model: Signal({N_sig:.0f}) + {model1_name}({N_bkg:.0f})")
    print(f"Fit Model: Signal + {model2_name}")
    print(f"\nFit Results:")
    print(f"  Signal: {Nsig.getVal():.2f} ± {Nsig.getError():.2f} (true: {N_sig:.0f})")
    print(f"  Background: {Nbkg.getVal():.2f} ± {Nbkg.getError():.2f} (true: {N_bkg:.0f})")
    print(f"\nBias Estimates:")
    sig_bias = (Nsig.getVal() - N_sig) / (Nsig.getError() if Nsig.getError() > 0 else 1)
    print(f"  Signal bias (pull): {sig_bias:.2f}σ")
    print(f"  Signal bias (absolute): {Nsig.getVal() - N_sig:.2f} events")
    print(f"{'='*70}\n")

def bias_study(tag_, binMVA, year, ntoys=100):

    """
    Perform bias study using toy MC generation
    
    FIXED: Use same parameter objects in gen and fit to avoid mismatches
    """
    if not do_bias_study:
        print("Bias study disabled (do_bias_study=False)")
        return

    print(f"\n{'='*70}")
    print(f"BIAS STUDY: {tag_}  {binMVA}  {year}")
    print(f"{'='*70}\n")

    if tag_ not in bkg_model_selection:
        print(f"ERROR: {tag_} not in bkg_model_selection")
        return

    model1_name = bkg_model_selection[tag_]["model1"]
    model2_name = bkg_model_selection[tag_]["model2"]

    print(f"Generation: Signal (FIXED) + {model1_name} (FIXED)")
    print(f"Fit: Signal (FLOATING) + {model2_name}\n")

    tag = tag_+"_"+binMVA+"_"+year
    x = setVar(tag_, xlowRange, xhighRange)

    pdfDef = PDFDefinitions()

    # ===== LOAD AND FIT SIGNAL MC =====

    print("Loading signal MC data...")

    sig_process_map = {
        "ggHcat": "ggH",
        "VBFcat": "qqH",
        "VHcat": "VH",
        "VLcat": "VH",
        "Zinvcat": "qqH",
        "TTHcat": "ttH",
        "TTLcat": "ttH",
    }

    sig_process = sig_process_map.get(tag_, "ggH")

    data_full_sig = getHisto(10*int(xhighRange - xlowRange), xlowRange, xhighRange, 
                             False, tag_, year, True, binMVA, sig_process)

    binLow = data_full_sig.GetBin(1)
    binUp = data_full_sig.GetNbinsX()
    N_sig = scaleSig*data_full_sig.Integral(binLow, binUp)

    print(f"Signal ({sig_process}) MC integral: {N_sig:.0f} events")

    if N_sig <= 0:
        print(f"ERROR: Signal MC is empty!")
        return

    data_sig = ROOT.RooDataHist(f'datahist_sig_{tag}', 'signal_data', 
                                ROOT.RooArgList(x), data_full_sig)

    # Create SINGLE set of signal PDFs (shared by gen and fit)
    print(f"Fitting signal PDF...")
    sig_pdfs = pdfDef.create_signal_pdfs(x, f'_{tag}')  # Use same tag, no _gen/_fit
    sig_pdf = sig_pdfs['pdf']

    try:
        fitresult_sig = sig_pdf.fitTo(data_sig, cmdList)
        print(f"✓ Signal PDF fitted")
    except Exception as e:
        print(f"ERROR fitting signal: {e}")
        return

    # ===== LOAD AND FIT BACKGROUND =====
    print(f"\nLoading background data...")
    nBins = 10*int(xhighRange-xlowRange)
    data_full_bkg = getHisto(nBins, xlowRange, xhighRange, False, tag_, year, False, binMVA)    
    data_bkg = ROOT.RooDataHist(f'datahist_bkg_{tag}', 'bkg_data', 
                                ROOT.RooArgList(x), data_full_bkg)

    N_bkg = data_full_bkg.Integral()
    if N_bkg <= 0:
        print(f"ERROR: Background data is empty!")
        return

    print(f"Background data integral: {N_bkg:.0f} events")

    # Create background PDFs
    bkg_pdfs_result = pdfDef.create_bkg_pdfs(x, tag)
    all_bkg_pdfs = bkg_pdfs_result['pdfs']

    gen_pdf_bkg = all_bkg_pdfs[model1_name]
    fit_pdf_bkg = all_bkg_pdfs[model2_name]

    # Fit background with both models
    print(f"Fitting background with {gen_pdf_bkg.GetName()}...")
    try:
        gen_pdf_bkg.fitTo(data_bkg, cmdList)
        print(f"✓ {gen_pdf_bkg.GetName()} fitted")
    except Exception as e:
        print(f"ERROR fitting {gen_pdf_bkg.GetName()}: {e}")
        return

    print(f"Fitting background with {fit_pdf_bkg.GetName()}...")
    try:
        fit_pdf_bkg.fitTo(data_bkg, cmdList)
        print(f"✓ {fit_pdf_bkg.GetName()} fitted")
    except Exception as e:
        print(f"ERROR fitting {fit_pdf_bkg.GetName()}: {e}")
        return

    # ===== PREPARE TOY MC STUDY =====
    print(f"\nPreparing toy MC study...")

    # Set signal parameters constant (same for both gen and fit!)
    for param in sig_pdfs['params'].values():
        param.setConstant()

    print(f"Signal yield: {N_sig:.0f} events")
    print(f"Background yield: {N_bkg:.0f} events")

    ntot = int(N_sig + N_bkg)
    print(f"Total events per toy: {ntot}\n")

    # ===== NORMALIZATION VARIABLES =====
    Nsig = ROOT.RooRealVar("Nsig", "signal_yield", N_sig, 0, N_sig*3)
    Nbkg = ROOT.RooRealVar("Nbkg", "bkg_yield", N_bkg)

    # ===== CREATE GENERATION PDF =====
    # Signal and background FIXED
    Nsig.setConstant(False)  # Will be set constant after building PDF
    Nbkg.setConstant(True)   # Keep background fixed

    genPDF = ROOT.RooAddPdf("totPDF_gen", "gen",
                           ROOT.RooArgList(sig_pdf, gen_pdf_bkg),
                           ROOT.RooArgList(Nsig, Nbkg))

    # NOW set signal constant for generation
    Nsig.setConstant(True)

    # ===== CREATE FIT PDF =====
    # Use SAME Nsig and Nbkg variables but allow Nsig to float during fit
    Nsig.setConstant(False)  # FLOAT for fitting
    Nbkg.setConstant(True)   # Keep background fixed

    fitPDF = ROOT.RooAddPdf("totPDF_fit", "fit",
                           ROOT.RooArgList(sig_pdf, fit_pdf_bkg),
                           ROOT.RooArgList(Nsig, Nbkg))

    # Run toy MC study
    print(f"Starting RooMCStudy...")
    print(f"  Number of toys: {ntoys}")
    print(f"  Events per toy: {ntot}")
    print(f"  Gen: Signal={N_sig:.0f} (FIXED) + {model1_name}={N_bkg:.0f} (FIXED)")
    print(f"  Fit: Signal (FLOATING) + {model2_name}={N_bkg:.0f} (FIXED)")
    print()

    try:
        mcstudy = ROOT.RooMCStudy(
            genPDF,
            ROOT.RooArgSet(x),
            ROOT.RooFit.Silence(),
            ROOT.RooFit.FitModel(fitPDF),
            ROOT.RooFit.Extended(),
            ROOT.RooFit.NumCPU(4),            
            ROOT.RooFit.FitOptions(
                ROOT.RooFit.Save(),
                ROOT.RooFit.Minimizer("Minuit2"),
                ROOT.RooFit.Strategy(2),
                ROOT.RooFit.Range("full"),
                ROOT.RooFit.PrintEvalErrors(0)
            )
        )

        start_time = time.time()
        print(f"Generating and fitting {ntoys} toys...")
        mcstudy.generateAndFit(ntoys, ntot, True)
        print(f"✓ Toy MC study completed\n")
        elapsed = time.time() - start_time
        avg_time_per_toy = elapsed / ntoys

        print(f"-----------------------------------------------") 
        print(f"✓ Toy MC study completed")
        print(f"  Total time: {elapsed:.1f}s")
        print(f"  Average per toy: {avg_time_per_toy:.2f}s")
        print(f"  Per 1000 toys: {avg_time_per_toy*1000:.0f}s\n")
        print(f"-----------------------------------------------") 
        
    except Exception as e:
        print(f"ERROR in RooMCStudy: {e}")
        import traceback
        traceback.print_exc()
        return

    try:
        # Plot 1: Signal yield distribution
        canvas1 = ROOT.TCanvas("sig_yield", "sig_yield", 800, 600)
        sig_frame = mcstudy.plotParam(Nsig, ROOT.RooFit.Bins(50))
        sig_frame.SetTitle(f"Signal Yield: {tag_}\nGen: {model1_name}, Fit: {model2_name}")
        sig_frame.GetXaxis().SetTitle("Signal Yield (events)")
        sig_frame.Draw()

        # Add line at true value
        line = ROOT.TLine(N_sig, 0, N_sig, sig_frame.GetMaximum())
        line.SetLineColor(ROOT.kRed)
        line.SetLineWidth(2)
        line.SetLineStyle(2)
        line.Draw()

        canvas1.SaveAs(f"{htmldir}/sig_yield_{tag}_{binMVA}_{year}_{model1_name}_vs_{model2_name}_scaleSig{scaleSig}.png")
        print(f"✓ Signal yield plot saved")

        # Plot 2: Pull distribution with Gaussian fit
        print("Creating pull distribution with Gaussian fit...")

        canvas2 = ROOT.TCanvas("sig_pull", "sig_pull", 1000, 700)
        canvas2.SetLeftMargin(0.12)
        canvas2.SetRightMargin(0.05)

        pull_frame = mcstudy.plotPull(Nsig,
                                      ROOT.RooFit.Range(-5., 5.),
                                      ROOT.RooFit.Bins(50),
                                      ROOT.RooFit.FitGauss())  # ← Gaussian fit added

        pull_frame.SetTitle(f"Signal Pull Distribution: {tag_}\n")
        pull_frame.GetYaxis().SetTitle("Number of Toys")
        pull_frame.GetYaxis().SetTitleSize(0.05)
        pull_frame.GetXaxis().SetTitle("(N_{sig}^{fit} - N_{sig}^{true}) / #sigma")
        pull_frame.GetXaxis().SetTitleSize(0.05)

        pull_frame.Draw()

        # Add information box with statistics
        info_box = ROOT.TPaveText(0.15, 0.75, 0.35, 0.92, "NDC")
        info_box.SetBorderSize(0)
        info_box.SetFillStyle(0)            
        info_box.SetTextAlign(12)
        info_box.SetTextSize(0.04)
        info_box.SetTextFont(42)

        info_box.AddText(f"Category: {tag_}")
        info_box.AddText(f"Toys: {ntoys}")
        info_box.AddText("")
        info_box.AddText(f"Generated: {sig_process}+{model1_name}")
        info_box.AddText(f"Fitted: {model2_name}")
        info_box.Draw()

        canvas2.SaveAs(f"~/public_html/HMUMU_FITS/MAY/BIAS/sig_pull_{tag}_{model1_name}_vs_{model2_name}_scaleSig{scaleSig}.png")
        print(f"✓ Signal pull plot saved with Gaussian fit\n")

    except Exception as e:
        print(f"ERROR plotting results: {e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":

    for year in ['Run3']:

        for binX in ["bdt0", "bdt1","bdt2"]:

            # Optional: Run bias studies
            if do_bias_study:
                print(f"\n{'='*80}")
                print(f"BIAS STUDY CAMPAIGN")
                print(f"{'='*80}")

                # First: visualize to debug
                visualize_single_toy('VBFcat', binX, year)
                visualize_single_toy('ggHcat', binX, year)
                visualize_single_toy('VLcat', binX, year)
                visualize_single_toy('VHcat', binX, year)
                visualize_single_toy('Zinvcat', binX, year)
                visualize_single_toy('TTLcat', binX, year)
                visualize_single_toy('TTHcat', binX, year)

                test_configs = [
                    # (tag, binX, year, BR, ntoys)
                    ('VBFcat', binX, year, 1000),    # Quick test
                    ('ggHcat', binX, year, 1000),    # Quick test
                    ('VLcat', binX, year, 1000),    # Quick test
                    ('VHcat', binX, year, 1000),    # Quick test
                    ('Zinvcat', binX, year, 1000),    # Quick test                    
                    ('TTLcat', binX, year, 1000),    # Quick test
                    ('TTHcat', binX, year, 1000),    # Quick test                    
               ]

                for tag_, bin_mva, yr, n_toys in test_configs:
                    try:
                        print(f"\nRunning: {tag_} {bin_mva} ntoys={n_toys}")
                        bias_study(tag_, bin_mva, yr, ntoys=n_toys)

                    except Exception as e:
                        print(f"ERROR: {e}")
                        import traceback
                        traceback.print_exc()
                        continue
