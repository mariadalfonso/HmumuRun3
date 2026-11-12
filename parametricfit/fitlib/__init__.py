"""Shared code of the parametric-fit chain.

  config          fit window, binning, modes per category, lumi, git commit
  style           CMS plot style and labels; sets ROOT up for batch plots
  signal          the signal fit (fit_one) and its plots (plot_fit)
  pdfDefinitions  signal and background PDFs, the candidate list
  dyCore          the Drell-Yan core families, dyexp and dypow
  bwz             the Breit-Wigner Z families of the Run 2 note (eqs. 5.1-5.4)
  fewz            FEWZ spline x Bernstein (Run 2 note eq. 5.5); fewz_spline.root
  agnostic        exp/pow + constant, Landau x Bernstein (Run 2 note eqs. 5.9-5.11)
  corePdf         the Run 2 core-PDF (note Sec. 6.4.2): a core shared by the
                  bins x a Bernstein shape modifier per bin
  prepareFits     input histograms, eras and luminosities
  makeWorkspaces  step of buildWorkspace.py: rebuild the models for Combine
  bwsHrare        step of buildWorkspace.py: write the datacard

Import from parametricfit/ (from fitlib import pdfDefinitions); run the two
steps as modules, e.g. python -m fitlib.makeWorkspaces, as
buildWorkspace.py does.
"""
