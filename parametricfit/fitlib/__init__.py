"""Shared code of the parametric-fit chain.

  pdfDefinitions  signal and background PDFs, the candidate list
  dyCore          the Drell-Yan core families, dyexp and dypow
  prepareFits     input histograms, eras and luminosities
  makeWorkspaces  step of buildWorkspace.py: rebuild the models for Combine
  bwsHrare        step of buildWorkspace.py: write the datacard

Import from parametricfit/ (from fitlib import pdfDefinitions); run the two
steps as modules, e.g. python -m fitlib.makeWorkspaces, as
buildWorkspace.py does.
"""
