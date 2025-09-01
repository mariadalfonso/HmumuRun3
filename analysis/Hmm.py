import ROOT
import json
import sys

from utilsAna import loadUserCode
from utilsAna import BuildDict, SwitchSample
from utilsAna import readDataQuality
from utilsAna import computeWeigths
from utilsAna import lumis,btagPNetBM,btagPNetBL
from utilsAna import loadCorrectionSet
from datetime import datetime

loadUserCode()

ROOT.gROOT.SetBatch()
ROOT.ROOT.EnableImplicitMT()

with open("/home/submit/mariadlf/Hrare/CMSSW_10_6_27_new/src/HmumuRun3/analysis/config/selection.json") as jsonFile:
    jsonObject = json.load(jsonFile)
    jsonFile.close()

selections = {key: jsonObject[key] for key in ["GOODMUON", "GOODJETSALL", "GOODelectrons", "LOOSEelectrons", "LOOSEmuons", "BJETS", "METFilters", "TRIGGER" ]}

JSON = "isGoodRunLS(isData, run, luminosityBlock)"

mode = sys.argv[2]  # expected: isVBF, isGGH, isVlep, isVhad, isTTlep, isTThad, isZinv

# --- place to write the files and strings  ---
myDir = '/work/submit/mariadlf/HmumuRun3/ROOTFILES/'

# --- Histogram output suffixes per mode ---
mode_map = {
    "isVBF":  "VBFcat",
    "isGGH":  "ggHcat",
    "isZinv": "Zinvcat",
    "isVlep": "VLcat",
    "isTTlep": "TTLcat",
}

def dfwithSYST(df,year):

    if year=='2024': dfFinal_withSF = df.Define("w_allSF", "w")

    else:

        dfFinal_withSF = (df
                          .Define("SFmuon1_ID_Nom",'corr_sf.eval_muonIDSF("{0}", "nominal", Muon1_eta, Muon1_pt, "M")'.format(year))
                          .Define("SFmuon1_ID_Up",'corr_sf.eval_muonIDSF("{0}", "systup", Muon1_eta, Muon1_pt, "M")'.format(year))
                          .Define("SFmuon1_ID_Dn",'corr_sf.eval_muonIDSF("{0}", "systdown", Muon1_eta, Muon1_pt, "M")'.format(year))
                          .Define("SFmuon2_ID_Nom",'corr_sf.eval_muonIDSF("{0}", "nominal", Muon2_eta, Muon2_pt, "M")'.format(year))
                          .Define("SFmuon2_ID_Up",'corr_sf.eval_muonIDSF("{0}", "systup", Muon2_eta, Muon2_pt, "M")'.format(year))
                          .Define("SFmuon2_ID_Dn",'corr_sf.eval_muonIDSF("{0}", "systdown", Muon2_eta, Muon2_pt, "M")'.format(year))
                          #
                          .Define("SFpu_Nom",'corr_sf.eval_puSF(Pileup_nTrueInt,"nominal")')
                          .Define("SFpu_Up",'corr_sf.eval_puSF(Pileup_nTrueInt,"up")')
                          .Define("SFpu_Dn",'corr_sf.eval_puSF(Pileup_nTrueInt,"down")')
                          #
                          #                      .Define("muoID_weights", "NomUpDownVar(SFmuon_ID_Nom, SFmuon_ID_Up, SFmuon_ID_Dn, w_allSF)")
                          #                      .Define("pu_weights", "NomUpDownVar(SFpu_Nom, SFpu_Up, SFpu_Dn, w_allSF)")
                          #                      .Define("idx_nom_up_down", "indices(3)")
                          .Define("w_allSF", "w*SFpu_Nom*SFmuon1_ID_Nom*SFmuon2_ID_Nom")
                          )

    return dfFinal_withSF



def doCategories(df,mc,year):

    #mode = sys.argv[2]  # expected: isVBF, isGGH, isW, isZinv

    # --- Common jet selections ---
    BJETSmedium = selections["BJETS"].format(btagPNetBM[year])
    BJETSloose  = selections["BJETS"].format(btagPNetBL[year])

    df = (df.Define("BJETS", f"{BJETSmedium}")
          .Define("BJETSloose", f"{BJETSloose}")
          )

    # --- Reusable blocks ---
    ele_veto = [
        ("FILTER", "Sum({})==0".format(selections["LOOSEelectrons"]), "no extra electrons"),
    ]

    bjet_veto = [
        ("FILTER", "Sum(BJETS)==0 or Sum(BJETSloose)<2", "events vetoed if at least 1bM or 2bL"),
    ]

    bjet = [
        ("FILTER", "Sum(BJETS)>0 or Sum(BJETSloose)>1", "events with at least 1bM or 2bL"),
    ]

    met_filters = [
        ("FILTER", "PuppiMET_pt>150", "MET > 150"),
        ("FILTER", selections["METFilters"], "apply MET filters"),
    ]

    ele = [
        ("DEFINE", "goodElectrons", selections["GOODelectrons"]),
        ("FILTER", "Sum(goodElectrons)>=1", "TT/W candidate with W to ele"),
    ]

    # --- Category recipes ---
    # PUT first the things that cut most of the events i.e. met in Zinv
    categories = {
        "isGGH":       ele_veto + bjet_veto, # TO add the met veto and VBFjet veto
        "isVBF":       ele_veto + bjet_veto, # TO add the met veto
        "isZinv":      ele_veto + bjet_veto + met_filters,
        "isWlep":      ele + bjet_veto,
#        "isWhad":      ele_vet + bjet_veto,
        "isTTlep":     ele + bjet,
        "isTTHhad":    ele_veto + bjet,
    }

   # --- Apply category ---
    if mode not in categories:
        raise ValueError(f"Unknown category mode: {mode}")

    for action, expr, *desc in categories[mode]:
        if action == "DEFINE":
            df = df.Define(expr, desc[0])   # expr = variable, desc[0] = RHS expression
        elif action == "FILTER":
            label = desc[0] if desc else ""
            df = df.Filter(expr, label)

    # ---------
    # ---  SPECIFIC THINGS for each MODE
    # ---------

    '''
    elif mode == "isGGH":

    elif mode == "isZinv":
        df= (df.Filter("PuppiMET_pt>150","MET > 150")
             .Define("METFilters","{}".format(selections["METFilters"]))
             .Filter("METFilters>0","Pass MET filters")
             )
    '''

    if mode == "isWlep":
        # 3 muons The charge of the three leptons must add up to ±1
        # At least one μ+μ− pair must have an invariant mass between 110 and 150 GeV
        # In 3μ events, the non-Higgs-candidate μ+μ− pair (μμOS) must not have an invariant mass between 81 and 101 GeV, to suppress WZ and Z+jets backgrounds
        df= (df.Define("Lepton_Pt","Electron_pt[goodElectrons][0]")
             .Define("looseEle","{}".format(selections["LOOSEelectrons"]))
             .Define("isMuorEle","(Sum(looseEle)==1 and Sum(goodElectrons)>=1)?1 : 0")
             )

    elif mode == "isTTlep":
        # 3 muons The charge of the three leptons must add up to ±1
        # At least one μ+μ− pair must have an invariant mass between 110 and 150 GeV
        # In 3μ events, the non-Higgs-candidate μ+μ− pair (μμOS) must not have an invariant mass between 81 and 101 GeV, to suppress WZ and Z+jets backgrounds
        df= (df.Define("Lepton_Pt","Electron_pt[goodElectrons][0]")
             .Define("looseEle","{}".format(selections["LOOSEelectrons"]))
             .Define("isMuorEle","(Sum(looseEle)==1 and Sum(goodElectrons)>=1)?1 : 0")
             )

    elif mode == "isttHhad":
        df= (df.Define("goodJetsAll","{}".format(selections["GOODJETSALL"]))
             .Define("nGoodJetsAll","Sum(goodJetsAll)*1.0f").Filter("Sum(goodJetsAll)>2","at least three jets")
             ## add leading jet pt>50 and triplet mass mjjj 100-300
             ## look the toponium and the ttHcc
             )

    elif mode == "isVBF":
        df= (df.Define("goodJetsAll","{}".format(selections["GOODJETSALL"]))
             .Define("jetAll_Pt","Jet_pt[goodJetsAll]")
             .Define("jetAll_Eta","Jet_eta[goodJetsAll]")
             .Define("jetAll_Phi","Jet_phi[goodJetsAll]")
             .Define("jetAll_Mass","Jet_mass[goodJetsAll]")
             #
             .Define("index_VBFJets","getVBFIndicies(jetAll_Pt, jetAll_Eta, jetAll_Phi, jetAll_Mass)")
             .Filter("index_VBFJets[0]!= -1", "both VBF jets")
             .Filter("index_VBFJets[1]!= -1", "both VBF jets")
             .Define("nGoodJetsAll","Sum(goodJetsAll)*1.0f").Filter("Sum(goodJetsAll)>1","two VBF jets")
             .Define("jetVBF1_Pt","Jet_pt[goodJetsAll][index_VBFJets[0]]")
             .Define("jetVBF1_Eta","Jet_eta[goodJetsAll][index_VBFJets[0]]")
             .Define("jetVBF1_Phi","Jet_phi[goodJetsAll][index_VBFJets[0]]")
             .Define("jetVBF1_Mass","Jet_mass[goodJetsAll][index_VBFJets[0]]")
             .Define("jetVBF2_Pt","Jet_pt[goodJetsAll][index_VBFJets[1]]")
             .Define("jetVBF2_Eta","Jet_eta[goodJetsAll][index_VBFJets[1]]")
             .Define("jetVBF2_Phi","Jet_phi[goodJetsAll][index_VBFJets[1]]")
             .Define("jetVBF2_Mass","Jet_mass[goodJetsAll][index_VBFJets[1]]")
             .Define("VBF1Vec", "MakeTLV(Jet_pt[goodJetsAll][index_VBFJets[0]],Jet_eta[goodJetsAll][index_VBFJets[0]], Jet_phi[goodJetsAll][index_VBFJets[0]],Jet_mass[goodJetsAll][index_VBFJets[0]])")
             .Define("VBF2Vec", "MakeTLV(Jet_pt[goodJetsAll][index_VBFJets[1]],Jet_eta[goodJetsAll][index_VBFJets[1]], Jet_phi[goodJetsAll][index_VBFJets[1]],Jet_mass[goodJetsAll][index_VBFJets[1]])")
             #
             .Define("jetVBF1_hfcentralEtaStripSize","Jet_hfcentralEtaStripSize[goodJetsAll][index_VBFJets[0]]")
             .Define("jetVBF1_hfadjacentEtaStripsSize","Jet_hfadjacentEtaStripsSize[goodJetsAll][index_VBFJets[0]]")
             .Define("jetVBF2_hfcentralEtaStripSize","Jet_hfcentralEtaStripSize[goodJetsAll][index_VBFJets[1]]")
             .Define("jetVBF2_hfadjacentEtaStripsSize","Jet_hfadjacentEtaStripsSize[goodJetsAll][index_VBFJets[1]]")
             #
             .Define("RPt","getRpt(Muon1Vec, Muon2Vec, VBF1Vec, VBF2Vec)")
             .Define("Mjj","Pair12Minv(jetVBF1_Pt, jetVBF1_Eta, jetVBF1_Phi, jetVBF1_Mass, jetVBF2_Pt, jetVBF2_Eta, jetVBF2_Phi, jetVBF2_Mass)")
             .Define("dEtaJJ","abs(jetVBF1_Eta-jetVBF2_Eta)")
             .Define("ZepVar","getZep(Muon1Vec, Muon2Vec, jetVBF1_Eta, jetVBF2_Eta)")
             .Define("minDetaDiMuVBF","minDeta(MinvCorr(Muon_bsConstrainedPt[goodMuons], Muon_eta[goodMuons], Muon_phi[goodMuons], Muon_mass[goodMuons], FsrPhoton_pt[goodFRSphoton], FsrPhoton_eta[goodFRSphoton],FsrPhoton_phi[goodFRSphoton],3), jetVBF1_Eta, jetVBF2_Eta)")
              )

        if mc>0:
            df= (df.Define("jetVBF1_partonFlavour","Jet_partonFlavour[goodJetsAll][index_VBFJets[0]]")
                 .Define("jetVBF2_partonFlavour","Jet_partonFlavour[goodJetsAll][index_VBFJets[1]]")
                 )

    return df


def analysis(files,year,mc,sumW):
    dfINI = ROOT.RDataFrame("Events", files)

    isData = "false"
    if mc<0: isData = "true"

    lumiIntegrated=0.
    if (isData == "false"):
        lumiIntegrated = lumis[year]
        print('lumiIntegrated = ',lumiIntegrated)

    weight = "{0}".format(1.)
    if mc>=0: weight = "{0}*genWeight*{1}".format(lumiIntegrated,sumW)

    rho = 'Rho_fixedGridRhoFastjetAll'
    if year=='12016' or year=='22016' or year=='2017' or  year=='2018':
        rho = 'fixedGridRhoFastjetAll'

    dfComm = (dfINI
              .Define("mc","{}".format(mc))
              .Define("isData","{}".format(isData))
              .Define("applyJson","{}".format(JSON)).Filter("applyJson","pass JSON")
              .Define("w","{}".format(weight))
              .Define("wraw","{}".format(weight))
              .Define("rho","{}".format(rho))
              .Define("lumiIntegrated","{}".format(lumiIntegrated))
              .Filter("PV_npvsGood>0","one good PV")
              .Filter("{}".format(selections["TRIGGER"]),"passing trigger")
              )


    # apply JEC
    dfComm = dfComm.Redefine("Jet_pt",'computeJECcorrection(corr_sf, Jet_pt, Jet_rawFactor, Jet_eta, Jet_phi, Jet_area, rho, run, isData, "{0}","{1}" )'.format(year,mc))

    # compute JetID (note v15 vs v12)
    if year=="2024":
        dfComm = dfComm.Define("jetID_mask", 'cleaningJetIDMask(Jet_eta, Jet_chHEF, Jet_neHEF, Jet_chEmEF, Jet_neEmEF, Jet_muEF, Jet_chMultiplicity, Jet_neMultiplicity, "{0}")'.format(year)) #  for nanov15 use the JetID
    else:
        dfComm = dfComm.Define("jetID_mask", "cleaningJetSelMask(0, Jet_eta, Jet_neHEF, Jet_chEmEF, Jet_muEF, Jet_neEmEF, Jet_jetId)") # redo JetID for nanoV12,

    # Zgamma has these: |dz| < 1.0 cm, |dxy| < 0.5 cm, SIP3D < 4
    df = (dfComm.Define("goodMuons","{}".format(selections["GOODMUON"]))
          .Filter("Sum(goodMuons)>=1 and Sum(Muon_charge[goodMuons])==0","at least two good muons OS") # fix the charge for WH and ZH          
          .Define("looseMu","{}".format(selections["LOOSEmuons"]))
          .Filter("Sum(looseMu)==2", "no extra muons")
          .Define("jetMuon_mask", "cleaningMask(Muon_jetIdx[goodMuons],nJet)")
          # Jet_puIdDisc only for nanov15
          .Define("jetVeto_mask",'cleaningJetVetoMapMask(Jet_eta,Jet_phi,"{0}")'.format(year))
          #
#          .Define("Muon1_pt","Muon_pt[goodMuons][0]")
#          .Define("Muon2_pt","Muon_pt[goodMuons][1]")
          .Define("Muon1_pt","Muon_bsConstrainedPt[goodMuons][0]")
          .Define("Muon2_pt","Muon_bsConstrainedPt[goodMuons][1]")
          .Define("Muon1_eta","Muon_eta[goodMuons][0]")
          .Define("Muon2_eta","Muon_eta[goodMuons][1]")
          .Define("Muon1Vec", "MakeTLV(Muon_bsConstrainedPt[goodMuons][0], Muon_eta[goodMuons][0], Muon_phi[goodMuons][0],Muon_mass[goodMuons][0])")
          .Define("Muon2Vec", "MakeTLV(Muon_bsConstrainedPt[goodMuons][1], Muon_eta[goodMuons][1], Muon_phi[goodMuons][1],Muon_mass[goodMuons][1])")
          ## end preselection
#          .Define("HiggsCandMass","Minv(Muon_pt[goodMuons], Muon_eta[goodMuons], Muon_phi[goodMuons], Muon_mass[goodMuons])")
          .Define("HiggsCandMass","Minv(Muon_bsConstrainedPt[goodMuons], Muon_eta[goodMuons], Muon_phi[goodMuons], Muon_mass[goodMuons])")
          .Define("goodFRSphoton", "fsrMask(Muon_fsrPhotonIdx[goodMuons],nFsrPhoton)")
          .Define("FsrPH_pt","FsrPhoton_pt[goodFRSphoton]")
          .Define("FsrPH_eta","FsrPhoton_eta[goodFRSphoton]")
          .Define("FsrPH_pt_ratio0","FsrPhoton_pt[goodFRSphoton][0]/Muon1_pt")
          .Define("FsrPH_pt_ratio1","FsrPhoton_pt[goodFRSphoton][1]/Muon2_pt")
          )

    df = doCategories(df,mc,year)

    if (isData == "false"):
        df = dfwithSYST(df,year)
    else:
        df = df.Define("w_allSF", "w")

    df = (df.Define("HiggsCandCorrMass","MinvCorr(Muon1Vec, Muon2Vec, FsrPhoton_pt[goodFRSphoton], FsrPhoton_eta[goodFRSphoton],FsrPhoton_phi[goodFRSphoton],0)")
          .Define("HiggsCandCorrPt","MinvCorr(Muon1Vec, Muon2Vec, FsrPhoton_pt[goodFRSphoton], FsrPhoton_eta[goodFRSphoton],FsrPhoton_phi[goodFRSphoton],1)")
          .Define("HiggsCandCorrRapidity","MinvCorr(Muon1Vec, Muon2Vec, FsrPhoton_pt[goodFRSphoton], FsrPhoton_eta[goodFRSphoton],FsrPhoton_phi[goodFRSphoton],2)")
          .Filter("HiggsCandCorrMass>50 and HiggsCandCorrMass<200","HiggsMass within reasonable range 50-200")
          )

    if False:
        print("---------------- SUMMARY -------------")
        ## this doens't work with the negative weights
        report = df.Report()
        report.Print()

    if False:
        print("writing histograms")
        hists = {
            "HiggsCandMass": {"name":"HiggsCandMass","title":"M(mu1,mu2); M(mu1,mu2) (GeV);N_{Events}","bin":250,"xmin":50.,"xmax":300.},
            "HiggsCandCorrMass": {"name":"HiggsCandCorrMass","title":"M(mu1,mu2,fsr); M(mu1,mu2,fsr) (GeV);N_{Events}","bin":250,"xmin":50.,"xmax":300.},
            "HiggsCandCorrPt": {"name":"HiggsCandCorrPt","title":"p_{T}(mu1,mu2,fsr); p_{T}(mu1,mu2,fsr) (GeV);N_{Events}","bin":250,"xmin":0.,"xmax":250.},
            "FsrPH_pt": {"name":"FsrPH_pt","title":"FSR0; FSR (GeV);N_{Events}","bin":200,"xmin":0.,"xmax":20.},
            "FsrPH_pt_ratio0": {"name":"FsrPH_pt_ratio0","title":"FSR_ratio0; FSR/p_{T}(mu1) (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":10.},
#            "PuppiMET_pt": {"name":"PuppiMET_pt","title":"PuppiMET_pt; PuppiMET (GeV);N_{Events}","bin":300,"xmin":0.,"xmax":300.},
        }

        histsVBF = {
            "Mjj": {"name":"Mjj","title":"M(VBFjet,VBFjet); M(VBFjet,VBFjet) (GeV);N_{Events}","bin":250,"xmin":0.,"xmax":1000},
            "RPt": {"name":"RPt","title":"RPt; RPt;N_{Events}","bin":100,"xmin":0.,"xmax":1},
            "jetVBF1_Pt":  {"name":"jetVBF1_Pt","title":"jetVBF1_Pt; p_{T}^{VBFjet} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":500.},
            "jetVBF2_Pt":  {"name":"jetVBF2_Pt","title":"jetVBF2_Pt; p_{T}^{VBFjet} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":500.},
            "ZepVar":  {"name":"ZepVar","title":"ZepVar; ZepVar;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
            ##
        }

        if mode == "isVBF": hists.update(histsVBF)

        histos = []
        for h in hists:
            model1d = (hists[h]["name"]+"_"+str(year), hists[h]["title"], hists[h]["bin"], hists[h]["xmin"], hists[h]["xmax"])
#            h1d = df.Histo1D(model1d, hists[h]["name"],"w")
            h1d = df.Histo1D(model1d, hists[h]["name"],"w_allSF")
            histos.append(h1d)

            for h in histos:
                canv = ROOT.TCanvas("stackcanvas","Stack canvas",800,800)
                canv.SetLogy()
                h.Draw("hist")
                canv.Draw()
                dirPlots="~/public_html/HMUMU_July/SignalPlots"
                if mc==10: canv.SaveAs(dirPlots+h.GetName()+"_VBF.png")
                if mc==11: canv.SaveAs(dirPlots+h.GetName()+"_ggH.png")
                if mc==12: canv.SaveAs(dirPlots+h.GetName()+"_Wm.png")
                if mc==13: canv.SaveAs(dirPlots+h.GetName()+"_Wp.png")
                if mc==14: canv.SaveAs(dirPlots+h.GetName()+"_ZH.png")
                if mc==15: canv.SaveAs(dirPlots+h.GetName()+"_TTH.png")
                if mc==100: canv.SaveAs(dirPlots+h.GetName()+"_DY.png")
                if mc==20: canv.SaveAs(dirPlots+h.GetName()+"_Zg_VBF.png")
                if mc==21: canv.SaveAs(dirPlots+h.GetName()+"_Zg_ggH.png")

        if mode not in mode_map:
            raise ValueError(f"Unknown mode: {mode}")

        outputFileHisto = f"{myDir}histoOUTname_{mc}_{year}_{mode_map[mode]}.root"
        print(outputFileHisto )

        myfile = ROOT.TFile(outputFileHisto,"RECREATE")

        myfile.ls()
        for h in histos:
            print(h.GetName())
            h.Write()

        myfile.Close()
        myfile.Write()

    if True:
        print("writing snapshot")

        branchList = ROOT.vector('string')()

        base_branches = [
                "mc",
                "w",
                "w_allSF",
                "lumiIntegrated",
                "PV_npvsGood",
                "run",
                "event",
                "luminosityBlock",
                #
                "HiggsCandCorrMass",
                "HiggsCandCorrPt",
                "HiggsCandCorrRapidity",
                "Muon1_pt",
                "Muon2_pt",
                "FsrPH_pt",
                "FsrPH_eta",
                "PuppiMET_pt",
                "PuppiMET_phi"
        ]

        # Mode-specific branches
        mode_branches = {
            "isVBF": [
                "Mjj",
                "RPt",
                "jetVBF2_Pt",
                "jetVBF1_Pt",
                "jetVBF2_Eta",
                "jetVBF1_Eta",
                "jetVBF2_Phi",
                "jetVBF1_Phi",
                "dEtaJJ",
                "ZepVar",
                "minDetaDiMuVBF",
                "jetVBF2_hfcentralEtaStripSize",
                "jetVBF1_hfcentralEtaStripSize",
                "jetVBF2_hfadjacentEtaStripsSize",
                "jetVBF1_hfadjacentEtaStripsSize",
            ],
            "isVlep": [
                "Lepton_Pt",
                "category",
            ],
        }

        # Add base branches
        for b in base_branches:
            branchList.push_back(b)

        # Add extra depending on mode
        for b in mode_branches.get(mode, []):
            branchList.push_back(b)

        if mode not in mode_map:
            raise ValueError(f"Unknown mode: {mode}")

        outputFile = f"{myDir}snapshot_mc_{mc}_{year}_{mode_map[mode]}.root"

        snapshotOptions = ROOT.RDF.RSnapshotOptions()
        snapshotOptions.fCompressionAlgorithm = ROOT.kLZ4
        snapshot_tdf = df.Snapshot("events", outputFile, branchList, snapshotOptions)
        print("snapshot_rdf DONE")
        print(outputFile)

        now = datetime.now()
        print('==> ends: ',now)


def loopOnDataset(year):

    thisdict = BuildDict(year)

    ## both data and MC
    loadCorrectionSet(int(year))

    mc = []
    mc.extend([10,11,12,13,14,15])
    if year in ["12022", "22022", "12023", "22023"]: mc.extend([20,21,22,23,24,25])
    if year=="2024": mc.extend([20,21,22,23,24,26])
    if year=="2024": mc.extend([103,104])
    else: mc.extend([100])
    mc.extend([102,105,106])
    mc.extend([201,202,203,204,205])
    mc.extend([211,212,213,214])

    print(mc)

    for sampleNOW in mc:
        files, xsec = SwitchSample(thisdict, sampleNOW)
        print(f"mc={mc}, outside the function: {len(files)}")
        rdf = ROOT.RDataFrame("Runs", files) # make sure this is not the distributed

        sumW = computeWeigths(rdf,xsec)
        analysis(files,year,sampleNOW,sumW)

    data_map = {
        "12022": [-11, -12, -13, -14],
        "22022": [-15, -16, -17],
        "12023": [-23, -24],
        "22023": [-31, -32],
        "2024":  list(range(-41, -55, -1)),  # generates -41 to -54
    }
    data = data_map.get(year, [])

    readDataQuality(year)

    for sampleNOW in data:
        files = SwitchSample(thisdict,sampleNOW)[0]
        print('outside the function: ', len(files))
        analysis(files,year,sampleNOW,1.)

if __name__ == "__main__":

    now = datetime.now()
    print('==> very beginning: ',now)

    year=sys.argv[1]
    print('year=',year)

    loopOnDataset(year)

