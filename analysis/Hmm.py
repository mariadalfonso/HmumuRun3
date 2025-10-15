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

mode = sys.argv[2]  # expected: isVBF, isGGH, isVlep, isWhad, isTTlep, isTThad, isZinv

# --- place to write the files and strings  ---
myDir = '/work/submit/mariadlf/HmumuRun3/ROOTFILES/'
# --- Histogram output suffixes per mode ---
mode_map = {
    "isVBF":   "VBFcat",
    "isGGH":   "ggHcat",
    "isZinv":  "Zinvcat",
    "isVlep":  "VLcat",
    "isTTlep": "TTLcat",
    "isTThad": "TTHcat",
}

def dfwithSYST(df,year):

    df = (df.Define("SFmuon1_ID_Nom",'corr_sf.eval_muonIDSF("{0}", "nominal", Muon1_eta, Muon1_pt, "M")'.format(year))
          .Define("SFmuon1_ID_Up",'corr_sf.eval_muonIDSF("{0}", "systup", Muon1_eta, Muon1_pt, "M")'.format(year))
          .Define("SFmuon1_ID_Dn",'corr_sf.eval_muonIDSF("{0}", "systdown", Muon1_eta, Muon1_pt, "M")'.format(year))
          .Define("SFmuon2_ID_Nom",'corr_sf.eval_muonIDSF("{0}", "nominal", Muon2_eta, Muon2_pt, "M")'.format(year))
          .Define("SFmuon2_ID_Up",'corr_sf.eval_muonIDSF("{0}", "systup", Muon2_eta, Muon2_pt, "M")'.format(year))
          .Define("SFmuon2_ID_Dn",'corr_sf.eval_muonIDSF("{0}", "systdown", Muon2_eta, Muon2_pt, "M")'.format(year))
          #
          .Define("SFmuon1_ISO_Nom",'corr_sf.eval_muonISOSF("{0}", "nominal", Muon1_eta, Muon1_pt, "T")'.format(year))
          .Define("SFmuon1_ISO_Up",'corr_sf.eval_muonISOSF("{0}", "systup", Muon1_eta, Muon1_pt, "T")'.format(year))
          .Define("SFmuon1_ISO_Dn",'corr_sf.eval_muonISOSF("{0}", "systdown", Muon1_eta, Muon1_pt, "T")'.format(year))
          .Define("SFmuon2_ISO_Nom",'corr_sf.eval_muonISOSF("{0}", "nominal", Muon2_eta, Muon2_pt, "T")'.format(year))
          .Define("SFmuon2_ISO_Up",'corr_sf.eval_muonISOSF("{0}", "systup", Muon2_eta, Muon2_pt, "T")'.format(year))
          .Define("SFmuon2_ISO_Dn",'corr_sf.eval_muonISOSF("{0}", "systdown", Muon2_eta, Muon2_pt, "T")'.format(year))
          )

    if year=='12022' or year=='22022' or year=='12023' or year=='22023':
        df = (df.Define("SFpu_Nom",'corr_sf.eval_puSF(Pileup_nTrueInt,"nominal")')
              .Define("SFpu_Up",'corr_sf.eval_puSF(Pileup_nTrueInt,"up")')
              .Define("SFpu_Dn",'corr_sf.eval_puSF(Pileup_nTrueInt,"down")')
        )


    # for now missing PU-2024
    if year=='2024':
        dfFinal_withSF = (df.Define("w_allSF", "w*SFmuon1_ID_Nom*SFmuon2_ID_Nom*SFmuon1_ISO_Nom*SFmuon2_ISO_Nom"))
#        dfFinal_withSF = (df.Define("w_allSF", "w*SFmuon1_ID_Nom*SFmuon2_ID_Nom"))
     else:
        dfFinal_withSF = (df.Define("w_allSF", "w*SFpu_Nom*SFmuon1_ID_Nom*SFmuon2_ID_Nom*SFmuon1_ISO_Nom*SFmuon2_ISO_Nom"))
#         dfFinal_withSF = (df.Define("w_allSF", "w*SFpu_Nom*SFmuon1_ID_Nom*SFmuon2_ID_Nom"))

    return dfFinal_withSF

def doCategories(df,mc,year):

    #mode = sys.argv[2]  # expected: isVBF, isGGH, isW, isZinv

    # --- Common jet selections ---
    BJETSmedium = selections["BJETS"].format(btagPNetBM[year])
    BJETSloose  = selections["BJETS"].format(btagPNetBL[year])

    muonSel = selections["GOODMUON"].format("")  
    if mode == "isVlep" or mode == "isTTlep": muonSel = selections["GOODMUON"].format(" and Muon_sip3d<5")
    print("GOODMUON = ",muonSel)

    #fix the muonPT it's 25 GeV for the one triggered
    df = (df.Define("goodMuons","{}".format(muonSel))
          .Filter("Sum(goodMuons)>=1","at least two good muons")
          .Define("index_Mu",'getMuonIndices(Muon_bsConstrainedPt[goodMuons], Muon_eta[goodMuons], Muon_phi[goodMuons], Muon_charge[goodMuons], "{}")'.format(mode))
          .Filter("index_Mu[0]!= -1", "both Muons OS")
          .Filter("index_Mu[1]!= -1", "both Muons OS")
          .Define("Muon1_pt","Muon_bsConstrainedPt[goodMuons][index_Mu[0]]")
          .Define("Muon2_pt","Muon_bsConstrainedPt[goodMuons][index_Mu[1]]")
          .Filter("(Muon1_pt > 26 || Muon2_pt > 26)","at least one >26 GeV to satisfy the trigger")
          .Define("Muon1_eta","Muon_eta[goodMuons][index_Mu[0]]")
          .Define("Muon2_eta","Muon_eta[goodMuons][index_Mu[1]]")
          .Define("Muon1_phi","Muon_phi[goodMuons][index_Mu[0]]")
          .Define("Muon2_phi","Muon_phi[goodMuons][index_Mu[1]]")
          .Define("Muon1_sip3d","Muon_sip3d[goodMuons][index_Mu[0]]")
          .Define("Muon2_sip3d","Muon_sip3d[goodMuons][index_Mu[1]]")
          .Define("classify","topology(Muon1_eta, Muon2_eta)")
          .Define("Muon1Vec", "MakeTLV(Muon1_pt, Muon1_eta, Muon1_phi, Muon_mass[goodMuons][index_Mu[0]])")
          .Define("Muon2Vec", "MakeTLV(Muon2_pt, Muon2_eta, Muon2_phi, Muon_mass[goodMuons][index_Mu[0]])")
          ###
          .Define("jetMuon_mask", "cleaningMask(Muon_jetIdx[goodMuons],nJet)")
          # Jet_puIdDisc only for nanov15
          .Define("jetVeto_mask",'cleaningJetVetoMapMask(Jet_eta,Jet_phi,"{0}")'.format(year))
          .Define("BJETS", f"{BJETSmedium}")
          .Define("BJETSloose", f"{BJETSloose}")
          )

    if mc>0:
        df = (df.Define("Muon1_genPartFlav","Muon_genPartFlav[goodMuons][index_Mu[0]]")
              .Define("Muon2_genPartFlav","Muon_genPartFlav[goodMuons][index_Mu[1]]")
              )

    # --- Reusable blocks ---

    n_looseEle = f"Sum({selections['LOOSEelectrons']})"
    n_goodEle  = "Sum(goodElectrons)"
    q_looseEle  = "Sum(Electron_charge[{}])".format(selections['LOOSEelectrons'])
    n_goodMu   = "Sum(goodMuons)"
    n_looseMu  = f"Sum({selections['LOOSEmuons']})"
    q_goodMu   = "Sum(Muon_charge[goodMuons])"
    freeOfZ    = "freeOfZ(Muon_bsConstrainedPt[goodMuons], Muon_eta[goodMuons], Muon_phi[goodMuons], Muon_charge[goodMuons])"
    n_allJets  = "Sum(goodJetsAll)*1.0f"

    mu_veto = [
        ("FILTER", f"({n_looseMu}==2 && {n_goodMu}==2 && {q_goodMu}==0)", "no extra muons, only 2 OS charge"),
        ]

    ele_veto = [
        ("FILTER", f"({n_looseEle}==0)", "no extra electrons"),
    ]

    bjet_veto = [
        ("FILTER", "Sum(BJETS)==0 && Sum(BJETSloose)<2", "events vetoed if at least 1bM and 2bL"),
    ]

    bjet = [
        ("FILTER", "Sum(BJETS)>0 or Sum(BJETSloose)>1", "events with at least 1bM or 2bL"),
    ]

    met_filters = [
        ("FILTER", "PuppiMET_pt>150", "MET > 150"),
        ("FILTER", selections["METFilters"], "apply MET filters"),
    ]

    # --- VBF-specific filters (add the ones you use in doCategories/isVBF) ---
#    vbf_filters = [
#        ("FILTER", "Sum(goodJetsAll)>1", "at least 2 jets"),
#        ("FILTER", "Mjj>400", "M_{jj} > 400 GeV"),
#        ("FILTER", "dEtaJJ>2.5", "|dEta_{jj}| > 2.5"),
#    ]

    # --- Category recipes ---
    # PUT first the things that cut most of the events i.e. met in Zinv
    categories = {
        "isGGH":       mu_veto + ele_veto + bjet_veto, # TO add the met veto and VBFjet veto
        "isVBF":       mu_veto + ele_veto + bjet_veto, # TO add the met veto
        "isZinv":      mu_veto + ele_veto + bjet_veto + met_filters,
        "isVlep":      bjet_veto,
#        "isWhad":      ele_vet + bjet_veto,
        "isTTlep":     bjet,
        "isTThad":    mu_veto + ele_veto + bjet,
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

    if mode == "isVlep":

        # Define leptonic category conditions
        category_map = {
            1: f"({n_looseEle}==1 && {n_goodEle}==1 && {n_goodMu}==2 && {q_goodMu}==0 && {n_looseMu}==2)",  # W→e  (WH)
            2: f"({n_looseEle}==2 && {q_looseEle}==0 && {q_goodMu}==0 && {n_looseMu}==2)",                  # Z→ee (ZH)
            3: f"({n_looseEle}==0 && {n_goodMu}==3 && abs({q_goodMu})==1 && {n_looseMu}==3 && {freeOfZ})",  # W→μ  (WH)
            4: f"({n_looseEle}==0 && {n_goodMu}==4 && {q_goodMu}==0 && {n_looseMu}==4)",                    # Z→μμ (ZH)
        }

        # Build a single expression with nested ternaries
        expr = " : ".join([f"({cond})?{cat}" for cat, cond in category_map.items()]) + " : 0"

        # 3 muons The charge of the three leptons must add up to ±1
        # At least one μ+μ− pair must have an invariant mass between 110 and 150 GeV
        # In 3μ events, the non-Higgs-candidate μ+μ− pair (μμOS) must not have an invariant mass between 81 and 101 GeV, to suppress WZ and Z+jets backgrounds
        # in 4μ events, if 2 inZ are found then discard
        df= (df.Define("goodElectrons","{}".format(selections["GOODelectrons"]))
             .Define("category", f"(int)({expr})")
             .Filter("category > 0", "Has valid category") #“discard events with no category”
             .Define("Lepton_Pt","(category ==1 or category ==2) ? Electron_pt[goodElectrons][0]:0")
             .Define("Lepton_sip3d","(category ==1 or category ==2) ? Electron_sip3d[goodElectrons][0]:0")
             )

    elif mode == "isTTlep":

        # Define leptonic category conditions
        category_map = {
            1: f"({n_allJets}>2 && {n_looseEle}==1 && {n_goodEle}==1 && {n_goodMu}==2 && {q_goodMu}==0 && {n_looseMu}==2)",  # W→e  (TTH semilep)
            2: f"({n_allJets}>0 && {n_looseEle}==2 && {q_looseEle}==0 && {q_goodMu}==0 && {n_looseMu}==2)",                  # 2W→e (TTH dilep)
            3: f"({n_allJets}>2 && {n_looseEle}==0 && {n_goodMu}==3 && abs({q_goodMu})==1 && {n_looseMu}==3)",               # W→μ  (TTH semilep)
            4: f"({n_allJets}>0 && {n_looseEle}==0 && {n_goodMu}==4 && {q_goodMu}==0 && {n_looseMu}==4)",                    # 2W→μ (TTH dilep)
        }

        JETSwithEleVeto = selections["GOODJETSALL"].format("and jetElectron_mask")
        print(JETSwithEleVeto)

        # Build a single expression with nested ternaries
        expr = " : ".join([f"({cond})?{cat}" for cat, cond in category_map.items()]) + " : 0"

        # 3 muons The charge of the three leptons must add up to ±1
        # At least one μ+μ− pair must have an invariant mass between 110 and 150 GeV
        # In 3μ events, the non-Higgs-candidate μ+μ− pair (μμOS) must not have an invariant mass between 81 and 101 GeV, to suppress WZ and Z+jets backgrounds
        df= (df.Define("goodElectrons","{}".format(selections["GOODelectrons"]))
             .Define("jetElectron_mask", "cleaningMask(Electron_jetIdx[goodElectrons],nJet)")
             .Define("goodJetsAll","{}".format(JETSwithEleVeto))
             .Define("category", f"(int)({expr})")
             .Filter("category > 0", "Has valid category")
             .Define("Lepton_Pt","(category ==1 or category ==2) ? Electron_pt[goodElectrons][0]:0")
             .Define("Lepton_sip3d","(category ==1 or category ==2) ? Electron_sip3d[goodElectrons][0]:0")
             .Define("hardestGoodJet_idx","hardest_pt_idx(Jet_pt[goodJetsAll])")
             .Define("Jet1_Pt","Jet_pt[goodJetsAll][hardestGoodJet_idx]")
             .Filter("Jet1_Pt>35")
             )

    elif mode == "isTThad":
        df= (df.Define("goodJetsAll","{}".format(selections["GOODJETSALL"].format("")))
             .Define("nGoodJetsAll","Sum(goodJetsAll)*1.0f").Filter("Sum(goodJetsAll)>2","at least three jets")
             ## add leading jet pt>50 and triplet mass mjjj 100-300
             ## look the toponium and the ttHcc
             )

    elif mode == "isVBF":
        df= (df.Define("goodJetsAll","{}".format(selections["GOODJETSALL"].format("")))
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
    dfComm = dfComm.Redefine("FatJet_pt",'computeJECcorrection(corr_sf, FatJet_pt, FatJet_rawFactor, FatJet_eta, FatJet_phi, FatJet_area, rho, run, isData, "{0}","{1}" )'.format(year,mc))
    # get the EGM scale
    dfComm = dfComm.Redefine("Electron_pt",'computeEleSSCorrection(corr_sf, Electron_pt, Electron_eta, Electron_r9, Electron_seedGain, event, run, isData, "{0}")'.format(year))
    # apply muonScale
    if year in ["12022", "22022", "12023", "22023"]:
        dfComm = dfComm.Redefine("Muon_bsConstrainedPt",'computeMUcorrection(Muon_bsConstrainedPt, Muon_eta, Muon_phi, Muon_charge, Muon_nTrackerLayers, isData, event, luminosityBlock)')

    # compute JetID (note v15 vs v12)
    if year=="2024":
        dfComm = dfComm.Define("jetID_mask", 'cleaningJetIDMask(Jet_eta, Jet_chHEF, Jet_neHEF, Jet_chEmEF, Jet_neEmEF, Jet_muEF, Jet_chMultiplicity, Jet_neMultiplicity, "{0}")'.format(year)) #  for nanov15 use the JetID
    else:
        dfComm = dfComm.Define("jetID_mask", "cleaningJetSelMask(0, Jet_eta, Jet_neHEF, Jet_chEmEF, Jet_muEF, Jet_neEmEF, Jet_jetId)") # redo JetID for nanoV12,

    df = doCategories(dfComm,mc,year)

    # Zgamma has these: |dz| < 1.0 cm, |dxy| < 0.5 cm, SIP3D < 4
    df = (df.Define("HiggsCandMass","Minv(Muon1Vec,Muon2Vec)")
          .Define("HiggsCandMassErr","MinvErr(Muon1_pt, Muon_bsConstrainedPtErr[goodMuons][index_Mu[0]], Muon2_pt, Muon_bsConstrainedPtErr[goodMuons][index_Mu[1]])")
          .Define("goodFRSphoton", "fsrMask(Muon_fsrPhotonIdx[goodMuons],nFsrPhoton)")
          .Define("FsrPH_pt","FsrPhoton_pt[goodFRSphoton]")
          .Define("FsrPH_eta","FsrPhoton_eta[goodFRSphoton]")
          .Define("FsrPH_pt_ratio0","FsrPhoton_pt[goodFRSphoton][0]/Muon1_pt")
          .Define("FsrPH_pt_ratio1","FsrPhoton_pt[goodFRSphoton][1]/Muon2_pt")
          .Define("HiggsCandCorrMass","MinvCorr(Muon1Vec,Muon2Vec,FsrPhoton_pt[goodFRSphoton], FsrPhoton_eta[goodFRSphoton],FsrPhoton_phi[goodFRSphoton],0)")
          .Define("HiggsCandCorrPt","MinvCorr(Muon1Vec,Muon2Vec,FsrPhoton_pt[goodFRSphoton], FsrPhoton_eta[goodFRSphoton],FsrPhoton_phi[goodFRSphoton],1)")
          .Define("HiggsCandCorrRapidity","MinvCorr(Muon1Vec,Muon2Vec,FsrPhoton_pt[goodFRSphoton], FsrPhoton_eta[goodFRSphoton],FsrPhoton_phi[goodFRSphoton],2)")
          .Filter("HiggsCandCorrMass>50 and HiggsCandCorrMass<200","HiggsMass within reasonable range 50-200")
          )

    if mode == "isVBF":
        df = (df.Define("minDetaDiMuVBF","minDeta(HiggsCandCorrRapidity, jetVBF1_Eta, jetVBF2_Eta)")
              )

    ## FSR and muonBeamSpot also for Z ?
    ## check the FSR for electrons

    if (isData == "false"):
        df = dfwithSYST(df,year)
    else:
        df = df.Define("w_allSF", "w")

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
                "HiggsCandMassErr",
                "classify",
                "Muon1_pt",
                "Muon2_pt",
                "Muon1_eta",
                "Muon2_eta",
                "Muon1_phi",
                "Muon2_phi",
                "Muon1_sip3d",
                "Muon2_sip3d",
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
                "Lepton_sip3d",
                "category",
            ],
            "isTTlep": [
                "Lepton_Pt",
                "Lepton_sip3d",
                "category",
                "Jet1_Pt"
            ],
        }

        mode_MConly_branches = {
            "isVBF": [
                "jetVBF1_partonFlavour",
                "jetVBF2_partonFlavour",
            ],
            "isVlep": [
                "Muon1_genPartFlav",
                "Muon2_genPartFlav",
            ],
            "isTTlep": [
                "Muon1_genPartFlav",
                "Muon2_genPartFlav",
            ],
        }

        # Add base branches
        for b in base_branches:
            branchList.push_back(b)

        # Add extra depending on mode
        for b in mode_branches.get(mode, []):
            branchList.push_back(b)

        # Add extra depending on mode
        if mc > 0:
            for b in mode_MConly_branches.get(mode, []):
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
    if year=="2024": mc.extend([103,104])   #DY
    else: mc.extend([100])    #DY
    mc.extend([101])    #DY EWK
    mc.extend([102])    #TT2l
    mc.extend([201,202,203,204,205]) #VV
    mc.extend([211,212,213,214])     #VVV
    if mode == "isVlep" or mode == "isTTlep" or mode == "isVhad" or mode == "isGGH":
        mc.extend([107,105,106]) # tt1l, tW
        if year in ["12022", "22022", "12023", "22023"]: mc.extend([221,222,223,224,225]) #ttV
        if year=="2024": mc.extend([222,223,224])    #ttV

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
#    data = []

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

