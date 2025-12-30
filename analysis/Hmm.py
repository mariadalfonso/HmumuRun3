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
from utilsAna import loadtmvahelper

loadUserCode()
loadtmvahelper()
import helper_tmva # need various definitions

ROOT.gROOT.SetBatch()
ROOT.ROOT.EnableImplicitMT()

with open("/home/submit/mariadlf/Hrare/CMSSW_10_6_27_new/src/HmumuRun3/analysis/config/selection.json") as jsonFile:
    jsonObject = json.load(jsonFile)
    jsonFile.close()

selections = {key: jsonObject[key] for key in ["GOODMUON", "GOODJETSALL", "GOODelectrons", "LOOSEelectrons", "LOOSEmuons", "BJETS", "METFilters", "TRIGGER", "GOODboson" ]}

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
    "isVhad":  "VHcat",
    "isTTlep": "TTLcat",
    "isTThad": "TTHcat",
}

MVA_map = {
    "isVBF":   "MVA/output/classification_model_VBFcat.root",
    "isGGH":   "MVA/output/classification_model_ggHcat.root",
}

def callMVAclassification(df):

    modelname = f"bdt_model_{mode_map[mode]}"
    fileName = MVA_map[mode]
    tmva_helper = helper_tmva.TMVAHelperXGB(fileName, modelname)
    print(tmva_helper.variables)
    dfWithMVA = tmva_helper.run_inference(df,"discrMVA")
    dfWithMVA = dfWithMVA.Define("discrMVA0","discrMVA[0]")

    return dfWithMVA

def dfwithSYST(df,year):

#    strWPID="L"
#    strWPISO="L"
    strWPID="M"
    strWPISO="T"

    df = (df.Define("SFmuon1_ID_Nom",'corr_sf.eval_muonIDSF("{0}", "nominal", Muon1_eta, Muon1_pt, {1})'.format(year,strWPID))
          .Define("SFmuon1_ID_Up",'corr_sf.eval_muonIDSF("{0}", "systup", Muon1_eta, Muon1_pt, {1})'.format(year,strWPID))
          .Define("SFmuon1_ID_Dn",'corr_sf.eval_muonIDSF("{0}", "systdown", Muon1_eta, Muon1_pt, {1})'.format(year,strWPID))
          .Define("SFmuon2_ID_Nom",'corr_sf.eval_muonIDSF("{0}", "nominal", Muon2_eta, Muon2_pt, {1})'.format(year,strWPID))
          .Define("SFmuon2_ID_Up",'corr_sf.eval_muonIDSF("{0}", "systup", Muon2_eta, Muon2_pt, {1})'.format(year,strWPID))
          .Define("SFmuon2_ID_Dn",'corr_sf.eval_muonIDSF("{0}", "systdown", Muon2_eta, Muon2_pt, {1})'.format(year,strWPID))
          #
          .Define("SFmuon1_ISO_Nom",'corr_sf.eval_muonISOSF("{0}", "nominal", Muon1_eta, Muon1_pt, {1})'.format(year,strWPISO))
          .Define("SFmuon1_ISO_Up",'corr_sf.eval_muonISOSF("{0}", "systup", Muon1_eta, Muon1_pt, {1})'.format(year,strWPISO))
          .Define("SFmuon1_ISO_Dn",'corr_sf.eval_muonISOSF("{0}", "systdown", Muon1_eta, Muon1_pt, {1})'.format(year,strWPISO))
          .Define("SFmuon2_ISO_Nom",'corr_sf.eval_muonISOSF("{0}", "nominal", Muon2_eta, Muon2_pt, {1})'.format(year,strWPISO))
          .Define("SFmuon2_ISO_Up",'corr_sf.eval_muonISOSF("{0}", "systup", Muon2_eta, Muon2_pt, {1})'.format(year,strWPISO))
          .Define("SFmuon2_ISO_Dn",'corr_sf.eval_muonISOSF("{0}", "systdown", Muon2_eta, Muon2_pt, {1})'.format(year,strWPISO))
        )

    if year=='12022' or year=='22022' or year=='12023' or year=='22023' or year=='2024':
        df = (df.Define("SFpu_Nom",'corr_sf.eval_puSF(Pileup_nTrueInt,"nominal")')
              .Define("SFpu_Up",'corr_sf.eval_puSF(Pileup_nTrueInt,"up")')
              .Define("SFpu_Dn",'corr_sf.eval_puSF(Pileup_nTrueInt,"down")')
        )

    dfFinal_withSF = (df.Define("w_allSF", "w*SFpu_Nom*SFmuon1_ID_Nom*SFmuon2_ID_Nom*SFmuon1_ISO_Nom*SFmuon2_ISO_Nom"))

    return dfFinal_withSF

def doCategories(df,mc,year):

    #mode = sys.argv[2]  # expected: isVBF, isGGH, isW, isZinv

    # --- Common jet selections ---
    if year == '202XXX':
        BJETSmedium = selections["BJETS"].format("UParTAK4B>"+str(btagPNetBM[year]))
        BJETSloose  = selections["BJETS"].format("UParTAK4B>"+str(btagPNetBL[year]))
    else:
        BJETSmedium = selections["BJETS"].format("PNetB>"+str(btagPNetBM[year]))  # year=='2024' in principle should be PNet
        BJETSloose  = selections["BJETS"].format("PNetB>"+str(btagPNetBL[year]))

    muonSel = selections["GOODMUON"].format("")  
    if mode == "isVlep" or mode == "isTTlep"or mode == "isZinv": muonSel = selections["GOODMUON"].format(" and Muon_sip3d<5")
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
        "isVhad":      mu_veto + ele_veto + bjet_veto,  # TO add the met veto
        "isTTlep":     bjet,
        "isTThad":     mu_veto + ele_veto + bjet,
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

    if mode == "isZinv":
        df= (df.Filter("PuppiMET_pt>150","MET > 150")
             .Define("METFilters","{}".format(selections["METFilters"]))
             .Filter("METFilters>0","Pass MET filters")
             )

    elif mode == "isVhad":

        df= (df.Define("fatjet_muon_clean","fatJetMask(FatJet_muonIdx3SJ, index_Mu[0], index_Mu[1])")
             .Define("dR2fatjetH_mask","deltaRMask((Muon1Vec+Muon2Vec).Eta(),(Muon1Vec+Muon2Vec).Phi(),FatJet_eta,FatJet_phi)")
             .Define("goodWjj","{}".format(selections["GOODboson"]))
             .Filter("Sum(goodWjj)>0","at least one W jet")
             # pick first W jet
             .Define("goodWjj_idx", "Nonzero(goodWjj)[0]")             
             .Define("goodWjj_mass", "FatJet_particleNet_massCorr[goodWjj_idx] * FatJet_mass[goodWjj_idx]")
             .Define("goodWjj_discr", "FatJet_particleNetWithMass_WvsQCD[goodWjj_idx]")
             .Define("goodWjj_pt", "FatJet_pt[goodWjj_idx]")
             .Define("goodWjj_eta", "FatJet_eta[goodWjj_idx]")
             .Define("goodWjj_phi", "FatJet_phi[goodWjj_idx]")
             .Define("dEtaWjjH","abs((Muon1Vec+Muon2Vec).Eta()-goodWjj_eta)")
             .Define("dPhiWjjH","deltaPhi((Muon1Vec+Muon2Vec).Phi(),goodWjj_phi)")
             )

    elif mode == "isVlep":

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

        eleSel = selections["GOODelectrons"].format("")
        if year=='2024' or year=='2025': eleSel = selections["GOODelectrons"].format(" and Electron_pfRelIso04_all < 0.15")
        else: eleSel = selections["GOODelectrons"].format(" and Electron_pfRelIso03_all < 0.15")

        df= (df.Define("goodElectrons","{}".format(eleSel))
             .Define("idx_goodEls", "Nonzero(goodElectrons)")
             .Define("category", f"(int)({expr})")
             .Filter("category > 0", "Has valid category") #“discard events with no category”
             .Define("Lepton_Pt","((category ==1 or category ==2) and idx_goodEls.size() > 0) ? Electron_pt[idx_goodEls[0]]:-1")
             .Define("Lepton_Phi","((category ==1 or category ==2) and idx_goodEls.size() > 0) ? Electron_phi[idx_goodEls[0]]:-999")
             .Define("mt","((category ==1) and idx_goodEls.size() > 0)? mt(Lepton_Pt, Lepton_Phi, PuppiMET_pt, PuppiMET_phi):-1")
             .Define("Lepton_sip3d","((category ==1 or category ==2) and idx_goodEls.size() > 0) ? Electron_sip3d[idx_goodEls[0]]:-999")
             )

    elif mode == "isTTlep":

        # Define leptonic category conditions
        category_map = {
            1: f"({n_allJets}>2 && {n_looseEle}==1 && {n_goodEle}==1 && {n_goodMu}==2 && {q_goodMu}==0 && {n_looseMu}==2)",  # W→e  (TTH semilep)
            2: f"({n_allJets}>0 && {n_looseEle}==2 && {q_looseEle}==0 && {q_goodMu}==0 && {n_looseMu}==2)",                  # 2W→e (TTH dilep)
            3: f"({n_allJets}>2 && {n_looseEle}==0 && {n_goodMu}==3 && abs({q_goodMu})==1 && {n_looseMu}==3 && {freeOfZ})",  # W→μ  (TTH semilep)
            4: f"({n_allJets}>0 && {n_looseEle}==0 && {n_goodMu}==4 && {q_goodMu}==0 && {n_looseMu}==4)",                    # 2W→μ (TTH dilep)
        }

        JETSwithEleVeto = selections["GOODJETSALL"].format("and jetElectron_mask")
        print(JETSwithEleVeto)

        # Build a single expression with nested ternaries
        expr = " : ".join([f"({cond})?{cat}" for cat, cond in category_map.items()]) + " : 0"

        # 3 muons The charge of the three leptons must add up to ±1
        # At least one μ+μ− pair must have an invariant mass between 110 and 150 GeV
        # In 3μ events, the non-Higgs-candidate μ+μ− pair (μμOS) must not have an invariant mass between 81 and 101 GeV, to suppress WZ and Z+jets backgrounds

        eleSel = selections["GOODelectrons"].format("")
        if year=='2024' or year=='2025': eleSel = selections["GOODelectrons"].format(" and Electron_pfRelIso04_all < 0.15")
        else: eleSel = selections["GOODelectrons"].format(" and Electron_pfRelIso03_all < 0.15")

        df= (df.Define("goodElectrons","{}".format(eleSel))
             .Define("idx_goodEls", "Nonzero(goodElectrons)")
             .Define("jetElectron_mask", "cleaningMask(Electron_jetIdx[goodElectrons],nJet)")
             .Define("goodJetsAll","{}".format(JETSwithEleVeto))
             .Define("nGoodJetsAll","Sum(goodJetsAll)*1.0f")
             .Define("category", f"(int)({expr})")
             .Filter("category > 0", "Has valid category")
             .Define("Lepton_Pt","((category ==1 or category ==2) and idx_goodEls.size() > 0) ? Electron_pt[idx_goodEls[0]]:-1")
             .Define("Lepton_Phi","((category ==1 or category ==2) and idx_goodEls.size() > 0) ? Electron_phi[idx_goodEls[0]]:-999")
             .Define("mt","((category ==1) and idx_goodEls.size() > 0)? mt(Lepton_Pt, Lepton_Phi, PuppiMET_pt, PuppiMET_phi):-1")
             .Define("Lepton_sip3d","((category ==1 or category ==2) and idx_goodEls.size() > 0) ? Electron_sip3d[idx_goodEls[0]]:-999")
             #
             .Define("hardestGoodJet_idx","hardest_pt_idx(Jet_pt[goodJetsAll])")
             .Define("Jet1_Pt","Jet_pt[goodJetsAll][hardestGoodJet_idx]")
             .Filter("Jet1_Pt>35")
             .Define("Jet1_Eta","Jet_eta[goodJetsAll][hardestGoodJet_idx]")
             .Define("HT", "ROOT::VecOps::Sum(Jet_pt[goodJetsAll])")
             )

    elif mode == "isTThad":

        myTOP = "FatJet_pt>400 && FatJet_particleNetWithMass_TvsQCD> 0.75 && FatJet_particleNet_massCorr*FatJet_mass>(175-80) && FatJet_particleNet_massCorr*FatJet_mass<(175+80) && fatjet_muon_clean"
        myW = "FatJet_pt>200 && FatJet_pt<400 && FatJet_particleNetWithMass_WvsQCD> 0.75 && FatJet_particleNet_massCorr*FatJet_mass>(80-40) && FatJet_particleNet_massCorr*FatJet_mass<(80+40) && fatjet_muon_clean"
        n_top = f"Sum({myTOP})"
        n_W = f"Sum({myW})"

        category_map = {
            1: f"({n_top}>=1)",     # with  TOP pt > 400 GeV
            2: f"({n_W}>=1 && {n_top}>=0)",       # with W pt > 200 GeV
            3: f"({n_top}==0 && {n_W}==0)",    # fully resolved   at least three jets
        }

        # Build a single expression with nested ternaries
        expr = f"({category_map[1]}) ? 1 : (({category_map[2]}) ? 2 : (({category_map[3]}) ? 3 : 0))"

        df= (df.Define("fatjet_muon_clean","fatJetMask(FatJet_muonIdx3SJ, index_Mu[0], index_Mu[1])")
             .Define("goodJetsAll","{}".format(selections["GOODJETSALL"].format(""))) # check eta distribution
             .Define("nGoodJetsAll","Sum(goodJetsAll)*1.0f")
             .Define("goodW","{}".format(myW)).Define("goodW_idx", "Nonzero(goodW)[0]")
             .Define("goodTOP","{}".format(myTOP)).Define("goodTOP_idx", "Nonzero(goodTOP)[0]")
             .Define("category", f"(int)({expr})")
             .Filter("category > 0", "Has valid category")
             .Define("hardestGoodJet_idx","hardest_pt_idx(Jet_pt[goodJetsAll])")
             .Define("Jet1_Pt","Jet_pt[goodJetsAll][hardestGoodJet_idx]")
             .Define("Jet1_Eta","Jet_eta[goodJetsAll][hardestGoodJet_idx]")
             .Define("HT", "ROOT::VecOps::Sum(Jet_pt[goodJetsAll])")
             .Filter("category==1 || category==2 || (category==3 && Jet1_Pt>50 && Sum(goodJetsAll)>3)")
             .Define("WTopJetMass","(category==1) ? (FatJet_particleNet_massCorr[goodTOP_idx]*FatJet_mass[goodTOP_idx]): (category==2) ? (FatJet_particleNet_massCorr[goodW_idx]*FatJet_mass[goodW_idx]): 0.f")
             .Define("WTopJetDiscr","(category==1) ? (FatJet_particleNetWithMass_TvsQCD[goodTOP_idx]): (category==2) ? (FatJet_particleNetWithMass_WvsQCD[goodW_idx]): 0.f")
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
             .Define("nGoodJetsAll","Sum(goodJetsAll)*1.0f").Filter("Sum(goodJetsAll)>1","two VBF jets") # this line is redundant
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
             .Define("Mjj","Pair12Minv(jetVBF1_Pt, jetVBF1_Eta, jetVBF1_Phi, jetVBF1_Mass, jetVBF2_Pt, jetVBF2_Eta, jetVBF2_Phi, jetVBF2_Mass)")
             .Define("dEtaJJ","abs(jetVBF1_Eta-jetVBF2_Eta)")
             .Define("dPhiJJ","deltaPhi(jetVBF1_Phi,jetVBF2_Phi)")
             .Filter("dEtaJJ>2.5","two jets with the Deta>2.5 ")
             .Define("minDR_jetVBF1_Mu","minDRmusJ(VBF1Vec, Muon1Vec, Muon2Vec)")
             .Define("minDR_jetVBF2_Mu","minDRmusJ(VBF2Vec, Muon1Vec, Muon2Vec)")
             .Define("RPt","getRpt(Muon1Vec, Muon2Vec, VBF1Vec, VBF2Vec)")
             .Define("ZepVar","getZep(Muon1Vec, Muon2Vec, jetVBF1_Eta, jetVBF2_Eta)") # should I pass the Jet.rapidity() and has very large values
             )

        if mc>0:
            df= (df.Define("jetVBF1_partonFlavour","Jet_partonFlavour[goodJetsAll][index_VBFJets[0]]")
                 .Define("jetVBF2_partonFlavour","Jet_partonFlavour[goodJetsAll][index_VBFJets[1]]")
                 .Define("jetAll_LHE","getLHEPart_match(jetAll_Eta, jetAll_Phi, LHEPart_eta, LHEPart_phi, LHEPart_status, LHEPart_pdgId, isData)")
                 .Define("jetVBF1_LHE","jetAll_LHE[index_VBFJets[0]]")
                 .Define("jetVBF2_LHE","jetAll_LHE[index_VBFJets[1]]")
                 )

    return df

def objScaleSmear(df, year, mc):

    if year in ["12022", "22022", "12023", "22023", "2024"]:
        # apply JEC AK4 and AK8
        df = df.Redefine("Jet_pt",'computeJECcorrection(corr_sf, Jet_pt, Jet_rawFactor, Jet_eta, Jet_phi, Jet_area, rho, run, isData, "{0}","{1}" )'.format(year,mc))
        df = df.Redefine("FatJet_pt",'computeJECcorrection(corr_sf, FatJet_pt, FatJet_rawFactor, FatJet_eta, FatJet_phi, FatJet_area, rho, run, isData, "{0}","{1}" )'.format(year,mc))

        # get the EGM scale
        df = df.Redefine("Electron_pt",'computeEleSSCorrection(corr_sf, Electron_pt, Electron_eta, Electron_r9, Electron_seedGain, event, run, isData, "{0}")'.format(year))

    # apply muonScale (2024 MC is still broken)
#    if year in ["12022", "22022", "12023", "22023", "2024", "2025"]:
    if year in ["12022", "22022", "12023", "22023"]:
        df = df.Redefine("Muon_bsConstrainedPt",'computeMUcorrection(Muon_bsConstrainedPt, Muon_eta, Muon_phi, Muon_charge, Muon_nTrackerLayers, isData, event, luminosityBlock)')

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

    dfComm = objScaleSmear(dfComm, year, mc)

    # compute JetID (note v15 vs v12)
    if year=="2024" or year=="2025":
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
          .Define("CSvars","CollinSopperAngles(Muon1Vec,Muon2Vec,Muon_charge[goodMuons][index_Mu[0]],Muon_charge[goodMuons][index_Mu[1]])")
          .Define("cosThetaCS","CSvars.first")
          .Define("phiStarCS","CSvars.second")
          .Filter("(HiggsCandCorrMass>70 and HiggsCandCorrMass<200)","HiggsMass within reasonable range 70-200")
          )

    if mode == "isVBF":
        df = (df.Define("minDetaDiMuVBF","minDeta(HiggsCandCorrRapidity, jetVBF1_Eta, jetVBF2_Eta)")
              )

    ## call MVA classification
    if mode == "isVBF" or mode == "isGGH": df = callMVAclassification(df)

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
#            "FsrPH_pt": {"name":"FsrPH_pt","title":"FSR0; FSR (GeV);N_{Events}","bin":200,"xmin":0.,"xmax":20.},
#            "FsrPH_pt_ratio0": {"name":"FsrPH_pt_ratio0","title":"FSR_ratio0; FSR/p_{T}(mu1) (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":10.},
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
                "cosThetaCS",
                "phiStarCS",
                "classify",
                "Muon1_pt",
                "Muon2_pt",
                "Muon1_eta",
                "Muon2_eta",
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
                "minDR_jetVBF1_Mu",
                "minDR_jetVBF2_Mu",
                "dEtaJJ",
                "dPhiJJ",
                "ZepVar",
                "minDetaDiMuVBF",
                "jetVBF2_hfcentralEtaStripSize",
                "jetVBF1_hfcentralEtaStripSize",
                "jetVBF2_hfadjacentEtaStripsSize",
                "jetVBF1_hfadjacentEtaStripsSize",
                "discrMVA0"
            ],
            "isGGH": [
                "discrMVA0"
            ],
            "isZinv": [
                "Muon1_phi",
                "Muon2_phi",
            ],
            "isVlep": [
                "Lepton_Pt",
                "Lepton_sip3d",
                "Muon1_sip3d",
                "Muon2_sip3d",
                "category",
                "mt"
            ],
            "isTTlep": [
                "Lepton_Pt",
                "Lepton_sip3d",
                "Muon1_sip3d",
                "Muon2_sip3d",
                "category",
                "Jet1_Pt",
                "Jet1_Eta",
                "HT",
                "nGoodJetsAll",
                "mt"
            ],
            "isTThad": [
                "category",
                "Muon1_sip3d",
                "Muon2_sip3d",
                "Jet1_Pt",
                "Jet1_Eta",
                "HT",
                "nGoodJetsAll",
                "WTopJetMass",
                "WTopJetDiscr"
            ],
            "isVhad": [
                "goodWjj_mass",
                "goodWjj_discr",
                "goodWjj_pt",
                "goodWjj_eta",
                "goodWjj_phi",
                "dEtaWjjH",
                "dPhiWjjH",
            ],
        }

        mode_MConly_branches = {
            "isVBF": [
                "jetVBF1_partonFlavour",
                "jetVBF2_partonFlavour",
                "jetVBF1_LHE",
                "jetVBF2_LHE",
                "jetAll_LHE"
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

        outputFile = f"{myDir}{mode_map[mode]}/snapshot_mc_{mc}_{year}_{mode_map[mode]}.root"

        snapshotOptions = ROOT.RDF.RSnapshotOptions()
        if mode == "isGGH":
            snapshotOptions.fCompressionAlgorithm = ROOT.kZLIB
            snapshotOptions.fCompressionLevel = 9
        else:
            snapshotOptions.fCompressionAlgorithm = ROOT.kLZ4
            snapshotOptions.fCompressionLevel = 9
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

    if mode == "isVlep" or mode == "isTTlep" or mode == "isVhad" or mode == "isGGH" or mode == "isTThad":
        mc.extend([107,105,106]) # tt1l, tW
        if year in ["12022", "22022", "12023", "22023"]: mc.extend([221,222,223,224,225]) #ttV
        if year=="2024": mc.extend([222,223,224,225])    #ttV

    # below for training
    if year in ["12022", "22022", "12023", "22023"]:
        if mode == "isVhad": mc.extend([111,112,113])  # extra DY jet binned for ML training
        if mode == "isVhad": mc.extend([114,115,116,117]) # extra DY pt binned for ML training
    if mode == "isVBF" or mode == "isGGH" or mode == "isVhad" or mode == "isTThad": mc.extend([108,109,110]) # extra DY jet mass binned

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
        "2025":  list(range(-61, -73, -1)),  # generates -61 to -70
    }

    data_map_jpsi = {
#        "12022": [-11, -12, -13, -14],
        "22022": [-116,-117], # missing (E -115)
        "12023": [-123],
        "22023": [-131],
        "2024":  list(range(-141, -148, -1)),
        "2025":  list(range(-161, -167, -1)),  # generates -61 to -70
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

