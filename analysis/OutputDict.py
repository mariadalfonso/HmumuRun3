import ROOT

def DefineBranchList(mode,mc):

    branchList = ROOT.vector('string')()
    branchListSynch = ROOT.vector('string')()        

    branches_synch = [
        "run",
        "event",
        "luminosityBlock",
        #            
        "HiggsCandCorrMass",
        "isZ",
        "isH",
        "isHSB",
        "isZandHSB",
        "Muon1_pt",
        "Muon2_pt",
        "Muon1_eta",
        "Muon2_eta",
        #
        "jetVBF2_Pt",
        "jetVBF1_Pt",
        "jetVBF2_Eta",
        "jetVBF1_Eta",
        "Mjj"            
    ]            
    
    # Add base branches
    for b in branches_synch:
        branchListSynch.push_back(b)
        
        base_branches = [
            "mc",
            "w",
            "w_allSF",
            "lumiIntegrated",
            "PV_npvsGood",
            "run",
            "event",
            "luminosityBlock",
            "boson_ptWeight", #DYturbo weights
            #
            "HiggsCandCorrMass",
            "HiggsCandCorrPt",
            "HiggsCandCorrRapidity",
            "HiggsCandMassErr",
            "cosThetaCS",
            "phiStarCS",
            "classify",
            "isZ",
            "isH",
            "isHSB",
            "isZandHSB",
            #
            "Muon1_pt",
            "Muon2_pt",
            "Muon1_eta",
            "Muon2_eta",
#                "Muon1_bsConstrainedChi2",
#                "Muon2_bsConstrainedChi2",
#                "Muon1_jetPtRel",
#                "Muon2_jetPtRel",
#                "Muon1_jetRelIso",
#                "Muon2_jetRelIso",
#                "FsrPH1_eta",
#                "FsrPH2_eta",
#                "FsrPH1_relIso03",
#                "FsrPH2_relIso03",
#                "FsrPH1_dROverEt2",
#                "FsrPH2_dROverEt2",
#                "FsrPH_pt_ratio0",
#                "FsrPH_pt_ratio1",
            "PuppiMET_pt",
            "PuppiMET_phi",
            "discrMVA"
        ]

        # Mode-specific branches
        mode_branches = {
            "isVBF": [
                "jetVBF2_Pt",
                "jetVBF1_Pt",
                "jetVBF2_Eta",
                "jetVBF1_Eta",
                "jetVBF2_Phi",
                "jetVBF1_Phi",
                "minDR_jetVBF1_Mu",
                "minDR_jetVBF2_Mu",
                #
                "Mjj",
                "dEtaJJ",
                "dPhiJJ",
                "ZepVar",
                "CenEta",
                "minDetaDiMuVBF",
                "minDphiDiMuVBF",
                "RPt",
                "CenPt",
                #
                "jetVBF2_hfcentralEtaStripSize",
                "jetVBF1_hfcentralEtaStripSize",
                "jetVBF2_hfadjacentEtaStripsSize",
                "jetVBF1_hfadjacentEtaStripsSize",
#                "jetVBF1_hfsigmaPhiPhi",
#                "jetVBF2_hfsigmaPhiPhi",
#                "jetVBF1_hfsigmaEtaEta",
#                "jetVBF2_hfsigmaEtaEta",
#                "jetVBF1_dPhiMET",
#                "jetVBF2_dPhiMET",
	    ],
            "isGGH": [
		"nGoodJetsAll",
	        "Jet1_Pt",
                "Jet1_Eta",
                "nGoodJetsTrk",
                "deltaRJet1H"                
            ],
            "isZinv": [
                "Muon1_phi",
                "Muon2_phi",
                "Muon1_sip3d",
                "Muon2_sip3d",
                "dPhiMETH",
                "RPt",
            ],
	    "isVlep": [
                "Muon1_promptMVA",
                "Muon2_promptMVA",
                #
                "category",
                "Muon1_sip3d",
                "Muon2_sip3d",
                "Lepton_promptMVA",
                "Lepton2_promptMVA",
                "Lepton_Eta",
                "Lepton2_Eta",
                "Lepton_Pt",
                "Lepton2_Pt",
                "Lepton_charge",
                "Lepton2_charge",
                "dEtaVH",
                "dPhiVH",
                "VMass",
                "ZMassPull",
                "WMassPull",
                "RPt",
            ],
            "isTTlep": [
                "Muon1_promptMVA",
                "Muon2_promptMVA",
                #
                "Lepton_promptMVA",
                "Lepton2_promptMVA",
                "Lepton_Pt",
                "Lepton2_Pt",
                "Lepton_Eta",
                "Lepton2_Eta",
                "Lepton_charge",
                "Lepton2_charge",
                #
                "dEtaLepH",
                "dPhiLepH",
                "dPhiMETH",
                "MetBisectorProj",
                "mbb",
                #
                "category",
                "Jet1_Pt",
                "Jet1_Eta",
                "HT",
                "ST",
                "dEta_j1j2",
                "Centrality",
                "mt"
            ],
            "isTThad": [
                "category",
                "Muon1_sip3d",
                "Muon2_sip3d",
                "Jet1_Pt",
                "Jet1_Eta",
                "JetAll_Eta",
                "HT",
                "WTopJetMass",
                "WTopJetDiscr",
                "nGoodJetsAll",
                "Centrality",
                "TopMassReco",
                "TopPairChi2",
                "dEta_j1j2",
                "MetBisectorProj",
                "dPhiMETH",
                "mindR_H_BJet", #used
                "nBMJets"
            ],
            "isVhad": [
                "goodWjj_mass",
                "goodWjj_discr",
#                "goodWjj_discr2",
                "goodWjj_pt",
                "goodWjj_eta",
                "goodWjj_phi",
                "dEtaWjjH",
                "dPhiWjjH",
                "RPt",
            ],
        }

        mode_MConly_branches = {
            "isGGH": [
                "Muon1_genPartFlav",
                "Muon2_genPartFlav",
#                "HTXS_Higgs_pt_new",
#                "HTXS_Higgs_pt",
#                "HTXS_njets30"
#                "boson_genMass",
#                "LHE_Vpt"
                "boson_genMassZ",
                "boson_genPtZ",
            ],
            "isVBF": [
                "Muon1_genPartFlav",
                "Muon2_genPartFlav",
                "jetVBF1_partonFlavour",
                "jetVBF2_partonFlavour",
                "jetVBF1_LHE",
                "jetVBF2_LHE",
                "jetAll_LHE",
                "boson_genMassZ",
                "boson_genPtZ",
            ],
            "isTTlep": [
                "Muon1_genPartFlav",
                "Muon2_genPartFlav",
                "Lepton_genPartFlav",
            ],
            "isTThad": [
                "Muon1_genPartFlav",
                "Muon2_genPartFlav",
            ],
            "isVlep": [
                "Muon1_genPartFlav",
                "Muon2_genPartFlav",
                "Lepton_genPartFlav",
#                "boson_genPtLHE",
#                "higgs_genPtLHE",
                "boson_genMassZ",
                "boson_genPtZ",
            ],
            "isZinv": [
                "Muon1_genPartFlav",
                "Muon2_genPartFlav",
#                "boson_genPtLHE",
#                "higgs_genPtLHE",
                "boson_genMassZ",
                "boson_genPtZ",
            ],
            "isVhad": [
                "Muon1_genPartFlav",
                "Muon2_genPartFlav",
#                "boson_genPtLHE",
#                "higgs_genPtLHE"
                "boson_genMassZ",
                "boson_genPtZ",
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

        return branchList
