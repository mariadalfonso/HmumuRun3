import ROOT
import json
import sys

from utilsAna import loadUserCode
from utilsAna import BuildDict, SwitchSample
from utilsAna import readDataQuality
from utilsAna import computeWeigths
from utilsAna import lumis
from utilsAna import loadCorrectionSet
from datetime import datetime

loadUserCode()

ROOT.gROOT.SetBatch()
ROOT.ROOT.EnableImplicitMT()

with open("/home/submit/mariadlf/Hrare/CMSSW_10_6_27_new/src/HmumuRun3/analysis/config/selection.json") as jsonFile:
    jsonObject = json.load(jsonFile)
    jsonFile.close()

GOODMUON = jsonObject['GOODMUON']
GOODJETSALL = jsonObject['GOODJETSALL']
LOOSEelectrons = jsonObject['LOOSEelectrons']
LOOSEmuons = jsonObject['LOOSEmuons']
BJETS = jsonObject['BJETS']
JSON = "isGoodRunLS(isData, run, luminosityBlock)"

doVBF = False
doGGH = False

if sys.argv[2]=='isVBF': doVBF = True
if sys.argv[2]=='isGGH': doGGH = True

def dfwithSYST(df,year):

    if year=='2024': dfFinal_withSF = df.Define("w_allSF", "w")

    else:

        dfFinal_withSF = (df
                          .Define("SFmuon1_ID_Nom",'corr_sf.eval_muonIDSF("{0}", "nominal", Muon_eta[goodMuons][0], Muon1_pt, "M")'.format(year))
                          .Define("SFmuon1_ID_Up",'corr_sf.eval_muonIDSF("{0}", "systup", Muon_eta[goodMuons][0], Muon1_pt, "M")'.format(year))
                          .Define("SFmuon1_ID_Dn",'corr_sf.eval_muonIDSF("{0}", "systdown", Muon_eta[goodMuons][0], Muon1_pt, "M")'.format(year))
                          .Define("SFmuon2_ID_Nom",'corr_sf.eval_muonIDSF("{0}", "nominal", Muon_eta[goodMuons][0], Muon2_pt, "M")'.format(year))
                          .Define("SFmuon2_ID_Up",'corr_sf.eval_muonIDSF("{0}", "systup", Muon_eta[goodMuons][0], Muon2_pt, "M")'.format(year))
                          .Define("SFmuon2_ID_Dn",'corr_sf.eval_muonIDSF("{0}", "systdown", Muon_eta[goodMuons][0], Muon2_pt, "M")'.format(year))
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

    TRIGGER="HLT_IsoMu24"
    print(TRIGGER)

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
              .Define("triggerAna","{}".format(TRIGGER))
              .Filter("triggerAna>0","passing trigger")
              )


    # apply JEC
    dfComm = dfComm.Redefine("Jet_pt",'computeJECcorrection(corr_sf, Jet_pt, Jet_rawFactor, Jet_eta, Jet_phi, Jet_area, rho, run, isData, "{0}","{1}" )'.format(year,mc))

    # compute JetID (note v15 vs v12)
    if year=="2024":
        dfComm = dfComm.Define("jetID_mask", 'cleaningJetIDMask(Jet_eta, Jet_chHEF, Jet_neHEF, Jet_chEmEF, Jet_neEmEF, Jet_muEF, Jet_chMultiplicity, Jet_neMultiplicity, "{0}")'.format(year)) #  for nanov15 use the JetID
    else:
        dfComm = dfComm.Define("jetID_mask", "cleaningJetSelMask(0, Jet_eta, Jet_neHEF, Jet_chEmEF, Jet_muEF, Jet_neEmEF, Jet_jetId)") # redo JetID for nanoV12,

    df = (dfComm.Define("goodMuons","{}".format(GOODMUON)+" and Muon_mediumId and Muon_pfRelIso04_all < 0.25") # add ID and ISO
          .Filter("Sum(goodMuons)>=1 and Sum(Muon_charge[goodMuons])==0","at least two good muons OS") # fix the charge for WH and ZH          
          .Define("looseMu","{}".format(LOOSEmuons))
          .Define("jetMuon_mask", "cleaningMask(Muon_jetIdx[goodMuons],nJet)")
          # Jet_puIdDisc only for nanov15
          .Define("jetVeto_mask",'cleaningJetVetoMapMask(Jet_eta,Jet_phi,"{0}")'.format(year))
          .Filter("Sum(looseMu)==2", "no extra muons")
          #
          .Define("looseEle","{}".format(LOOSEelectrons))
          .Filter("Sum(looseEle)==0", "no extra electrons")
          .Define("BJETS","{}".format(BJETS))
          .Filter("Sum(BJETS)==0","no bjets")
          .Define("Muon1_pt","Muon_pt[goodMuons][0]")
          .Define("Muon2_pt","Muon_pt[goodMuons][1]")
          ## end preselection
          .Define("HiggsCandMass","Minv(Muon_pt[goodMuons], Muon_eta[goodMuons], Muon_phi[goodMuons], Muon_mass[goodMuons])")
          .Define("goodFRSphoton", "fsrMask(Muon_fsrPhotonIdx[goodMuons],nFsrPhoton)")
          .Define("FsrPH_pt","FsrPhoton_pt[goodFRSphoton]")
          .Define("FsrPH_eta","FsrPhoton_eta[goodFRSphoton]")
          .Define("FsrPH_pt_ratio0","FsrPhoton_pt[goodFRSphoton][0]/Muon1_pt")
          .Define("FsrPH_pt_ratio1","FsrPhoton_pt[goodFRSphoton][1]/Muon2_pt")
          .Define("HiggsCandCorrMass","MinvCorr(Muon_pt[goodMuons], Muon_eta[goodMuons], Muon_phi[goodMuons], Muon_mass[goodMuons], FsrPhoton_pt[goodFRSphoton], FsrPhoton_eta[goodFRSphoton],FsrPhoton_phi[goodFRSphoton],0)")
          .Define("HiggsCandCorrPt","MinvCorr(Muon_pt[goodMuons], Muon_eta[goodMuons], Muon_phi[goodMuons], Muon_mass[goodMuons], FsrPhoton_pt[goodFRSphoton], FsrPhoton_eta[goodFRSphoton],FsrPhoton_phi[goodFRSphoton],1)")
          .Define("HiggsCandCorrRapidity","MinvCorr(Muon_pt[goodMuons], Muon_eta[goodMuons], Muon_phi[goodMuons], Muon_mass[goodMuons], FsrPhoton_pt[goodFRSphoton], FsrPhoton_eta[goodFRSphoton],FsrPhoton_phi[goodFRSphoton],2)")
          .Filter("HiggsCandCorrMass>50 and HiggsCandCorrMass<200","HiggsMass within reasonable range 50-200")
          )

    if (isData == "false"):
        df = dfwithSYST(df,year)
    else:
        df = df.Define("w_allSF", "w")
    
    if doVBF:
    
        print('GOODJETS',GOODJETSALL)
        df= (df.Define("goodJetsAll","{}".format(GOODJETSALL))
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
	     .Define("RPt","getRpt(Muon1_pt,Muon_eta[goodMuons][0],Muon_phi[goodMuons][0],Muon_mass[goodMuons][0],Muon2_pt,Muon_eta[goodMuons][1],Muon_phi[goodMuons][1],Muon_mass[goodMuons][1],jetVBF1_Pt, jetVBF1_Eta, jetVBF1_Phi, jetVBF1_Mass, jetVBF2_Pt, jetVBF2_Eta, jetVBF2_Phi, jetVBF2_Mass)")
             .Define("Mjj","Pair12Minv(jetVBF1_Pt, jetVBF1_Eta, jetVBF1_Phi, jetVBF1_Mass, jetVBF2_Pt, jetVBF2_Eta, jetVBF2_Phi, jetVBF2_Mass)")
             .Define("dEtaJJ","abs(jetVBF1_Eta-jetVBF2_Eta)")
             .Define("ZepVar","getZep(Muon1_pt,Muon_eta[goodMuons][0],Muon_phi[goodMuons][0],Muon_mass[goodMuons][0],Muon2_pt,Muon_eta[goodMuons][1],Muon_phi[goodMuons][1],Muon_mass[goodMuons][1], jetVBF1_Eta, jetVBF2_Eta)")
             .Define("minDetaDiMuVBF","minDeta(MinvCorr(Muon_pt[goodMuons], Muon_eta[goodMuons], Muon_phi[goodMuons], Muon_mass[goodMuons], FsrPhoton_pt[goodFRSphoton], FsrPhoton_eta[goodFRSphoton],FsrPhoton_phi[goodFRSphoton],3), jetVBF1_Eta, jetVBF2_Eta)")             
              )

        if mc>0:
            df= (df.Define("jetVBF1_partonFlavour","Jet_partonFlavour[goodJetsAll][index_VBFJets[0]]")
                 .Define("jetVBF2_partonFlavour","Jet_partonFlavour[goodJetsAll][index_VBFJets[1]]")
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
        }

        histsVBF = {
            "Mjj": {"name":"Mjj","title":"M(VBFjet,VBFjet); M(VBFjet,VBFjet) (GeV);N_{Events}","bin":250,"xmin":0.,"xmax":1000},
            "RPt": {"name":"RPt","title":"RPt; RPt;N_{Events}","bin":100,"xmin":0.,"xmax":1},
            "jetVBF1_Pt":  {"name":"jetVBF1_Pt","title":"jetVBF1_Pt; p_{T}^{VBFjet} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":500.},
            "jetVBF2_Pt":  {"name":"jetVBF2_Pt","title":"jetVBF2_Pt; p_{T}^{VBFjet} (GeV);N_{Events}","bin":100,"xmin":0.,"xmax":500.},
            "ZepVar":  {"name":"ZepVar","title":"ZepVar; ZepVar;N_{Events}","bin":100,"xmin":0.,"xmax":2.},
            ##
        }

        if doVBF: hists.update(histsVBF)

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
                if mc==10: canv.SaveAs("~/public_html/HMUMU_July/"+h.GetName()+"_VBF.png")
                if mc==11: canv.SaveAs("~/public_html/HMUMU_July/"+h.GetName()+"_ggH.png")
                if mc==12: canv.SaveAs("~/public_html/HMUMU_July/"+h.GetName()+"_Wm.png")
                if mc==13: canv.SaveAs("~/public_html/HMUMU_July/"+h.GetName()+"_Wp.png")
                if mc==14: canv.SaveAs("~/public_html/HMUMU_July/"+h.GetName()+"_ZH.png")
                if mc==15: canv.SaveAs("~/public_html/HMUMU_July/"+h.GetName()+"_TTH.png")
                if mc==100: canv.SaveAs("~/public_html/HMUMU_July/"+h.GetName()+"_DY.png")
                if mc==20: canv.SaveAs("~/public_html/HMUMU_July/"+h.GetName()+"_Zg_VBF.png")
                if mc==21: canv.SaveAs("~/public_html/HMUMU_July/"+h.GetName()+"_Zg_ggH.png")


        myDir='/work/submit/mariadlf/HmumuRun3/ROOTFILES/'
        if doVBF: outputFileHisto = myDir+"histoOUTname_"+str(mc)+"_"+str(year)+"_VBFcat.root"
        else : outputFileHisto = myDir+"histoOUTname_"+str(mc)+"_"+str(year)+"_ggHcat.root"

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
        for branchName in [
                "mc",
                "w",
                "w_allSF",
                "lumiIntegrated",
                "PV_npvsGood",
                "triggerAna",
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
                "FsrPH_eta"
        ]:
            branchList.push_back(branchName)

        if doVBF:
            for branchName in [
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
                    "minDetaDiMuVBF"
            ]:
                branchList.push_back(branchName)

        myDir='/work/submit/mariadlf/HmumuRun3/ROOTFILES/'
        if doVBF: outputFile = myDir+"snapshot_mc"+str(mc)+"_"+str(year)+"_VBFcat.root"
        else : outputFile = myDir+"snapshot_mc"+str(mc)+"_"+str(year)+"_ggHcat.root"

        print(outputFile)
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
    if year=="12022" or year=="22022" or year=="12023" or year=="22023": mc.extend([20,21,22,23,24,25])
    if year=="2024": mc.extend([20,21,22,23,24,26])
    if year=="2024": mc.extend([103,104])
    else: mc.extend([100])
    mc.extend([102])
    mc.extend([201,202,203,204,205,211,212,213,214])

    print(mc)

    for sampleNOW in mc:
        print('mc=',mc)
        files = SwitchSample(thisdict,sampleNOW)[0]
        print('outside the function: ', len(files))
        rdf = ROOT.RDataFrame("Runs", files) # make sure this is not the distributed

        sumW = computeWeigths(rdf,SwitchSample(thisdict,sampleNOW)[1])
        analysis(files,year,sampleNOW,sumW)

    data = []
    if year=="12022": data = [-11,-12,-13,-14]
    if year=="22022": data = [-15,-16,-17]
    if year=="12023": data = [-23,-24]
    if year=="22023": data = [-31,-32]
    if year=="2024": data = [-41,-42,-43,-44,-45,-46,-47,-48,-49,-50,-51,-52,-53,-54]

    readDataQuality(year)

    for sampleNOW in data:
        files = SwitchSample(thisdict,sampleNOW)[0]
        print('outside the function: ', len(files))
        analysis(files,year,sampleNOW,1.)



if __name__ == "__main__":

    now = datetime.now()
    print('==> very beginning: ',now)

    print('year=',sys.argv[1])
    year=sys.argv[1]

    loopOnDataset(year)

