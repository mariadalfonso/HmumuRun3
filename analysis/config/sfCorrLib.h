//#include "mysf.h"
//#include <ostream>
//#include <iostream>

#ifndef sfDistr_h
#define sfDistr_h

//#include "correction.h"
#include <stdio.h>
#include <string.h>
#include <ostream>
#include <iostream>

class MyCorrections {

public:
  MyCorrections(int year);

  double eval_jetCORR   (double area, double eta, double phi, double pt, double rho, bool isData, int run, std::string year, std::string mc);
  double eval_fatJetCORR(double area, double eta, double phi, double pt, double rho, bool isData, int run, std::string year, std::string mc);
  double eval_jesUnc    (double eta, double pt, int type);
  double eval_jer       (double pt, double eta, double rho, double area);
  double eval_jetID     (float eta, float chHEF, float neHEF, float chEmEF, float neEmEF, float muEF, int chMultiplicity, int neMultiplicity);
  double eval_jetVeto   (double eta, double phi);
  double eval_puSF      (double NumTrueInteractions, std::string weights);
  double eval_photonSF  (std::string year, std::string valType, std::string workingPoint, double eta, double pt);
  double eval_photonPixVetoSF  (std::string year, std::string valType, std::string workingPoint, double eta, double pt);
  double eval_electronSF(std::string year, std::string valType, std::string workingPoint, double eta, double pt, double minVal);
  double eval_electronScaleData(std::string year, const char *valType, int run, const double eta, const double r9, double pt, const double gain);
  double eval_electronSmearingSystMC(const char *valType, double pt, const double eta, const double r9);
  double eval_muonTRKSF (std::string year, std::string valType, double eta, double pt);
  double eval_muonIDSF  (std::string year, std::string valType, double eta, double pt, std::string workingPoint);
  double eval_muonISOSF (std::string year, std::string valType, double eta, double pt, std::string workingPoint);
  
private:
  correction::Correction::Ref puSF_;
  correction::Correction::Ref photonSF_;
  correction::Correction::Ref photonPixVetoSF_;
  correction::Correction::Ref electronIDSF_;
  correction::CompoundCorrection::Ref electronScaleData_;
  correction::Correction::Ref electronSmearAndSystMC_;

  correction::Correction::Ref muonTRKSF_;
  correction::Correction::Ref muonIDMSF_;
  correction::Correction::Ref muonISOLSF_;
  correction::Correction::Ref muonIDTSF_;
  correction::Correction::Ref muonISOTSF_;
  correction::Correction::Ref muonTIso_;
  correction::Correction::Ref muonLIso_;
  correction::Correction::Ref muonTrgSF_;

  correction::Correction::Ref JER_;
  correction::Correction::Ref JERsf_;
  correction::CompoundCorrection::Ref JEC_;
  correction::CompoundCorrection::Ref JECdata_;
  correction::CompoundCorrection::Ref JECdata22E_;
  correction::CompoundCorrection::Ref JECdata22F_;
  correction::CompoundCorrection::Ref JECdata22G_;
  correction::CompoundCorrection::Ref fatJEC_;
  correction::CompoundCorrection::Ref fatJECdata_;
  correction::CompoundCorrection::Ref fatJECdata22E_;
  correction::CompoundCorrection::Ref fatJECdata22F_;
  correction::CompoundCorrection::Ref fatJECdata22G_;
  correction::Correction::Ref jesUnc_;
  correction::Correction::Ref fatjesUnc_;
  correction::Correction::Ref vetoMaps_;
  correction::Correction::Ref jetTightID_;

};

MyCorrections::MyCorrections(int year) {

  std::string dirName = "/cvmfs/cms-griddata.cern.ch/cat/metadata/";

  std::string subDirName = "";
  std::string dataName = "";

  if(year == 2018)  { subDirName += "Run2-2018-UL-NanoAODv9/latest/"; dataName = "UL18"; }
  if(year == 2017)  { subDirName += "Run2-2017-UL-NanoAODv9/latest/"; dataName = "UL17"; }
  if(year == 22016)  { subDirName += "Run2-2016postVFP-UL-NanoAODv9/latest/"; dataName = "UL16"; }
  if(year == 12016)  { subDirName += "Run2-2016preVFP-UL-NanoAODv9/latest/"; dataName = "UL16APV"; }
  if(year == 12022)  { subDirName += "Run3-22CDSep23-Summer22-NanoAODv12/latest/"; }
  if(year == 22022)  { subDirName += "Run3-22EFGSep23-Summer22EE-NanoAODv12/latest/"; }
  if(year == 12023)  { subDirName += "Run3-23CSep23-Summer23-NanoAODv12/latest/"; }
  if(year == 22023)  { subDirName += "Run3-23DSep23-Summer23BPix-NanoAODv12/latest/"; }
  if(year == 2024)   { subDirName += "Run3-24CDEReprocessingFGHIPrompt-Summer24-NanoAODv15/latest/"; }

  /*
  if(year == 2018)  { subDirName += "2018_UL/"; dataName = "UL18"; }
  if(year == 2017)  { subDirName += "2017_UL/"; dataName = "UL17"; }
  if(year == 22016)  { subDirName += "2016postVFP_UL/"; dataName = "UL16"; }
  if(year == 12016)  { subDirName += "2016preVFP_UL/"; dataName = "UL16APV"; }
  if(year == 12022)  { subDirName += "2022_Summer22/"; }
  if(year == 22022)  { subDirName += "2022_Summer22EE/"; }
  if(year == 12023)  { subDirName += "2023_Summer23/"; }
  if(year == 22023)  { subDirName += "2023_Summer23BPix/"; }
  if(year == 2024)   { subDirName += "2024_Summer24/"; }
  */

  std::string fileNameLUM = dirName+"LUM/"+subDirName+"puWeights.json.gz";
  if(year == 2024) fileNameLUM = dirName+"LUM/"+subDirName+"puWeights_BCDEFGHI.json.gz";

  std::string corrNameLUM = "";  
  if(year == 2018) corrNameLUM = "Collisions18_UltraLegacy_goldenJSON";
  if(year == 2017) corrNameLUM = "Collisions17_UltraLegacy_goldenJSON";
  if(year == 22016 or year == 12016) corrNameLUM = "Collisions16_UltraLegacy_goldenJSON";
  if(year == 12022) corrNameLUM = "Collisions2022_355100_357900_eraBCD_GoldenJson";
  if(year == 22022) corrNameLUM = "Collisions2022_359022_362760_eraEFG_GoldenJson";
  if(year == 12023) corrNameLUM = "Collisions2023_366403_369802_eraBC_GoldenJson";
  if(year == 22023) corrNameLUM = "Collisions2023_369803_370790_eraD_GoldenJson";
  if(year == 2024) corrNameLUM = "Collisions24_BCDEFGHI_goldenJSON";

  auto csetPU = correction::CorrectionSet::from_file(fileNameLUM);
  puSF_ = csetPU->at(corrNameLUM);

  std::string fileNameMU = dirName+"MUO/"+subDirName+"muon_Z.json.gz";

  //// muon WP https://muon-wiki.docs.cern.ch/guidelines/recommendations/#particle-flow-isolation

  if(year == 2018 or year == 2017 or year == 12016 or year == 22016) {
    auto csetMu = correction::CorrectionSet::from_file(fileNameMU);

    muonTRKSF_ = csetMu->at("NUM_TrackerMuons_DEN_genTracks");
    muonIDTSF_ = csetMu->at("NUM_TightID_DEN_genTracks");
    muonIDMSF_ = csetMu->at("NUM_MediumID_DEN_genTracks");
    muonISOTSF_ = csetMu->at("NUM_TightRelIso_DEN_TightIDandIPCut");
    muonISOLSF_ = csetMu->at("NUM_LooseRelIso_DEN_MediumID");
  }

  if(year == 2024 or year == 22023 or year == 12023 or year == 22022 or year == 12022) {
    auto csetMu = correction::CorrectionSet::from_file(fileNameMU);

    muonIDMSF_ = csetMu->at("NUM_MediumID_DEN_TrackerMuons");
    muonTIso_ = csetMu->at("NUM_TightPFIso_DEN_MediumID"); //NUM_TightPFIso_DEN_MediumID
    muonLIso_ = csetMu->at("NUM_LoosePFIso_DEN_MediumID"); //NUM_LoosePFIso_DEN_MediumID
    if(year == 22023 or year == 12023 or year == 22022 or year == 12022) {
      // 2024 missing
      muonTrgSF_ = csetMu->at("NUM_IsoMu24_DEN_CutBasedIdMedium_and_PFIsoMedium");
    }
  }

  //////////////
  ////  JME
  //////////////

  const std::string fileNameJEC = dirName+"JME/"+subDirName+"jet_jerc.json.gz";
  auto csetJEC = correction::CorrectionSet::from_file(fileNameJEC);

  if(year == 12022 or year == 22022 or year == 12023 or year == 22023 or year == 2024) {
    const std::string jetType="AK4PFPuppi";
    const std::string suffix = "_DATA_L1L2L3Res_";

    std::string tagName = "";
    if(year == 12022)  tagName = "Summer22_22Sep2023_RunCD_V3";
    if(year == 12023)  tagName = "Summer23Prompt23_V2";
    if(year == 22023)  tagName = "Summer23BPixPrompt23_V3";
    if(year == 2024)   tagName = "Summer24Prompt24_V2";
    if(year == 22022)  {
      tagName = "Summer22EE_22Sep2023_RunF_V3";
      JECdata22E_ = csetJEC->compound().at("Summer22EE_22Sep2023_RunE_V3"+suffix+jetType);
      JECdata22F_ = csetJEC->compound().at("Summer22EE_22Sep2023_RunF_V3"+suffix+jetType);
      JECdata22G_ = csetJEC->compound().at("Summer22EE_22Sep2023_RunG_V3"+suffix+jetType);
    }
    JECdata_ = csetJEC->compound().at(tagName+suffix+jetType);

    std::string tagNameMC = "";
    if(year == 12022)  tagNameMC = "Summer22_22Sep2023_V3";
    if(year == 22022)  tagNameMC = "Summer22EE_22Sep2023_V3";
    if(year == 12023)  tagNameMC = "Summer23Prompt23_V2";
    if(year == 22023)  tagNameMC = "Summer23BPixPrompt23_V3";
    if(year == 2024)   tagNameMC = "Summer24Prompt24_V2";
    JEC_ = csetJEC->compound().at(tagNameMC+"_MC_L1L2L3Res_"+jetType);
    jesUnc_ = csetJEC->at(tagNameMC+"_MC_Total_"+jetType);

  }

  // 2024 not yet there
  if(year == 12022 or year == 22022 or year == 12023 or year == 22023 or year == 2024) {

    const std::string fileNameFatJEC = dirName+"JME/"+subDirName+"fatJet_jerc.json.gz";
    auto csetFatJEC = correction::CorrectionSet::from_file(fileNameFatJEC);

    const std::string jetType="AK8PFPuppi";
    const std::string suffix = "_DATA_L1L2L3Res_";

    std::string tagName = "";
    if(year == 12022)  tagName = "Summer22_22Sep2023_RunCD_V3";
    if(year == 12023)  tagName = "Summer23Prompt23_V2";
    if(year == 22023)  tagName = "Summer23BPixPrompt23_V3";
    if(year == 2024)   tagName = "Summer24Prompt24_V2";
    if(year == 22022)  {
      tagName = "Summer22EE_22Sep2023_RunF_V3";
      fatJECdata22E_ = csetFatJEC->compound().at("Summer22EE_22Sep2023_RunE_V3"+suffix+jetType);
      fatJECdata22F_ = csetFatJEC->compound().at("Summer22EE_22Sep2023_RunF_V3"+suffix+jetType);
      fatJECdata22G_ = csetFatJEC->compound().at("Summer22EE_22Sep2023_RunG_V3"+suffix+jetType);
    }
    fatJECdata_ = csetFatJEC->compound().at(tagName+suffix+jetType);

    std::string tagNameMC = "";
    if(year == 12022)  tagNameMC = "Summer22_22Sep2023_V3";
    if(year == 22022)  tagNameMC = "Summer22EE_22Sep2023_V3";
    if(year == 12023)  tagNameMC = "Summer23Prompt23_V2";
    if(year == 22023)  tagNameMC = "Summer23BPixPrompt23_V3";
    if(year == 2024)   tagNameMC = "Summer24Prompt24_V2";
    fatJEC_ = csetFatJEC->compound().at(tagNameMC+"_MC_L1L2L3Res_"+jetType);
    fatjesUnc_ = csetFatJEC->at(tagNameMC+"_MC_Total_"+jetType);

  }

  // jetID
  std::string fileNameJetID = dirName+"JME/"+subDirName+"jetid.json.gz";
  auto csetJetID = correction::CorrectionSet::from_file(fileNameJetID);
  const std::string tagNameJetID = "AK4PUPPI_Tight";
  //  std::string tagNameJetID = "AK4PUPPI_TightLeptonVeto";
  jetTightID_           = csetJetID->at(tagNameJetID);

  // veto Map the jet
  std::string fileNameJetVeto = dirName+"JME/"+subDirName+"jetvetomaps.json.gz";
  auto csetVeto = correction::CorrectionSet::from_file(fileNameJetVeto);
  std::string tagNameVeto = "";
  if(year == 22016 or year == 2017 or year == 2018) tagNameVeto = "Summer19"+dataName+"_V1";
  if(year == 12016) tagNameVeto = "Summer19UL16_V1";
  if(year == 12022) tagNameVeto = "Summer22_23Sep2023_RunCD_V1";
  if(year == 22022) tagNameVeto = "Summer22EE_23Sep2023_RunEFG_V1";
  if(year == 12023) tagNameVeto = "Summer23Prompt23_RunC_V1";
  if(year == 22023) tagNameVeto = "Summer23BPixPrompt23_RunD_V1";
  if(year == 2024) tagNameVeto = "Summer24Prompt24_RunBCDEFGHI_V1";
  vetoMaps_ = csetVeto->at(tagNameVeto);

  //////////////
  ////  EGM https://twiki.cern.ch/twiki/bin/view/CMS/EgammSFandSSRun3
  //////////////

  std::string fileNameIDELE   = dirName+"EGM/"+subDirName+"electron.json.gz";
  auto csetIDELE = correction::CorrectionSet::from_file(fileNameIDELE);

  const std::string tagNameEleID = "Electron-ID-SF"; // both reco and ID
  csetIDELE->at(tagNameEleID);

  std::string fileNameEnergyEtDependentELE = dirName+"EGM/"+subDirName+"electronSS_EtDependent.json.gz";
  auto csetEnergyEtDependentELE = correction::CorrectionSet::from_file(fileNameEnergyEtDependentELE);

  std::string suffix = "";
  if     (year == 12022) suffix = "2022preEE";
  else if(year == 22022) suffix = "2022postEE";
  else if(year == 12023) suffix = "2023preBPIX";
  else if(year == 22023) suffix = "2023postBPIX";
  else if(year == 2024) suffix = "2024";

  electronScaleData_ = csetEnergyEtDependentELE->compound().at("Scale");
  electronSmearAndSystMC_ = csetEnergyEtDependentELE->at("SmearAndSyst");

};

double MyCorrections::eval_jetCORR(double area, double eta, double phi, double pt, double rho, bool isData, int run, std::string year, std::string mc) {

  if (year == "2024" or year == "22023") {
    if(isData) return JECdata_->evaluate({area, eta,  pt, rho, phi, (float) run});
    else JEC_->evaluate({area, eta, pt, rho, phi});
  } else if (year == "12023") {
    if(isData) return JECdata_->evaluate({area, eta, pt, rho, (float) run});
    else return JEC_->evaluate({area, eta, pt, rho});
  } else if (year == "22022") {
    if(isData and mc == "-15") return JECdata22E_->evaluate({area, eta, pt, rho});
    if(isData and mc == "-16") return JECdata22F_->evaluate({area, eta, pt, rho});
    if(isData and mc == "-17") return JECdata22G_->evaluate({area, eta, pt, rho});
    else return JEC_->evaluate({area, eta, pt, rho});
  } else {
    if(isData) return JECdata_->evaluate({area, eta, pt, rho});
    else return JEC_->evaluate({area, eta, pt, rho});
  }

  return 1.0;

};

double MyCorrections::eval_fatJetCORR(double area, double eta, double phi, double pt, double rho, bool isData, int run, std::string year, std::string mc) {

  if (year == "2024" or year == "22023") {
    if(isData) return fatJECdata_->evaluate({area, eta,  pt, rho, phi, (float) run});
    else fatJEC_->evaluate({area, eta, pt, rho, phi});
  } else if (year == "12023") {
    if(isData) return fatJECdata_->evaluate({area, eta, pt, rho, (float) run});
    else return fatJEC_->evaluate({area, eta, pt, rho});
  } else if (year == "22022") {
    if(isData and mc == "-15") return fatJECdata22E_->evaluate({area, eta, pt, rho});
    if(isData and mc == "-16") return fatJECdata22F_->evaluate({area, eta, pt, rho});
    if(isData and mc == "-17") return fatJECdata22G_->evaluate({area, eta, pt, rho});
    else return fatJEC_->evaluate({area, eta, pt, rho});
  } else {
    if(isData) return fatJECdata_->evaluate({area, eta, pt, rho});
    else return fatJEC_->evaluate({area, eta, pt, rho});
  }

  return 1.0;

};

double MyCorrections::eval_jesUnc(double eta, double pt, int type) {
  if(type == 0) return jesUnc_->evaluate({eta, pt});
  return 0.0;
};

double MyCorrections::eval_jer(double double1, double double2, double double3, double double4) {
  return JER_->evaluate({double1, double2, double3, double4});
};

double MyCorrections::eval_jetID(float eta, float chHEF, float neHEF, float chEmEF, float neEmEF, float muEF, int chMultiplicity, int neMultiplicity) {

  eta = fabs(eta);
  int multiplicity = chMultiplicity + neMultiplicity;
  //"chEmEF" and "muEF" unused in 2024 but still needed
  return jetTightID_->evaluate({eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity});

};

double MyCorrections::eval_jetVeto(double eta, double phi) {
  std::string typeMaps = "jetvetomap";
  return vetoMaps_->evaluate({typeMaps,eta, phi});
};

double MyCorrections::eval_electronSF(std::string year, std::string valType,  std::string workingPoint, double eta, double pt, double minVal) {
  //  pt = std::max(pt,10.001); for wp80 is 10 while for the above10 is 20
  pt = std::max(pt,minVal);
  return electronIDSF_->evaluate({year, valType, workingPoint, eta, pt});
};

double MyCorrections::eval_electronScaleData(std::string year, const char *valType, int run, const double eta, const double r9, double pt, const double gain) {
  // "scale"
  return electronScaleData_->evaluate({valType, (float) run, eta, r9, pt, gain});
};

double MyCorrections::eval_electronSmearingSystMC(const char *valType, double pt, const double eta, const double r9) {
  // "smear"
  pt = std::max(pt,15.000);
  return electronSmearAndSystMC_->evaluate({valType,pt,r9,abs(eta)});
};

double MyCorrections::eval_photonSF(std::string year, std::string valType,  std::string workingPoint, double eta, double pt) {
  pt = std::max(pt,20.001);
  return photonSF_->evaluate({year, valType, workingPoint, eta, pt});
};

double MyCorrections::eval_photonPixVetoSF(std::string year, std::string valType,  std::string workingPoint, double eta, double pt) {
  pt = std::max(pt,20.001);
  return photonPixVetoSF_->evaluate({year, valType, workingPoint, eta, pt});
};

double MyCorrections::eval_muonTRKSF(std::string year, std::string valType, double eta, double pt) {
  eta = std::min(std::abs(eta),2.399);
  pt = std::max(pt,20.001);
  return muonTRKSF_->evaluate({year, eta, pt, valType});
};

double MyCorrections::eval_muonIDSF(std::string year, std::string valType, double eta, double pt, std::string workingPoint) {
  if (year == "22023" or year == "2024") eta = eta;
  else eta = std::min(std::abs(eta),2.399);

  pt = std::max(pt,15.001);

  if (workingPoint=="T") {
    return muonIDTSF_->evaluate({eta, pt, valType});
  } else if (workingPoint=="M") {
    return muonIDMSF_->evaluate({eta, pt, valType});
  }
  return 1.;
};

double MyCorrections::eval_muonISOSF(std::string year, std::string valType, double eta, double pt, std::string workingPoint) {

  if (year == "22023" or year == "2024") eta = eta;
  else eta = std::min(std::abs(eta),2.399);

  pt = std::max(pt,15.001);

  if (workingPoint=="T") {
    if (year == "12022" or year == "22022" or year == "12023" or year == "22023" or year == "2024") return muonTIso_->evaluate({ eta, pt, valType});
    return muonISOTSF_->evaluate({year, eta, pt, valType});
  } else if (workingPoint=="L") {
    if (year == "12022" or year == "22022" or year == "12023" or year == "22023" or year == "2024") return muonLIso_->evaluate({eta, pt, valType});
    return muonISOLSF_->evaluate({year, eta, pt, valType});
  }
  return 1.;
};
  
double MyCorrections::eval_puSF(double int1, std::string str1) {
  return puSF_->evaluate({int1, str1});
};

#endif


