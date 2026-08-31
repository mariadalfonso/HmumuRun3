#ifndef MYFUN                                                                                                                                                                       
#define MYFUN                                                                                                                                                                       

using Vec_b = ROOT::VecOps::RVec<bool>;
using Vec_f = ROOT::VecOps::RVec<float>;
using Vec_i = ROOT::VecOps::RVec<int>;
using Vec_ui = ROOT::VecOps::RVec<unsigned int>;

/*
Vec_f computeMUcorrection(Vec_f mu_pt, Vec_f mu_eta, Vec_f mu_phi, Vec_i mu_charge, Vec_ui mu_nTrackerLayers, bool isData, float event, float luminosityBlock ){

  Vec_f new_mu(mu_pt.size());

  for (unsigned int idx = 0; idx < mu_pt.size(); ++idx) {

    // Data apply scale correction
    // MC apply scale shift and resolution smearing

    if (isData) new_mu[idx] = pt_scale(isData, mu_pt[idx], mu_eta[idx], mu_phi[idx], mu_charge[idx]);
    else {
      auto pt_scale_corr = pt_scale(isData, mu_pt[idx], mu_eta[idx], mu_phi[idx], mu_charge[idx]);
      new_mu[idx] = pt_resol(pt_scale_corr, mu_eta[idx], mu_phi[idx], float(mu_nTrackerLayers[idx]), event, luminosityBlock);
    }

  }

  return new_mu;

}
*/

//on the Z
float computeDYturbo(float boson_genPt, float boson_genMass){

  float weight = 1.;

  float pt_eval = std::min(boson_genPt, 50.f);

  if ( boson_genMass > 70 && boson_genMass<=106 ) weight = (csetZ1->at("qt_norm_reweight")->evaluate({pt_eval}));
  else if ( boson_genMass > 106 && boson_genMass<160 ) weight = csetZ2->at("qt_norm_reweight")->evaluate({pt_eval});

  // protect against bad values
  if (weight > 1.5 || weight < 0.01 || weight < -998.)
    weight = 1.;

  return weight;

}

Vec_f getLeptonJetBTag(const Vec_f &Jet_btagDeepFlavB,
		       const Vec_s &Lepton_jetIdx)
{
    Vec_f out(Lepton_jetIdx.size(), 0.f);

    for (size_t i = 0; i < Lepton_jetIdx.size(); ++i) {
        if (Lepton_jetIdx[i] >= 0)
            out[i] = Jet_btagDeepFlavB[Lepton_jetIdx[i]];
    }

    return out;
}

Vec_f computeLeptonJetPtRatio(
    const Vec_f& relIso)
{
    Vec_f out(relIso.size());

    for (size_t i = 0; i < relIso.size(); ++i) {
        out[i] = std::min(1.f/(1.f + relIso[i]), 1.5f);
    }

    return out;
}

Vec_f reshapeBTVdiscr(MyCorrections corrSFs, Vec_f Jet_btagRobustParTAK4B, Vec_f Jet_partonFlavour, Vec_f Jet_eta, Vec_f Jet_pt, bool isData) {

  if (isData) return Jet_btagRobustParTAK4B; // do nothing

  Vec_f new_discr(Jet_btagRobustParTAK4B.size());

  for (unsigned int idx = 0; idx < Jet_btagRobustParTAK4B.size(); ++idx) {
    new_discr[idx] = Jet_btagRobustParTAK4B[idx] * corrSFs.eval_btv("central",abs(Jet_partonFlavour[idx])==5?5:abs(Jet_partonFlavour[idx])==4?4:0, Jet_eta[idx], Jet_pt[idx], Jet_btagRobustParTAK4B[idx], "L");
    std::cout << " new_discr[idx] = " << new_discr[idx] << " org Jet_btagRobustParTAK4B[idx] = " << Jet_btagRobustParTAK4B[idx] << std::endl;
  }

  return new_discr;

}

Vec_f computeEleSSCorrection(MyCorrections corrSFs, Vec_f ele_pt, Vec_f ele_eta, Vec_f ele_r9, Vec_f ele_gain, float event, float run, bool isData, string year){

  std::vector<double> random_numbers(ele_pt.size(), 0.0);

  // this is what's implemented: each event get the same reproducible random_numbers; each electron in the list get a different value
  // to recheck what should be
  if (!isData) {
    unsigned int seed = static_cast<unsigned int>(event + 0.5);
    std::mt19937 gen(seed);
    std::normal_distribution<> dist(0.0, 1.0);
    for (auto &x : random_numbers) x = dist(gen);
  }

  Vec_f new_ele; new_ele.resize(ele_pt.size());
  for (unsigned int idx = 0; idx < ele_pt.size(); ++idx) {
    if (isData) new_ele[idx] = ele_pt[idx] * corrSFs.eval_electronScaleData(year, "scale", run, ele_eta[idx], ele_r9[idx], ele_pt[idx], ele_gain[idx]);
    else new_ele[idx] = ele_pt[idx] * (1 + random_numbers[idx] * corrSFs.eval_electronSmearingSystMC("smear", ele_pt[idx], ele_eta[idx], ele_r9[idx]));
  }
  return new_ele;

}

Vec_f computeJECcorrection(MyCorrections corrSFs, Vec_f jet_pt, Vec_f jet_rawFactor, Vec_f jet_eta, Vec_f jet_phi, Vec_f jet_area, float rho, float run, bool isData, string year, string mc){

  Vec_f new_jet; new_jet.resize(jet_pt.size());
  Vec_f raw_jet; raw_jet.resize(jet_pt.size());
  for (unsigned int idx = 0; idx < jet_pt.size(); ++idx) {                                                                                                                            
    raw_jet[idx] = jet_pt[idx] * (1.0 - jet_rawFactor[idx]);                                                                                                                       
    new_jet[idx] = raw_jet[idx] * corrSFs.eval_jetCORR(jet_area[idx], jet_eta[idx], jet_phi[idx], raw_jet[idx], rho, isData, run, year, mc );                                      
  }                                                                                                                                                                                   
  return new_jet;                                                                                                                                                                     
}                                                                                                                                                                                   

Vec_f computeJECuncertainties(MyCorrections corrSFs, Vec_f jet_pt, Vec_f jet_eta){                                                                                                  
  Vec_f new_jet_delta; new_jet_delta.resize(jet_pt.size());
  int type = 0;                                                                                                                                                                       
  for (unsigned int idx = 0; idx < jet_pt.size(); ++idx) new_jet_delta[idx] = corrSFs.eval_jesUnc(jet_eta[idx], jet_pt[idx], type );                                                  
  return new_jet_delta;                                                                                                                                                               
}                                                                                                                                                                                   

Vec_b cleaningJetVetoMapMask(const Vec_f& jet_eta, const Vec_f& jet_phi, const string year) {                                                                                       
  Vec_b jet_vetoMap_mask(jet_eta.size(), true);                                                                                                                                       
  for (unsigned int idx = 0; idx < jet_eta.size(); ++idx) {                                                                                                                           
    double jetVetoMap = corr_sf.eval_jetVeto(jet_eta[idx], jet_phi[idx]);                                                                                                               
    if(jetVetoMap > 0) jet_vetoMap_mask[idx] = false;                                                                                                                                   
  }                                                                                                                                                                                   
  return jet_vetoMap_mask;                                                                                                                                                            
}                                                                                                                                                                                   

Vec_b cleaningJetIDMask(Vec_f jet_eta, Vec_f jet_chHEF, Vec_f jet_neHEF, Vec_f jet_chEmEF, Vec_f jet_neEmEF, Vec_f jet_muEF, Vec_f jet_chMultiplicity, Vec_f jet_neMultiplicity, string year) {                                                                                                                                                                                
  Vec_b jetID_mask(jet_eta.size(), true);                                                                                                                                             
  for (unsigned int idx = 0; idx < jet_eta.size(); ++idx) {                                                                                                                           
    double jetID = corr_sf.eval_jetID(jet_eta[idx], jet_chHEF[idx], jet_neHEF[idx], jet_chEmEF[idx], jet_neEmEF[idx], jet_muEF[idx], jet_chMultiplicity[idx], jet_neMultiplicity[idx]); 
    if (jetID < 0) jetID_mask[idx] = false;                                                                                                                                             
  }                                                                                                                                                                                   
  return jetID_mask;                                                                                                                                                                  
}

Vec_b cleaningJetSelMask(int sel, Vec_f jet_eta, Vec_f jet_neHEF, Vec_f jet_chEmEF, Vec_f jet_muEF, Vec_f jet_neEmEF, Vec_ui jet_jetId){

  Vec_b jet_Sel_mask(jet_eta.size(), true);
  bool debug = false;
  float result = 0.0;

  for(unsigned int i=0;i<jet_eta.size();i++) {

    result = 0.0;
    if(sel >= 0){
      bool jet_passJetIdTight = false;
      if      (abs(jet_eta[i]) <= 2.7) jet_passJetIdTight = jet_jetId[i] & (1 << 1);
      else if (abs(jet_eta[i]) > 2.7 && abs(jet_eta[i]) <= 3.0) jet_passJetIdTight = (jet_jetId[i] & (1 << 1)) && (jet_neHEF[i] < 0.99);
      else if (abs(jet_eta[i]) > 3.0) jet_passJetIdTight = (jet_jetId[i] & (1 << 1)) && (jet_neEmEF[i] < 0.4);

      if(sel == 0 && jet_passJetIdTight == true) result = 1.0;

      if(sel == 1) {
        bool jet_passJetIdTightLepVeto = false;
        if (abs(jet_eta[i]) <= 2.7) jet_passJetIdTightLepVeto = jet_passJetIdTight && (jet_muEF[i] < 0.8) && (jet_chEmEF[i] < 0.8);
        else jet_passJetIdTightLepVeto = jet_passJetIdTight;

        if(jet_passJetIdTightLepVeto == true) result = 1.0;
      } // sel == 1
    } // sel == 0 / 1

    if(result == 0) jet_Sel_mask[i] = false;
  }

  return jet_Sel_mask;
}


#endif                 
