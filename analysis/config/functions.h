#ifndef FUNCTIONS_H
#define FUNCTIONS_H

using Vec_b = ROOT::VecOps::RVec<bool>;
using Vec_d = ROOT::VecOps::RVec<double>;
using Vec_f = ROOT::VecOps::RVec<float>;
using Vec_i = ROOT::VecOps::RVec<int>;
using Vec_ui = ROOT::VecOps::RVec<unsigned int>;

using stdVec_i = std::vector<int>;
using stdVec_b = std::vector<bool>;
using stdVec_f = std::vector<float>;


typedef ROOT::Math::LorentzVector<ROOT::Math::PtEtaPhiM4D<float> > PtEtaPhiMVector;

std::unordered_map< UInt_t, std::vector< std::pair<UInt_t,UInt_t> > > jsonMap;

bool isGoodRunLS(const bool isData, const UInt_t run, const UInt_t lumi) {

  if(not isData) return true;
  if(jsonMap.find(run) == jsonMap.end()) return false; // run not found

  auto& validlumis = jsonMap.at(run);
  auto match = std::lower_bound(std::begin(validlumis), std::end(validlumis), lumi,
                                [](std::pair<unsigned int, unsigned int>& range, unsigned int val) { return range.second < val; });
  return match->first <= lumi && match->second >= lumi;
}

Vec_i indices(const int& size, const int& start = 0) {
    Vec_i res(size, 0);
    std::iota(std::begin(res), std::end(res), start);
    return res;
}

Vec_f NomUpDownVar(const float nom, const float up, const float down, float weight) {

  Vec_f res(3, 1);
  res[0] = weight;  // nom - already mutliplied for the Nom
  res[1] = (nom!=0) ? weight*up/nom : weight;  // up
  res[2] = (nom!=0) ? weight*down/nom : weight;  // down

  return res;
}

Vec_b cleaningMask(Vec_i indices, int size) {

  Vec_b mask(size, true);
  for (int idx : indices) {
    if(idx < 0) continue;
    mask[idx] = false;
  }
  return mask;
}

Vec_b fsrMask(Vec_i indices, int size) {

  Vec_b mask(size, false);
  for (int idx : indices) {
    if(idx < 0) continue;
    mask[idx] = true;
  }
  return mask;
}

TLorentzVector MakeTLV(const float pt, const float eta, const float phi, const float mass) {
    TLorentzVector v;
    v.SetPtEtaPhiM(pt, eta, phi, mass);
    return v;
}

// https://twiki.cern.ch/twiki/bin/view/CMS/JetID13p6TeV
// bug fix for v12 nano
// going to use the json for v15

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

float Minv(const Vec_f& pt, const Vec_f& eta, const Vec_f& phi, const Vec_f& m) {
  PtEtaPhiMVector p1(pt[0], eta[0], phi[0], m[0]);
  PtEtaPhiMVector p2(pt[1], eta[1], phi[1], m[1]);
  return (p1 + p2).M();
}


float minDeta(const float& etaDiMu, const float& jetEta1, const float& jetEta2) {

  return   std::min(abs(etaDiMu-jetEta1),abs(etaDiMu-jetEta2));
  
}  

float MinvCorr(const TLorentzVector& p1, const TLorentzVector& p2,
               const Vec_f& fsr_pt, const Vec_f& fsr_eta, const Vec_f& fsr_phi,
               const int& typeVar) {

  // Index of the lowest-dR/ET2 among associated FSR photons
  TLorentzVector p4 = p1 + p2;

  for (int i = 0; i < 2; i++){

    if (i==0 and fsr_pt[i]/p1.Pt() > 0.4) continue;
    if (i==1 and fsr_pt[i]/p2.Pt() > 0.4) continue;

    TLorentzVector ph;
    ph.SetPtEtaPhiM(fsr_pt[i], fsr_eta[i], fsr_phi[i], 0);

    p4 = p4 + ph;

  }

  if (typeVar==0) { return p4.M(); }
  if (typeVar==1) { return p4.Pt(); }
  if (typeVar==2) { return p4.Rapidity(); }
  if (typeVar==3) { return p4.Eta(); }

  return 0;

}

stdVec_i getVBFIndicies(const Vec_f& pts, const Vec_f& etas, const Vec_f& phis, const Vec_f& masses){

  stdVec_i idx_(2, -1);

  int n = etas.size(), index0 = -1, index1 = -1;
  float pt_limit = 20;

  float max = 0;

  for (int i = 0; i < n; i++){
    if (pts[i] < pt_limit) continue;

    for (int j = i+1; j < n; j++){
      if (pts[j] < pt_limit) continue;

      if(std::max(pts[i],pts[j]) < 35 ) continue;
      if(etas[i]*etas[j] > 0) continue;

      PtEtaPhiMVector pi(pts[i], etas[i], phis[i], masses[i]);
      PtEtaPhiMVector pj(pts[j], etas[j], phis[j], masses[j]);
      float MJJ = (pi+pj).M();

      //      if (max < abs(etas[i] - etas[j]) && etas[i]*etas[j] < 0) {
      if (max < MJJ ) {
        max = MJJ;
        index0 = i;
        index1 = j;
      }
    }
  }

  idx_[0] = index0;
  idx_[1] = index1;

  return idx_;

}

float getZep(const TLorentzVector& mu1, const TLorentzVector& mu2,
             const float& VBF1_eta, const float& VBF2_eta
             ){
  TLorentzVector HiggsCand = mu1 + mu2;

  float ZepVar = ( HiggsCand.Rapidity() - (VBF1_eta+VBF2_eta)/2.0)/abs(VBF1_eta - VBF2_eta);

  return ZepVar;

}

float getRpt(const TLorentzVector& mu1, const TLorentzVector& mu2,
	     const TLorentzVector& VBF1Cand, const TLorentzVector& VBF2Cand ){

  TLorentzVector HiggsCand = mu1 + mu2;
  TLorentzVector All_Vec = VBF1Cand + VBF2Cand + HiggsCand;

  float Rpt = abs(All_Vec.Pt())/( abs(VBF1Cand.Pt()) + abs(VBF2Cand.Pt()) + abs(HiggsCand.Pt()) );

  return Rpt;

}

float Pair12Minv(const float& v1_Pt, const float& v1_Eta, const float& v1_Phi, const float& v1_Mass,
                 const float& v2_Pt, const float& v2_Eta, const float& v2_Phi, const float& v2_Mass) {
  PtEtaPhiMVector p1(v1_Pt, v1_Eta, v1_Phi, v1_Mass);
  PtEtaPhiMVector p2(v2_Pt, v2_Eta, v2_Phi, v2_Mass);
  return (p1 + p2).M();

}

#endif

