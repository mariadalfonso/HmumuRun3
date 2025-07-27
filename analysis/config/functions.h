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

float Minv(const Vec_f& pt, const Vec_f& eta, const Vec_f& phi, const Vec_f& m) {
  PtEtaPhiMVector p1(pt[0], eta[0], phi[0], m[0]);
  PtEtaPhiMVector p2(pt[1], eta[1], phi[1], m[1]);
  return (p1 + p2).M();
}


float minDeta(const float& etaDiMu, const float& jetEta1, const float& jetEta2) {

  return   std::min(abs(etaDiMu-jetEta1),abs(etaDiMu-jetEta2));
  
}  

float MinvCorr(const Vec_f& pt, const Vec_f& eta, const Vec_f& phi, const Vec_f& m,
	       const Vec_f& fsr_pt, const Vec_f& fsr_eta, const Vec_f& fsr_phi,
	       const int& typeVar) {

  // Index of the lowest-dR/ET2 among associated FSR photons
  PtEtaPhiMVector p1(pt[0], eta[0], phi[0], m[0]);
  PtEtaPhiMVector p2(pt[1], eta[1], phi[1], m[1]);
  PtEtaPhiMVector ph(fsr_pt[0], fsr_eta[0], fsr_phi[0], 0);
  // fixME add also the second FSR photon if exist

  if (typeVar==0) { return (p1 + p2 + ph).M(); }
  if (typeVar==1) { return (p1 + p2 + ph).Pt(); }
  if (typeVar==2) { return (p1 + p2 + ph).Rapidity(); }
  if (typeVar==3) { return (p1 + p2 + ph).Eta(); }
  return 0;
  
}

stdVec_i getVBFIndicies(const Vec_f& etas, const Vec_f& phis, const Vec_f& pts){

  stdVec_i idx_(2, -1);

  int n = etas.size(), index0 = -1, index1 = -1;
  float pt_limit = 20;

  float max = 0;

  for (int i = 0; i < n; i++){
    if (pts[i] < pt_limit) continue;

    for (int j = i+1; j < n; j++){
      if (pts[j] < pt_limit) continue;

      if (max < abs(etas[i] - etas[j]) && etas[i]*etas[j] < 0) {
        max = abs(etas[i] - etas[j]);
        index0 = i;
        index1 = j;
      }
    }
  }

  idx_[0] = index0;
  idx_[1] = index1;

  return idx_;

}


float getZep(const float& mu1_pt, const float& mu1_eta, const float& mu1_phi, const float& mu1_m,
             const float& mu2_pt, const float& mu2_eta, const float& mu2_phi, const float& mu2_m,
             const float& VBF1_eta, const float& VBF2_eta
             ){

  PtEtaPhiMVector mu1(mu1_pt, mu1_eta, mu1_phi, mu1_m);
  PtEtaPhiMVector mu2(mu2_pt, mu2_eta, mu2_phi, mu2_m);

  PtEtaPhiMVector HiggsCand = mu1 + mu2;

  float ZepVar = ( HiggsCand.Rapidity() - (VBF1_eta+VBF2_eta)/2.0)/abs(VBF1_eta - VBF2_eta);

  return ZepVar;

}

float getRpt(const float& mu1_pt, const float& mu1_eta, const float& mu1_phi, const float& mu1_m,
             const float& mu2_pt, const float& mu2_eta, const float& mu2_phi, const float& mu2_m,
             const float& VBF1_pt, const float& VBF1_eta, const float& VBF1_phi, const float& VBF1_mass,
             const float& VBF2_pt, const float& VBF2_eta, const float& VBF2_phi, const float& VBF2_mass
             ){

  PtEtaPhiMVector mu1(mu1_pt, mu1_eta, mu1_phi, mu1_m);
  PtEtaPhiMVector mu2(mu2_pt, mu2_eta, mu2_phi, mu2_m);  

  PtEtaPhiMVector HiggsCand = mu1 + mu2;

  PtEtaPhiMVector VBF1Cand(VBF1_pt,VBF1_eta,VBF1_phi,VBF1_mass);
  PtEtaPhiMVector VBF2Cand(VBF2_pt,VBF2_eta,VBF2_phi,VBF2_mass);

  PtEtaPhiMVector All_Vec = VBF1Cand + VBF2Cand + HiggsCand;

  float Rpt = abs(All_Vec.Pt())/( abs(VBF1_pt) + abs(VBF2_pt) + abs(HiggsCand.Pt()) );
  return Rpt;

}

float Pair12Minv(const float& v1_Pt, const float& v1_Eta, const float& v1_Phi, const float& v1_Mass,
                 const float& v2_Pt, const float& v2_Eta, const float& v2_Phi, const float& v2_Mass) {
  PtEtaPhiMVector p1(v1_Pt, v1_Eta, v1_Phi, v1_Mass);
  PtEtaPhiMVector p2(v2_Pt, v2_Eta, v2_Phi, v2_Mass);
  return (p1 + p2).M();

}

#endif

