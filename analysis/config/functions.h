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

const float muon_mass_ = 0.10566;
const float Z_mass_ = 91.1880; // GeV
const float H_mass_ = 125.0; // GeV

// TO DO: implement the correction
// https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopPtReweighting

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


// $$$$$$$$$$
// above generic RDF analysis functinos
// below analysis oriented functions


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

bool freeOfZ(const Vec_f& pts, const Vec_f& etas, const Vec_f& phis, const Vec_f& charges){
// useful for the Wmunu + Hmm 3mu final state

  bool freeOfZ=true;

  int n = etas.size();

  for (int i = 0; i < n; i++){

    for (int j = i+1; j < n; j++){

      if (charges[i]*charges[j]>0) continue; // looking for OS

      PtEtaPhiMVector pi(pts[i], etas[i], phis[i], muon_mass_);
      PtEtaPhiMVector pj(pts[j], etas[j], phis[j], muon_mass_);

      float delta = abs((pi + pj).M() - Z_mass_);

      if ( delta <10 ) freeOfZ=false;

    }
  }

  return freeOfZ;

}

struct Pair {
    int i, j;
    float mass;
};

struct PairingResult {
    Pair Zpair;
    Pair Hpair;
    float score; // total deviation from MZ + MH
    bool valid;
};


PairingResult findBestZHCombo(const Vec_f& pts, const Vec_f& etas, const Vec_f& phis, const Vec_f& charges) {
// useful for the Zmumu + Hmm 4mu final state

  auto makePair = [&](int i, int j) {
    PtEtaPhiMVector pi(pts[i], etas[i], phis[i], muon_mass_);
    PtEtaPhiMVector pj(pts[j], etas[j], phis[j], muon_mass_);
    return Pair{i, j, (pi+pj).M()};
  };

  // possible disjoint pairings
  std::vector<std::pair<std::pair<int,int>, std::pair<int,int>>> partitions = {
    {{0,1},{2,3}}, {{0,2},{1,3}}, {{0,3},{1,2}}
  };

  PairingResult best;
  best.valid = false;
  best.score = std::numeric_limits<float>::max();

  for (auto &p : partitions) {
    int i1 = p.first.first, j1 = p.first.second;
    int i2 = p.second.first, j2 = p.second.second;

    // enforce OS requirement
    if (charges[i1] * charges[j1] != -1) continue;
    if (charges[i2] * charges[j2] != -1) continue;

    Pair p1 = makePair(i1, j1);
    Pair p2 = makePair(i2, j2);

    // assign which pair is Z vs H
    float d1Z = std::abs(p1.mass - Z_mass_);
    float d1H = std::abs(p1.mass - H_mass_);
    float d2Z = std::abs(p2.mass - Z_mass_);
    float d2H = std::abs(p2.mass - H_mass_);

    PairingResult candidate;
    candidate.valid = true;

    if (d1Z + d2H < d1H + d2Z) {
      candidate.Zpair = p1;
      candidate.Hpair = p2;
      candidate.score = d1Z + d2H;
    } else {
      candidate.Zpair = p2;
      candidate.Hpair = p1;
      candidate.score = d2Z + d1H;
    }

    if (candidate.score < best.score) {
      best = candidate;
    }
  }

  return best;

}

stdVec_i getMuonIndices(const Vec_f& pts, const Vec_f& etas, const Vec_f& phis, const Vec_f& charges, const std::string mode){

  stdVec_i idx_(2, -1);

  int n = etas.size();

  if (n <2 )  { return idx_; }
  else if (n == 2)  {
    // target hadronic final state
    idx_[0] = 0; idx_[1] = 1;

  } else if (n == 4 and mode=="isVlep") {
    // target Z-->mumu and H--> mumu

    auto result = findBestZHCombo(pts, etas, phis, charges);

    if (result.valid) {
      idx_[0] = result.Hpair.i;
      idx_[1] = result.Hpair.j;
    } // check if not valid if can reuse those with high HPT

  } else {

    float max = 0;

    int n = etas.size(), index0 = -1, index1 = -1;

    for (int i = 0; i < n; i++){

      for (int j = i+1; j < n; j++){

	if (charges[i]*charges[j]>0) continue; // looking for OS

	PtEtaPhiMVector pi(pts[i], etas[i], phis[i], muon_mass_);
	PtEtaPhiMVector pj(pts[j], etas[j], phis[j], muon_mass_);

	float PT = (pi+pj).Pt();
	if (max < PT ) {
	  max = PT;
	  index0 = i;
	  index1 = j;
	}
      }
    }

    idx_[0] = index0;
    idx_[1] = index1;

  }

  return idx_;

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

