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

Vec_b fatJetMask(Vec_i FatJet_muonIdx3SJ, int muon1_idx, int muon2_idx) {

  Vec_b mask(FatJet_muonIdx3SJ.size(), true);
  for (int idx : FatJet_muonIdx3SJ) {
    if(idx < 0) continue;
    if(idx == muon1_idx or idx == muon2_idx) mask[idx] = false;
  }
  return mask;
}


TLorentzVector MakeTLV(const float pt, const float eta, const float phi, const float mass) {
    TLorentzVector v;
    v.SetPtEtaPhiM(pt, eta, phi, mass);
    return v;
}

float deltaPhi(float phi1, float phi2) {
  float result = phi1 - phi2;
  while (result > float(M_PI)) result -= float(2*M_PI);
  while (result <= -float(M_PI)) result += float(2*M_PI);
  return result;
}

float deltaR2(float eta1, float phi1, float eta2, float phi2) {
  float deta = eta1-eta2;
  float dphi = deltaPhi(phi1,phi2);
  return deta*deta + dphi*dphi;
}

float deltaR(float eta1, float phi1, float eta2, float phi2) {
  return std::sqrt(deltaR2(eta1,phi1,eta2,phi2));
}

Vec_b deltaRMask(float eta1, float phi1,
		 const Vec_f & eta2,
		 const Vec_f & phi2,
		 float minDR = 0.8)
{
    Vec_b mask(eta2.size(), true);

    for (size_t i = 0; i < eta2.size(); i++) {
        float dr = deltaR(eta1, phi1, eta2[i], phi2[i]);
        mask[i] = (dr < minDR) ? false : true;
    }
    return mask;
}


float minDRmusJ(const TLorentzVector& VBFCand, const TLorentzVector& mu1, const TLorentzVector& mu2 ){

  float dr1 = deltaR2(VBFCand.Eta(), VBFCand.Phi(), mu1.Eta(), mu1.Phi());
  float dr2 = deltaR2(VBFCand.Eta(), VBFCand.Phi(), mu2.Eta(), mu2.Phi());

  return std::min(sqrt(dr1),sqrt(dr2));
}

// $$$$$$$$$$
// here the angles i.e. Collin Sopper will move these to another file

std::pair<float, float> CollinSopperAngles(const TLorentzVector& mu1, const TLorentzVector& mu2, float mu1_charge, float mu2_charge) {

 // Set references correctly at initialization
  const TLorentzVector& pMinus = (mu1_charge < 0 ? mu1 : mu2);
  const TLorentzVector& pPlus  = (mu1_charge < 0 ? mu2 : mu1);

  // Dilepton 4-vector and boost vector to go to its rest frame
  TLorentzVector Higgs = pPlus + pMinus;
  TVector3 beta = Higgs.BoostVector(); // boost to lab->dilepton rest: Boost(-beta)

  // Copy the mu- and boost into the dilepton rest frame
  TLorentzVector pMinus_rf = pMinus;
  pMinus_rf.Boost(-beta);

  // Construct (approximate) beam four-vectors in lab frame.
  // Only directions matter; energies can be arbitrary positive numbers.
  // We use unit z directions for the two beams.
  TLorentzVector beam1_lab; // +z
  TLorentzVector beam2_lab; // -z
  beam1_lab.SetPxPyPzE(0.0, 0.0, 1.0, 1.0);
  beam2_lab.SetPxPyPzE(0.0, 0.0,-1.0, 1.0);

  // Boost beams into dilepton rest frame
  TLorentzVector beam1_rf = beam1_lab;
  TLorentzVector beam2_rf = beam2_lab;
  beam1_rf.Boost(-beta);
  beam2_rf.Boost(-beta);

  // Use spatial parts to make unit vectors
  TVector3 b1 = beam1_rf.Vect();
  TVector3 b2 = beam2_rf.Vect();
  // Protect against degenerate cases
  if (b1.Mag() == 0 || b2.Mag() == 0) {
    return std::make_pair(std::nanf(""), std::nanf(""));
  }
  TVector3 b1u = b1.Unit();
  TVector3 b2u = b2.Unit();

  // Collins-Soper axes:
  TVector3 zCS = (b1u - b2u).Unit(); // z axis (bisector)
  TVector3 xCS = (b1u + b2u).Unit(); // x axis (sum)
  TVector3 yCS = zCS.Cross(xCS).Unit();

  // Momentum of the selected lepton (mu-) in CS frame
  TVector3 p = pMinus_rf.Vect();
  double pMag = p.Mag();
  if (pMag == 0) {
    return std::make_pair(std::nanf(""), std::nanf(""));
  }

  // cos(theta_CS)
  double cosTheta = p.Dot(zCS) / pMag;
  // phi_CS from x,y projections
  double projX = p.Dot(xCS);
  double projY = p.Dot(yCS);
  double phi = std::atan2(projY, projX); // range (-pi, pi)

  return std::make_pair(static_cast<float>(cosTheta),
			static_cast<float>(phi));

}

// $$$$$$$$$$
// above generic RDF analysis functinos
// below analysis oriented functions

Vec_i getLHEPart_match(Vec_f & Jeta, Vec_f & Jphi, Vec_f & LHEPart_eta, Vec_f & LHEPart_phi, Vec_i & LHEPart_status, Vec_i & LHEPart_pdgId, bool isData) {

  Vec_i idxMatch_(Jeta.size(), -1);
  if(isData) return idxMatch_;

  for (unsigned int idx = 0; idx < Jeta.size(); ++idx) {
    for (unsigned int j = 0; j < LHEPart_eta.size(); ++j) {
      if(LHEPart_status[j]!=1) continue; // outgoing
      float dr2 = deltaR2(Jeta[idx], Jphi[idx], LHEPart_eta[j], LHEPart_phi[j]);
      if(dr2>0.5*0.5) continue;
      idxMatch_[idx] = j;
    }
  }
  return idxMatch_;
}

float mt(float pt1, float phi1, float pt2, float phi2) {
  return std::sqrt(2*pt1*pt2*(1-std::cos(phi1-phi2)));
}

int topology(float eta1, float eta2) {

  int topology = 0;
  if (abs(eta1)<1.4 and abs(eta2)<1.4) topology = 1; // this is BB
  if (abs(eta1)<1.4 and abs(eta2)>1.4) topology = 2; // this is BE
  if (abs(eta1)>1.4 and abs(eta2)<1.4) topology = 3; // this is EB
  if (abs(eta1)>1.4 and abs(eta2)>1.4) topology = 4; // this is EE

  return topology;
}



int hardest_pt_idx(const Vec_f pts) {

  int idx = -1;
  float maxpt = -1.0f;
  for (size_t i=0; i<pts.size(); ++i) {
    if (pts[i] > maxpt) {
      maxpt = pts[i];
      idx = i;
    }
  }
  return idx;

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

float Minv(const TLorentzVector& p1, const TLorentzVector& p2) {
  return (p1 + p2).M();
}

float MinvErr(const float pt1, const float err1, const float pt2, const float err2) {
  return sqrt((err1*err1)/(pt1*pt1) + (err2*err2)/(pt2*pt2));
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

float wrongOSSFmass(const ROOT::VecOps::RVec<float>& pts,
		    const ROOT::VecOps::RVec<float>& etas,
		    const ROOT::VecOps::RVec<float>& phis,
		    const ROOT::VecOps::RVec<int>& charges,
		    int mu1_idx, int mu2_idx)
{
    const float Z_MASS = 91.1876;

    float bestMass = -1.f;
    float bestScore = 1e9;  // distance from Z or low-mass threshold

    int n = pts.size();

    // Loop over all other muons
    for (int k = 0; k < n; k++) {
        if (k == mu1_idx || k == mu2_idx) continue;

        // Check pair (mu1, k)
        if (charges[mu1_idx] * charges[k] < 0) {
            TLorentzVector a, b;
            a.SetPtEtaPhiM(pts[mu1_idx], etas[mu1_idx], phis[mu1_idx], muon_mass_);
            b.SetPtEtaPhiM(pts[k],       etas[k],       phis[k],       muon_mass_);
            float M = (a + b).M();

            bool nearZ   = fabs(M - Z_MASS) < 15.;
            bool lowMass = M < 12.;

            if (nearZ || lowMass) {
                float score = nearZ ? fabs(M - Z_MASS) : (12. - M);
                if (score < bestScore) {
                    bestScore = score;
                    bestMass = M;
                }
            }
        }

        // Check pair (mu2, k)
        if (charges[mu2_idx] * charges[k] < 0) {
            TLorentzVector a, b;
            a.SetPtEtaPhiM(pts[mu2_idx], etas[mu2_idx], phis[mu2_idx], muon_mass_);
            b.SetPtEtaPhiM(pts[k],       etas[k],       phis[k],       muon_mass_);
            float M = (a + b).M();

            bool nearZ   = fabs(M - Z_MASS) < 15.;
            bool lowMass = M < 12.;

            if (nearZ || lowMass) {
                float score = nearZ ? fabs(M - Z_MASS) : (12. - M);
                if (score < bestScore) {
                    bestScore = score;
                    bestMass = M;
                }
            }
        }
    }

    return bestMass;   // -1 if no alternative OS pair found in signal regions
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

    if (charges[0]*charges[1]>0) {idx_[0] = -1; idx_[1] = -1; } // looking for OS
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
      float DetaJJ = abs(etas[i] - etas[j]);

      //      if (max < abs(etas[i] - etas[j]) && etas[i]*etas[j] < 0) {
      if (max < MJJ ) {
	max = MJJ;
	//      if (max < DetaJJ ) {
	//        max = DetaJJ;
        index0 = i;
        index1 = j;
      }
    }
  }

  // If a valid pair was found, reorder as leading/subleading pT
  if (index0 >= 0 && index1 >= 0){
    int leading  = (pts[index0] >= pts[index1]) ? index0 : index1;
    int sublead  = (pts[index0] >= pts[index1]) ? index1 : index0;

    idx_[0] = leading;
    idx_[1] = sublead;
  }

  //  idx_[0] = index0;
  //  idx_[1] = index1;

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

