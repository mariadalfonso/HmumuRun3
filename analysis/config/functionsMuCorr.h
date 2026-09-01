#ifndef MYFUN2                                                                                                                                                                       
#define MYFUN2                                                                                                                                                                      

using Vec_b = ROOT::VecOps::RVec<bool>;
using Vec_f = ROOT::VecOps::RVec<float>;
using Vec_i = ROOT::VecOps::RVec<int>;
using Vec_ui = ROOT::VecOps::RVec<unsigned int>;

Vec_f computeMUcorrection(Vec_f mu_pt, Vec_f mu_eta, Vec_f mu_phi, Vec_i mu_charge, Vec_ui mu_nTrackerLayers, bool isData, float event, float luminosityBlock ){

  Vec_f new_mu(mu_pt.size());

  for (unsigned int idx = 0; idx < mu_pt.size(); ++idx) {

    // Data apply scale correction
    // MC apply scale shift and resolution smearing

    if (isData) {
      new_mu[idx] = pt_scale(isData, mu_pt[idx], mu_eta[idx], mu_phi[idx], mu_charge[idx]);
      //      std::cout << " applying to the Data: org = " << mu_pt[idx] << " new = " << new_mu[idx] << std::endl;
    }
    else {
      auto pt_scale_corr = pt_scale(isData, mu_pt[idx], mu_eta[idx], mu_phi[idx], mu_charge[idx]); // this does not correct for below 26 GeV
      new_mu[idx] = pt_resol(pt_scale_corr, mu_eta[idx], mu_phi[idx], float(mu_nTrackerLayers[idx]), event, luminosityBlock); // this need to be debugged !!!!
      //      new_mu[idx] = pt_scale_corr;

      /*
      auto delta = new_mu[idx] - pt_scale_corr;
      std::cout << " applying to the MC: org = " << mu_pt[idx] << " new = " << new_mu[idx]
	        << " scale (before sm) = " << pt_scale_corr
	        << " delta res = " << delta 
		<< " mu_eta[idx] = " << mu_eta[idx]
		<< " mu_phi[idx] = " << mu_phi[idx]
	        << " mu_charge[idx] = "  << mu_charge[idx]
		<< " mu_nTrackerLayers[idx] = " << mu_nTrackerLayers[idx]
		<< std::endl;      
      */
    }
  }

  return new_mu;

}

#endif
