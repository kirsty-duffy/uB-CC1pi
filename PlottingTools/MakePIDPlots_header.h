#include "../Algorithms/TopologyEnums.h"
#include "TTree.h"
#include "TH1.h"
#include "TColor.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TStyle.h"
#include "TFile.h"
#include <vector>

struct treevars{
  // These are the variables that are filled directly from the tree
  std::vector<int> *TPCObj_PFP_truePDG = nullptr;
  std::vector<bool> *TPCObj_PFP_PandoraClassedAsTrack = nullptr;
  std::vector<std::vector<double>> *TPCObj_PFP_n2LLH_fwd_mu = nullptr;
  std::vector<std::vector<double>> *TPCObj_PFP_n2LLH_fwd_p = nullptr;
  std::vector<std::vector<double>> *TPCObj_PFP_n2LLH_fwd_pi = nullptr;
  std::vector<std::vector<double>> *TPCObj_PFP_n2LLH_bwd_mu = nullptr;
  std::vector<std::vector<double>> *TPCObj_PFP_n2LLH_bwd_p = nullptr;
  std::vector<std::vector<double>> *TPCObj_PFP_n2LLH_bwd_pi = nullptr;
  std::vector<std::vector<double>> *TPCObj_PFP_n2LLH_MIP = nullptr;
  std::vector<std::vector<double>> *TPCObj_PFP_PIDA = nullptr;
  std::vector<std::vector<double>> *TPCObj_PFP_track_depE = nullptr;
  std::vector<std::vector<double>> *TPCObj_PFP_track_Chi2Proton = nullptr;
  std::vector<std::vector<double>> *TPCObj_PFP_track_Chi2Muon = nullptr;
  std::vector<std::vector<double>> *TPCObj_PFP_track_Chi2Pion = nullptr;
  std::vector<double> *TPCObj_PFP_track_rangeE_mu = nullptr;
  std::vector<double> *TPCObj_PFP_track_rangeE_p = nullptr;

  std::vector<double> *TPCObj_PFP_trueEndP = nullptr;
  NuIntTopology Truth_topology = kUnknown;

  // These are derived quantities - derived from the values above in CalcPIDvars
  std::vector<std::vector<double>> *TPCObj_PFP_n2LLH_p;
  std::vector<std::vector<double>> *TPCObj_PFP_n2LLH_mu;
  std::vector<std::vector<double>> *TPCObj_PFP_n2LLH_pi;
  std::vector<std::vector<double>> *TPCObj_PFP_n2LLH_mip;
  std::vector<std::vector<double>> *TPCObj_PFP_n2LLH_minmumip;
  std::vector<std::vector<double>> *TPCObj_PFP_n2LLH_muminusp;
  std::vector<std::vector<double>> *TPCObj_PFP_n2LLH_mipminusp;
  std::vector<std::vector<double>> *TPCObj_PFP_n2LLH_minmumipminusp;
  std::vector<std::vector<double>> *TPCObj_PFP_depE_minus_rangeE_mu;
  std::vector<std::vector<double>> *TPCObj_PFP_depE_minus_rangeE_p;
  std::vector<std::vector<double>> *TPCObj_PFP_chi2_muminusp;

  // These variables are for when "neglogl" is actually not a log, but just a likelihood (then we want the maximum likelihood, and likelihood ratio instead of subtraction)
  std::vector<std::vector<double>> *TPCObj_PFP_max_Lp;
  std::vector<std::vector<double>> *TPCObj_PFP_max_Lmu;
  std::vector<std::vector<double>> *TPCObj_PFP_max_Lpi;
  std::vector<std::vector<double>> *TPCObj_PFP_max_Lk;
  std::vector<std::vector<double>> *TPCObj_PFP_max_Lmip;
  std::vector<std::vector<double>> *TPCObj_PFP_max_Lmaxmumip;
  std::vector<std::vector<double>> *TPCObj_PFP_Lmuoverp;
  std::vector<std::vector<double>> *TPCObj_PFP_Lmipoverp;
  std::vector<std::vector<double>> *TPCObj_PFP_Lmaxmumipoverp;
  std::vector<std::vector<double>> *TPCObj_PFP_Lmu_0to1;
  std::vector<std::vector<double>> *TPCObj_PFP_Lmip_0to1;
  std::vector<std::vector<double>> *TPCObj_PFP_Lpi_0to1;
  std::vector<std::vector<double>> *TPCObj_PFP_Lp_0to1;
  std::vector<std::vector<double>> *TPCObj_PFP_Lmumip_0to1;
  std::vector<std::vector<double>> *TPCObj_PFP_Lmumippi_0to1;
};

void settreevars(TTree *intree, treevars *varstoset){
  intree->SetBranchStatus("*",0);
  intree->SetBranchStatus("TPCObj_PFP_truePDG",1);
  intree->SetBranchAddress("TPCObj_PFP_truePDG", &(varstoset->TPCObj_PFP_truePDG));
  intree->SetBranchStatus("TPCObj_PFP_PandoraClassedAsTrack",1);
  intree->SetBranchAddress("TPCObj_PFP_PandoraClassedAsTrack", &(varstoset->TPCObj_PFP_PandoraClassedAsTrack));
  intree->SetBranchStatus("TPCObj_PFP_n2LLH_fwd_mu",1);
  intree->SetBranchAddress("TPCObj_PFP_n2LLH_fwd_mu", &(varstoset->TPCObj_PFP_n2LLH_fwd_mu));
  intree->SetBranchStatus("TPCObj_PFP_n2LLH_fwd_p",1);
  intree->SetBranchAddress("TPCObj_PFP_n2LLH_fwd_p", &(varstoset->TPCObj_PFP_n2LLH_fwd_p));
  intree->SetBranchStatus("TPCObj_PFP_n2LLH_fwd_pi",1);
  intree->SetBranchAddress("TPCObj_PFP_n2LLH_fwd_pi", &(varstoset->TPCObj_PFP_n2LLH_fwd_pi));
  intree->SetBranchStatus("TPCObj_PFP_n2LLH_bwd_mu",1);
  intree->SetBranchAddress("TPCObj_PFP_n2LLH_bwd_mu", &(varstoset->TPCObj_PFP_n2LLH_bwd_mu));
  intree->SetBranchStatus("TPCObj_PFP_n2LLH_bwd_p",1);
  intree->SetBranchAddress("TPCObj_PFP_n2LLH_bwd_p", &(varstoset->TPCObj_PFP_n2LLH_bwd_p));
  intree->SetBranchStatus("TPCObj_PFP_n2LLH_bwd_pi",1);
  intree->SetBranchAddress("TPCObj_PFP_n2LLH_bwd_pi", &(varstoset->TPCObj_PFP_n2LLH_bwd_pi));
  intree->SetBranchStatus("TPCObj_PFP_n2LLH_MIP",1);
  intree->SetBranchAddress("TPCObj_PFP_n2LLH_MIP", &(varstoset->TPCObj_PFP_n2LLH_MIP));
  intree->SetBranchStatus("TPCObj_PFP_PIDA",1);
  intree->SetBranchAddress("TPCObj_PFP_PIDA", &(varstoset->TPCObj_PFP_PIDA));
  intree->SetBranchStatus("TPCObj_PFP_track_depE",1);
  intree->SetBranchAddress("TPCObj_PFP_track_depE", &(varstoset->TPCObj_PFP_track_depE));
  intree->SetBranchStatus("TPCObj_PFP_track_Chi2Proton",1);
  intree->SetBranchAddress("TPCObj_PFP_track_Chi2Proton", &(varstoset->TPCObj_PFP_track_Chi2Proton));
  intree->SetBranchStatus("TPCObj_PFP_track_Chi2Muon",1);
  intree->SetBranchAddress("TPCObj_PFP_track_Chi2Muon", &(varstoset->TPCObj_PFP_track_Chi2Muon));
  intree->SetBranchStatus("TPCObj_PFP_track_Chi2Pion",1);
  intree->SetBranchAddress("TPCObj_PFP_track_Chi2Pion", &(varstoset->TPCObj_PFP_track_Chi2Pion));
  intree->SetBranchStatus("TPCObj_PFP_track_rangeE_mu",1);
  intree->SetBranchAddress("TPCObj_PFP_track_rangeE_mu", &(varstoset->TPCObj_PFP_track_rangeE_mu));
  intree->SetBranchStatus("TPCObj_PFP_track_rangeE_p",1);
  intree->SetBranchAddress("TPCObj_PFP_track_rangeE_p", &(varstoset->TPCObj_PFP_track_rangeE_p));
  intree->SetBranchStatus("TPCObj_PFP_trueEndP",1);
  intree->SetBranchAddress("TPCObj_PFP_trueEndP", &(varstoset->TPCObj_PFP_trueEndP));
  intree->SetBranchStatus("Truth_topology",1);
  intree->SetBranchAddress("Truth_topology", &(varstoset->Truth_topology));

  // intree->GetEntry(0);
  // size_t ntracks = varstoset->TPCObj_PFP_n2LLH_fwd_p->size();
  // size_t nplanes = varstoset->TPCObj_PFP_n2LLH_fwd_p->at(0)->size();

  // Assume 3 planes

  // varstoset->TPCObj_PFP_n2LLH_p = new std::vector<std::vector<double>>;
  // varstoset->TPCObj_PFP_n2LLH_mu = new std::vector<std::vector<double>>;
  // varstoset->TPCObj_PFP_n2LLH_pi = new std::vector<std::vector<double>>;
  // varstoset->TPCObj_PFP_n2LLH_mip = new std::vector<std::vector<double>>;
  // varstoset->TPCObj_PFP_n2LLH_minmumip = new std::vector<std::vector<double>>;
  // varstoset->TPCObj_PFP_n2LLH_muminusp = new std::vector<std::vector<double>>;
  // varstoset->TPCObj_PFP_n2LLH_mipminusp = new std::vector<std::vector<double>>;
  // varstoset->TPCObj_PFP_n2LLH_minmumipminusp = new std::vector<std::vector<double>>;
  // varstoset->TPCObj_PFP_depE_minus_rangeE_mu = new std::vector<std::vector<double>>;
  // varstoset->TPCObj_PFP_depE_minus_rangeE_p = new std::vector<std::vector<double>>;
  // varstoset->TPCObj_PFP_chi2_muminusp = new std::vector<std::vector<double>>;
  //
  // varstoset->TPCObj_PFP_n2LLH_p.at(0) = new std::vector<double>(3);
  // varstoset->TPCObj_PFP_n2LLH_mu.at(0) = new std::vector<double>(3);
  // varstoset->TPCObj_PFP_n2LLH_pi.at(0) = new std::vector<double>(3);
  // varstoset->TPCObj_PFP_n2LLH_mip.at(0) = new std::vector<double>(3);
  // varstoset->TPCObj_PFP_n2LLH_minmumip.at(0) = new std::vector<double>(3);
  // varstoset->TPCObj_PFP_n2LLH_muminusp.at(0) = new std::vector<double>(3);
  // varstoset->TPCObj_PFP_n2LLH_mipminusp.at(0) = new std::vector<double>(3);
  // varstoset->TPCObj_PFP_n2LLH_minmumipminusp.at(0) = new std::vector<double>(3);
  // varstoset->TPCObj_PFP_depE_minus_rangeE_mu.at(0) = new std::vector<double>(3);
  // varstoset->TPCObj_PFP_depE_minus_rangeE_p.at(0) = new std::vector<double>(3);
  // varstoset->TPCObj_PFP_chi2_muminusp.at(0) = new std::vector<double>(3);
  //
  // varstoset->TPCObj_PFP_max_Lp = new std::vector<double>(3);
  // varstoset->TPCObj_PFP_max_Lmu = new std::vector<double>(3);
  // varstoset->TPCObj_PFP_max_Lpi = new std::vector<double>(3);
  // varstoset->TPCObj_PFP_max_Lmip = new std::vector<double>(3);
  // varstoset->TPCObj_PFP_max_Lmaxmumip = new std::vector<double>(3);
  // varstoset->TPCObj_PFP_Lmuoverp = new std::vector<double>(3);
  // varstoset->TPCObj_PFP_Lmipoverp = new std::vector<double>(3);
  // varstoset->TPCObj_PFP_Lminmumipoverp = new std::vector<double>(3);
  // varstoset->TPCObj_PFP_Lmu_0to1 = new std::vector<double>(3);
  // varstoset->TPCObj_PFP_Lmip_0to1 = new std::vector<double>(3);
  // varstoset->TPCObj_PFP_Lpi_0to1 = new std::vector<double>(3);
  // varstoset->TPCObj_PFP_Lp_0to1 = new std::vector<double>(3);
  // varstoset->TPCObj_PFP_Lmumip_0to1 = new std::vector<double>(3);
  // varstoset->TPCObj_PFP_Lmumippi_0to1 = new std::vector<double>(3);
}

void CalcPIDvars(treevars *vars){

  // Initialise vectors that we are going to fill with calculated values
  int vecsize = vars->TPCObj_PFP_n2LLH_fwd_p->size();
  int nplanes = 3;

  vars->TPCObj_PFP_n2LLH_p = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));
  vars->TPCObj_PFP_n2LLH_mu = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));
  vars->TPCObj_PFP_n2LLH_pi = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));
  vars->TPCObj_PFP_n2LLH_mip = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));
  vars->TPCObj_PFP_n2LLH_minmumip = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));
  vars->TPCObj_PFP_n2LLH_muminusp = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));
  vars->TPCObj_PFP_n2LLH_mipminusp = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));
  vars->TPCObj_PFP_n2LLH_minmumipminusp = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));
  vars->TPCObj_PFP_chi2_muminusp = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));
  vars->TPCObj_PFP_depE_minus_rangeE_mu = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));
  vars->TPCObj_PFP_depE_minus_rangeE_p = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));
  vars->TPCObj_PFP_max_Lp = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));
  vars->TPCObj_PFP_max_Lmu = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));
  vars->TPCObj_PFP_max_Lpi = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));
  vars->TPCObj_PFP_max_Lmip = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));
  vars->TPCObj_PFP_max_Lmaxmumip = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));
  vars->TPCObj_PFP_Lmuoverp = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));
  vars->TPCObj_PFP_Lmipoverp = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));
  vars->TPCObj_PFP_Lmaxmumipoverp = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));
  vars->TPCObj_PFP_Lmu_0to1 = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));
  vars->TPCObj_PFP_Lmip_0to1 = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));
  vars->TPCObj_PFP_Lpi_0to1 = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));
  vars->TPCObj_PFP_Lp_0to1 = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));
  vars->TPCObj_PFP_Lmumip_0to1 = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));
  vars->TPCObj_PFP_Lmumippi_0to1 = new std::vector<std::vector<double>>(vecsize, std::vector<double>(nplanes));

  // Now calculate the values for all variables
  for (int i_track=0; i_track < vecsize; i_track++){
    for (int i_pl=0; i_pl < nplanes; i_pl++){
      vars->TPCObj_PFP_n2LLH_p->at(i_track).at(i_pl)     = std::min(vars->TPCObj_PFP_n2LLH_fwd_p->at(i_track).at(i_pl)     , vars->TPCObj_PFP_n2LLH_bwd_p->at(i_track).at(i_pl));
      vars->TPCObj_PFP_n2LLH_mu->at(i_track).at(i_pl)    = std::min(vars->TPCObj_PFP_n2LLH_fwd_mu->at(i_track).at(i_pl)    , vars->TPCObj_PFP_n2LLH_bwd_mu->at(i_track).at(i_pl));
      vars->TPCObj_PFP_n2LLH_pi->at(i_track).at(i_pl)    = std::min(vars->TPCObj_PFP_n2LLH_fwd_pi->at(i_track).at(i_pl)    , vars->TPCObj_PFP_n2LLH_bwd_pi->at(i_track).at(i_pl));
      vars->TPCObj_PFP_n2LLH_mip->at(i_track).at(i_pl)   = vars->TPCObj_PFP_n2LLH_MIP->at(i_track).at(i_pl);

      vars->TPCObj_PFP_n2LLH_minmumip->at(i_track).at(i_pl) = std::min(vars->TPCObj_PFP_n2LLH_mu->at(i_track).at(i_pl), vars->TPCObj_PFP_n2LLH_mip->at(i_track).at(i_pl));

      vars->TPCObj_PFP_n2LLH_muminusp->at(i_track).at(i_pl) = vars->TPCObj_PFP_n2LLH_mu->at(i_track).at(i_pl) - vars->TPCObj_PFP_n2LLH_p->at(i_track).at(i_pl);
      vars->TPCObj_PFP_n2LLH_mipminusp->at(i_track).at(i_pl) = vars->TPCObj_PFP_n2LLH_mip->at(i_track).at(i_pl) - vars->TPCObj_PFP_n2LLH_p->at(i_track).at(i_pl);
      vars->TPCObj_PFP_n2LLH_minmumipminusp->at(i_track).at(i_pl) = vars->TPCObj_PFP_n2LLH_minmumip->at(i_track).at(i_pl) - vars->TPCObj_PFP_n2LLH_p->at(i_track).at(i_pl);

      vars->TPCObj_PFP_chi2_muminusp->at(i_track).at(i_pl) = vars->TPCObj_PFP_track_Chi2Muon->at(i_track).at(i_pl) - vars->TPCObj_PFP_track_Chi2Proton->at(i_track).at(i_pl);

      vars->TPCObj_PFP_depE_minus_rangeE_mu->at(i_track).at(i_pl) = vars->TPCObj_PFP_track_depE->at(i_track).at(i_pl) - vars->TPCObj_PFP_track_rangeE_mu->at(i_track);
      vars->TPCObj_PFP_depE_minus_rangeE_p->at(i_track).at(i_pl) = vars->TPCObj_PFP_track_depE->at(i_track).at(i_pl) - vars->TPCObj_PFP_track_rangeE_p->at(i_track);

      // These variables are for when "neglogl" is actually not a log, but just a likelihood (then we want the maximum likelihood, and likelihood ratio instead of subtraction)
      vars->TPCObj_PFP_max_Lp->at(i_track).at(i_pl) = std::max(vars->TPCObj_PFP_n2LLH_fwd_p->at(i_track).at(i_pl)     , vars->TPCObj_PFP_n2LLH_bwd_p->at(i_track).at(i_pl));
      vars->TPCObj_PFP_max_Lmu->at(i_track).at(i_pl) = std::max(vars->TPCObj_PFP_n2LLH_fwd_mu->at(i_track).at(i_pl)    , vars->TPCObj_PFP_n2LLH_bwd_mu->at(i_track).at(i_pl));
      vars->TPCObj_PFP_max_Lpi->at(i_track).at(i_pl) = std::max(vars->TPCObj_PFP_n2LLH_fwd_pi->at(i_track).at(i_pl)    , vars->TPCObj_PFP_n2LLH_bwd_pi->at(i_track).at(i_pl));
      vars->TPCObj_PFP_max_Lmip->at(i_track).at(i_pl) = vars->TPCObj_PFP_n2LLH_MIP->at(i_track).at(i_pl);
      vars->TPCObj_PFP_max_Lmaxmumip->at(i_track).at(i_pl) = std::max(vars->TPCObj_PFP_n2LLH_mu->at(i_track).at(i_pl), vars->TPCObj_PFP_n2LLH_mip->at(i_track).at(i_pl));
      vars->TPCObj_PFP_Lmuoverp->at(i_track).at(i_pl) = vars->TPCObj_PFP_n2LLH_mu->at(i_track).at(i_pl) / vars->TPCObj_PFP_n2LLH_p->at(i_track).at(i_pl);
      vars->TPCObj_PFP_Lmipoverp->at(i_track).at(i_pl) = vars->TPCObj_PFP_n2LLH_mip->at(i_track).at(i_pl) / vars->TPCObj_PFP_n2LLH_p->at(i_track).at(i_pl);
      vars->TPCObj_PFP_Lmaxmumipoverp->at(i_track).at(i_pl) = vars->TPCObj_PFP_n2LLH_minmumip->at(i_track).at(i_pl) / vars->TPCObj_PFP_n2LLH_p->at(i_track).at(i_pl);

      vars->TPCObj_PFP_Lmu_0to1->at(i_track).at(i_pl) = vars->TPCObj_PFP_n2LLH_mu->at(i_track).at(i_pl)/(vars->TPCObj_PFP_n2LLH_p->at(i_track).at(i_pl)+vars->TPCObj_PFP_n2LLH_mu->at(i_track).at(i_pl)+vars->TPCObj_PFP_n2LLH_pi->at(i_track).at(i_pl)+vars->TPCObj_PFP_n2LLH_mip->at(i_track).at(i_pl));
      vars->TPCObj_PFP_Lmip_0to1->at(i_track).at(i_pl) = vars->TPCObj_PFP_n2LLH_mip->at(i_track).at(i_pl)/(vars->TPCObj_PFP_n2LLH_p->at(i_track).at(i_pl)+vars->TPCObj_PFP_n2LLH_mu->at(i_track).at(i_pl)+vars->TPCObj_PFP_n2LLH_pi->at(i_track).at(i_pl)+vars->TPCObj_PFP_n2LLH_mip->at(i_track).at(i_pl));
      vars->TPCObj_PFP_Lpi_0to1->at(i_track).at(i_pl) = vars->TPCObj_PFP_n2LLH_pi->at(i_track).at(i_pl)/(vars->TPCObj_PFP_n2LLH_p->at(i_track).at(i_pl)+vars->TPCObj_PFP_n2LLH_mu->at(i_track).at(i_pl)+vars->TPCObj_PFP_n2LLH_pi->at(i_track).at(i_pl)+vars->TPCObj_PFP_n2LLH_mip->at(i_track).at(i_pl));
      vars->TPCObj_PFP_Lp_0to1->at(i_track).at(i_pl) = vars->TPCObj_PFP_n2LLH_p->at(i_track).at(i_pl)/(vars->TPCObj_PFP_n2LLH_p->at(i_track).at(i_pl)+vars->TPCObj_PFP_n2LLH_mu->at(i_track).at(i_pl)+vars->TPCObj_PFP_n2LLH_pi->at(i_track).at(i_pl)+vars->TPCObj_PFP_n2LLH_mip->at(i_track).at(i_pl));

      vars->TPCObj_PFP_Lmumip_0to1->at(i_track).at(i_pl) = (vars->TPCObj_PFP_n2LLH_mu->at(i_track).at(i_pl)+vars->TPCObj_PFP_n2LLH_mip->at(i_track).at(i_pl))/(vars->TPCObj_PFP_n2LLH_p->at(i_track).at(i_pl)+vars->TPCObj_PFP_n2LLH_mu->at(i_track).at(i_pl)+vars->TPCObj_PFP_n2LLH_pi->at(i_track).at(i_pl)+vars->TPCObj_PFP_n2LLH_mip->at(i_track).at(i_pl));
      vars->TPCObj_PFP_Lmumippi_0to1->at(i_track).at(i_pl) = (vars->TPCObj_PFP_n2LLH_mu->at(i_track).at(i_pl)+vars->TPCObj_PFP_n2LLH_mip->at(i_track).at(i_pl)+vars->TPCObj_PFP_n2LLH_pi->at(i_track).at(i_pl))/(vars->TPCObj_PFP_n2LLH_p->at(i_track).at(i_pl)+vars->TPCObj_PFP_n2LLH_mu->at(i_track).at(i_pl)+vars->TPCObj_PFP_n2LLH_pi->at(i_track).at(i_pl)+vars->TPCObj_PFP_n2LLH_mip->at(i_track).at(i_pl));
    }
  }
}


// --------------------------------------------------- //
// Function to clear memory from calculated PID variables


void ClearPIDvars(treevars *vars){

  delete vars->TPCObj_PFP_n2LLH_p;
  delete vars->TPCObj_PFP_n2LLH_mu;
  delete vars->TPCObj_PFP_n2LLH_pi;
  delete vars->TPCObj_PFP_n2LLH_mip;
  delete vars->TPCObj_PFP_n2LLH_minmumip;
  delete vars->TPCObj_PFP_n2LLH_muminusp;
  delete vars->TPCObj_PFP_n2LLH_mipminusp;
  delete vars->TPCObj_PFP_n2LLH_minmumipminusp;
  delete vars->TPCObj_PFP_chi2_muminusp;
  delete vars->TPCObj_PFP_depE_minus_rangeE_mu;
  delete vars->TPCObj_PFP_depE_minus_rangeE_p;
  delete vars->TPCObj_PFP_max_Lp;
  delete vars->TPCObj_PFP_max_Lmu;
  delete vars->TPCObj_PFP_max_Lpi;
  delete vars->TPCObj_PFP_max_Lmip;
  delete vars->TPCObj_PFP_max_Lmaxmumip;
  delete vars->TPCObj_PFP_Lmuoverp;
  delete vars->TPCObj_PFP_Lmipoverp;
  delete vars->TPCObj_PFP_Lmaxmumipoverp;
  delete vars->TPCObj_PFP_Lmu_0to1;
  delete vars->TPCObj_PFP_Lmip_0to1;
  delete vars->TPCObj_PFP_Lpi_0to1;
  delete vars->TPCObj_PFP_Lp_0to1;
  delete vars->TPCObj_PFP_Lmumip_0to1;
  delete vars->TPCObj_PFP_Lmumippi_0to1;
}

// --------------------------------------------------- //
// This struct contains the histograms and all functions related to them

struct hist1D{
  TH1D *h_mu;
  TH1D *h_p;
  TH1D *h_pi;
  TH1D *h_k;
  TH1D *h_other;
  TH1D *h_all;

  TLegend *l;

  // Constructor for this struct of hists
  hist1D(std::string name, std::string title, double nbins, double binlow, double binhigh){
    h_mu = new TH1D(std::string(name+"_mu").c_str(),title.c_str(),nbins,binlow,binhigh);
    h_p = new TH1D(std::string(name+"_p").c_str(),title.c_str(),nbins,binlow,binhigh);
    h_pi = new TH1D(std::string(name+"_pi").c_str(),title.c_str(),nbins,binlow,binhigh);
    h_k = new TH1D(std::string(name+"_k").c_str(),title.c_str(),nbins,binlow,binhigh);
    h_other = new TH1D(std::string(name+"_other").c_str(),title.c_str(),nbins,binlow,binhigh);

    h_all = new TH1D(std::string(name+"_all").c_str(),title.c_str(),nbins,binlow,binhigh);

    h_mu->SetFillColor(TColor::GetColor(8,64,129));
    h_p->SetFillColor(TColor::GetColor(215, 48, 39));
    h_pi->SetFillColor(TColor::GetColor(166,217,106));
    h_k->SetFillColor(TColor::GetColor(133,1,98));
    h_other->SetFillColor(TColor::GetColor(197,197,197));
    // "All" styling for MC (data styling is set in DrawData)
    h_all->SetFillStyle(3345);
    h_all->SetFillColor(kGray+2);
    h_all->SetMarkerSize(0.); // bad hack because root keeps drawing markers and I can't make it stop

    l = new TLegend(0.59,0.64,0.81,0.87);
    l->AddEntry(h_p,"True proton","f");
    l->AddEntry(h_mu,"True muon","f");
    l->AddEntry(h_pi,"True pion","f");
    l->AddEntry(h_k,"True kaon","f");
    l->AddEntry(h_other,"True other","f");
  }
};

  void FillHist(hist1D *hists, std::vector<std::vector<double>> value_vec, int plane, std::vector<bool> *pandoraclassedastrack, std::vector<int> *pdg_vec=nullptr){
    // value_vec is a vector corresponding to the TPCObject. It has one entry per track in the TPCObject, and each of those entries has three entries (one for each plane)
    // Loop through value_vec and fill histogram for the correct plane for every track
    // Remember that plane=3 is code for (plane0+plane1+plane2)/3
    for (size_t i_track=0; i_track<value_vec.size(); i_track++){
      double value = -9999;
      if (plane==3){
        value = (value_vec.at(i_track).at(0)+value_vec.at(i_track).at(1)+value_vec.at(i_track).at(2))/3.;
      }
      else{
        value = value_vec.at(i_track).at(plane);
      }

      if (value == -9999 || value == -999) continue;
      if (!pandoraclassedastrack->at(i_track)) continue;

      int pdg = 0;
      if (pdg_vec!=nullptr) pdg = pdg_vec->at(i_track);

      // Fill "all" histogram for every entry
      hists->h_all->Fill(value);

      // Now fill histograms by particle type
      if (TMath::Abs(pdg)==13){ // muon
        hists->h_mu->Fill(value);
      }
      else if (TMath::Abs(pdg)==2212){ // proton
        hists->h_p->Fill(value);
      }
      else if (TMath::Abs(pdg)==211){ // pion
        hists->h_pi->Fill(value);
      }
      else if (TMath::Abs(pdg)==321){ // kaon
        hists->h_k->Fill(value);
      }
      else{ // other
        hists->h_other->Fill(value);
      }
    } // end loop over tracks in TPCObject
  }

  void DrawMC(hist1D *hists, double POTScaling){
    if (POTScaling == 0.){ // area normalise
      POTScaling = 1./hists->h_all->Integral();
      hists->h_all->GetYaxis()->SetTitle("No. tracks (area normalised)");
    }
    else hists->h_all->GetYaxis()->SetTitle("No. tracks (POT normalised)");

    hists->h_mu->Sumw2();
    hists->h_p->Sumw2();
    hists->h_pi->Sumw2();
    hists->h_k->Sumw2();
    hists->h_other->Sumw2();
    hists->h_all->Sumw2();

    hists->h_mu->Scale(POTScaling);
    hists->h_p->Scale(POTScaling);
    hists->h_pi->Scale(POTScaling);
    hists->h_k->Scale(POTScaling);
    hists->h_other->Scale(POTScaling);
    hists->h_all->Scale(POTScaling);

    std::cout << "h_all MC->Integral() = " << hists->h_all->Integral() << std::endl;

    THStack *hs = new THStack("hs","hs");
    hs->Add(hists->h_p);
    hs->Add(hists->h_mu);
    hs->Add(hists->h_pi);
    hs->Add(hists->h_k);
    hs->Add(hists->h_other);

    hists->h_all->SetMaximum(hists->h_all->GetMaximum()*1.2);
    hists->h_all->SetMinimum(0);
    hists->h_all->Draw("hist"); // Draw this one first because it knows about the axis titles
    hs->Draw("same hist");
    hists->h_all->Draw("same E2"); // Draw it again so errors are on top
    hists->l->Draw();
}

void DrawMCPlusOffbeam(hist1D *hists, hist1D *offbeam, double POTScaling, double OffBeamScaling){
  // Note that there are no area-normalised options here because I'm not sure that makes sense
  hists->h_all->GetYaxis()->SetTitle("No. tracks (POT normalised)");

  hists->h_mu->Sumw2();
  hists->h_p->Sumw2();
  hists->h_pi->Sumw2();
  hists->h_k->Sumw2();
  hists->h_other->Sumw2();
  hists->h_all->Sumw2();
  offbeam->h_all->Sumw2();

  hists->h_mu->Scale(POTScaling);
  hists->h_p->Scale(POTScaling);
  hists->h_pi->Scale(POTScaling);
  hists->h_k->Scale(POTScaling);
  hists->h_other->Scale(POTScaling);
  hists->h_all->Scale(POTScaling);

  offbeam->h_all->Scale(OffBeamScaling);
  offbeam->h_all->SetFillColor(kWhite);
  offbeam->h_all->SetLineColor(kBlack);

  THStack *hs = new THStack("hs","hs");
  hs->Add(hists->h_p);
  hs->Add(hists->h_mu);
  hs->Add(hists->h_pi);
  hs->Add(hists->h_k);
  hs->Add(hists->h_other);
  hs->Add(offbeam->h_all);

  hists->h_all->SetMaximum((hists->h_all->GetMaximum()+offbeam->h_all->GetMaximum())*1.2);
  hists->h_all->SetMinimum(0);
  hists->h_all->Draw("hist"); // Draw this one first because it knows about the axis titles
  hs->Draw("same hist");
  hists->h_all->Draw("same E2"); // Draw it again so errors are on top
  hists->l->AddEntry(offbeam->h_all,"Data (off-beam)","f");
  hists->l->Draw();
}

void OverlayOnMinusOffData(TCanvas *c, hist1D *onbeam, hist1D *offbeam, double OffBeamScaling, double POTScaling){
  TH1D *h_onminusoff = (TH1D*)onbeam->h_all->Clone("h_onminusoff");

  h_onminusoff->Sumw2();
  offbeam->h_all->Sumw2();

  h_onminusoff->Add(offbeam->h_all,-1.0*OffBeamScaling);
  if (POTScaling==0){
    h_onminusoff->Scale(1.0/h_onminusoff->Integral());
  }

  std::cout << "h_onminusoff->Integral() = " << h_onminusoff->Integral() << std::endl;

  h_onminusoff->SetMarkerStyle(20);
  h_onminusoff->SetMarkerSize(0.6);

  c->cd();
  h_onminusoff->Draw("same p E1");

  TLegend *l = (TLegend*)c->GetPrimitive("TPave");
  l->AddEntry(h_onminusoff,"Data (on-off beam)","lp");
}

void OverlayOnBeamData(TCanvas *c, hist1D *onbeam){

  onbeam->h_all->SetMarkerStyle(20);
  onbeam->h_all->SetMarkerSize(0.6);

  c->cd();
  onbeam->h_all->Draw("same p E1");

  TLegend *l = (TLegend*)c->GetPrimitive("TPave");
  l->AddEntry(onbeam->h_all,"Data (on-beam)","lp");
}


void DrawMCEffPur(TCanvas *c, hist1D *hists, bool MIPlow){
  std::vector<TH1D*> histstoeval = {
                    hists->h_mu,
                    hists->h_pi,
                    hists->h_p
                    };

  std::vector<std::string> histtitles = {
                    "True muons",
                    "True pions",
                    "True protons"
                    };

  c->Divide(2,2,0.0005,0.0005);

  TLegend *l = new TLegend(0.2,0.2,0.8,0.8);
  l->SetTextFont(132);
  l->SetLineColor(kWhite);
  l->SetFillColor(kWhite);

  for (size_t i_h=0; i_h<histstoeval.size(); i_h++){
    TH1D *heff = (TH1D*)hists->h_all->Clone("heff");
    TH1D *hpur = (TH1D*)hists->h_all->Clone("hpur");
    TH1D *heffpur = (TH1D*)hists->h_all->Clone("heffpur");

    heff->Clear();
    hpur->Clear();
    heffpur->Clear();

    heff->SetTitle(histtitles.at(i_h).c_str());

    for (int i_bin=1; i_bin <= histstoeval.at(i_h)->GetXaxis()->GetNbins(); i_bin++){

      double eff, pur;//, efferr, purerr;
      if ((MIPlow && i_h < 2) || (!MIPlow && i_h ==2)){ // integrate from the bottom
        double selected_i = histstoeval.at(i_h)->Integral(1,i_bin);
        double total_i = histstoeval.at(i_h)->Integral();
        double selected_all = hists->h_all->Integral(1,i_bin);

        eff = selected_i/total_i;
        pur = selected_i/selected_all;

        if (selected_i==0 && selected_all==0) pur = 0;
      }
      else{ // integrate up to the top
        double selected_i = histstoeval.at(i_h)->Integral(i_bin,histstoeval.at(i_h)->GetXaxis()->GetNbins());
        double total_i = histstoeval.at(i_h)->Integral();
        double selected_all = hists->h_all->Integral(i_bin,histstoeval.at(i_h)->GetXaxis()->GetNbins());

        eff = selected_i/total_i;
        pur = selected_i/selected_all;

        if (selected_i==0 && selected_all==0) pur = 0;
      }

      double effpur = eff*pur;
      heff->SetBinContent(i_bin,eff);
      hpur->SetBinContent(i_bin,pur);
      heffpur->SetBinContent(i_bin,effpur);
    }

    heff->SetLineColor(kRed);
    heff->SetMarkerColor(kRed);
    heff->SetMarkerStyle(20);
    heff->SetMarkerSize(.3
    );
    hpur->SetLineColor(kBlue);
    hpur->SetMarkerColor(kBlue);
    hpur->SetMarkerStyle(20);
    hpur->SetMarkerSize(.3
    );
    heffpur->SetLineColor(kBlack);
    heffpur->SetMarkerColor(kBlack);
    heffpur->SetMarkerStyle(20);
    heffpur->SetMarkerSize(.3
    );

    heff->GetYaxis()->SetRangeUser(0,1.1);

    c->cd(i_h+1);
    heff->Draw("p");
    hpur->Draw("same p");
    heffpur->Draw("same p");

    if (i_h==0){
      l->AddEntry(heff,"Efficiency","p");
      l->AddEntry(hpur,"Purity","p");
      l->AddEntry(heffpur,"Efficiency #times Purity","p");
    }

    c->cd(histstoeval.size()+1);
    l->Draw();
  }

}





// --------------------------------------------------- //
// This struct contains signal vs background histograms and all functions related to them

struct histCC1piselEffPur{
  TH1D *h_cc1pi_sel;
  TH1D *h_bg_sel;
  TH1D *h_cc1pi_notsel;
  // TH1D *h_bg_notsel;

  TLegend *l;

  // Constructor for this struct of hists
  histCC1piselEffPur(std::string name, std::string title, double nbins, double binlow, double binhigh){
    h_cc1pi_sel = new TH1D(std::string(name+"_cc1pi_sel").c_str(),title.c_str(),nbins,binlow,binhigh);
    h_bg_sel = new TH1D(std::string(name+"_bg_sel").c_str(),title.c_str(),nbins,binlow,binhigh);
    h_cc1pi_notsel = new TH1D(std::string(name+"_cc1pi_notsel").c_str(),title.c_str(),nbins,binlow,binhigh);
    // h_bg_notsel = new TH1D(std::string(name+"_bg_notsel").c_str(),title.c_str(),nbins,binlow,binhigh);


    // h_cc1pi_sel->SetFillColor(TColor::GetColor(8,64,129));
    // h_cc1pi_all->SetFillColor(TColor::GetColor(8,64,129));
    // h_bg_sel->SetFillColor(TColor::GetColor(197,197,197));
    // h_bg_all->SetFillColor(TColor::GetColor(197,197,197));
    //
    // l = new TLegend(0.59,0.64,0.81,0.87);
    // l->AddEntry(h_p,"True cc1#pi^{+}","f");
    // l->AddEntry(h_mu,"True other","f");
  }
};

  void FillCC1piEffPurHist(histCC1piselEffPur *hists, std::vector<std::vector<double>> value_vec, int plane, std::vector<bool> *pandoraclassedastrack, NuIntTopology topology, bool MIPlow){
    // Loop through all bins in the histograms, and evaluate a cut value at the centre of each bin
    for (int i_bin=1; i_bin < hists->h_cc1pi_sel->GetXaxis()->GetNbins()+1; i_bin++){
      double cutval = hists->h_cc1pi_sel->GetXaxis()->GetBinCenter(i_bin);

      // value_vec is a vector corresponding to the TPCObject. It has one entry per track in the TPCObject, and each of those entries has three entries (one for each plane)
      // Loop through value_vec and calculate number of MIPs for the correct plane for every track
      // Remember that plane=3 is code for (plane0+plane1+plane2)/3
      int n_mips = 0;
      for (size_t i_track=0; i_track<value_vec.size(); i_track++){
        double value = -9999;
        if (plane==3){
          value = (value_vec.at(i_track).at(0)+value_vec.at(i_track).at(1)+value_vec.at(i_track).at(2))/3.;
        }
        else{
          value = value_vec.at(i_track).at(plane);
        }

        // if (value == -9999 || value == -999) continue;
         if (!pandoraclassedastrack->at(i_track)) {};

        if (MIPlow && value < cutval){
          n_mips++;
        }
        else if (!MIPlow && value > cutval){
          n_mips++;
        }
      } // end loop over tracks in TPCObject

      // If event has at least 2 MIPs, it is selected. Fill that into the relevant histograms
      if (n_mips >= 2){
        if (topology == kCC1piplus0p || topology == kCC1piplus1p || topology == kCC1piplusNp){
          hists->h_cc1pi_sel->Fill(cutval);
        }
        else {
          hists->h_bg_sel->Fill(cutval);
        }
      }
      else { // if not at least 2 MIPs
        if (topology == kCC1piplus0p || topology == kCC1piplus1p || topology == kCC1piplusNp){
          hists->h_cc1pi_notsel->Fill(cutval);
        }
        // else {
        //   hists->h_bg_notsel->Fill(cutval);
        // }
      }
    } // End loop over bins in the histograms
  }

void DrawCC1piMCEffPur(TCanvas *c, histCC1piselEffPur *hists){
  TH1D *heff = (TH1D*)hists->h_cc1pi_sel->Clone("heff");
  TH1D *hpur = (TH1D*)hists->h_cc1pi_sel->Clone("hpur");
  TH1D *heffpur = (TH1D*)hists->h_cc1pi_sel->Clone("heffpur");

  heff->Clear();
  hpur->Clear();
  heffpur->Clear();

  TLegend *l = new TLegend(0.59,0.64,0.81,0.87);
  l->SetTextFont(132);
  l->SetLineColor(kWhite);
  l->SetFillColor(kWhite);

  for (int i_bin=1; i_bin < heff->GetXaxis()->GetNbins()+1; i_bin++){
    double selected_cc1pi = hists->h_cc1pi_sel->GetBinContent(i_bin);
    double total_cc1pi = hists->h_cc1pi_sel->GetBinContent(i_bin)+hists->h_cc1pi_notsel->GetBinContent(i_bin);
    double selected_all = hists->h_cc1pi_sel->GetBinContent(i_bin)+hists->h_bg_sel->GetBinContent(i_bin);

    double eff = selected_cc1pi/total_cc1pi;
    double pur = selected_cc1pi/selected_all;

    heff->SetBinContent(i_bin,eff);
    hpur->SetBinContent(i_bin,pur);
    heffpur->SetBinContent(i_bin,eff*pur);
  }

  heff->SetLineColor(kRed);
  heff->SetMarkerColor(kRed);
  heff->SetMarkerStyle(20);
  heff->SetMarkerSize(.3
  );
  hpur->SetLineColor(kBlue);
  hpur->SetMarkerColor(kBlue);
  hpur->SetMarkerStyle(20);
  hpur->SetMarkerSize(.3
  );
  heffpur->SetLineColor(kBlack);
  heffpur->SetMarkerColor(kBlack);
  heffpur->SetMarkerStyle(20);
  heffpur->SetMarkerSize(.3
  );

  heff->GetYaxis()->SetRangeUser(0,1.1);

  c->cd();
  heff->Draw("p");
  hpur->Draw("same p");
  heffpur->Draw("same p");

  l->AddEntry(heff,"Efficiency","p");
  l->AddEntry(hpur,"Purity","p");
  l->AddEntry(heffpur,"Efficiency #times Purity","p");
  l->Draw();

  c->Draw();
  c->Modified();
  c->Update();

}
