#include "StackedHistPDGCode.h"
#include "StackedHistTopology.h"
#include "../Algorithms/TopologyEnums.h"
#include "../Algorithms/PDGEnums.h"
#include "CC1pi_treevars.h"
#include "CC1pi_treevars.cxx"
#include "CC1pi_plotvars_def.h"
#include "CC1pi_cuts.h"
#include "CC1pi_cuts.cxx"
#include "ProtondEdxrrPlots.h"
#include "getSimPot.C"
#include "MakeLocalLinearityPlots.C"

#include <boost/algorithm/string.hpp>

// What variables do we want these plots as a function of?
std::vector<CC1piPlotVars> GetVarstoplot(treevars *vars){
   std::vector<CC1piPlotVars> varstoplot = {
      Var_TPCObj_PFP_lnLmipoverp(vars)
      ,Var_TPCObj_PFP_lnLmipoverp_SecondMIP(vars)
      ,Var_TPCObj_PFP_Lmumipovermumipp(vars)
      ,Var_TPCObj_PFP_lnLmipoverp_SecondMIPcont(vars)
      // ,Var_TPCObj_PFP_Lmumipovermumipp_cont(vars)
      // ,Var_TPCObj_PFP_BrokenTrackAngle(vars)
      // ,Var_TPCObj_PFP_track_residual_mean_low(vars)
      // ,Var_TPCObj_PFP_track_residual_mean_high(vars)
      // ,Var_TPCObj_PFP_track_residual_std_low(vars)
      // ,Var_TPCObj_PFP_track_residual_std_high(vars)
      ,Var_TPCObj_PFP_track_Chi2Proton_SecondMIP(vars)
      ,Var_TPCObj_PFP_track_Chi2Proton_SecondMIPcont(vars)
      ,Var_TPCObj_AngleBetweenMIPs(vars)
      ,Var_TPCObj_PFP_track_perc_used_hits(vars)
      ,Var_TPCObj_PFP_track_perc_used_hits_SecondMIP(vars)
      ,Var_TPCObj_PFP_VtxTrackDist(vars)
      ,Var_TPCObj_PFP_VtxTrackDist_SecondMIP(vars)
      ,Var_TPCObj_AngleBetweenMIPs_mumupairs(vars)
      ,Var_TPCObj_PFP_VtxTrackDist_mumupairs(vars)
      ,Var_TPCObj_PFP_VtxTrackDist_SecondMIP_mumupairs(vars)
      ,Var_TPCObj_PFP_track_length_LeadingMIP_mumupairs(vars)
      ,Var_TPCObj_PFP_track_length_SecondMIP_mumupairs(vars)
      ,Var_TPCObj_PFP_track_theta_SecondMIP(vars)
      ,Var_TPCObj_PFP_track_phi_SecondMIP(vars)
      ,Var_TPCObj_PFP_track_dEdx_stddev_SecondMIP(vars)
      ,Var_TPCObj_PFP_track_dEdx_stddev_start_SecondMIP(vars)
      ,Var_TPCObj_PFP_track_dEdx_mean_start_SecondMIP(vars)
      ,Var_TPCObj_PFP_track_dEdx_truncmoverm_start_SecondMIP(vars)
      // ,Var_TPCObj_PFP_isContained(vars)
      ,Var_TPCObj_PFP_track_dEdx_truncmean_start(vars)
      // ,Var_TPCObj_DaughterTracks_Order_dEdxtr(vars)
      // ,Var_TPCObj_DaughterTracks_Order_dEdxtr_selMIPs(vars)
      // ,Var_TPCObj_DaughterTracks_Order_trklen_selMIPs(vars)

      // Var_TPCObj_PFP_MCSLLpiMinusLLp_LeadingMIP(vars)
      ,Var_TPCObj_PFP_MCSLLpiMinusLLp_SecondMIP(vars)
      // ,Var_TPCObj_PFP_MCSLLpiMinusLLp_NotPion(vars)
      // ,Var_TPCObj_PFP_MCSLLpiMinusLLp_LeadingMIPcont(vars)
      ,Var_TPCObj_PFP_MCSLLpiMinusLLp_SecondMIPcont(vars)
      // ,Var_TPCObj_PFP_MCSLLpiMinusLLp_NotPioncont(vars)
      // ,Var_TPCObj_PFP_track_MomRangeMinusMCS_p_LeadingMIPcont(vars)
      ,Var_TPCObj_PFP_track_MomRangeMinusMCS_p_SecondMIPcont(vars)
      // ,Var_TPCObj_PFP_track_MomRangeMinusMCS_p_NotPioncont(vars)
      // ,Var_TPCObj_PFP_track_MomRangeMinusMCS_mu_LeadingMIPcont(vars)
      // ,Var_TPCObj_PFP_track_MomRangeMinusMCS_mu_SecondMIPcont(vars)
      // ,Var_TPCObj_PFP_track_MomRangeMinusMCS_mu_SecondMIPcont_mumupairs(vars)
      // ,Var_TPCObj_PFP_track_MomRangeMinusMCS_mu_NotPioncont(vars)
      // ,Var_TPCObj_PFP_track_nhits_LeadingMIP(vars)
      // ,Var_TPCObj_PFP_track_nhits_SecondMIP(vars)
      ,Var_TPCObj_PFP_track_MCSpi_maxScatter_SecondMIP(vars)
      ,Var_TPCObj_PFP_track_MCSp_maxScatter_SecondMIP(vars)
      ,Var_TPCObj_PFP_track_MCSpi_meanScatter_SecondMIP(vars)
      ,Var_TPCObj_PFP_track_MCSp_meanScatter_SecondMIP(vars)
      ,Var_TPCObj_PFP_track_MCSpi_maxScatter_SecondMIP_mumupairs(vars)
      ,Var_TPCObj_PFP_track_MCSp_maxScatter_SecondMIP_mumupairs(vars)
      ,Var_TPCObj_PFP_track_MCSpi_meanScatter_SecondMIP_mumupairs(vars)
      ,Var_TPCObj_PFP_track_MCSp_meanScatter_SecondMIP_mumupairs(vars)
      // ,Var_TPCObj_LeadingMIP_PandoraClassedAsTrack(vars)
      // ,Var_TPCObj_PFP_trueKE_selMIPs(vars)
      // ,Var_TPCObj_PFP_trueKE_SecondMIP(vars)
      ,Var_TPCObj_PFP_trueEndP_SecondMIP(vars)
      // ,Var_TPCObj_NotMIPs_isContained(vars)
      // ,Var_TPCObj_PFP_ndaughters_SecondMIP(vars)
      // ,Var_TPCObj_PFP_ndaughters_SecondMIP_mumupairs(vars)
      // ,Var_TPCObj_PFP_MCP_PDG_mTruePDG(vars)
      // ,Var_TPCObj_PFP_MCP_numdaughters_SecondMIP(vars)
      ,Var_TPCObj_PFP_MCP_motherIDeq0_SecondMIP(vars)
      ,Var_TPCObj_dEdx_truncmean_MIPdiff(vars)
      ,Var_TPCObj_dEdx_truncmean_MIPdiff_mupi(vars)
      ,Var_TPCObj_dEdx_truncmean_MIPdiff_muproton(vars)
      ,Var_TPCObj_dEdx_truncmean_MIPdiff_other(vars)
      ,Var_TPCObj_PFP_track_dEdx_nhits(vars)
      //,Var_TPCObj_PFP_track_BDTscore_contained(vars)
      //,Var_TPCObj_PFP_track_BDTscore_uncontained(vars)
      ,Var_TPCObj_PFP_track_passesMIPcut(vars)
      //,Var_TPCObj_NDaughterPFPs(vars)
   };
   return varstoplot;
}


// What variables do we want these 2D plots as a function of?
std::vector<std::pair<CC1piPlotVars,CC1piPlotVars>> GetVarstoplot2D(treevars *vars){
   std::vector<std::pair<CC1piPlotVars,CC1piPlotVars>> varstoplot2D = {
      // {Var_TPCObj_PFP_track_dEdx_stddev_SecondMIP(vars),Var_TPCObj_PFP_track_dEdx_truncmean_start(vars)}
      // ,{Var_TPCObj_PFP_track_dEdx_mean_start_SecondMIP(vars),Var_TPCObj_PFP_track_dEdx_truncmean_start(vars)}
      // ,{Var_TPCObj_PFP_track_dEdx_stddev_start_SecondMIP(vars),Var_TPCObj_PFP_track_dEdx_truncmean_start(vars)}
      // ,{Var_TPCObj_PFP_track_theta_SecondMIP(vars),Var_TPCObj_PFP_track_dEdx_truncmean_start(vars)}
      {Var_TPCObj_PFP_track_theta(vars),Var_TPCObj_PFP_track_dEdx_truncmean_start(vars)}
      // {Var_TPCObj_PFP_VtxTrackDist_mumupairs(vars),Var_TPCObj_PFP_VtxTrackDist_SecondMIP_mumupairs(vars)}
      // {Var_TPCObj_PFP_MCSLLpiMinusLLp_NotPion(vars),Var_TPCObj_PFP_track_length(vars)}
      // ,{Var_TPCObj_PFP_MCSLLpiMinusLLp_NotPioncont(vars),Var_TPCObj_PFP_track_length(vars)}
      // ,{Var_TPCObj_PFP_track_MomRangeMinusMCS_p_NotPioncont(vars),Var_TPCObj_PFP_track_length(vars)}
      // ,{Var_TPCObj_PFP_track_MomRangeMinusMCS_mu_NotPioncont(vars),Var_TPCObj_PFP_track_length(vars)}
      // ,{Var_TPCObj_PFP_track_SimpleCluster_hitLinearity_minimum_NotPioncont(vars),Var_TPCObj_PFP_track_length(vars)}
      // ,{Var_TPCObj_PFP_track_SimpleCluster_hitLinearity_mean_NotPioncont(vars),Var_TPCObj_PFP_track_length(vars)}
      // ,{Var_TPCObj_PFP_track_SimpleCluster_hitLinearity_truncated_mean_NotPioncont(vars),Var_TPCObj_PFP_track_length(vars)}
      // ,{Var_TPCObj_PFP_track_SimpleCluster_hitLinearity_mean_NotPioncont(vars),Var_TPCObj_PFP_track_SimpleCluster_hitLinearity_minimum_NotPioncont(vars)}
   };
   return varstoplot2D;
}



// ---------------------------------------------------- //
//  Now the function starts
// ---------------------------------------------------- //

void MakeCC1piPlots(std::string mcfile, double POTscaling=0., std::string onbeamdatafile="", std::string offbeamdatafile="", double offbeamscaling=0., bool onminusoffbeam=false, bool MakeKinkFindingPlots=false, bool MakeProtonPlots=false){

  // Note: MVA trees are made for MC only
   TFile f_MVA("MVA_Trees.root", "RECREATE");
   MVAvars *MVA_vars = new MVAvars();
   TTree *muon = new TTree("muon","muon");
   TTree *pion = new TTree("pion","pion");
   TTree *background = new TTree("background","background");
   muon -> Branch("dEdx_truncmean_start", &(MVA_vars->dEdx_truncmean_start));
   muon -> Branch("VtxTrackDist", &(MVA_vars->VtxTrackDist));
   muon -> Branch("nhits", &(MVA_vars->nhits));
   muon -> Branch("lnLmipoverp", &(MVA_vars->lnLmipoverp));
   pion -> Branch("dEdx_truncmean_start", &(MVA_vars->dEdx_truncmean_start));
   pion -> Branch("VtxTrackDist", &(MVA_vars->VtxTrackDist));
   pion -> Branch("nhits", &(MVA_vars->nhits));
   pion -> Branch("lnLmipoverp", &(MVA_vars->lnLmipoverp));
   background -> Branch("dEdx_truncmean_start", &(MVA_vars->dEdx_truncmean_start));
   background -> Branch("VtxTrackDist", &(MVA_vars->VtxTrackDist));
   background -> Branch("nhits", &(MVA_vars->nhits));
   background -> Branch("lnLmipoverp", &(MVA_vars->lnLmipoverp));

   gStyle->SetTitleX(0.5);
   gStyle->SetTitleAlign(23);
   gStyle->SetOptTitle(1);
   gStyle->SetTitleBorderSize(0.);
   gStyle->SetOptStat(0);

   TFile *f_bnbcos = new TFile(mcfile.c_str(), "read");
   TTree *t_bnbcos = (TTree*)f_bnbcos->Get("cc1piselec/outtree");

    if (!t_bnbcos){
      std::cout << "Error: did not find tree cc1piselec/outtree in file" << std::endl;
      return;
    }

   treevars mc_vars;
   settreevars(t_bnbcos,&mc_vars);

   std::string containedBookMVAType = "BDTG";
   std::string containedBookMVALoc = "/uboone/app/users/ddevitt/LArSoft_v06_26_01_14_uboonecode_v06_26_01_22/srcs/uboonecode/uboone/CC1pi/MVA/dataset_contained/weights/TMVAClassification_BDTG.weights.xml";
   std::string uncontainedBookMVAType = "BDTG";
   std::string uncontainedBookMVALoc = "/uboone/app/users/ddevitt/LArSoft_v06_26_01_14_uboonecode_v06_26_01_22/srcs/uboonecode/uboone/CC1pi/MVA/dataset_uncontained/weights/TMVAClassification_BDTG.weights.xml";

   TMVA::Reader fReader_contained("");
   fReader_contained.AddVariable("dEdx_truncmean_start", &(mc_vars.float_dEdx_truncmean_start));
   fReader_contained.AddVariable("VtxTrackDist", &(mc_vars.float_VtxTrackDist));
   fReader_contained.AddVariable("nhits", &(mc_vars.float_nhits));
   fReader_contained.AddVariable("lnLmipoverp", &(mc_vars.float_lnLmipoverp));
   fReader_contained.BookMVA(containedBookMVAType.c_str(), containedBookMVALoc.c_str());

   TMVA::Reader fReader_uncontained("");
   fReader_uncontained.AddVariable("dEdx_truncmean_start", &(mc_vars.float_dEdx_truncmean_start));
   fReader_uncontained.AddVariable("VtxTrackDist", &(mc_vars.float_VtxTrackDist));
   fReader_uncontained.AddVariable("nhits", &(mc_vars.float_nhits));
   fReader_uncontained.AddVariable("lnLmipoverp", &(mc_vars.float_lnLmipoverp));
   fReader_uncontained.BookMVA(uncontainedBookMVAType.c_str(), uncontainedBookMVALoc.c_str());

   TFile *f_onbeam=nullptr;
   TTree *t_onbeam=nullptr;
   treevars onbeam_vars;
   TMVA::Reader fReader_onbeam_contained("");
   TMVA::Reader fReader_onbeam_uncontained("");
   if (onbeamdatafile!=""){
     std::cout << "Making data-MC comparisons" << std::endl;
     f_onbeam = new TFile(onbeamdatafile.c_str(), "read");
     t_onbeam = (TTree*)f_onbeam->Get("cc1piselec/outtree");
     settreevars(t_onbeam,&onbeam_vars);

     fReader_onbeam_contained.AddVariable("dEdx_truncmean_start", &(onbeam_vars.float_dEdx_truncmean_start));
     fReader_onbeam_contained.AddVariable("VtxTrackDist", &(onbeam_vars.float_VtxTrackDist));
     fReader_onbeam_contained.AddVariable("nhits", &(onbeam_vars.float_nhits));
     fReader_onbeam_contained.AddVariable("lnLmipoverp", &(onbeam_vars.float_lnLmipoverp));
     fReader_onbeam_contained.BookMVA(containedBookMVAType.c_str(), containedBookMVALoc.c_str());

     fReader_onbeam_uncontained.AddVariable("dEdx_truncmean_start", &(onbeam_vars.float_dEdx_truncmean_start));
     fReader_onbeam_uncontained.AddVariable("VtxTrackDist", &(onbeam_vars.float_VtxTrackDist));
     fReader_onbeam_uncontained.AddVariable("nhits", &(onbeam_vars.float_nhits));
     fReader_onbeam_uncontained.AddVariable("lnLmipoverp", &(onbeam_vars.float_lnLmipoverp));
     fReader_onbeam_uncontained.BookMVA(uncontainedBookMVAType.c_str(), uncontainedBookMVALoc.c_str());
   }

   TFile *f_offbeam=nullptr;
   TTree *t_offbeam=nullptr;
   treevars offbeam_vars;
   TMVA::Reader fReader_offbeam_contained("");
   TMVA::Reader fReader_offbeam_uncontained("");
   if (offbeamdatafile!=""){
     std::cout << "Making data-MC comparisons" << std::endl;
     f_offbeam = new TFile(offbeamdatafile.c_str(), "read");
     t_offbeam = (TTree*)f_offbeam->Get("cc1piselec/outtree");
     settreevars(t_offbeam,&offbeam_vars);

     fReader_offbeam_contained.AddVariable("dEdx_truncmean_start", &(offbeam_vars.float_dEdx_truncmean_start));
     fReader_offbeam_contained.AddVariable("VtxTrackDist", &(offbeam_vars.float_VtxTrackDist));
     fReader_offbeam_contained.AddVariable("nhits", &(offbeam_vars.float_nhits));
     fReader_offbeam_contained.AddVariable("lnLmipoverp", &(offbeam_vars.float_lnLmipoverp));
     fReader_offbeam_contained.BookMVA(containedBookMVAType.c_str(), containedBookMVALoc.c_str());

     fReader_offbeam_uncontained.AddVariable("dEdx_truncmean_start", &(offbeam_vars.float_dEdx_truncmean_start));
     fReader_offbeam_uncontained.AddVariable("VtxTrackDist", &(offbeam_vars.float_VtxTrackDist));
     fReader_offbeam_uncontained.AddVariable("nhits", &(offbeam_vars.float_nhits));
     fReader_offbeam_uncontained.AddVariable("lnLmipoverp", &(offbeam_vars.float_lnLmipoverp));
     fReader_offbeam_uncontained.BookMVA(uncontainedBookMVAType.c_str(), uncontainedBookMVALoc.c_str());
   }

   ofstream evdinfo;
   evdinfo.open("evdinfo.txt");

   // Sanity check: the plot vectors should be the same size
   t_bnbcos->GetEntry(0);
   Calcvars(&mc_vars, &fReader_contained, &fReader_uncontained);
   std::vector<CC1piPlotVars> varstoplot_dummy = GetVarstoplot(&mc_vars);
   std::vector<std::pair<CC1piPlotVars,CC1piPlotVars>> varstoplot2D_dummy = GetVarstoplot2D(&mc_vars);
   std::vector<CC1piPlotVars> cutvars_dummy = GetCutVars(&mc_vars);

   // ----------------- MC

   // Make histograms to fill
   const size_t nplots = varstoplot_dummy.size();
   StackedHistPDGCode *mc_hists_cc1pi_pdg_beforecuts[nplots];
   StackedHistTopology *mc_hists_cc1pi_top_beforecuts[nplots];
   StackedHistPDGCode *mc_hists_cc1pi_pdg_aftercuts[nplots];
   StackedHistTopology *mc_hists_cc1pi_top_aftercuts[nplots];
   for (size_t i_h=0; i_h<nplots; i_h++){
     std::string histtitle_i = varstoplot_dummy.at(i_h).histtitle;
     std::string histname_i = varstoplot_dummy.at(i_h).histname;
     std::vector<double> bins_i = varstoplot_dummy.at(i_h).bins;

      mc_hists_cc1pi_pdg_beforecuts[i_h] = new StackedHistPDGCode(std::string("hCC1pi_PDG_beforecuts_")+histname_i,histtitle_i,bins_i.at(0),bins_i.at(1),bins_i.at(2));
      mc_hists_cc1pi_top_beforecuts[i_h] = new StackedHistTopology(std::string("hCC1pi_Top_beforecuts_")+histname_i,histtitle_i,bins_i.at(0),bins_i.at(1),bins_i.at(2));

      mc_hists_cc1pi_pdg_aftercuts[i_h] = new StackedHistPDGCode(std::string("hCC1pi_PDG_aftercuts_")+histname_i,histtitle_i,bins_i.at(0),bins_i.at(1),bins_i.at(2));
      mc_hists_cc1pi_top_aftercuts[i_h] = new StackedHistTopology(std::string("hCC1pi_Top_aftercuts_")+histname_i,histtitle_i,bins_i.at(0),bins_i.at(1),bins_i.at(2));
  }

  // Make 2Dhistograms to fill
  const size_t nplots2D = varstoplot2D_dummy.size();
  StackedHistPDGCode *mc_hists2D_cc1pi_pdg_beforecuts[nplots2D];
  StackedHistTopology *mc_hists2D_cc1pi_top_beforecuts[nplots2D];
  StackedHistPDGCode *mc_hists2D_cc1pi_pdg_aftercuts[nplots2D];
  StackedHistTopology *mc_hists2D_cc1pi_top_aftercuts[nplots2D];
  for (size_t i_h=0; i_h<nplots2D; i_h++){
    std::string histname_x = varstoplot2D_dummy.at(i_h).first.histname;
    std::vector<double> bins_x = varstoplot2D_dummy.at(i_h).first.bins;
    std::string histname_y = varstoplot2D_dummy.at(i_h).second.histname;
    std::vector<double> bins_y = varstoplot2D_dummy.at(i_h).second.bins;

    std::vector<std::string> scratch_x;
    boost::split(scratch_x,varstoplot2D_dummy.at(i_h).first.histtitle,[](char c){return c == ';';});
    std::vector<std::string> scratch_y;
    boost::split(scratch_y,varstoplot2D_dummy.at(i_h).second.histtitle,[](char c){return c == ';';});

    mc_hists2D_cc1pi_pdg_beforecuts[i_h] = new StackedHistPDGCode(std::string("hCC1pi2D_PDG_beforecuts_")+histname_x+histname_y,std::string(";")+scratch_x.at(1)+std::string(";")+scratch_y.at(1),bins_x.at(0),bins_x.at(1),bins_x.at(2),bins_y.at(0),bins_y.at(1),bins_y.at(2));
    // mc_hists2D_cc1pi_top_beforecuts[i_h] = new StackedHistTopology(std::string("hCC1pi2D_Top_beforecuts_")+histname_x+histname_y,std::string(";")+scratch_x.at(1)+std::string(";")+scratch_y.at(1),bins_x.at(0),bins_x.at(1),bins_x.at(2),bins_y.at(0),bins_y.at(1),bins_y.at(2));

    mc_hists2D_cc1pi_pdg_aftercuts[i_h] = new StackedHistPDGCode(std::string("hCC1pi2D_PDG_aftercuts_")+histname_x+histname_y,std::string(";")+scratch_x.at(1)+std::string(";")+scratch_y.at(1),bins_x.at(0),bins_x.at(1),bins_x.at(2),bins_y.at(0),bins_y.at(1),bins_y.at(2));
    // mc_hists2D_cc1pi_top_aftercuts[i_h] = new StackedHistTopology(std::string("hCC1pi2D_Top_aftercuts_")+histname_x+histname_y,std::string(";")+scratch_x.at(1)+std::string(";")+scratch_y.at(1),bins_x.at(0),bins_x.at(1),bins_x.at(2),bins_y.at(0),bins_y.at(1),bins_y.at(2));
}

// Custom 2D histograms: leading vs second MIP true PDG
std::vector<PDGCode> MIPpdgs = {kMuMinus, kMuPlus, kPiMinus, kPiPlus, kProton};
int nMIPpdgs = MIPpdgs.size()+1;
TH2D *selMIPs2D = new TH2D("selMIPs2D","Selected MIP-like tracks;Leading MIP; Second MIP",nMIPpdgs,0,nMIPpdgs,nMIPpdgs,0,nMIPpdgs);
for (size_t i_bin=1; i_bin<nMIPpdgs+1; i_bin++){
   if (i_bin==nMIPpdgs){
      selMIPs2D->GetXaxis()->SetBinLabel(i_bin,"Other");
      selMIPs2D->GetYaxis()->SetBinLabel(i_bin,"Other");
      continue;
   }
   selMIPs2D->GetXaxis()->SetBinLabel(i_bin,PDGenum2str(MIPpdgs.at(i_bin-1)).c_str());
   selMIPs2D->GetYaxis()->SetBinLabel(i_bin,PDGenum2str(MIPpdgs.at(i_bin-1)).c_str());
}

   std::cout << std::endl << "--- Considering the following 1D variables --- " << std::endl;
   for (size_t i_var=0; i_var < nplots; i_var++){
         std::cout << varstoplot_dummy.at(i_var).histtitle.substr(1,varstoplot_dummy.at(i_var).histtitle.size()-2) << std::endl;
   }
   std::cout << "------" << std::endl << std::endl;

   std::cout << "--- Considering the following 2D variables --- " << std::endl;
   for (size_t i_var=0; i_var < nplots2D; i_var++){

      std::vector<std::string> scratch_x;
      boost::split(scratch_x,varstoplot2D_dummy.at(i_var).first.histtitle,[](char c){return c == ';';});
      std::vector<std::string> scratch_y;
      boost::split(scratch_y,varstoplot2D_dummy.at(i_var).second.histtitle,[](char c){return c == ';';});

      std::cout << scratch_x.at(1) << " vs " << scratch_y.at(1) << std::endl;
   }
   std::cout << "------" << std::endl << std::endl;

    std::cout << "--- Applying the following cuts --- " << std::endl;
    for (size_t i_var=0; i_var < cutvars_dummy.size(); i_var++){
          std::cout << cutvars_dummy.at(i_var).histtitle.substr(1,cutvars_dummy.at(i_var).histtitle.size()-2);
          if (cutvars_dummy.at(i_var).KeepBelowCut) std::cout << " < ";
          else std::cout << " > ";
          std::cout << cutvars_dummy.at(i_var).CutValue << std::endl;
          if (cutvars_dummy.at(i_var).histname == "isMIP"){
            std::vector<CC1piPlotVars> cutvars_dummy_MIPs = GetMIPCutVars(&mc_vars);
            for (size_t i_mipcut=0; i_mipcut<cutvars_dummy_MIPs.size(); i_mipcut++){
              std::cout << "  -> " << cutvars_dummy_MIPs.at(i_mipcut).histtitle.substr(1,cutvars_dummy_MIPs.at(i_mipcut).histtitle.size()-2);
              if (cutvars_dummy_MIPs.at(i_mipcut).KeepBelowCut) std::cout << " < ";
              else std::cout << " > ";
              std::cout << cutvars_dummy_MIPs.at(i_mipcut).CutValue << std::endl;
            }
          }
    }
    std::cout << "------" << std::endl << std::endl;

   // Dummy plot in order to get percentage of different topologies that are selected
   StackedHistTopology *SelectedEvents = new StackedHistTopology("SelectedEvents",";;",1,0,1);
   StackedHistTopology *AllEvents = new StackedHistTopology("AllEvents",";;",1,0,1);
   StackedHistTopology *NDaughters = new StackedHistTopology("NDaughters",";Number of #nu daughter PFPs in selected TPCObject;",9,1,10);

   int protonplotsmade=0;

   // Loop through MC tree and fill plots
   for (int i = 0; i < /*t_bnbcos->GetEntries()*/100; i++){
     if (i%1000==0) std::cout << "MC: " << i << "/" << t_bnbcos->GetEntries() << std::endl;
      t_bnbcos->GetEntry(i);

      Calcvars(&mc_vars, &fReader_contained, &fReader_uncontained);
      MakeMVATrees(muon, pion, background, MVA_vars, &mc_vars);
      std::vector<CC1piPlotVars> Varstoplot = GetVarstoplot(&mc_vars);
      std::vector<std::pair<CC1piPlotVars,CC1piPlotVars>> Varstoplot2D = GetVarstoplot2D(&mc_vars);

      AllEvents->Fill((NuIntTopology)mc_vars.Truth_topology,0.5);
      NDaughters->Fill((NuIntTopology)mc_vars.Truth_topology,mc_vars.TPCObj_NDaughterPFPs->at(0));

      bool isSelected = IsEventSelected(&mc_vars);


      // For all events, make hit dQds and local linearity plots for true pi+ and mu tracks
      for (int i_tr=0; i_tr<mc_vars.TPCObj_PFP_isDaughter->size(); i_tr++){
         // Do this for the first 100 events only because otherwise we'll just have a ridiculous number of plots
         if (MakeKinkFindingPlots && i<50){
           if (mc_vars.TPCObj_PFP_track_SpacepointsXYZ->at(i_tr).size()==0) continue;
           if (!(mc_vars.TPCObj_PFP_truePDG->at(i_tr)==13 ||mc_vars.TPCObj_PFP_truePDG->at(i_tr)==211)) continue;

           TCanvas *c1 = new TCanvas("c1","",500,500);
           bool saveplot = PlotLocalLinearityDetails(10,&mc_vars,c1,i_tr,isSelected);

           if (!saveplot) break;

           TString savename = TString::Format("./linearity_evt%d_pfp%d.pdf",i,i_tr);
           c1->Print(savename.Data());

           delete c1;
         }
      }

      if(isSelected){
         SelectedEvents->Fill((NuIntTopology)mc_vars.Truth_topology,0.5);

         PDGCode LeadingMIPpdg = PDGCode(mc_vars.TPCObj_PFP_truePDG->at(mc_vars.TPCObj_LeadingMIPtrackIndex));
         PDGCode SecondMIPpdg = PDGCode(mc_vars.TPCObj_PFP_truePDG->at(mc_vars.TPCObj_SecondMIPtrackIndex));

         int LeadingMIPbin = nMIPpdgs;
         int SecondMIPbin = nMIPpdgs;
         for (int i_bin=0; i_bin<MIPpdgs.size(); i_bin++){
            if (MIPpdgs.at(i_bin) == LeadingMIPpdg) LeadingMIPbin = i_bin+1;
            if (MIPpdgs.at(i_bin) == SecondMIPpdg) SecondMIPbin = i_bin+1;
         }
         selMIPs2D->Fill((double)LeadingMIPbin-0.5,(double)SecondMIPbin-0.5);



         // For selected events, make proton plots
         for (int i_tr=0; i_tr<mc_vars.TPCObj_PFP_isDaughter->size(); i_tr++){

            if (MakeProtonPlots && protonplotsmade<30 && i_tr == mc_vars.TPCObj_SecondMIPtrackIndex && mc_vars.TPCObj_PFP_truePDG->at(i_tr) == 2212){
               TCanvas *c0 = new TCanvas();
               MakeProtondEdxrrPlot_SingleTrack(c0, &mc_vars, i_tr);
               c0->Print(std::string(std::string("ProtonPlots_TPCObj")+std::to_string(i)+std::string("_track")+std::to_string(i_tr)+std::string(".png")).c_str());
               protonplotsmade++;
            }
         }


         //Record info for event displays
         if(SecondMIPpdg==kProton) {
            evdinfo << "Run/Subrun/Event: " << mc_vars.run_num << " " << mc_vars.subrun_num << " " << mc_vars.event_num << std::endl;
            evdinfo << "True topology: " << topologyenum2str(mc_vars.Truth_topology) << std::endl;
            evdinfo << "Reco nu vertex position (x,y,z): " << mc_vars.TPCObj_reco_vtx->at(0) << " " << mc_vars.TPCObj_reco_vtx->at(1) << " " << mc_vars.TPCObj_reco_vtx->at(2) << std::endl;
            for(int PFP = 0; PFP < mc_vars.TPCObj_PFP_truePDG->size(); PFP++) {
               evdinfo << "PDG: " << mc_vars.TPCObj_PFP_truePDG->at(PFP) << " – ";
               evdinfo << "Start (z,x): " << mc_vars.TPCObj_PFP_track_start->at(PFP).at(2) << " " << mc_vars.TPCObj_PFP_track_start->at(PFP).at(0) << " – ";
               evdinfo << "End (z,x): " << mc_vars.TPCObj_PFP_track_end->at(PFP).at(2) << " " << mc_vars.TPCObj_PFP_track_end->at(PFP).at(0);
               evdinfo << std::endl;
            }
            evdinfo << std::endl;
         }

      } // end if (isSelected)

      // Fill 1D plots
      // Loop over Varstoplot
      for (size_t i_h = 0; i_h < nplots; i_h++){
        std::vector<double> vartoplot = *(Varstoplot.at(i_h).Var);
         // Loop over tracks
         for (size_t i_tr = 0; i_tr < vartoplot.size(); i_tr++){

            // if (!(mc_vars.TPCObj_PFP_isDaughter->at(i_tr) && bool(mc_vars.TPCObj_PFP_track_passesMIPcut->at(i_tr)))) continue;

            if (!(FillPlotForTrack(&(Varstoplot.at(i_h)), &mc_vars, i_tr))) continue;

            mc_hists_cc1pi_pdg_beforecuts[i_h]->Fill((PDGCode)mc_vars.TPCObj_PFP_truePDG->at(i_tr),vartoplot.at(i_tr));
            mc_hists_cc1pi_top_beforecuts[i_h]->Fill((NuIntTopology)mc_vars.Truth_topology,vartoplot.at(i_tr),1.0/vartoplot.size());

            if(isSelected) {
               mc_hists_cc1pi_pdg_aftercuts[i_h]->Fill((PDGCode)mc_vars.TPCObj_PFP_truePDG->at(i_tr),vartoplot.at(i_tr));
               mc_hists_cc1pi_top_aftercuts[i_h]->Fill((NuIntTopology)mc_vars.Truth_topology,vartoplot.at(i_tr),1.0/vartoplot.size());
            } // end if(isSelected)
         } // end loop over tracks
      } // end loop over 1D Varstoplot

      // Fill 2D plots
      // Loop over Varstoplot2D
      for (size_t i_h = 0; i_h < nplots2D; i_h++){
        std::vector<double> vartoplot_x = *(Varstoplot2D.at(i_h).first.Var);
         std::vector<double> vartoplot_y = *(Varstoplot2D.at(i_h).second.Var);
         // Loop over tracks
         for (size_t i_tr = 0; i_tr < vartoplot_x.size(); i_tr++){
            // if (!(mc_vars.TPCObj_PFP_isDaughter->at(i_tr) && bool(mc_vars.TPCObj_PFP_track_passesMIPcut->at(i_tr)))) continue;

            if (!(FillPlotForTrack(&(Varstoplot2D.at(i_h).first), &mc_vars, i_tr) && FillPlotForTrack(&(Varstoplot2D.at(i_h).second), &mc_vars, i_tr))) continue;

            mc_hists2D_cc1pi_pdg_beforecuts[i_h]->Fill2D((PDGCode)mc_vars.TPCObj_PFP_truePDG->at(i_tr),vartoplot_x.at(i_tr),vartoplot_y.at(i_tr));
            // mc_hists2D_cc1pi_top_beforecuts[i_h]->Fill2D((NuIntTopology)mc_vars.Truth_topology,vartoplot_x.at(i_tr),vartoplot_y.at(i_tr));

            if(isSelected) {
               mc_hists2D_cc1pi_pdg_aftercuts[i_h]->Fill2D((PDGCode)mc_vars.TPCObj_PFP_truePDG->at(i_tr),vartoplot_x.at(i_tr),vartoplot_y.at(i_tr));
               // mc_hists2D_cc1pi_top_aftercuts[i_h]->Fill2D((NuIntTopology)mc_vars.Truth_topology,vartoplot_x.at(i_tr),vartoplot_y.at(i_tr));
            } // end if(isSelected)
         } // end loop over tracks
      } // end loop over 2D Varstoplot
   } // end loop over entries in tree

   evdinfo.close();

   // ----------------- Data: on-beam

   int nsel_onbeam=0;
   // Make histograms to fill
   TH1F *onb_hists_cc1pi_pdg_beforecuts[nplots];
   TH1F *onb_hists_cc1pi_top_beforecuts[nplots];
   TH1F *onb_hists_cc1pi_pdg_aftercuts[nplots];
   TH1F *onb_hists_cc1pi_top_aftercuts[nplots];
   if (t_onbeam){
      for (size_t i_h=0; i_h<nplots; i_h++){
        std::string histtitle_i = varstoplot_dummy.at(i_h).histtitle;
        std::string histname_i = varstoplot_dummy.at(i_h).histname;
        std::vector<double> bins_i = varstoplot_dummy.at(i_h).bins;

         onb_hists_cc1pi_pdg_beforecuts[i_h] = new TH1F(TString::Format("hCC1pi_onbeam_PDG_beforecuts_%s",histname_i.c_str()).Data(),histtitle_i.c_str(),bins_i.at(0),bins_i.at(1),bins_i.at(2));
         onb_hists_cc1pi_top_beforecuts[i_h] = new TH1F(TString::Format("hCC1pi_onbeam_Top_beforecuts_%s",histname_i.c_str()).Data(),histtitle_i.c_str(),bins_i.at(0),bins_i.at(1),bins_i.at(2));

         onb_hists_cc1pi_pdg_aftercuts[i_h] = new TH1F(TString::Format("hCC1pi_onbeam_PDG_aftercuts_%s",histname_i.c_str()).Data(),histtitle_i.c_str(),bins_i.at(0),bins_i.at(1),bins_i.at(2));
         onb_hists_cc1pi_top_aftercuts[i_h] = new TH1F(TString::Format("hCC1pi_onbeam_Top_aftercuts_%s",histname_i.c_str()).Data(),histtitle_i.c_str(),bins_i.at(0),bins_i.at(1),bins_i.at(2));
     }

     // Loop through on-beam data tree and fill plots
     for (int i = 0; i < t_onbeam->GetEntries(); i++){
       if (i%1000==0) std::cout << "Onbeam: " << i << "/" << t_onbeam->GetEntries() << std::endl;
        t_onbeam->GetEntry(i);
        Calcvars(&onbeam_vars, &fReader_onbeam_contained, &fReader_onbeam_uncontained);
        std::vector<CC1piPlotVars> Varstoplot = GetVarstoplot(&onbeam_vars);

        bool isSelected = IsEventSelected(&onbeam_vars);

        if (isSelected) nsel_onbeam++;

        // For all events, make hit dQds and local linearity plots for true pi+ and mu tracks
        /*for (int i_tr=0; i_tr<onbeam_vars.TPCObj_PFP_isDaughter->size(); i_tr++){
           // Do this for the first 100 events only because otherwise we'll just have a ridiculous number of plots
           if (MakeKinkFindingPlots && i<50 && isSelected){
             if (onbeam_vars.TPCObj_PFP_track_SpacepointsXYZ->at(i_tr).size()==0) continue;

             TCanvas *c1 = new TCanvas("c1","",400,500);
             bool saveplot = PlotLocalLinearityDetails(10,&onbeam_vars,c1,i_tr,isSelected);

             if (!saveplot) break;

             TString savename = TString::Format("./linearity_bnbdata_evt%d_pfp%d.pdf",i,i_tr);
             c1->Print(savename.Data());

             delete c1;
           }
        }*/


        // Fill 1D plots
        // Loop over Varstoplot
        for (size_t i_h = 0; i_h < nplots; i_h++){
          std::vector<double> vartoplot = *(Varstoplot.at(i_h).Var);
           // Loop over tracks
           for (size_t i_tr = 0; i_tr < vartoplot.size(); i_tr++){

              // if (!(onbeam_vars.TPCObj_PFP_isDaughter->at(i_tr) && bool(onbeam_vars.TPCObj_PFP_track_passesMIPcut->at(i_tr)))) continue;

              if (!(FillPlotForTrack(&(Varstoplot.at(i_h)), &onbeam_vars, i_tr))) continue;

              onb_hists_cc1pi_pdg_beforecuts[i_h]->Fill(vartoplot.at(i_tr));
              onb_hists_cc1pi_top_beforecuts[i_h]->Fill(vartoplot.at(i_tr),1.0/vartoplot.size());

              if(isSelected) {
                 onb_hists_cc1pi_pdg_aftercuts[i_h]->Fill(vartoplot.at(i_tr));
                 onb_hists_cc1pi_top_aftercuts[i_h]->Fill(vartoplot.at(i_tr),1.0/vartoplot.size());
              } // end if(isSelected)
           } // end loop over tracks
        } // end loop over 1D Varstoplot

     } // end loop over entries in tree

  } // end if onbeam data tree exists
  else{ // this is important for plotting (to make sure DrawStack doesn't segfault when the histograms are given to it)
     for (size_t i_h=0; i_h<nplots; i_h++){
        onb_hists_cc1pi_pdg_beforecuts[i_h] = nullptr;
        onb_hists_cc1pi_top_beforecuts[i_h] = nullptr;
        onb_hists_cc1pi_pdg_aftercuts[i_h] = nullptr;
        onb_hists_cc1pi_top_aftercuts[i_h] = nullptr;
    }
  }

  // ----------------- Data: EXT/off-beam

  int nsel_offbeam=0;
  // Make histograms to fill
  TH1F *offb_hists_cc1pi_pdg_beforecuts[nplots];
  TH1F *offb_hists_cc1pi_top_beforecuts[nplots];
  TH1F *offb_hists_cc1pi_pdg_aftercuts[nplots];
  TH1F *offb_hists_cc1pi_top_aftercuts[nplots];
  if (t_offbeam){
     for (size_t i_h=0; i_h<nplots; i_h++){
      std::string histtitle_i = varstoplot_dummy.at(i_h).histtitle;
      std::string histname_i = varstoplot_dummy.at(i_h).histname;
      std::vector<double> bins_i = varstoplot_dummy.at(i_h).bins;

        offb_hists_cc1pi_pdg_beforecuts[i_h] = new TH1F(TString::Format("hCC1pi_offbeam_PDG_beforecuts_%s",histname_i.c_str()).Data(),histtitle_i.c_str(),bins_i.at(0),bins_i.at(1),bins_i.at(2));
        offb_hists_cc1pi_top_beforecuts[i_h] = new TH1F(TString::Format("hCC1pi_offbeam_Top_beforecuts_%s",histname_i.c_str()).Data(),histtitle_i.c_str(),bins_i.at(0),bins_i.at(1),bins_i.at(2));

        offb_hists_cc1pi_pdg_aftercuts[i_h] = new TH1F(TString::Format("hCC1pi_offbeam_PDG_aftercuts_%s",histname_i.c_str()).Data(),histtitle_i.c_str(),bins_i.at(0),bins_i.at(1),bins_i.at(2));
        offb_hists_cc1pi_top_aftercuts[i_h] = new TH1F(TString::Format("hCC1pi_offbeam_Top_aftercuts_%s",histname_i.c_str()).Data(),histtitle_i.c_str(),bins_i.at(0),bins_i.at(1),bins_i.at(2));
    }

    // Loop through on-beam data tree and fill plots
    for (int i = 0; i < t_offbeam->GetEntries(); i++){
      if (i%1000==0) std::cout << "Offbeam: " << i << "/" << t_offbeam->GetEntries() << std::endl;
      t_offbeam->GetEntry(i);
      Calcvars(&offbeam_vars, &fReader_offbeam_contained, &fReader_offbeam_uncontained);
      std::vector<CC1piPlotVars> Varstoplot = GetVarstoplot(&offbeam_vars);

      bool isSelected = IsEventSelected(&offbeam_vars);

      if (isSelected) nsel_offbeam++;

      // Fill 1D plots
      // Loop over Varstoplot
      for (size_t i_h = 0; i_h < nplots; i_h++){
         std::vector<double> vartoplot = *(Varstoplot.at(i_h).Var);
          // Loop over tracks
          for (size_t i_tr = 0; i_tr < vartoplot.size(); i_tr++){
             // if (!(onbeam_vars.TPCObj_PFP_isDaughter->at(i_tr) && bool(onbeam_vars.TPCObj_PFP_track_passesMIPcut->at(i_tr)))) continue;

             if (!(FillPlotForTrack(&(Varstoplot.at(i_h)), &offbeam_vars, i_tr))) continue;

             offb_hists_cc1pi_pdg_beforecuts[i_h]->Fill(vartoplot.at(i_tr));
             offb_hists_cc1pi_top_beforecuts[i_h]->Fill(vartoplot.at(i_tr),1.0/vartoplot.size());

             if(isSelected) {
                offb_hists_cc1pi_pdg_aftercuts[i_h]->Fill(vartoplot.at(i_tr));
                offb_hists_cc1pi_top_aftercuts[i_h]->Fill(vartoplot.at(i_tr),1.0/vartoplot.size());
             } // end if(isSelected)
          } // end loop over tracks
      } // end loop over 1D Varstoplot

    } // end loop over entries in tree

   } // end if onbeam data tree exists
   else{ // this is important for plotting (to make sure DrawStack doesn't segfault when the histograms are given to it)
      for (size_t i_h=0; i_h<nplots; i_h++){
         offb_hists_cc1pi_pdg_beforecuts[i_h] = nullptr;
         offb_hists_cc1pi_top_beforecuts[i_h] = nullptr;
         offb_hists_cc1pi_pdg_aftercuts[i_h] = nullptr;
         offb_hists_cc1pi_top_aftercuts[i_h] = nullptr;
     }
   }



   // -------------------- Now make all the plots

   // Make 1D plots
   for (size_t i_h=0; i_h < nplots; i_h++){
      TCanvas *c1 = new TCanvas("c1","c1");
      std::string printname = std::string(varstoplot_dummy.at(i_h).histname+".png");

      mc_hists_cc1pi_pdg_beforecuts[i_h]->DrawStack(POTscaling,c1,"",onb_hists_cc1pi_pdg_beforecuts[i_h],offb_hists_cc1pi_pdg_beforecuts[i_h],offbeamscaling,onminusoffbeam);
      c1->Print(std::string(std::string("CC1pi_pdg_beforecuts")+printname).c_str());
      c1->Clear();

      mc_hists_cc1pi_top_beforecuts[i_h]->DrawStack(POTscaling,c1,true,onb_hists_cc1pi_top_beforecuts[i_h],offb_hists_cc1pi_top_beforecuts[i_h],offbeamscaling,onminusoffbeam);
      c1->Print(std::string(std::string("CC1pi_top_beforecuts")+printname).c_str());
      c1->Clear();

      mc_hists_cc1pi_pdg_aftercuts[i_h]->DrawStack(POTscaling,c1,"",onb_hists_cc1pi_pdg_aftercuts[i_h],offb_hists_cc1pi_pdg_aftercuts[i_h],offbeamscaling,onminusoffbeam);
      c1->Print(std::string(std::string("CC1pi_pdg_aftercuts")+printname).c_str());
      c1->Clear();

      mc_hists_cc1pi_top_aftercuts[i_h]->DrawStack(POTscaling,c1,true,onb_hists_cc1pi_top_aftercuts[i_h],offb_hists_cc1pi_top_aftercuts[i_h],offbeamscaling,onminusoffbeam);
      c1->Print(std::string(std::string("CC1pi_top_aftercuts")+printname).c_str());

      delete c1;
   }

   // Make 2D plots
   for (size_t i_h=0; i_h < nplots2D; i_h++){
      TCanvas *c1 = new TCanvas("c1","c1");
      std::string printname = std::string(varstoplot2D_dummy.at(i_h).first.histname+"_vs_"+varstoplot2D_dummy.at(i_h).second.histname+".png");

      mc_hists2D_cc1pi_pdg_beforecuts[i_h]->Draw2D(c1);
      c1->Print(std::string(std::string("CC1pi2D_pdg_beforecuts")+printname).c_str());
      c1->Clear();

      // mc_hists2D_cc1pi_top_beforecuts[i_h]->Draw2D(c1,true);
      // c1->Print(std::string(std::string("CC1pi2D_top_beforecuts")+printname).c_str());
      // c1->Clear();

      mc_hists2D_cc1pi_pdg_aftercuts[i_h]->Draw2D(c1);
      c1->Print(std::string(std::string("CC1pi2D_pdg_aftercuts")+printname).c_str());
      c1->Clear();

      // mc_hists2D_cc1pi_top_aftercuts[i_h]->Draw2D(c1,true);
      // c1->Print(std::string(std::string("CC1pi2D_top_aftercuts")+printname).c_str());

      delete c1;
   }

   // Make custom plot: leading vs second MIP
   TCanvas *c1 = new TCanvas();
   selMIPs2D->Scale(1.0/selMIPs2D->Integral());
   selMIPs2D->Draw("colz text");
   c1->Print("CC1pi_selMIPs2D.png");

   SelectedEvents->DrawStack(1.,c1);
   c1->Print("CC1pi_SelectedEvents.png");

   AllEvents->DrawStack(1.,c1);
   c1->Print("CC1pi_AllEvents.png");

   NDaughters->DrawStack(1.,c1);
   c1->Print("CC1pi_NDaughters.png");

   delete c1;

   // Print out integrals
   SelectedEvents->PrintHistIntegrals();
   double CC1pi_selected = SelectedEvents->GetCC1piIntegral();
   double CC1pi_all = AllEvents->GetCC1piIntegral();
   std::cout << std::endl << "Total number of selected events: " << SelectedEvents->GetTotalIntegral() << " MC, " << nsel_onbeam << " beam-on data, " << nsel_offbeam << " ext data" << std::endl;
   std::cout << "CC1pi+ selection efficiency: " << CC1pi_selected << "/" << CC1pi_all << " = " << CC1pi_selected/CC1pi_all << std::endl;
   std::cout << "Background rejection: 1 - " << SelectedEvents->GetTotalIntegral() - CC1pi_selected << "/" << AllEvents->GetTotalIntegral() << " = " << 1-(SelectedEvents->GetTotalIntegral() - CC1pi_selected)/(1.0*AllEvents->GetTotalIntegral()) << std::endl;
   getSimPot(mcfile);

   f_MVA.Write();
}
