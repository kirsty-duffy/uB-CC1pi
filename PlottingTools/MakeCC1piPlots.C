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
#include "CC1pi_EffPurHists.h"

#include <boost/algorithm/string.hpp>

// What variables do we want these plots as a function of?
std::vector<CC1piPlotVars> GetVarstoplot(treevars *vars){
   std::vector<CC1piPlotVars> varstoplot = {
      Var_TPCObj_PFP_track_BDTscore(vars)
      ,Var_TPCObj_PFP_track_BDTscore_IsMIP(vars)
      ,Var_TPCObj_PFP_track_BDTscore_NotMIP(vars)
      ,Var_TPCObj_PFP_track_BDTscore_MuMuPairs(vars)
      ,Var_TPCObj_AngleBetweenMIPs(vars)
      ,Var_TPCObj_AngleBetweenMIPs_mumupairs(vars)
      ,Var_TPCObj_PFP_track_dEdx_truncmean_start_IsMIP(vars)
      ,Var_TPCObj_PFP_track_dEdx_truncmean_start_NotMIP(vars)
      ,Var_TPCObj_PFP_track_dEdx_truncmean_start_mumupairs(vars)
      ,Var_TPCObj_PFP_VtxTrackDist_IsMIP(vars)
      ,Var_TPCObj_PFP_VtxTrackDist_NotMIP(vars)
      ,Var_TPCObj_PFP_VtxTrackDist_mumupairs(vars)
      ,Var_TPCObj_PFP_track_dedx_grminhits_IsMIP(vars)
      ,Var_TPCObj_PFP_track_dedx_grminhits_NotMIP(vars)
      ,Var_TPCObj_PFP_track_dedx_grminhits_mumupairs(vars)
      ,Var_TPCObj_PFP_lnLmipoverp_IsMIP(vars)
      ,Var_TPCObj_PFP_lnLmipoverp_NotMIP(vars)
      ,Var_TPCObj_PFP_lnLmipoverp_mumupairs(vars)
      ,Var_TPCObj_BDTscore_MIPdiv(vars)
      ,Var_TPCObj_BDTscore_MIPdiff(vars)
      ,Var_TPCObj_PFP_trueKE_selMIPs(vars)
      ,Var_TPCObj_PFP_trueEndP_selMIPs(vars)
      ,Var_TPCObj_PFP_track_theta_selMIPs(vars)
      ,Var_TPCObj_PFP_track_phi_selMIPs(vars)
      ,Var_TPCObj_PFP_track_nhits_LeadingMIP(vars)
      ,Var_TPCObj_PFP_track_nhits_SecondMIP(vars)
      ,Var_TPCObj_PFP_track_length_LeadingMIP_mumupairs(vars)
      ,Var_TPCObj_PFP_track_length_SecondMIP_mumupairs(vars)
      ,Var_TPCObj_PFP_track_length_LeadingMIP(vars)
      ,Var_TPCObj_PFP_track_length_SecondMIP(vars)
      ,Var_TPCObj_PFP_track_length(vars)
      ,Var_TPCObj_PFP_lnLmipovermu(vars)
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
      ,{Var_TPCObj_PFP_track_theta(vars),Var_TPCObj_PFP_track_BDTscore(vars)}
      // ,{Var_TPCObj_PFP_track_phi(vars),Var_TPCObj_PFP_track_BDTscore(vars)}
      ,{Var_TPCObj_AngleBetweenMIPs(vars),Var_TPCObj_PFP_track_BDTscore(vars)}
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
//   muon -> Branch("nhits", &(MVA_vars->nhits));
   muon -> Branch("lnLmipoverp", &(MVA_vars->lnLmipoverp));
   pion -> Branch("dEdx_truncmean_start", &(MVA_vars->dEdx_truncmean_start));
   pion -> Branch("VtxTrackDist", &(MVA_vars->VtxTrackDist));
//   pion -> Branch("nhits", &(MVA_vars->nhits));
   pion -> Branch("lnLmipoverp", &(MVA_vars->lnLmipoverp));
   background -> Branch("dEdx_truncmean_start", &(MVA_vars->dEdx_truncmean_start));
   background -> Branch("VtxTrackDist", &(MVA_vars->VtxTrackDist));
//   background -> Branch("nhits", &(MVA_vars->nhits));
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

   std::string BookMVAType = "BDTG";
   std::string BookMVALoc = "/uboone/app/users/ddevitt/LArSoft_v06_26_01_14_uboonecode_v06_26_01_22/srcs/uboonecode/uboone/CC1pi/MVA/dataset_noNhits/weights/TMVAClassification_BDTG.weights.xml";

   TMVA::Reader fReader("");
   fReader.AddVariable("dEdx_truncmean_start", &(mc_vars.float_dEdx_truncmean_start));
//   fReader.AddVariable("VtxTrackDist", &(mc_vars.float_VtxTrackDist));
//   fReader.AddVariable("nhits", &(mc_vars.float_nhits));
   fReader.AddVariable("lnLmipoverp", &(mc_vars.float_lnLmipoverp));
   fReader.BookMVA(BookMVAType.c_str(), BookMVALoc.c_str());

   TFile *f_onbeam=nullptr;
   TTree *t_onbeam=nullptr;
   treevars onbeam_vars;
   TMVA::Reader fReader_onbeam("");
   if (onbeamdatafile!=""){
     std::cout << "Making data-MC comparisons" << std::endl;
     f_onbeam = new TFile(onbeamdatafile.c_str(), "read");
     t_onbeam = (TTree*)f_onbeam->Get("cc1piselec/outtree");
     settreevars(t_onbeam,&onbeam_vars);

     fReader_onbeam.AddVariable("dEdx_truncmean_start", &(onbeam_vars.float_dEdx_truncmean_start));
     fReader_onbeam.AddVariable("VtxTrackDist", &(onbeam_vars.float_VtxTrackDist));
     fReader_onbeam.AddVariable("nhits", &(onbeam_vars.float_nhits));
     fReader_onbeam.AddVariable("lnLmipoverp", &(onbeam_vars.float_lnLmipoverp));
     fReader_onbeam.BookMVA(BookMVAType.c_str(), BookMVALoc.c_str());
   }

   TFile *f_offbeam=nullptr;
   TTree *t_offbeam=nullptr;
   treevars offbeam_vars;
   TMVA::Reader fReader_offbeam("");
   if (offbeamdatafile!=""){
     std::cout << "Making data-MC comparisons" << std::endl;
     f_offbeam = new TFile(offbeamdatafile.c_str(), "read");
     t_offbeam = (TTree*)f_offbeam->Get("cc1piselec/outtree");
     settreevars(t_offbeam,&offbeam_vars);

     fReader_offbeam.AddVariable("dEdx_truncmean_start", &(offbeam_vars.float_dEdx_truncmean_start));
     fReader_offbeam.AddVariable("VtxTrackDist", &(offbeam_vars.float_VtxTrackDist));
     fReader_offbeam.AddVariable("nhits", &(offbeam_vars.float_nhits));
     fReader_offbeam.AddVariable("lnLmipoverp", &(offbeam_vars.float_lnLmipoverp));
     fReader_offbeam.BookMVA(BookMVAType.c_str(), BookMVALoc.c_str());
   }

   ofstream evdinfo;
   evdinfo.open("evdinfo.txt");

   // Sanity check: the plot vectors should be the same size
   t_bnbcos->GetEntry(0);
   Calcvars(&mc_vars, &fReader);
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

   // Custom 2D histograms: muon candidate vs pion candidate true PDG
   TH2D *selMIPsPID2D = new TH2D("selMIPsPID2D","Selected MIP-like tracks;Muon Candidate; Pion Candidate",nMIPpdgs,0,nMIPpdgs,nMIPpdgs,0,nMIPpdgs);
   for (size_t i_bin=1; i_bin<nMIPpdgs+1; i_bin++){
      if (i_bin==nMIPpdgs){
         selMIPsPID2D->GetXaxis()->SetBinLabel(i_bin,"Other");
         selMIPsPID2D->GetYaxis()->SetBinLabel(i_bin,"Other");
         continue;
      }
      selMIPsPID2D->GetXaxis()->SetBinLabel(i_bin,PDGenum2str(MIPpdgs.at(i_bin-1)).c_str());
      selMIPsPID2D->GetYaxis()->SetBinLabel(i_bin,PDGenum2str(MIPpdgs.at(i_bin-1)).c_str());
   }

   // Custom 2D histograms: leading vs second MIP true PDG for signal events only
   TH2D *selMIPs2D_signal = new TH2D("selMIPs2D_signal","Selected MIP-like tracks from true signal events;Muon Candidate; Pion Candidate",nMIPpdgs+1,0,nMIPpdgs+1,nMIPpdgs+1,0,nMIPpdgs+1);
   for (size_t i_bin=1; i_bin<nMIPpdgs+1; i_bin++){
      if (i_bin==nMIPpdgs){
         selMIPs2D_signal->GetXaxis()->SetBinLabel(i_bin,"Other");
         selMIPs2D_signal->GetYaxis()->SetBinLabel(i_bin,"Other");
         selMIPs2D_signal->GetXaxis()->SetBinLabel(i_bin+1,"Row Sums");
         selMIPs2D_signal->GetYaxis()->SetBinLabel(i_bin+1,"Col Sums");
         continue;
      }
      selMIPs2D_signal->GetXaxis()->SetBinLabel(i_bin,PDGenum2str(MIPpdgs.at(i_bin-1)).c_str());
      selMIPs2D_signal->GetYaxis()->SetBinLabel(i_bin,PDGenum2str(MIPpdgs.at(i_bin-1)).c_str());
   }

   // Custom 2D histogram: containment vs length for muons and pions from true signal events
   TH2D *ContainmentLength2D = new TH2D("ContainmentLength2D","True CC1#pi^{+} events;Contained;Longer",4,0,4,2,0,2);
   ContainmentLength2D->GetXaxis()->SetBinLabel(1,"neither");
   ContainmentLength2D->GetXaxis()->SetBinLabel(2,"muon only");
   ContainmentLength2D->GetXaxis()->SetBinLabel(3,"pion only");
   ContainmentLength2D->GetXaxis()->SetBinLabel(4,"both");
   ContainmentLength2D->GetYaxis()->SetBinLabel(1,"#pi^{+}");
   ContainmentLength2D->GetYaxis()->SetBinLabel(2,"#mu^{-}");

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

   // Standard efficiency plots as a function of kinematic variables
   // True mu P, true mu theta, true mu phi
   // True pi P, true pi theta, true pi phi
   // True mu-pi opening angle
   histCC1piselEffPur *mc_effpur_truemuP = new histCC1piselEffPur("mc_effpur_truemuP","CC1pi selection efficiency;True muon momentum (GeV/c);CC1pi Selection Efficiency",50,0,1);
   histCC1piselEffPur *mc_effpur_truemuTheta = new histCC1piselEffPur("mc_effpur_truemuTheta","CC1pi selection efficiency;True muon theta;CC1pi Selection Efficiency",32,0,3.2);
   histCC1piselEffPur *mc_effpur_truemuPhi = new histCC1piselEffPur("mc_effpur_truemuPhi","CC1pi selection efficiency;True muon phi;CC1pi Selection Efficiency",32,-3.2,3.2);
   histCC1piselEffPur *mc_effpur_truepiP = new histCC1piselEffPur("mc_effpur_truepiP","CC1pi selection efficiency;True pion momentum (GeV/c);CC1pi Selection Efficiency",50,0,1);
   histCC1piselEffPur *mc_effpur_truepiTheta = new histCC1piselEffPur("mc_effpur_truepiTheta","CC1pi selection efficiency;True pion theta;CC1pi Selection Efficiency",32,0,3.2);
   histCC1piselEffPur *mc_effpur_truepiPhi = new histCC1piselEffPur("mc_effpur_truepiPhi","CC1pi selection efficiency;True pion phi;CC1pi Selection Efficiency",32,-3.2,3.2);
   histCC1piselEffPur *mc_effpur_truemupiOpeningAngle = new histCC1piselEffPur("mc_effpur_truemupiOpeningAngle","CC1pi selection efficiency;True muon-pion opening angle;CC1pi Selection Efficiency",16,0,3.2);
   histCC1piselEffPur *mc_effpur_trueQ2 = new histCC1piselEffPur("mc_effpur_trueQ2","CC1pi selection efficiency;True Q^{2} (GeV^{2});CC1pi Selection Efficiency",50,0,2);
   histCC1piselEffPur *mc_effpur_trueW = new histCC1piselEffPur("mc_effpur_trueW","CC1pi selection efficiency;True W (GeV);CC1pi Selection Efficiency",100,0,4);


   int protonplotsmade=0;

   // Loop through MC tree and fill plots
   for (int i = 0; i < t_bnbcos->GetEntries(); i++){
     if (i%1000==0) std::cout << "MC: " << i << "/" << t_bnbcos->GetEntries() << std::endl;

      t_bnbcos->GetEntry(i);

      Calcvars(&mc_vars, &fReader);
      MakeMVATrees(muon, pion, background, MVA_vars, &mc_vars);
      std::vector<CC1piPlotVars> Varstoplot = GetVarstoplot(&mc_vars);
      std::vector<std::pair<CC1piPlotVars,CC1piPlotVars>> Varstoplot2D = GetVarstoplot2D(&mc_vars);

      AllEvents->Fill((NuIntTopology)mc_vars.Truth_topology,0.5);
      NDaughters->Fill((NuIntTopology)mc_vars.Truth_topology,mc_vars.TPCObj_NDaughterPFPs->at(0));

      bool isSelected = IsEventSelected(&mc_vars);

      // Fill efficiency plots
      if ((NuIntTopology)mc_vars.Truth_topology == kCC1piplus0p || (NuIntTopology)mc_vars.Truth_topology == kCC1piplus1p || (NuIntTopology)mc_vars.Truth_topology == kCC1piplusNp){
         double true_mu_P=-9999, true_mu_theta=-9999, true_mu_phi=-9999, true_mu_E=-9999, true_pi_P=-9999, true_pi_theta=-9999, true_pi_phi=-9999, true_opening_angle=-9999;
         TVector3 true_mu_startdir(1,1,1);
         TVector3 true_pi_startdir(1,1,1);

         for (size_t i_mcp=0; i_mcp<mc_vars.MCP_ID->size(); i_mcp++){
            // Skip not-neutrino-daughters
            if (mc_vars.MCP_MotherID->at(i_mcp)!=0) continue;

            // True muon
            if ((PDGCode)mc_vars.MCP_PDG->at(i_mcp)==kMuMinus){
               true_mu_E = mc_vars.MCP_E->at(i_mcp);
               true_mu_P = mc_vars.MCP_P->at(i_mcp);
               true_mu_startdir.SetXYZ(mc_vars.MCP_Px->at(i_mcp),mc_vars.MCP_Py->at(i_mcp),mc_vars.MCP_Pz->at(i_mcp));
               true_mu_theta = true_mu_startdir.Theta();
               true_mu_phi = true_mu_startdir.Phi();
            }

            // True pion
            if ((PDGCode)mc_vars.MCP_PDG->at(i_mcp)==kPiPlus){
               true_pi_P = mc_vars.MCP_P->at(i_mcp);
               true_pi_startdir.SetXYZ(mc_vars.MCP_Px->at(i_mcp),mc_vars.MCP_Py->at(i_mcp),mc_vars.MCP_Pz->at(i_mcp));
               true_pi_theta = true_pi_startdir.Theta();
               true_pi_phi = true_pi_startdir.Phi();
            }
         }

         true_mu_startdir.SetMag(1);
         true_pi_startdir.SetMag(1);
         true_opening_angle = TMath::ACos(true_mu_startdir.Dot(true_pi_startdir));

         double true_nu_E=mc_vars.nu_E;
         // Calculate true W and Q2 as in MINERvA: https://arxiv.org/pdf/1406.6415.pdf
         // Q^2 = 2Enu(Emu − |pmu| cos(thetamu)) − m_mu^2
         // W^2 = M_p^2 − Q^2 + 2M_p*E_recoil
         // where E_recoil = Enu − Emu
         double Q2 = 2*true_nu_E*(true_mu_E - TMath::Abs(true_mu_P)*TMath::Cos(true_mu_theta)) - 0.105658*0.105658;
         double E_recoil = true_nu_E - true_mu_E;
         double W2 = 0.938272*0.938272 - Q2 + 2*0.938272*E_recoil;

         // std::cout << "true_nu_E = " << true_nu_E << ", true_mu_E = " << true_mu_E << ", true_mu_P = " << true_mu_P << ", Q2 = " << Q2 << ", Erecoil = " << E_recoil << ", W = " << TMath::Sqrt(W2) << std::endl;

         if (isSelected){
            mc_effpur_truemuP->h_cc1pi_sel->Fill(true_mu_P);
            mc_effpur_truemuTheta->h_cc1pi_sel->Fill(true_mu_theta);
            mc_effpur_truemuPhi->h_cc1pi_sel->Fill(true_mu_phi);
            mc_effpur_truepiP->h_cc1pi_sel->Fill(true_pi_P);
            mc_effpur_truepiTheta->h_cc1pi_sel->Fill(true_pi_theta);
            mc_effpur_truepiPhi->h_cc1pi_sel->Fill(true_pi_phi);
            mc_effpur_truemupiOpeningAngle->h_cc1pi_sel->Fill(true_opening_angle);
            mc_effpur_trueQ2->h_cc1pi_sel->Fill(Q2);
            mc_effpur_trueW->h_cc1pi_sel->Fill(TMath::Sqrt(W2));
         }
         else{
            mc_effpur_truemuP->h_cc1pi_notsel->Fill(true_mu_P);
            mc_effpur_truemuTheta->h_cc1pi_notsel->Fill(true_mu_theta);
            mc_effpur_truemuPhi->h_cc1pi_notsel->Fill(true_mu_phi);
            mc_effpur_truepiP->h_cc1pi_notsel->Fill(true_pi_P);
            mc_effpur_truepiTheta->h_cc1pi_notsel->Fill(true_pi_theta);
            mc_effpur_truepiPhi->h_cc1pi_notsel->Fill(true_pi_phi);
            mc_effpur_truemupiOpeningAngle->h_cc1pi_notsel->Fill(true_opening_angle);
            mc_effpur_trueQ2->h_cc1pi_notsel->Fill(Q2);
            mc_effpur_trueW->h_cc1pi_notsel->Fill(TMath::Sqrt(W2));
         }
      }


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
         // std::cout << i << " " << mc_vars.TPCObj_LeadingMIPtrackIndex << " " << mc_vars.TPCObj_SecondMIPtrackIndex << std::endl;

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
         if ((NuIntTopology)mc_vars.Truth_topology == kCC1piplus0p || (NuIntTopology)mc_vars.Truth_topology == kCC1piplus1p || (NuIntTopology)mc_vars.Truth_topology == kCC1piplusNp)
            selMIPs2D_signal->Fill((double)LeadingMIPbin-0.5,(double)SecondMIPbin-0.5);


         //PDGCode MuonCandpdg = PDGCode(mc_vars.muoncandidatePDG);
         //PDGCode PionCandpdg = PDGCode(mc_vars.pioncandidatePDG);
         PDGCode MuonCandpdg = PDGCode(mc_vars.BDT_muoncandidatePDG);
         PDGCode PionCandpdg = PDGCode(mc_vars.BDT_pioncandidatePDG);

         int MuonCandbin = nMIPpdgs;
         int PionCandbin = nMIPpdgs;
         for (int i_bin=0; i_bin<MIPpdgs.size(); i_bin++){
            if (MIPpdgs.at(i_bin) == MuonCandpdg) MuonCandbin = i_bin+1;
            if (MIPpdgs.at(i_bin) == PionCandpdg) PionCandbin = i_bin+1;
         }
         selMIPsPID2D->Fill((double)MuonCandbin-0.5,(double)PionCandbin-0.5);



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

            // Hack: fill plots for selected MIPs only (for diagnosing issues with MIP selection)
            // if (!(mc_vars.TPCObj_PFP_isDaughter->at(i_tr))) continue;
            // if (EvalMIPCut(&mc_vars,i_tr,nullptr)<1) continue;

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

      // Fill custom 2D plot
      ContainmentLength2D->Fill(mc_vars.TPCObj_mupiContained->at(0),mc_vars.TPCObj_muonLonger->at(0));

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
        Calcvars(&onbeam_vars, &fReader_onbeam);
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
      Calcvars(&offbeam_vars, &fReader_offbeam);
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

   selMIPs2D_signal->Scale(1.0/selMIPs2D_signal->Integral());
   for (size_t x_bin=1; x_bin<nMIPpdgs+1; x_bin++){
      float total = 0;
      for (size_t y_bin=1; y_bin<nMIPpdgs+1; y_bin++){
         total+=selMIPs2D_signal->GetBinContent(x_bin,y_bin);
      }
      selMIPs2D_signal->SetBinContent(x_bin,nMIPpdgs+1,total);
   }
   for (size_t y_bin=1; y_bin<nMIPpdgs+1; y_bin++){
      float total = 0;
      for (size_t x_bin=1; x_bin<nMIPpdgs+1; x_bin++){
         total+=selMIPs2D_signal->GetBinContent(x_bin,y_bin);
      }
      selMIPs2D_signal->SetBinContent(nMIPpdgs+1,y_bin,total);
   }
   selMIPs2D_signal->Draw("colz text");
   c1->Print("CC1pi_selMIPs2D_signal.png");

   selMIPsPID2D->Scale(1.0/selMIPsPID2D->Integral());
   selMIPsPID2D->Draw("colz text");
   c1->Print("CC1pi_selMIPsPID2D.png");

   SelectedEvents->DrawStack(1.,c1);
   c1->Print("CC1pi_SelectedEvents.png");

   AllEvents->DrawStack(1.,c1);
   c1->Print("CC1pi_AllEvents.png");

   NDaughters->DrawStack(1.,c1);
   c1->Print("CC1pi_NDaughters.png");

   ContainmentLength2D->Scale(1.0/ContainmentLength2D->Integral());
   ContainmentLength2D->Draw("colz text");
   c1->Print("ContainmentLength2D.png");

   DrawCC1piMCEffOnly(c1,mc_effpur_truemuP);
   c1->Print("Efficiency_truemuP.png");
   DrawCC1piMCEffOnly(c1,mc_effpur_truemuTheta);
   c1->Print("Efficiency_truemuTheta.png");
   DrawCC1piMCEffOnly(c1,mc_effpur_truemuPhi);
   c1->Print("Efficiency_truemuPhi.png");
   DrawCC1piMCEffOnly(c1,mc_effpur_truepiP);
   c1->Print("Efficiency_truepiP.png");
   DrawCC1piMCEffOnly(c1,mc_effpur_truepiTheta);
   c1->Print("Efficiency_truepiTheta.png");
   DrawCC1piMCEffOnly(c1,mc_effpur_truepiPhi);
   c1->Print("Efficiency_truepiPhi.png");
   DrawCC1piMCEffOnly(c1,mc_effpur_truemupiOpeningAngle);
   c1->Print("Efficiency_trueOpeningAngle.png");
   DrawCC1piMCEffOnly(c1,mc_effpur_trueQ2);
   c1->Print("Efficiency_trueQ2.png");
   DrawCC1piMCEffOnly(c1,mc_effpur_trueW);
   c1->Print("Efficiency_trueW.png");

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
