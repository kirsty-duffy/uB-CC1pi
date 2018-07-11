#include "MakeCutEfficiencyPlots_header.h"

// What variables do we want these plots as a function of?
std::vector<std::vector<double>> GetVarstoplot(treevars *vars){
   std::vector<std::vector<double>> varstoplot = {
      *(vars->TPCObj_PFP_Lmipoverp),
      *(vars->TPCObj_PFP_Lmumipovermumipp),
      *(vars->TPCObj_PFP_BrokenTrackAngle),
      *(vars->TPCObj_PFP_track_residual_mean),
      *(vars->TPCObj_PFP_track_residual_mean),
      *(vars->TPCObj_PFP_track_residual_std),
      *(vars->TPCObj_PFP_track_residual_std),
      *(vars->TPCObj_PFP_track_perc_used_hits),
      *(vars->TPCObj_PFP_VtxTrackDist),
      *(vars->TPCObj_PFP_isContained_double),
      *(vars->TPCObj_PFP_track_dEdx_truncmean_start)
   };
   return varstoplot;
};


// Binning (nbins, binlow, binhigh) in the same order as the vector above
std::vector<std::vector<double>> bins = {
   {25,0.3,3},    // Lmipoverp
   {25,0.5,0.9},  // Lmumipovermumipp
   {25,2.8,3.15}, // BrokenTrackAngle
   {25,0,2.8},    // residual_mean_up
   {25,-2.8,0},   // residual_mean_down
   {25,0,4},      // residual_std_up
   {25,0,2},      // residual_std_down
   {25,0,1},      // perc_used_hits
   {25,0,20},     // VtxTrackDist
   {2,0,2},        // isContained
   {50,0,10}      // dEdx truncmean at start
};

// Histogram titles in the same order as the vector above
std::vector<std::string> histtitles = {
   ";(L_{MIP})/(L_{p});",
   ";(L_{#mu}+L_{MIP})/(L_{#mu}+L_{MIP}+L_{p});",
   ";Angle [rad];",
   ";<r_{i}>;",
   ";<r_{i}>;",
   ";#sigma_{r_{i}};",
   ";#sigma_{r_{i}};",
   ";Fraction of used hits in cluster;",
   ";Distance from reconstructed vertex;",
   ";isContained;",
   ";Truncated Mean dE/dx at start of track;"
};

// What to call saved plots in the same order as the vector above
std::vector<std::string> histnames = {
   "Lmipoverp",
   "Lmumipovermumipp",
   "BrokenTrackAngle",
   "residual_mean_up",
   "residual_mean_down",
   "residual_std_up",
   "residual_std_down",
   "perc_used_hits",
   "VtxTrackDist",
   "isContained",
   "dEdx_truncmean_atstart"
};

// ---------------------------------------------------- //
//  Now the function starts
// ---------------------------------------------------- //

void MakeCutPlots(std::string mcfile){

   gStyle->SetTitleX(0.5);
   gStyle->SetTitleAlign(23);
   gStyle->SetOptTitle(1);
   gStyle->SetTitleBorderSize(0.);

   TFile *f_bnbcos = new TFile(mcfile.c_str(), "read");
   TTree *t_bnbcos = (TTree*)f_bnbcos->Get("cc1piselec/outtree");
   treevars mc_vars;
   settreevars(t_bnbcos,&mc_vars);

   // Sanity check: the plot vectors should be the same size
   t_bnbcos->GetEntry(0);
   Calcvars(&mc_vars);
   std::vector<std::vector<double>> varstoplot_dummy = GetVarstoplot(&mc_vars);
   std::vector<std::vector<double>> cutvarstoplot_dummy = GetCutvarstoplot(&mc_vars);

   // if (varstoplot_dummy.size() != bins.size()) std::cout << "WARNING varstoplot_dummy.size() = " << varstoplot_dummy.size() << "and bins.size() = " << bins.size() << ". This is going to cause you problems!" << std::endl;
   std::cout << "varstoplot_dummy.size() = " << varstoplot_dummy.size() << std::endl;
   std::cout << "bins.size() = " << bins.size() << std::endl;
   std::cout << "histtitles.size() = " << histtitles.size() << std::endl;
   std::cout << "histnames.size() = " << histnames.size() << std::endl;
   std::cout << "cutvarstoplot_dummy.size() = " << cutvarstoplot_dummy.size() << std::endl;
   std::cout << "KeepBelowCut.size() = " << KeepBelowCut.size() << std::endl;
   std::cout << "OnlyDaughters.size() = " << OnlyDaughters.size() << std::endl;
   std::cout << "TracksNeeded.size() = " << TracksNeeded.size() << std::endl;

   // ----------------- MC

   // Make histograms to fill
   const size_t nplots = varstoplot_dummy.size();
   StackedHistPDGCode *mc_hists_cc1pi_pdg_beforecuts[nplots];
   StackedHistTopology *mc_hists_cc1pi_top_beforecuts[nplots];
   StackedHistPDGCode *mc_hists_cc1pi_pdg_aftercuts[nplots];
   StackedHistTopology *mc_hists_cc1pi_top_aftercuts[nplots];
   for (size_t i_h=0; i_h<nplots; i_h++){
      mc_hists_cc1pi_pdg_beforecuts[i_h] = new StackedHistPDGCode(std::string("hCC1pi_PDG_beforecuts_")+histnames.at(i_h),histtitles.at(i_h),bins.at(i_h).at(0),bins.at(i_h).at(1),bins.at(i_h).at(2));
      mc_hists_cc1pi_top_beforecuts[i_h] = new StackedHistTopology(std::string("hCC1pi_Top_beforecuts_")+histnames.at(i_h),histtitles.at(i_h),bins.at(i_h).at(0),bins.at(i_h).at(1),bins.at(i_h).at(2));

      mc_hists_cc1pi_pdg_aftercuts[i_h] = new StackedHistPDGCode(std::string("hCC1pi_PDG_aftercuts_")+histnames.at(i_h),histtitles.at(i_h),bins.at(i_h).at(0),bins.at(i_h).at(1),bins.at(i_h).at(2));
      mc_hists_cc1pi_top_aftercuts[i_h] = new StackedHistTopology(std::string("hCC1pi_Top_aftercuts_")+histnames.at(i_h),histtitles.at(i_h),bins.at(i_h).at(0),bins.at(i_h).at(1),bins.at(i_h).at(2));
  }

   std::cout << "--- Considering the following variables --- " << std::endl;
   for (size_t i_var=1; i_var < nplots+1; i_var++){
         std::cout << histtitles.at(i_var-1).substr(1,histtitles.at(i_var-1).size()-2) << std::endl;
   }
   std::cout << "------" << std::endl;

   // Dummy plot in order to get percentage of different topologies that are selected
   StackedHistTopology *SelectedEvents = new StackedHistTopology("SelectedEvents",";;",1,0,1);

   // Loop through MC tree and fill plots
   for (int i = 0; i < t_bnbcos->GetEntries(); i++){
      t_bnbcos->GetEntry(i);
      Calcvars(&mc_vars);
      std::vector<std::vector<double>> Varstoplot = GetVarstoplot(&mc_vars);
      std::vector<std::vector<double>> Cutvarstoplot = GetCutvarstoplot(&mc_vars);

      
      // Determine if event is selected
      bool isSelected = true;
      for (size_t i_cut = 0; i_cut < Cutvarstoplot.size(); i_cut++){
         if(!IsEventSelected(CutValues.at(i_cut), Cutvarstoplot.at(i_cut), KeepBelowCut.at(i_cut), mc_vars.Truth_topology, *(mc_vars.TPCObj_PFP_isDaughter), OnlyDaughters.at(i_cut), TracksNeeded.at(i_cut))) {
            isSelected = false;
            break;
         }
      }

      if(isSelected) SelectedEvents->Fill((NuIntTopology)mc_vars.Truth_topology,0);

      // Loop over Varstoplot
      for (size_t i_h = 0; i_h < nplots; i_h++){

         // Loop over tracks
         for (size_t i_tr = 0; i_tr < Varstoplot.at(0).size(); i_tr++){
            mc_hists_cc1pi_pdg_beforecuts[i_h]->Fill((PDGCode)mc_vars.TPCObj_PFP_truePDG->at(i_tr),Varstoplot.at(i_h).at(i_tr));
            mc_hists_cc1pi_top_beforecuts[i_h]->Fill((NuIntTopology)mc_vars.Truth_topology,Varstoplot.at(i_h).at(i_tr),1.0/Varstoplot.at(0).size());

            if(isSelected) {
               mc_hists_cc1pi_pdg_aftercuts[i_h]->Fill((PDGCode)mc_vars.TPCObj_PFP_truePDG->at(i_tr),Varstoplot.at(i_h).at(i_tr));
               mc_hists_cc1pi_top_aftercuts[i_h]->Fill((NuIntTopology)mc_vars.Truth_topology,Varstoplot.at(i_h).at(i_tr),1.0/Varstoplot.at(0).size());
              
            } // end if(isSelected)

         } // end loop over tracks

      } // end loop over Varstoplot

   } // end loop over entries in tree

   // -------------------- Now make all the plots

   for (size_t i_h=0; i_h < nplots; i_h++){
      TCanvas *c1 = new TCanvas("c1","c1");
      std::string printname = std::string(histnames[i_h]+".png");

      mc_hists_cc1pi_pdg_beforecuts[i_h]->DrawStack(1.,c1);
      c1->Print(std::string(std::string("CC1pi_pdg_beforecuts")+printname).c_str());
      c1->Clear();

      mc_hists_cc1pi_top_beforecuts[i_h]->DrawStack(1.,c1,true);
      c1->Print(std::string(std::string("CC1pi_top_beforecuts")+printname).c_str());
      c1->Clear();

      mc_hists_cc1pi_pdg_aftercuts[i_h]->DrawStack(1.,c1);
      c1->Print(std::string(std::string("CC1pi_pdg_aftercuts")+printname).c_str());
      c1->Clear();

      mc_hists_cc1pi_top_aftercuts[i_h]->DrawStack(1.,c1,true);
      c1->Print(std::string(std::string("CC1pi_top_aftercuts")+printname).c_str());

      delete c1;
   }

   // Print out integrals
   SelectedEvents->GetHistIntegrals();
}
