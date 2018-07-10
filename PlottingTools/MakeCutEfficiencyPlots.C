#include "MakeCutEfficiencyPlots_header.h"

// What variables do we want these plots as a function of?
std::vector<std::vector<double>> GetCutvarstoplot(treevars *vars){
   std::vector<std::vector<double>> varstoplot = {
      *(vars->TPCObj_PFP_Lmipoverp),
      *(vars->TPCObj_PFP_Lmumipovermumipp),
      *(vars->TPCObj_PFP_BrokenTrackAngle),
      *(vars->TPCObj_PFP_track_residual_mean),
      *(vars->TPCObj_PFP_track_residual_mean),
      *(vars->TPCObj_PFP_track_residual_std),
      *(vars->TPCObj_PFP_track_residual_std),
      *(vars->TPCObj_PFP_track_perc_used_hits),
      *(vars->TPCObj_PFP_VtxTrackDist)
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
   {25,0,2},      // residual_std_up
   {25,0,3},      // residual_std_down
   {25,0,1},      // perc_used_hits
   {25,0,20}      // VtxTrackDist
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
   ";Distance from reconstructed vertex;"
};

// What to call saved plots in the same order as the vector above
std::vector<std::string> histnames = {
   "effpur_Lmipoverp",
   "effpur_Lmumipovermumipp",
   "effpur_BrokenTrackAngle",
   "effpur_residual_mean_up",
   "effpur_residual_mean_down",
   "effpur_residual_std_up",
   "effpur_residual_std_down",
   "effpur_perc_used_hits",
   "effpur_VtxTrackDist"
};

// For efficiency/purity we need to know whether we want to be keep tracks that have values above or below the cut value
std::vector<bool> KeepBelowCut = {
   false, // Lmipoverp
   false, // Lmumipovermumipp
   true,  // BrokenTrackAngle
   true,  // residual_mean_up
   false, // residual_mean_down
   true,  // residual_std_up
   false, // residual_std_down
   false, // perc_used_hits
   true   // VtxTrackDist
};

// Do we want to consider just the direct daughters of the neutrino?
std::vector<bool> OnlyDaughters = {
   true,  // Lmipoverp
   true,  // Lmumipovermumipp
   false, // BrokenTrackAngle
   true,  // residual_mean_up
   true,  // residual_mean_down
   true,  // residual_std_up
   true,  // residual_std_down
   true,  // perc_used_hits
   false  // VtxTrackDist
};

// How many tracks do we want to pass the cut? (Options are atleasttwo, exactlytwo, all)
std::vector<std::string> TracksNeeded = {
   "exactlytwo",  // Lmipoverp
   "exactlytwo",  // Lmumipovermumipp
   "all",         // BrokenTrackAngle
   "atleasttwo",  // residual_mean_up
   "atleasttwo",  // residual_mean_down
   "atleasttwo",  // residual_std_up
   "atleasttwo",  // residual_std_down
   "atleasttwo",  // perc_used_hits
   "atleasttwo"   // VtxTrackDist

};

// Cut values for N-1 plot
std::vector<double> CutValues = {
   1.,  // Lmipoverp
   0.67,  // Lmumipovermumipp
   3.05, // BrokenTrackAngle
   0.7,  // residual_mean_up
   -0.7,  // residual_mean_down
   2.5,  // residual_std_up
   0.,  // residual_std_down
   0.7,   // perc_used_hits
   15.  // VtxTrackDist
};

// ---------------------------------------------------- //
//  Now the function starts
// ---------------------------------------------------- //

void MakeCutEfficiencyPlots(std::string mcfile){

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
   std::vector<std::vector<double>> varstoplot_dummy = GetCutvarstoplot(&mc_vars);
   // if (varstoplot_dummy.size() != bins.size()) std::cout << "WARNING varstoplot_dummy.size() = " << varstoplot_dummy.size() << "and bins.size() = " << bins.size() << ". This is going to cause you problems!" << std::endl;
   std::cout << "varstoplot_dummy.size() = " << varstoplot_dummy.size() << std::endl;
   std::cout << "bins.size() = " << bins.size() << std::endl;
   std::cout << "histtitles.size() = " << histtitles.size() << std::endl;
   std::cout << "histnames.size() = " << histnames.size() << std::endl;
   std::cout << "KeepBelowCut.size() = " << KeepBelowCut.size() << std::endl;
   std::cout << "OnlyDaughters.size() = " << OnlyDaughters.size() << std::endl;
   std::cout << "TracksNeeded.size() = " << TracksNeeded.size() << std::endl;
   std::cout << "CutValues.size() = " << CutValues.size() << std::endl;

   // ----------------- MC

   // Make histograms to fill
   const size_t nplots = varstoplot_dummy.size();
   histCC1piselEffPur *mc_hists_cc1pieffpur[nplots];
   for (size_t i_h=0; i_h<nplots; i_h++){
      mc_hists_cc1pieffpur[i_h] = new histCC1piselEffPur(std::string("hCC1pieffpur_")+histnames.at(i_h),histtitles.at(i_h),bins.at(i_h).at(0),bins.at(i_h).at(1),bins.at(i_h).at(2));
   }

   histCC1piselEffPur *Nminus1plots = new histCC1piselEffPur("hCC1pi_nminus1","N-1 plot;;Efficiency, purity",nplots+2,0,nplots+2);

   std::cout << "Considering the following cuts --- " << std::endl;
   for (size_t i_bin=1; i_bin < nplots+3; i_bin++){
      if (i_bin==1) Nminus1plots->SetBinLabel(i_bin, "CC incl. 2 tracks");
      else if (i_bin==2) Nminus1plots->SetBinLabel(i_bin,"All cuts");
      else{
         std::string tmpname = histtitles.at(i_bin-3).substr(1,histtitles.at(i_bin-3).size()-2);
         std::cout << tmpname << std::endl; Nminus1plots->SetBinLabel(i_bin,tmpname);
      }
   }
   std::cout << "---" << std::endl;

   histCC1piselEffPur2D *hCC1pieffpur_residual_mean = new histCC1piselEffPur2D("hCC1pieffpur_residual_mean",";residual_mean_up;residual_mean_down",bins.at(3).at(0),bins.at(3).at(1),bins.at(3).at(2),bins.at(4).at(0),bins.at(4).at(1),bins.at(4).at(2));

   // Loop through MC tree and fill plots
   for (int i = 0; i < t_bnbcos->GetEntries(); i++){
      t_bnbcos->GetEntry(i);
      Calcvars(&mc_vars);
      std::vector<std::vector<double>> Cutvarstoplot = GetCutvarstoplot(&mc_vars);

      for (size_t i_h = 0; i_h < nplots; i_h++){
         FillCC1piEffPurHist(mc_hists_cc1pieffpur[i_h],Cutvarstoplot.at(i_h),mc_vars.Truth_topology,KeepBelowCut.at(i_h),mc_vars.Marco_selected,*mc_vars.TPCObj_PFP_isDaughter, OnlyDaughters.at(i_h), TracksNeeded.at(i_h));
      }

      // Now fill histograms for N-1 plots
      FillNminus1EffPurHist(Nminus1plots,Cutvarstoplot,mc_vars.Truth_topology,KeepBelowCut,mc_vars.Marco_selected,*(mc_vars.TPCObj_PFP_isDaughter),OnlyDaughters,TracksNeeded,CutValues);

      // Fill 2D histograms
      // Make sure the indices actually point to the correct points in the other vectors. For instance, here we have "3" and "4" corresponding to residual_mean_up and residual_mean_down.
      // Also do the same in the histCC1piselEffPur2D delcaration above.
      std::vector<std::vector<double>> value_vec_residual_mean = {Cutvarstoplot.at(3),Cutvarstoplot.at(4)};
      std::vector<bool> KeepBelowCut_residual_mean = {KeepBelowCut.at(3),KeepBelowCut.at(4)};
      std::vector<bool> OnlyDaughters_residual_mean = {OnlyDaughters.at(3),OnlyDaughters.at(4)};
      std::vector<std::string> TracksNeeded_residual_mean = {TracksNeeded.at(3),TracksNeeded.at(4)};
      FillCC1piEffPurHist2D(hCC1pieffpur_residual_mean,value_vec_residual_mean,mc_vars.Truth_topology,KeepBelowCut_residual_mean,mc_vars.Marco_selected,*(mc_vars.TPCObj_PFP_isDaughter),OnlyDaughters_residual_mean,TracksNeeded_residual_mean);

   } // end loop over entries in tree

   // -------------------- Now make all the plots

   for (size_t i_h=0; i_h < nplots; i_h++){
      TCanvas *c1 = new TCanvas("c1","c1");

      DrawCC1piMCEffPur(c1, mc_hists_cc1pieffpur[i_h]);
      std::string printname = std::string(histnames[i_h]+".png");
      c1->Print(std::string(std::string("CC1pi_")+printname).c_str());
/*
      mc_hists_cc1pieffpur[i_h]->h_cc1pi_sel->Draw();
      c1->Print(std::string(std::string("CC1pi_sel_")+printname).c_str());
      mc_hists_cc1pieffpur[i_h]->h_cc1pi_notsel->Draw();
      c1->Print(std::string(std::string("CC1pi_notsel_")+printname).c_str());
      mc_hists_cc1pieffpur[i_h]->h_bg_sel->Draw();
      c1->Print(std::string(std::string("Background_sel_")+printname).c_str());
*/
      delete c1;
   }
   TCanvas *c1 = new TCanvas("c1","c1");
   DrawCC1piMCEffPur(c1, Nminus1plots,"hist",true);
   c1->Print("Nminus1plots.png");

   DrawCC1piMCEffPur2D(c1, hCC1pieffpur_residual_mean);
   c1->Print("hCC1pieffpur_residual_mean.png");
   delete c1;
}
