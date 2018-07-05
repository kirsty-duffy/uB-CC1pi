#include "MakeCutEfficiencyPlots_header.h"

// What variables do we want these plots as a function of?
std::vector<std::vector<double>> GetCutvarstoplot(treevars *vars){
   std::vector<std::vector<double>> varstoplot = {
      *(vars->TPCObj_PFP_Lmipoverp),
      *(vars->TPCObj_PFP_Lmumipovermumipp),
      *(vars->TPCObj_PFP_BrokenTrackAngle)
   };
   return varstoplot;
};


// Binning (nbins, binlow, binhigh) in the same order as the vector above
std::vector<std::vector<double>> bins = {
   {20,0.3,3},    // Lmipoverp
   {20,0.5,0.9},  // Lmumipovermumipp
   {20,1,3.15} // BrokenTrackAngle
};

// Histogram titles in the same order as the vector above
std::vector<std::string> histtitles = {
   ";(L_{MIP})/(L_{p});",
   ";(L_{#mu}+L_{MIP})/(L_{#mu}+L_{MIP}+L_{p});",
   ";Angle [rad];"
};

// What to call saved plots in the same order as the vector above
std::vector<std::string> histnames = {
   "effpur_Lmipoverp",
   "effpur_Lmumipovermumipp",
   "effpur_BrokenTrackAngle"
};

// For efficiency/purity we need to know whether we want to be higher or lower than the variable
std::vector<bool> Cutlow = {
   false, // Lmipoverp
   false, // Lmumipovermumipp
   true   // BrokenTrackAngle
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
   std::cout << "Cutlow.size() = " << Cutlow.size() << std::endl;

   // ----------------- MC

   // Make histograms to fill
   const size_t nplots = varstoplot_dummy.size();
   histCC1piselEffPur *mc_hists_cc1pieffpur[nplots];
   for (size_t i_h=0; i_h<nplots; i_h++){
      mc_hists_cc1pieffpur[i_h] = new histCC1piselEffPur(std::string("hCC1pieffpur_")+histnames.at(i_h),histtitles.at(i_h),bins.at(i_h).at(0),bins.at(i_h).at(1),bins.at(i_h).at(2));
   }

   // Loop through MC tree and fill plots
   for (int i = 0; i < t_bnbcos->GetEntries(); i++){
      t_bnbcos->GetEntry(i);
      Calcvars(&mc_vars);
      std::vector<std::vector<double>> Cutvarstoplot = GetCutvarstoplot(&mc_vars);

      for (size_t i_h = 0; i_h < nplots; i_h++){
         FillCC1piEffPurHist(mc_hists_cc1pieffpur[i_h],Cutvarstoplot.at(i_h),mc_vars.Truth_topology,Cutlow.at(i_h),mc_vars.Marco_selected,*mc_vars.TPCObj_PFP_isDaughter);
      }


   } // end loop over entries in tree


   // -------------------- Now make all the plots

   for (size_t i_h=0; i_h < nplots; i_h++){
      TCanvas *c1 = new TCanvas("c1","c1");

      DrawCC1piMCEffPur(c1, mc_hists_cc1pieffpur[i_h]);
      std::string printname = std::string(histnames[i_h]+".png");
      c1->Print(std::string(std::string("CC1pi_")+printname).c_str());

      mc_hists_cc1pieffpur[i_h]->h_cc1pi_sel->Draw();
      c1->Print(std::string(std::string("CC1pi_sel_")+printname).c_str());
      mc_hists_cc1pieffpur[i_h]->h_cc1pi_notsel->Draw();
      c1->Print(std::string(std::string("CC1pi_notsel_")+printname).c_str());
      mc_hists_cc1pieffpur[i_h]->h_bg_sel->Draw();
      c1->Print(std::string(std::string("Background_sel_")+printname).c_str());

      delete c1;
   }

}
