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
      *(vars->TPCObj_PFP_isContained_double)
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
   {2,0,2}        // isContained
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
   ";isContained;"
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
   "effpur_VtxTrackDist",
   "effpur_isContained"
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

   // EndProcess Truth Study
   StackedHistPDGCode *EndProcess_selected = new StackedHistPDGCode("EndProcess_selected",";;",32,0,32);
   StackedHistPDGCode *EndProcess_CCincl = new StackedHistPDGCode("EndProcess_CCincl",";;",32,0,32);
   std::vector<std::string> EndProcess_vect;
   TH2D *EndProcess_selected_2D = new TH2D("EndProcess_selected_2D",";Muons;Pions",4,0,4,4,0,4);
   TH2D *EndProcess_CCincl_2D = new TH2D("EndProcess_CCincl_2D",";Muons;Pions",4,0,4,4,0,4);
   // Set bins to be in a sensible order
   EndProcess_selected_2D->Fill("Exiting","Exiting",0);
   EndProcess_selected_2D->Fill("FastScintillation","FastScintillation",0);
   EndProcess_selected_2D->Fill("Decay","Decay",0);
   EndProcess_selected_2D->Fill("Other","Other",0);
   EndProcess_CCincl_2D->Fill("Exiting","Exiting",0);
   EndProcess_CCincl_2D->Fill("FastScintillation","FastScintillation",0);
   EndProcess_CCincl_2D->Fill("Decay","Decay",0);
   EndProcess_CCincl_2D->Fill("Other","Other",0);

   // Loop through MC tree and fill plots
   for (int i = 0; i < t_bnbcos->GetEntries(); i++){
      t_bnbcos->GetEntry(i);
      Calcvars(&mc_vars);
      std::vector<std::vector<double>> Varstoplot = GetVarstoplot(&mc_vars);
      std::vector<std::vector<double>> Cutvarstoplot = GetCutvarstoplot(&mc_vars);

      
      // Determine if event is selected
      bool isSelected = true;
      for (size_t i_cut = 0; i_cut < Cutvarstoplot.size(); i_cut++){
         if(!IsEventSelected(CutValues.at(i_cut), Cutvarstoplot.at(i_cut), KeepBelowCut.at(i_cut), mc_vars.Marco_selected, *(mc_vars.TPCObj_PFP_isDaughter), OnlyDaughters.at(i_cut), TracksNeeded.at(i_cut))) {
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

      // EndProcess Truth Study
      std::string EndProcess_string;
      double EndProcess_double;
      std::string EndProcess_muon;
      std::string EndProcess_pion;
      for (size_t i_mcp = 0; i_mcp < mc_vars.MCP_endprocess->size(); i_mcp++) {
         if(!mc_vars.MCP_isContained->at(i_mcp)) EndProcess_string = "Exiting";
         else EndProcess_string = mc_vars.MCP_endprocess->at(i_mcp);

         auto EndProcess_iter = std::find(EndProcess_vect.begin(),EndProcess_vect.end(),EndProcess_string);
         if (EndProcess_iter == EndProcess_vect.end()) {
            EndProcess_vect.emplace_back(EndProcess_string);
            EndProcess_double = (double)(EndProcess_vect.size()-1);
         }
         else EndProcess_double = 2*(double)(EndProcess_iter - EndProcess_vect.begin());

         if(isSelected) EndProcess_selected -> Fill((PDGCode)mc_vars.MCP_PDG->at(i_mcp), EndProcess_double);
         if(mc_vars.Marco_selected) EndProcess_CCincl -> Fill((PDGCode)mc_vars.MCP_PDG->at(i_mcp), EndProcess_double);

         if (mc_vars.Truth_topology == kCC1piplus0p || mc_vars.Truth_topology == kCC1piplus1p || mc_vars.Truth_topology == kCC1piplusNp){
            if(mc_vars.MCP_PDG->at(i_mcp) == kMuMinus) {
               if(EndProcess_muon.empty()) EndProcess_muon = EndProcess_string;
               else std::cout << "More than one muon in true CC1pi+ event" << std::endl;
            }
            else if(mc_vars.MCP_PDG->at(i_mcp) == kPiPlus) {
               if(EndProcess_pion.empty()) EndProcess_pion = EndProcess_string;
               else std::cout << "More than one pion in true CC1pi+ event" << std::endl;
            }
         }

      }

      if (!EndProcess_muon.empty() && !EndProcess_pion.empty()) {
         if(!(EndProcess_muon=="Exiting" || EndProcess_muon=="FastScintillation" || EndProcess_muon=="Decay")) EndProcess_muon = "Other";
         if(!(EndProcess_pion=="Exiting" || EndProcess_pion=="FastScintillation" || EndProcess_pion=="Decay")) EndProcess_pion = "Other";

         if(isSelected) EndProcess_selected_2D->Fill(EndProcess_muon.c_str(),EndProcess_pion.c_str(),1);
         if(mc_vars.Marco_selected) EndProcess_CCincl_2D->Fill(EndProcess_muon.c_str(),EndProcess_pion.c_str(),1);
      }


   } // end loop over entries in tree

   // -------------------- Now make all the plots

   gStyle->SetOptStat(0);

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
      c1->Clear();

      delete c1;
   }
   
   TCanvas *c1 = new TCanvas("c1","c1");
   EndProcess_selected->DrawOverlayMuPi(0.,c1);
   c1->Print(std::string("EndProcess_selected.png").c_str());
   c1->Clear();

   EndProcess_CCincl->DrawOverlayMuPi(0.,c1);
   c1->Print(std::string("EndProcess_CCincl.png").c_str());
   c1->Clear();

   EndProcess_selected_2D->Draw("TEXT COLZ");
   c1->Print(std::string("EndProcess_selected_2D.png").c_str());
   c1->Clear();

   EndProcess_CCincl_2D->Draw("TEXT COLZ");
   c1->Print(std::string("EndProcess_CCincl_2D.png").c_str());
   c1->Clear();

   delete c1;

   // Print out integrals
   SelectedEvents->GetHistIntegrals();

   // Print out the "map" of doubles to EndProcesses
   for(size_t i_end = 0; i_end < EndProcess_vect.size(); i_end++) {
      std::cout << 2*i_end << " is EndProcess " << EndProcess_vect.at(i_end) << std::endl;
   }
}
