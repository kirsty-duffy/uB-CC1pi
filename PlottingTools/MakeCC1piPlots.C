#include "StackedHistPDGCode.h"
#include "StackedHistTopology.h"
#include "../Algorithms/TopologyEnums.h"
#include "CC1pi_treevars.h"
#include "CC1pi_treevars.cxx"
#include "CC1pi_plotvars_def.h"
#include "CC1pi_cuts.h"
#include "CC1pi_cuts.cxx"
#include "StoppingParticlePlots.h"

// What variables do we want these plots as a function of?
std::vector<CC1piPlotVars> GetVarstoplot(treevars *vars){
   std::vector<CC1piPlotVars> varstoplot = {
      Var_TPCObj_PFP_lnLmipoverp(vars)
      ,Var_TPCObj_PFP_Lmumipovermumipp(vars)
      ,Var_TPCObj_PFP_BrokenTrackAngle(vars)
      ,Var_TPCObj_PFP_track_residual_mean_low(vars)
      ,Var_TPCObj_PFP_track_residual_mean_high(vars)
      ,Var_TPCObj_PFP_track_residual_std_low(vars)
      ,Var_TPCObj_PFP_track_residual_std_high(vars)
      ,Var_TPCObj_PFP_track_perc_used_hits(vars)
      ,Var_TPCObj_PFP_VtxTrackDist(vars)
      ,Var_TPCObj_PFP_isContained(vars)
      ,Var_TPCObj_PFP_track_dEdx_truncmean_start(vars)
      ,Var_TPCObj_DaughterTracks_Order_dEdxtr(vars)
      ,Var_TPCObj_DaughterTracks_Order_dEdxtr_selMIPs(vars)
      ,Var_TPCObj_DaughterTracks_Order_trklen_selMIPs(vars)
   };
   return varstoplot;
}

// ---------------------------------------------------- //
//  Now the function starts
// ---------------------------------------------------- //

void MakeCC1piPlots(std::string mcfile, bool MakeKinkFindingPlots=false){

   gStyle->SetTitleX(0.5);
   gStyle->SetTitleAlign(23);
   gStyle->SetOptTitle(1);
   gStyle->SetTitleBorderSize(0.);

   TFile *f_bnbcos = new TFile(mcfile.c_str(), "read");
   TTree *t_bnbcos = (TTree*)f_bnbcos->Get("cc1piselec/outtree");
   TTree *t_bnbcos_friend = (TTree*)f_bnbcos->Get("StoppingParticleTagger/StopStoppingTaggerTree");

   if (t_bnbcos->GetEntries() != t_bnbcos_friend->GetEntries()){
     std::cout << "ERROR: cc1piselec/outtree has " << t_bnbcos->GetEntries() << " entries, and StoppingParticleTagger/StoppingTaggerTree has " << t_bnbcos_friend->GetEntries() << " entries. Cannot make friend tree, exiting." << std::endl;
     return;
   }
   t_bnbcos->AddFriend(t_bnbcos_friend);

   treevars mc_vars;
   settreevars(t_bnbcos,&mc_vars);

   // Sanity check: the plot vectors should be the same size
   t_bnbcos->GetEntry(0);
   Calcvars(&mc_vars);
   std::vector<CC1piPlotVars> varstoplot_dummy = GetVarstoplot(&mc_vars);
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

   std::cout << std::endl << "--- Considering the following variables --- " << std::endl;
   for (size_t i_var=0; i_var < nplots; i_var++){
         std::cout << varstoplot_dummy.at(i_var).histtitle.substr(1,varstoplot_dummy.at(i_var).histtitle.size()-2) << std::endl;
   }
   std::cout << "------" << std::endl << std::endl;

    std::cout << "--- Applying the following cuts --- " << std::endl;
    for (size_t i_var=0; i_var < cutvars_dummy.size(); i_var++){
          std::cout << cutvars_dummy.at(i_var).histtitle.substr(1,cutvars_dummy.at(i_var).histtitle.size()-2) << std::endl;
    }
    std::cout << "------" << std::endl;

   // Dummy plot in order to get percentage of different topologies that are selected
   StackedHistTopology *SelectedEvents = new StackedHistTopology("SelectedEvents",";;",1,0,1);

   // Loop through MC tree and fill plots
   for (int i = 0; i < t_bnbcos->GetEntries(); i++){
      t_bnbcos->GetEntry(i);
      Calcvars(&mc_vars);
      std::vector<CC1piPlotVars> Varstoplot = GetVarstoplot(&mc_vars);

      bool isSelected = IsEventSelected(&mc_vars);

      if(isSelected) SelectedEvents->Fill((NuIntTopology)mc_vars.Truth_topology,0);

      // Loop over Varstoplot
      for (size_t i_h = 0; i_h < nplots; i_h++){
        std::vector<double> vartoplot = *(Varstoplot.at(i_h).Var);
         // Loop over tracks
         for (size_t i_tr = 0; i_tr < vartoplot.size(); i_tr++){
            // if (!(mc_vars.TPCObj_PFP_isDaughter->at(i_tr) && bool(mc_vars.TPCObj_PFP_track_passesMIPcut->at(i_tr)))) continue;

            mc_hists_cc1pi_pdg_beforecuts[i_h]->Fill((PDGCode)mc_vars.TPCObj_PFP_truePDG->at(i_tr),vartoplot.at(i_tr));
            mc_hists_cc1pi_top_beforecuts[i_h]->Fill((NuIntTopology)mc_vars.Truth_topology,vartoplot.at(i_tr),1.0/vartoplot.size());

            if(isSelected) {
               mc_hists_cc1pi_pdg_aftercuts[i_h]->Fill((PDGCode)mc_vars.TPCObj_PFP_truePDG->at(i_tr),vartoplot.at(i_tr));
               mc_hists_cc1pi_top_aftercuts[i_h]->Fill((NuIntTopology)mc_vars.Truth_topology,vartoplot.at(i_tr),1.0/vartoplot.size());

               // For selected events, make hit dQds and local linearity plots for the two MIP-like tracks.
               // Do this for the first 100 events only because otherwise we'll just have a ridiculous number of plots
               if (MakeKinkFindingPlots && i<100){
                 if (mc_vars.TPCObj_PFP_isDaughter->at(i_tr) && bool(mc_vars.TPCObj_PFP_track_passesMIPcut->at(i_tr))){
                   TCanvas *c0 = new TCanvas("","",400,500);
                   MakeStoppingParticlePlots_SingleTrack(c0, &mc_vars, i_tr);
                   c0->Print(std::string(std::string("StoppingParticlePlots_TPCObj")+i+std::string("_track")+i_tr+std::string(".pdf")).c_str());
                 }
               }

            } // end if(isSelected)

         } // end loop over tracks

      } // end loop over Varstoplot

   } // end loop over entries in tree

   // -------------------- Now make all the plots

   for (size_t i_h=0; i_h < nplots; i_h++){
      TCanvas *c1 = new TCanvas("c1","c1");
      std::string printname = std::string(varstoplot_dummy.at(i_h).histname+".png");

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
