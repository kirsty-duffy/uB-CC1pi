
#include "../Algorithms/TopologyEnums.h"
#include "CC1pi_treevars.h"
#include "CC1pi_treevars.cxx"
#include "CC1pi_plotvars_def.h"
#include "CC1pi_cuts.h"
#include "CC1pi_cuts.cxx"
#include "CC1pi_EffPurHists.h"
#include "CC1pi_EffPurHists2D.h"

#include "TFile.h"
#include "TTree.h"


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
   TTree *t_bnbcos_friend = (TTree*)f_bnbcos->Get("StoppingParticleTagger/StoppingTaggerTree");

    if (!t_bnbcos){
      std::cout << "Error: did not find tree cc1piselec/outtree in file" << std::endl;
      return;
    }
    if (!t_bnbcos_friend){
      std::cout << "Error: did not find tree StoppingParticleTagger/StoppingTaggerTree in file" << std::endl;
      return;
    }
   if (t_bnbcos->GetEntries() != t_bnbcos_friend->GetEntries()){
     std::cout << "ERROR: cc1piselec/outtree has " << t_bnbcos->GetEntries() << " entries, and StoppingParticleTagger/StoppingTaggerTree has " << t_bnbcos_friend->GetEntries() << " entries. Cannot make friend tree, exiting." << std::endl;
     return;
   }
   t_bnbcos->AddFriend(t_bnbcos_friend);

   treevars mc_vars;
   settreevars(t_bnbcos,&mc_vars);

   TMVA::Reader fReader("");
   fReader.AddVariable("dEdx_truncmean_start", &(mc_vars.float_dEdx_truncmean_start));
   fReader.AddVariable("VtxTrackDist", &(mc_vars.float_VtxTrackDist));
   fReader.AddVariable("nhits", &(mc_vars.float_nhits));
   fReader.AddVariable("lnLmipoverp", &(mc_vars.float_lnLmipoverp));
   fReader.BookMVA("BDT", "/uboone/app/users/ddevitt/LArSoft_v06_26_01_10/srcs/uboonecode/uboone/CC1pi/MVA/dataset_mupi_4MIPvars_equalweighting_binsadjusted/weights/TMVAClassification_BDT.weights.xml");

   // Get vector of cuts we want to plot
   // We want both the "normal" cut vars and the ones that define a MIP
   t_bnbcos->GetEntry(0);
   Calcvars(&mc_vars, &fReader);
   std::vector<CC1piPlotVars> cutvars_notMIPs = GetCutVars(&mc_vars);
   std::vector<CC1piPlotVars> cutvars_MIPs = GetMIPCutVars(&mc_vars);
   std::vector<CC1piPlotVars> varstoplot_dummy;
   for (size_t i_cut=0; i_cut<cutvars_notMIPs.size(); i_cut++){
      varstoplot_dummy.push_back(cutvars_notMIPs.at(i_cut));
   }
   for (size_t i_cut=0; i_cut<cutvars_MIPs.size(); i_cut++){
      varstoplot_dummy.push_back(cutvars_MIPs.at(i_cut));
   }


   // ----------------- MC

   // Make histograms to fill
   const size_t nplots = varstoplot_dummy.size();

   histCC1piselEffPur *mc_hists_cc1pieffpur[nplots];
   for (size_t i_h=0; i_h<nplots; i_h++){
      std::string histtitle_i = varstoplot_dummy.at(i_h).histtitle;
      std::string histname_i = varstoplot_dummy.at(i_h).histname;
      std::vector<double> bins_i = varstoplot_dummy.at(i_h).bins;

      mc_hists_cc1pieffpur[i_h] = new histCC1piselEffPur(std::string("hCC1pieffpur_")+histname_i,histtitle_i,bins_i.at(0),bins_i.at(1),bins_i.at(2));
   }

   histCC1piselEffPur *Nminus1plots = new histCC1piselEffPur("hCC1pi_nminus1","N-1 plot;;Efficiency, purity",nplots+3,0,nplots+3);

   std::cout << "--- Considering the following cuts --- " << std::endl;
   for (size_t i_bin=1; i_bin < nplots+4; i_bin++){
      if (i_bin<nplots+1){
         std::string tmpname = varstoplot_dummy.at(i_bin-1).histtitle.substr(1,varstoplot_dummy.at(i_bin-1).histtitle.size()-2);
         std::cout << tmpname << std::endl; Nminus1plots->SetBinLabel(i_bin,tmpname);
      }
      if (i_bin==nplots+1) Nminus1plots->SetBinLabel(i_bin,"All cuts");
      if (i_bin==nplots+2) Nminus1plots->SetBinLabel(i_bin, "CC incl. 2 tracks");
      if (i_bin==nplots+3) Nminus1plots->SetBinLabel(i_bin, "CC incl.");
   }
   std::cout << "------" << std::endl;

   // Loop through MC tree and fill plots
   for (int i = 0; i < t_bnbcos->GetEntries(); i++){
//   for (int i = 0; i < 1000; i++){
      if (i%1000==0 || i==t_bnbcos->GetEntries()-1) std::cout << i << "/" << t_bnbcos->GetEntries() << std::endl;
      t_bnbcos->GetEntry(i);
      Calcvars(&mc_vars, &fReader);

       for (size_t i_h = 2; i_h < 3; i_h++){
          // std::cout << i_h << std::endl;
          FillCC1piEffPurHist(mc_hists_cc1pieffpur[i_h], &mc_vars, &fReader, (int)i_h);
       }

      // Now fill histograms for N-1 plots
      FillNminus1EffPurHist(Nminus1plots,&mc_vars,&fReader);

   } // end loop over entries in tree

   // -------------------- Now make all the plots
   std::cout << "Drawing plots..." << std::endl;
    for (size_t i_h=2; i_h < 3; i_h++){
       TCanvas *c1 = new TCanvas("c1","c1");
   
       DrawCC1piMCEffPur(c1, mc_hists_cc1pieffpur[i_h]);
       std::string printname = std::string(varstoplot_dummy.at(i_h).histname+".png");
       c1->Print(std::string(std::string("CC1pi_effpur_")+printname).c_str());
   
       delete c1;
    }
   TCanvas *c1 = new TCanvas("c1","c1");
   DrawCC1piMCEffPur(c1, Nminus1plots,"hist",true);
   c1->Print("Nminus1plots.png");

   delete c1;
}
