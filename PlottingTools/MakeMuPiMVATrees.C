#include "../Algorithms/TopologyEnums.h"
#include "../Algorithms/PDGEnums.h"
#include "CC1pi_treevars.h"
#include "CC1pi_treevars.cxx"

#include <boost/algorithm/string.hpp>

void MakeMVATrees_mupi(TTree *muon_cont, TTree *pion_cont, TTree *muon_exit, TTree *pion_exit, treevars *vars) {

   for(int i_track=0; i_track < vars->TPCObj_PFP_truePDG->size(); i_track++){

      //Only daughters get considered as possible MIP candidates, so also only train with them
      if(!(vars->TPCObj_PFP_isDaughter->at(i_track))) continue;

      // Only train on neutrino-induced particles
      if(vars->Truth_topology == kCosmic || vars->Truth_topology == kMixed || vars->Truth_topology == kUnknown) continue;

      // Quality pre-cut: tracks parallel to collection plane have bad dE/dx
      if(vars->TPCObj_PFP_track_dEdx_truncmean_start->at(i_track) < 1) continue;

      // Quality pre-cut: tracks further than 5 cm away from the vertex are typically misreconstructed (or background)
      if(vars->TPCObj_PFP_VtxTrackDist->at(i_track) > 5) continue;

      // Now set variables
      SetMVAvars_mupiBDT(vars,i_track);


      //Don't train on bogus values (though also, we should probably try to understand why there are bogus values)
      if(vars->mupiBDT.CheckMVAvars()==false) continue;



      // Separate trees (and BDTs) for contained and uncontained tracks
      if (vars->TPCObj_PFP_track_isContained->at(i_track)==1){
        if(vars->TPCObj_PFP_truePDG->at(i_track)==13) {
          muon_cont->Fill();
        }
        else if(vars->TPCObj_PFP_truePDG->at(i_track)==211) {
          pion_cont->Fill();
        }
      }
      else if (vars->TPCObj_PFP_track_isContained->at(i_track)==0){
        if(vars->TPCObj_PFP_truePDG->at(i_track)==13) {
          muon_exit->Fill();
        }
        else if(vars->TPCObj_PFP_truePDG->at(i_track)==211) {
          pion_exit->Fill();
        }
      } // end if contained/exiting
   } // end loop over tracks
}


void MakeMuPiMVATrees(std::string mcfile){

  TFile *f_bnbcos = new TFile(mcfile.c_str(), "read");
  TTree *t_bnbcos = (TTree*)f_bnbcos->Get("cc1piselec/outtree");

   if (!t_bnbcos){
     std::cout << "Error: did not find tree cc1piselec/outtree in file" << std::endl;
     return;
   }

  treevars mc_vars;
  settreevars(t_bnbcos,&mc_vars);

  // Note: MVA trees are made for MC only
   TFile f_MVA("MVA_Trees_mupi.root", "RECREATE");

   TTree *muon_cont = new TTree("muon_cont","muon_cont");
   InitialiseTree(muon_cont,&mc_vars);
   TTree *pion_cont = new  TTree("pion_cont","pion_cont");
   InitialiseTree(pion_cont,&mc_vars);
   TTree *muon_exit = new TTree("muon_exit","muon_exit");
   InitialiseTree(muon_exit,&mc_vars);
   TTree *pion_exit = new  TTree("pion_exit","pion_exit");
   InitialiseTree(pion_exit,&mc_vars);


   // Setup for MIP-proton BDT
   std::string BookMVAType = "BDTG";
   std::string BookMVALoc = "/uboone/app/users/ddevitt/LArSoft_v06_26_01_14_uboonecode_v06_26_01_22/srcs/uboonecode/uboone/CC1pi/MVA/dataset_vtxtrackprecut/weights/TMVAClassification_BDTG.weights.xml";

   TMVA::Reader fReader("");
   fReader.AddVariable("dEdx_truncmean_start", &(mc_vars.float_dEdx_truncmean_start));
   fReader.AddVariable("nhits", &(mc_vars.float_nhits));
   fReader.AddVariable("lnLmipoverp", &(mc_vars.float_lnLmipoverp));
   fReader.BookMVA(BookMVAType.c_str(), BookMVALoc.c_str());


   // ----------------- MC

   // Loop through MC tree and fill plots
   for (int i = 0; i < t_bnbcos->GetEntries(); i++){
     if (i%1000==0) std::cout << "MC: " << i << "/" << t_bnbcos->GetEntries() << std::endl;

      t_bnbcos->GetEntry(i);

      Calcvars(&mc_vars, &fReader);
      MakeMVATrees_mupi(muon_cont, pion_cont, muon_exit, pion_exit, &mc_vars);

   } // end loop over entries in tree

   std::cout << "Made trees: " << std::endl;
   std::cout << "   " << muon_cont->GetEntries() << " contained muons" << std::endl;
   std::cout << "   " << pion_cont->GetEntries() << " contained pions" << std::endl;
   std::cout << "   " << muon_exit->GetEntries() << " exiting muons" << std::endl;
   std::cout << "   " << pion_exit->GetEntries() << " exiting pions" << std::endl;

   f_MVA.Write();
}
