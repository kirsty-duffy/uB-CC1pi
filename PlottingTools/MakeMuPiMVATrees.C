#include "../Algorithms/TopologyEnums.h"
#include "../Algorithms/PDGEnums.h"
#include "CC1pi_treevars.h"
#include "CC1pi_treevars.cxx"

#include <boost/algorithm/string.hpp>


struct MVAvars_mupi {
   double length_over_startend; // Track length vs start-end distance (division)
   double track_length; // Track length
   double track_length_over_longest; // Track length divided by length of longest MIP track in event
   double ndaughters; // No. reco daughters
   double lnLmipovermu; // MIP vs mu LLH (is there a Bragg peak?)
   double lnLmipoverpi; // MIP vs pi LLH (is there a Bragg peak? Is this duplicating above?)
   double pandoraclassedastrack; // Pandora track/shower classification (watch out for problems with DIC and the fact that this is really a bool)
   double perc_used_hits; // Percentage of used hits in track
   double residuals_mean; // Mean of track residuals (watch out for problems with DIC)
   double residuals_stddev; // Std dev of track residuals (watch out for problems with DIC)
   double MCS_pi_maxScatter; // Largest MCS scattering angle
   double MCS_pi_meanScatter; // Mean MCS scattering angle
   double n_unused_hits_nearend; // No. hits near the end of the track not matched to any other track
   double unmatched_charge_nearend_plane2; // Total charge of hits on plane 2 near end of track not matched to any other track
   // TODO Energy near end of track not associated with a track
};

bool CheckMVAvars_mupi(MVAvars_mupi *MVA_vars){
  bool returnval = true;

  if (MVA_vars->length_over_startend<0) returnval=false;
  if (MVA_vars->track_length==-9999) returnval=false;
  if (MVA_vars->track_length_over_longest==-9999) returnval=false;
  if (MVA_vars->ndaughters==-9999) returnval=false;
  if (MVA_vars->lnLmipovermu==-9999) returnval=false;
  if (MVA_vars->lnLmipoverpi==-9999) returnval=false;
  if (MVA_vars->pandoraclassedastrack==-9999) returnval=false;
  if (MVA_vars->perc_used_hits==-9999) returnval=false;
  if (MVA_vars->residuals_mean==-9999) returnval=false;
  if (MVA_vars->residuals_stddev==-9999) returnval=false;
  if (MVA_vars->MCS_pi_maxScatter==-9999) returnval=false;
  if (MVA_vars->MCS_pi_meanScatter==-9999) returnval=false;

  return returnval;
}

void MakeMVATrees_mupi(TTree *muon_tree, TTree *pion_tree, MVAvars_mupi *MVA_vars, treevars *vars) {

   for(int i_track=0; i_track < vars->TPCObj_PFP_truePDG->size(); i_track++){

      //Only daughters get considered as possible MIP candidates, so also only train with them
      if(!(vars->TPCObj_PFP_isDaughter->at(i_track))) continue;

      // Only train on neutrino-induced particles
      if(vars->Truth_topology == kCosmic || vars->Truth_topology == kMixed || vars->Truth_topology == kUnknown) continue;

      // Only train on (and apply to) contained tracks
      if (vars->TPCObj_PFP_track_isContained->at(i_track)!=1) continue;

      // Quality pre-cut: tracks parallel to collection plane have bad dE/dx
      if(vars->TPCObj_PFP_track_dEdx_truncmean_start->at(i_track) < 1) continue;

      // Quality pre-cut: tracks further than 5 cm away from the vertex are typically misreconstructed (or background)
      if(vars->TPCObj_PFP_VtxTrackDist->at(i_track) > 5) continue;

      // Now set variables
      TVector3 start(vars->TPCObj_PFP_track_start->at(i_track).at(0),vars->TPCObj_PFP_track_start->at(i_track).at(1),vars->TPCObj_PFP_track_start->at(i_track).at(2));
      TVector3 end(vars->TPCObj_PFP_track_end->at(i_track).at(0),vars->TPCObj_PFP_track_end->at(i_track).at(1),vars->TPCObj_PFP_track_end->at(i_track).at(2));
      double startenddist = (end-start).Mag();
      MVA_vars->length_over_startend = vars->TPCObj_PFP_track_length->at(i_track)/startenddist;

      MVA_vars->track_length = vars->TPCObj_PFP_track_length->at(i_track);

      if (vars->TPCObj_LeadingMIPtrackIndex>=0){
        MVA_vars->track_length_over_longest=vars->TPCObj_PFP_track_length->at(i_track)/vars->TPCObj_PFP_track_length->at(vars->TPCObj_LeadingMIPtrackIndex);
      }
      else{
        MVA_vars->track_length_over_longest=-9999;
      }

      MVA_vars->ndaughters = vars->TPCObj_PFP_ndaughters->at(i_track);
      MVA_vars->lnLmipovermu = vars->TPCObj_PFP_lnLmipovermu->at(i_track);
      MVA_vars->lnLmipoverpi = vars->TPCObj_PFP_lnLmipoverpi->at(i_track);
      MVA_vars->pandoraclassedastrack = vars->TPCObj_PFP_PandoraClassedAsTrack_double->at(i_track);
      MVA_vars->perc_used_hits = vars->TPCObj_PFP_track_perc_used_hits->at(i_track);
      MVA_vars->residuals_mean = vars->TPCObj_PFP_track_residual_mean->at(i_track);
      MVA_vars->residuals_stddev = vars->TPCObj_PFP_track_residual_std->at(i_track);
      MVA_vars->MCS_pi_maxScatter = vars->TPCObj_PFP_track_MCS_pi_maxScatter->at(i_track);
      MVA_vars->MCS_pi_meanScatter = vars->TPCObj_PFP_track_MCS_pi_meanScatter->at(i_track);

      double n_unused_hits = 0;
      double unused_charge = 0;
      // Plane 0: nhits only
      for (size_t i_h=0; i_h<vars->TPCObj_PFP_track_unusedhits_endtimedist_plane0->at(i_track).size(); i_h++){
        double dt = (double)vars->TPCObj_PFP_track_unusedhits_endtimedist_plane0->at(i_track).at(i_h);
        double dwire = vars->TPCObj_PFP_track_unusedhits_endwiredist_plane0->at(i_track).at(i_h);

        // Fairly loose cut: look at 5 time ticks and 10 wires in each direction
        if (TMath::Abs(dt)<5 && TMath::Abs(dwire)<10){
          n_unused_hits++;
        }
      }
      // Plane 1: nhits only
      for (size_t i_h=0; i_h<vars->TPCObj_PFP_track_unusedhits_endtimedist_plane1->at(i_track).size(); i_h++){
        double dt = (double)vars->TPCObj_PFP_track_unusedhits_endtimedist_plane1->at(i_track).at(i_h);
        double dwire = vars->TPCObj_PFP_track_unusedhits_endwiredist_plane1->at(i_track).at(i_h);

        // Fairly loose cut: look at 5 time ticks and 10 wires in each direction
        if (TMath::Abs(dt)<5 && TMath::Abs(dwire)<10){
          n_unused_hits++;
        }
      }
      // Plane 2: nhits and charge
      for (size_t i_h=0; i_h<vars->TPCObj_PFP_track_unusedhits_endtimedist_plane2->at(i_track).size(); i_h++){
        double dt = (double)vars->TPCObj_PFP_track_unusedhits_endtimedist_plane2->at(i_track).at(i_h);
        double dwire = vars->TPCObj_PFP_track_unusedhits_endwiredist_plane2->at(i_track).at(i_h);

        // Fairly loose cut: look at 5 time ticks and 10 wires in each direction
        if (TMath::Abs(dt)<5 && TMath::Abs(dwire)<10){
          n_unused_hits++;
          unused_charge+=vars->TPCObj_PFP_track_unusedhits_charge_plane2->at(i_track).at(i_h);
        }
      }
      MVA_vars->n_unused_hits_nearend = n_unused_hits;
      MVA_vars->unmatched_charge_nearend_plane2 = unused_charge;

      //Don't train on bogus values (though also, we should probably try to understand why there are bogus values)
      if(CheckMVAvars_mupi(MVA_vars)==false) continue;

      if(vars->TPCObj_PFP_truePDG->at(i_track)==13) {
         muon_tree->Fill();
      }
      else if(vars->TPCObj_PFP_truePDG->at(i_track)==211) {
         pion_tree->Fill();
      }
   }
}


void MakeMuPiMVATrees(std::string mcfile){

  // Note: MVA trees are made for MC only
   TFile f_MVA("MVA_Trees_mupi.root", "RECREATE");
   MVAvars_mupi *MVA_vars = new MVAvars_mupi();
   TTree *muon = new TTree("muon","muon");
   TTree *pion = new TTree("pion","pion");

   muon -> Branch("length_over_startend", &(MVA_vars->length_over_startend));
   muon -> Branch("track_length", &(MVA_vars->track_length));
   muon -> Branch("track_length_over_longest", &(MVA_vars->track_length_over_longest));
   muon -> Branch("ndaughters", &(MVA_vars->ndaughters));
   muon -> Branch("lnLmipovermu", &(MVA_vars->lnLmipovermu));
   muon -> Branch("lnLmipoverpi", &(MVA_vars->lnLmipoverpi));
   muon -> Branch("pandoraclassedastrack", &(MVA_vars->pandoraclassedastrack));
   muon -> Branch("perc_used_hits", &(MVA_vars->perc_used_hits));
   muon -> Branch("residuals_mean", &(MVA_vars->residuals_mean));
   muon -> Branch("residuals_stddev", &(MVA_vars->residuals_stddev));
   muon -> Branch("MCS_pi_maxScatter", &(MVA_vars->MCS_pi_maxScatter));
   muon -> Branch("MCS_pi_meanScatter", &(MVA_vars->MCS_pi_meanScatter));
   muon -> Branch("n_unused_hits_nearend", &(MVA_vars->n_unused_hits_nearend));
   muon -> Branch("unmatched_charge_nearend_plane2", &(MVA_vars->unmatched_charge_nearend_plane2));

   pion -> Branch("length_over_startend", &(MVA_vars->length_over_startend));
   pion -> Branch("track_length", &(MVA_vars->track_length));
   pion -> Branch("track_length_over_longest", &(MVA_vars->track_length_over_longest));
   pion -> Branch("ndaughters", &(MVA_vars->ndaughters));
   pion -> Branch("lnLmipovermu", &(MVA_vars->lnLmipovermu));
   pion -> Branch("lnLmipoverpi", &(MVA_vars->lnLmipoverpi));
   pion -> Branch("pandoraclassedastrack", &(MVA_vars->pandoraclassedastrack));
   pion -> Branch("perc_used_hits", &(MVA_vars->perc_used_hits));
   pion -> Branch("residuals_mean", &(MVA_vars->residuals_mean));
   pion -> Branch("residuals_stddev", &(MVA_vars->residuals_stddev));
   pion -> Branch("MCS_pi_maxScatter", &(MVA_vars->MCS_pi_maxScatter));
   pion -> Branch("MCS_pi_meanScatter", &(MVA_vars->MCS_pi_meanScatter));
   pion -> Branch("n_unused_hits_nearend", &(MVA_vars->n_unused_hits_nearend));
   pion -> Branch("unmatched_charge_nearend_plane2", &(MVA_vars->unmatched_charge_nearend_plane2));


   TFile *f_bnbcos = new TFile(mcfile.c_str(), "read");
   TTree *t_bnbcos = (TTree*)f_bnbcos->Get("cc1piselec/outtree");

    if (!t_bnbcos){
      std::cout << "Error: did not find tree cc1piselec/outtree in file" << std::endl;
      return;
    }

   treevars mc_vars;
   settreevars(t_bnbcos,&mc_vars);

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
      MakeMVATrees_mupi(muon, pion, MVA_vars, &mc_vars);

   } // end loop over entries in tree

   f_MVA.Write();
}
