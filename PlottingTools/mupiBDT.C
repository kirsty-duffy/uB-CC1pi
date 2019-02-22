#include "CC1pi_treevars.h"
#include "TMVA/Reader.h"

 bool mupiBDT::CheckMVAvars(){
   bool returnval = true;

   // std::cout << "length_over_startend = " << length_over_startend << std::endl
   // << "track_length = " << track_length << std::endl
   // << "track_length_over_longest = " << track_length_over_longest << std::endl
   // << "ndaughters = " << ndaughters << std::endl
   // << "lnLmipovermu = " << lnLmipovermu << std::endl
   // << "lnLmipoverpi = " << lnLmipoverpi << std::endl
   // << "pandoraclassedastrack = " << pandoraclassedastrack << std::endl
   // << "perc_used_hits = " << perc_used_hits << std::endl
   // << "residuals_mean = " << residuals_mean << std::endl
   // << "residuals_stddev = " << residuals_stddev << std::endl
   // << "MCS_pi_maxScatter = " << MCS_pi_maxScatter << std::endl
   // << "MCS_pi_meanScatter = " << MCS_pi_meanScatter << std::endl << std::endl;

   if (length_over_startend<0 || isnan(length_over_startend)) returnval=false;
   if (track_length==-9999 || isnan(track_length)) returnval=false;
   if (track_length_over_longest<0 || isnan(track_length_over_longest)) returnval=false;
   if (ndaughters==-9999 || isnan(ndaughters)) returnval=false;
   if (lnLmipovermu==-9999 || isnan(lnLmipovermu)) returnval=false;
   if (lnLmipoverpi==-9999 || isnan(lnLmipoverpi)) returnval=false;
   if (pandoraclassedastrack==-9999 || isnan(pandoraclassedastrack)) returnval=false;
   if (perc_used_hits==-9999 || isnan(perc_used_hits)) returnval=false;
   if (residuals_mean==-9999 || isnan(residuals_mean)) returnval=false;
   if (residuals_stddev==-9999 || isnan(residuals_stddev)) returnval=false;
   if (MCS_pi_maxScatter==-9999 || isnan(MCS_pi_maxScatter)) returnval=false;
   if (MCS_pi_meanScatter==-9999 || isnan(MCS_pi_meanScatter)) returnval=false;

   return returnval;
}

void mupiBDT::initialise_BDT_contained(){
   fReader_contained = new TMVA::Reader("");
   fReader_contained->AddVariable("length_over_startend",&length_over_startend);
   fReader_contained->AddVariable("track_length_over_longest",&track_length_over_longest);
   fReader_contained->AddVariable("ndaughters",&ndaughters);
   fReader_contained->AddVariable("lnLmipovermu",&lnLmipovermu);
   fReader_contained->AddVariable("perc_used_hits",&perc_used_hits);
   // fReader_contained->AddVariable("residuals_mean",&residuals_mean);
   // fReader_contained->AddVariable("residuals_stddev",&residuals_stddev);
   fReader_contained->AddVariable("MCS_pi_maxScatter",&MCS_pi_maxScatter);
   fReader_contained->AddVariable("MCS_pi_meanScatter",&MCS_pi_meanScatter);
   fReader_contained->AddVariable("n_unused_hits_nearend",&n_unused_hits_nearend);
   fReader_contained->AddVariable("unmatched_charge_nearend_plane2",&unmatched_charge_nearend_plane2);

   fReader_contained->AddSpectator("track_length", &track_length);
   fReader_contained->AddSpectator("pandoraclassedastrack", &pandoraclassedastrack);
   fReader_contained->AddSpectator("lnLmipoverpi", &lnLmipoverpi);

   fReader_contained->BookMVA(BookMVAType_contained.c_str(), BookMVALoc_contained.c_str());
}

void mupiBDT::initialise_BDT_exiting(){
   fReader_exiting = new TMVA::Reader("");
   fReader_exiting->AddVariable("length_over_startend",&length_over_startend);
   fReader_exiting->AddVariable("track_length_over_longest",&track_length_over_longest);
   fReader_exiting->AddVariable("perc_used_hits",&perc_used_hits);
   // fReader_exiting->AddVariable("residuals_mean",&residuals_mean);
   // fReader_exiting->AddVariable("residuals_stddev",&residuals_stddev);
   fReader_exiting->AddVariable("MCS_pi_maxScatter",&MCS_pi_maxScatter);
   fReader_exiting->AddVariable("MCS_pi_meanScatter",&MCS_pi_meanScatter);

   fReader_exiting->AddSpectator("track_length", &track_length); fReader_exiting->AddSpectator("pandoraclassedastrack", &pandoraclassedastrack);
   fReader_exiting->AddSpectator("lnLmipovermu", &lnLmipovermu);
   fReader_exiting->AddSpectator("lnLmipoverpi", &lnLmipoverpi);
   fReader_exiting->AddSpectator("ndaughters", &ndaughters);
   fReader_exiting->AddSpectator("n_unused_hits_nearend", &n_unused_hits_nearend);
   fReader_exiting->AddSpectator("unmatched_charge_nearend_plane2", &unmatched_charge_nearend_plane2);

   fReader_exiting->BookMVA(BookMVAType_exiting.c_str(), BookMVALoc_exiting.c_str());
}



void SetMVAvars_mupiBDT(treevars *vars, int i_track){
 TVector3 start(vars->TPCObj_PFP_track_start->at(i_track).at(0),vars->TPCObj_PFP_track_start->at(i_track).at(1),vars->TPCObj_PFP_track_start->at(i_track).at(2));
 TVector3 end(vars->TPCObj_PFP_track_end->at(i_track).at(0),vars->TPCObj_PFP_track_end->at(i_track).at(1),vars->TPCObj_PFP_track_end->at(i_track).at(2));
 double startenddist = (end-start).Mag();
 vars->mupiBDT.length_over_startend = (float)(vars->TPCObj_PFP_track_length->at(i_track)/startenddist);

 vars->mupiBDT.track_length = (float)vars->TPCObj_PFP_track_length->at(i_track);

 if (vars->TPCObj_LeadingMIPtrackIndex>=0){
   vars->mupiBDT.track_length_over_longest=(float)vars->TPCObj_PFP_track_length->at(i_track)/(float)vars->TPCObj_PFP_track_length->at(vars->TPCObj_LeadingMIPtrackIndex);
 }
 else{
   vars->mupiBDT.track_length_over_longest=-9999;
 }


 vars->mupiBDT.ndaughters = (float)vars->TPCObj_PFP_ndaughters->at(i_track);
 vars->mupiBDT.lnLmipovermu = (float)vars->TPCObj_PFP_lnLmipovermu->at(i_track);
 vars->mupiBDT.lnLmipoverpi = (float)vars->TPCObj_PFP_lnLmipoverpi->at(i_track);
 vars->mupiBDT.pandoraclassedastrack = (float)vars->TPCObj_PFP_PandoraClassedAsTrack_double->at(i_track);
 vars->mupiBDT.perc_used_hits = (float)vars->TPCObj_PFP_track_perc_used_hits->at(i_track);
 vars->mupiBDT.residuals_mean = (float)vars->TPCObj_PFP_track_residual_mean->at(i_track);
 vars->mupiBDT.residuals_stddev = (float)vars->TPCObj_PFP_track_residual_std->at(i_track);
 vars->mupiBDT.MCS_pi_maxScatter = (float)vars->TPCObj_PFP_track_MCS_pi_maxScatter->at(i_track);
 vars->mupiBDT.MCS_pi_meanScatter = (float)vars->TPCObj_PFP_track_MCS_pi_meanScatter->at(i_track);

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
 vars->mupiBDT.n_unused_hits_nearend = (float)n_unused_hits;
 vars->mupiBDT.unmatched_charge_nearend_plane2 = (float)unused_charge;
}

void InitialiseTree(TTree *tree, treevars *vars){
  tree -> Branch("length_over_startend", &(vars->mupiBDT.length_over_startend));
  tree -> Branch("track_length", &(vars->mupiBDT.track_length));
  tree -> Branch("track_length_over_longest", &(vars->mupiBDT.track_length_over_longest));
  tree -> Branch("ndaughters", &(vars->mupiBDT.ndaughters));
  tree -> Branch("lnLmipovermu", &(vars->mupiBDT.lnLmipovermu));
  tree -> Branch("lnLmipoverpi", &(vars->mupiBDT.lnLmipoverpi));
  tree -> Branch("pandoraclassedastrack", &(vars->mupiBDT.pandoraclassedastrack));
  tree -> Branch("perc_used_hits", &(vars->mupiBDT.perc_used_hits));
  tree -> Branch("residuals_mean", &(vars->mupiBDT.residuals_mean));
  tree -> Branch("residuals_stddev", &(vars->mupiBDT.residuals_stddev));
  tree -> Branch("MCS_pi_maxScatter", &(vars->mupiBDT.MCS_pi_maxScatter));
  tree -> Branch("MCS_pi_meanScatter", &(vars->mupiBDT.MCS_pi_meanScatter));
  tree -> Branch("n_unused_hits_nearend", &(vars->mupiBDT.n_unused_hits_nearend));
  tree -> Branch("unmatched_charge_nearend_plane2", &(vars->mupiBDT.unmatched_charge_nearend_plane2));
}
