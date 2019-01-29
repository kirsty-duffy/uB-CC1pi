#include "CC1pi_treevars.h"
#include "TMVA/Reader.h"

 bool mupiBDT::CheckMVAvars(){
   bool returnval = true;

   if (length_over_startend<0 || isnan(length_over_startend)) returnval=false;
   if (track_length==-9999 || isnan(track_length)) returnval=false;
   if (ndaughters==-9999 || isnan(ndaughters)) returnval=false;
   if (lnLmipovermu==-9999 || isnan(lnLmipovermu)) returnval=false;
   if (lnLmipoverpi==-9999 || isnan(lnLmipoverpi)) returnval=false;
   //if (pandoraclassedastrack==-9999 || isnan(pandoraclassedastrack)) returnval=false;
   if (perc_used_hits==-9999 || isnan(perc_used_hits)) returnval=false;
   if (residuals_mean==-9999 || isnan(residuals_mean)) returnval=false;
   if (residuals_stddev==-9999 || isnan(residuals_stddev)) returnval=false;
   if (MCS_pi_maxScatter==-9999 || isnan(MCS_pi_maxScatter)) returnval=false;
   if (MCS_pi_meanScatter==-9999 || isnan(MCS_pi_meanScatter)) returnval=false;

   return returnval;
}

void mupiBDT::initialise_BDT(){
   fReader = new TMVA::Reader("");
   fReader->AddVariable("length_over_startend",&length_over_startend);
   fReader->AddVariable("track_length",&track_length);
   fReader->AddVariable("ndaughters",&ndaughters);
   fReader->AddVariable("lnLmipovermu",&lnLmipovermu);
   fReader->AddVariable("lnLmipoverpi",&lnLmipoverpi);
   //fReader->AddVariable("pandoraclassedastrack",&pandoraclassedastrack);
   fReader->AddVariable("perc_used_hits",&perc_used_hits);
   fReader->AddVariable("residuals_mean",&residuals_mean);
   fReader->AddVariable("residuals_stddev",&residuals_stddev);
   fReader->AddVariable("MCS_pi_maxScatter",&MCS_pi_maxScatter);
   fReader->AddVariable("MCS_pi_meanScatter",&MCS_pi_meanScatter);
   fReader->BookMVA(BookMVAType.c_str(), BookMVALoc.c_str());
}



void SetMVAvars_mupiBDT(treevars *vars, int i_track){
 TVector3 start(vars->TPCObj_PFP_track_start->at(i_track).at(0),vars->TPCObj_PFP_track_start->at(i_track).at(1),vars->TPCObj_PFP_track_start->at(i_track).at(2));
 TVector3 end(vars->TPCObj_PFP_track_end->at(i_track).at(0),vars->TPCObj_PFP_track_end->at(i_track).at(1),vars->TPCObj_PFP_track_end->at(i_track).at(2));
 double startenddist = (end-start).Mag();
 vars->mupiBDT.length_over_startend = (float)(vars->TPCObj_PFP_track_length->at(i_track)/startenddist);

 vars->mupiBDT.track_length = (float)vars->TPCObj_PFP_track_length->at(i_track);
 vars->mupiBDT.ndaughters = (float)vars->TPCObj_PFP_ndaughters->at(i_track);
 vars->mupiBDT.lnLmipovermu = (float)vars->TPCObj_PFP_lnLmipovermu->at(i_track);
 vars->mupiBDT.lnLmipoverpi = (float)vars->TPCObj_PFP_lnLmipoverpi->at(i_track);
 // vars->mupiBDT.pandoraclassedastrack = (float)vars->TPCObj_PFP_PandoraClassedAsTrack_double->at(i_track);
 vars->mupiBDT.perc_used_hits = (float)vars->TPCObj_PFP_track_perc_used_hits->at(i_track);
 vars->mupiBDT.residuals_mean = (float)vars->TPCObj_PFP_track_residual_mean->at(i_track);
 vars->mupiBDT.residuals_stddev = (float)vars->TPCObj_PFP_track_residual_std->at(i_track);
 vars->mupiBDT.MCS_pi_maxScatter = (float)vars->TPCObj_PFP_track_MCS_pi_maxScatter->at(i_track);
 vars->mupiBDT.MCS_pi_meanScatter = (float)vars->TPCObj_PFP_track_MCS_pi_meanScatter->at(i_track);
}
