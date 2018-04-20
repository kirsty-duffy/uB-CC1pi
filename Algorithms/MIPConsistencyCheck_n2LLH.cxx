// MIP Consistency Check copied from Marco's code on 25th Jan 2018
// This is not intended to be the final MIP consistency cut, just a placeholder
// until we can put in something more concrete.

// It uses a hardcoded array for the cut in length-dqdx space. I have copied
// this from Marco's github (fcl file UBXSec/job/muoncandidatefinder.fcl) but
// put the hardcoded array directly in the header file

// -- Kirsty Duffy, Fermilab, 25/01/2018

// Check if we've already defined a MIP consistency function (don't redefine!)
#ifndef MIPCONSISTENCY_CXX
#define MIPCONSISTENCY_CXX

// Don't need any include statements other than to include MIPConsistencyCheck_PIDA.h
// All the other necessary includes should be in there
#include "MIPConsistencyCheck_n2LLH.h"

// ------------------------------------------------------------------------------- //

bool IsMIP(art::FindManyP<anab::ParticleID> trackPIDAssn, int TrackID)
{
  // Return false for tracks that we have no PID information about
  if (!trackPIDAssn.isValid()){
    return false;
  }

   std::vector<art::Ptr<anab::ParticleID>> trackPID = trackPIDAssn.at(TrackID);

   if (trackPID.size() == 0){
     return false;
   }

   // Only fill variables if the track-PID association exists
   std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();

   //double Bragg_fwd_mu = 99999;
   double Bragg_fwd_p  = 99999;
   //double Bragg_bwd_mu = 99999;
   double Bragg_bwd_p  = 99999;
   //double noBragg_MIP  = 99999;

   // Loop through AlgScoresVec and find the variables we want
   for (size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){

     anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);

     if (AlgScore.fAlgName == "BraggPeakLLH"){
       if (anab::kVariableType(AlgScore.fVariableType) == anab::kLogL_fwd){
         //if (AlgScore.fAssumedPdg == 13)   Bragg_fwd_mu = AlgScore.fValue;
         if (AlgScore.fAssumedPdg == 2212) Bragg_fwd_p =  AlgScore.fValue;
         //if (AlgScore.fAssumedPdg == 0)    noBragg_MIP = AlgScore.fValue;
       }// if fVariableType == anab::kLogL_fwd
       else if (anab::kVariableType(AlgScore.fVariableType) == anab::kLogL_bwd){
         //if (AlgScore.fAssumedPdg == 13)   Bragg_bwd_mu = AlgScore.fValue;
         if (AlgScore.fAssumedPdg == 2212) Bragg_bwd_p =  AlgScore.fValue;
       } // if fVariableType == anab::kLogL_bwd
     } // if fAlName = BraggPeakLLH
   } // end loop through AlgScoresVec


  // Now evaluate whether it passes the cut!
  //double Bragg_mu = std::min({Bragg_fwd_mu, Bragg_bwd_mu, noBragg_MIP});
  double Bragg_p = std::min(Bragg_fwd_p,Bragg_bwd_p);

  if (Bragg_p > 100){ return true;}
  else{ return false;}

}

#endif
