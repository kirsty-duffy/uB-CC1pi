#ifndef __CC1PI_CUTS_CXX__
#define __CC1PI_CUTS_CXX__

// #include <vector>

#include "CC1pi_treevars.h"
// #include "CC1pi_treevars.cxx"
#include "CC1pi_cuts.h"
#include "CC1pi_MIPcut.cxx"


bool IsEventSelected_SingleCut(double cutval, treevars *vars, int i_cut){

   bool Marco_selected = vars->Marco_selected;
   std::vector<bool> TPCObj_PFP_isDaughter = *(vars->TPCObj_PFP_isDaughter);

   CC1piPlotVars cutvar = GetCutVars(vars).at(i_cut);
   // cutvars is a vector corresponding to the TPCObject. It has one entry per track in the TPCObject.

   bool KeepBelowCut_i = cutvar.KeepBelowCut;
   bool OnlyDaughters_i = cutvar.OnlyDaughters;
   std::string TracksNeeded_i = cutvar.TracksNeeded;

   // std::cout << "i_cut = " << i_cut << ", cutval = " << cutval << ", KeepBelowCut = " << KeepBelowCut_i << ", OnlyDaughters = " << OnlyDaughters_i << ", TracksNeeded = " << TracksNeeded_i << std::endl;

   int tracks_in_event = cutvar.Var->size();
   int n_tracks = 0;
   int n_failed = 0;
   bool passes_LeadingMIP = false;
   bool passes_SecondMIP = false;
   for (int i_track=0; i_track<tracks_in_event; i_track++){
      double value = cutvar.Var->at(i_track);
      // std::cout << " --- Track " << i_track << " value " << value << std::endl;

      // Ignore PFPs that failed to reco as tracks. The neutrino PFP will also have a bogus value.
      if (value == -9999 || value == -999) {
         n_failed++;
         continue;
      }

      bool track_passes_cut=false;
      // (Optionally) only conisder direct daughters of the neutrino
      if (OnlyDaughters_i && !TPCObj_PFP_isDaughter.at(i_track)){
         track_passes_cut=false;
         continue;
      }

      if (KeepBelowCut_i && value < cutval){
         track_passes_cut=true;
         if (i_track == vars->TPCObj_LeadingMIPtrackIndex) passes_LeadingMIP = true;
         if (i_track == vars->TPCObj_SecondMIPtrackIndex) passes_SecondMIP = true;
      }
      else if (!KeepBelowCut_i && value > cutval){
         track_passes_cut=true;
         if (i_track == vars->TPCObj_LeadingMIPtrackIndex) passes_LeadingMIP = true;
         if (i_track == vars->TPCObj_SecondMIPtrackIndex) passes_SecondMIP = true;
      }

      if (track_passes_cut){
         n_tracks++;
      }

   } // end loop over tracks in TPCObject

   // Throw away events that failed the CC inclusive selection or that have PFPs that failed to reco as tracks.
   // Then check if the correct number of tracks pass the cut, given the chosen options
   bool isSelected = false;
   if(n_failed < 2 && Marco_selected) { // The neutrino PFP will look "failed", so we expect every event to have 1
      if(TracksNeeded_i == "atleasttwo" && n_tracks >= 2) isSelected = true;

      else if (TracksNeeded_i == "exactlytwo" && n_tracks == 2) isSelected = true;

      else if (TracksNeeded_i == "exactlyone" && n_tracks == 1) isSelected = true;

      else if (TracksNeeded_i == "all") {
         if (OnlyDaughters_i && n_tracks == std::count(TPCObj_PFP_isDaughter.begin(),TPCObj_PFP_isDaughter.end(),true)) isSelected = true;
         else if(!OnlyDaughters_i && n_tracks == tracks_in_event-1) isSelected = true; // Again, 1 less because of the neutrino
      }

      else if (TracksNeeded_i == "NA" && n_tracks == tracks_in_event) isSelected = true;

      else if (TracksNeeded_i == "LeadingMIP" && passes_LeadingMIP) isSelected = true;

      else if (TracksNeeded_i == "SecondMIP" && passes_SecondMIP) isSelected = true;

      else if (TracksNeeded_i == "BothMIPs" && passes_LeadingMIP && passes_SecondMIP) isSelected = true;

      else if (TracksNeeded_i == "allExceptLeadingMIP"){
        if (OnlyDaughters_i && passes_LeadingMIP && n_tracks == std::count(TPCObj_PFP_isDaughter.begin(),TPCObj_PFP_isDaughter.end(),true)) isSelected = true;
        else if (OnlyDaughters_i && !passes_LeadingMIP && n_tracks == std::count(TPCObj_PFP_isDaughter.begin(),TPCObj_PFP_isDaughter.end(),true)-1) isSelected = true;
        else if (!OnlyDaughters_i && passes_LeadingMIP && n_tracks == tracks_in_event-1) isSelected = true;
        else if (!OnlyDaughters_i && !passes_LeadingMIP && n_tracks == tracks_in_event-2) isSelected = true;
      }
   }
   // std::cout << "isSelected_i = " << isSelected << std::endl;
   return isSelected;
}


// --------------------------------------------

bool IsEventSelected(treevars *vars){

   bool isSelected = true;
   std::vector<CC1piPlotVars> Cuts = GetCutVars(vars);
   for (size_t i_cut = 0; i_cut < Cuts.size(); i_cut++){
      if(!IsEventSelected_SingleCut(Cuts.at(i_cut).CutValue, vars, i_cut)){
         isSelected = false;
         break;
      }
   }

   return isSelected;
}

#endif
