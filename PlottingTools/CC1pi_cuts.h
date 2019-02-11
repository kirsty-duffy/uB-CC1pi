#ifndef __CC1PI_CUTS_H__
#define __CC1PI_CUTS_H__


// This code defines the cuts used in the CC1pi analysis. These cuts will automatically be applied in MakeCutPlots and MakeCutEfficiencyPlots.

// #include <vector>
#include "CC1pi_treevars.h"
#include "CC1pi_plotvars_def.h"

// --------------------------------------------------- //
// First: define cuts that we want to apply. Anything uncommented will be applied, and you need to make sure the vectors are in the same order! Note that we do *NOT* want to include individual cuts that define a MIP here (e.g. <dEdx> or Lmipoverp) - those will be defined in the next section and applied in EvalMIPCut. In this section, they are all just combined into one variable, vars->TPCObj_PFP_track_passesMIPcut.

// What variables do we want to cut on?
std::vector<CC1piPlotVars> GetCutVars(treevars *vars) {
   std::vector<CC1piPlotVars> cut_vars = {
   //Var_TPCObj_AngleBetweenMIPs(vars)
   Var_TPCObj_AngleBetweenMIPs_high(vars)
   // ,Var_TPCObj_LeadingMIP_PandoraClassedAsTrack(vars)
   ,Var_TPCObj_PFP_track_passesMIPcut(vars)
   ,Var_MIP_containment(vars)
   // ,Var_TPCObj_PFP_track_passesPioncut(vars)

   // ,Var_TPCObj_SecondMIP_isContained(vars)
   // ,Var_TPCObj_FirstMIP_isContained(vars)
   // ,Var_TPCObj_AllDaughters_isContained(vars)
   // ,Var_TPCObj_AllTracks_isContained(vars)
   // ,Var_TPCObj_AllDaughtersExceptLeadingMIP_isContained(vars)
   // ,Var_TPCObj_PFP_track_nhits_LeadingMIP(vars) // Doesn't matter that it's got LeadingMIP in name, applies to all tracks

   //,Var_TPCObj_dEdx_truncmean_MIPdiff(vars)
   };
   return cut_vars;
};

// --------------------------------------------------- //
// Now define the vectors that go into the MIP cut. In CC1pi_MIPcut.cxx we define EvalMIPCut - that function will loop over all of these variables and evaluate whether a track passes a cut in them. Only if a track passes all cuts will it be defined as a MIP.

// What variables do we want to cut on?
std::vector<CC1piPlotVars> GetMIPCutVars(treevars *vars) {
   std::vector<CC1piPlotVars> cut_vars = {
   //Var_TPCObj_PFP_track_dEdx_truncmean_start(vars)
   //Var_TPCObj_PFP_track_dEdx_truncmean_start_lowcut(vars)
   Var_TPCObj_PFP_VtxTrackDist(vars)
   //,Var_TPCObj_PFP_track_dedx_grminhits(vars)
   //,Var_TPCObj_PFP_lnLmipoverp(vars)
   ,Var_TPCObj_PFP_track_BDTscore(vars)
   ,Var_TPCObj_PFP_track_theta_parallel(vars)
   };
   return cut_vars;
};

// --------------------------------------------------- //
// Now define the vectors that go into the pion cut. In CC1pi_Pioncut.cxx we define EvalPionCut - that function will loop over all of these variables and evaluate whether a track passes a cut in them. Only if a track passes all cuts will it be defined as a pion.

// What variables do we want to cut on?
// std::vector<CC1piPlotVars> GetPionCutVars(treevars *vars) {
//    std::vector<CC1piPlotVars> cut_vars = {
//    Var_TPCObj_PFP_track_passesMIPcut(vars)
//    ,Var_TPCObj_PFP_isDaughter(vars)
//    ,
//    };
//    return cut_vars;
// };


#endif
