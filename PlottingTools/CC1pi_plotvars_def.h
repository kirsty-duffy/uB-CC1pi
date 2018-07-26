#ifndef __CC1PI_PLOTVARS_DEF_H__
#define __CC1PI_PLOTVARS_DEF_H__

#include "CC1pi_treevars.h"

struct CC1piPlotVars{
  std::vector<double> const *Var;

  // For cut variables
  bool KeepBelowCut;
  bool OnlyDaughters;
  std::string TracksNeeded;
  double CutValue;

  // For plotting variables
  std::vector<double> bins;
  std::string histtitle;
  std::string histname;
};

// TPCObj_PFP_BrokenTrackAngle
CC1piPlotVars Var_TPCObj_PFP_BrokenTrackAngle(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_BrokenTrackAngle;
  tmp.KeepBelowCut = true;
  tmp.OnlyDaughters = false;
  tmp.TracksNeeded = "all";
  tmp.CutValue = 3.05;
  tmp.bins = {25,2.8,3.15};
  tmp.histtitle = ";Angle [rad];";
  tmp.histname = "BrokenTrackAngle";
  return tmp;
}

// TPCObj_PFP_track_perc_used_hits
CC1piPlotVars Var_TPCObj_PFP_track_perc_used_hits(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_perc_used_hits;
  tmp.KeepBelowCut = false;
  tmp.OnlyDaughters = true;
  tmp.TracksNeeded = "atleasttwo";
  tmp.CutValue = 0.7;
  tmp.bins = {25,0,1};
  tmp.histtitle = ";Fraction of used hits in cluster;";
  tmp.histname = "perc_used_hits";
  return tmp;
}

// TPCObj_PFP_VtxTrackDist
CC1piPlotVars Var_TPCObj_PFP_VtxTrackDist(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_VtxTrackDist;
  tmp.KeepBelowCut = true;
  tmp.OnlyDaughters = false;
  tmp.TracksNeeded = "atleasttwo";
  tmp.CutValue = 15.;
  tmp.bins = {25,0,20};
  tmp.histtitle = ";Distance from reconstructed vertex;";
  tmp.histname = "VtxTrackDist";
  return tmp;
}

// TPCObj_PFP_track_dEdx_truncmean_start
CC1piPlotVars Var_TPCObj_PFP_track_dEdx_truncmean_start(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_dEdx_truncmean_start;
  tmp.KeepBelowCut = true;
  tmp.OnlyDaughters = true;
  tmp.TracksNeeded = "exactlytwo";
  tmp.CutValue = 2.4;
  tmp.bins = {50,0,10};
  tmp.histtitle = ";Truncated Mean dE/dx at start of track;";
  tmp.histname = "dEdx_truncmean_atstart";
  return tmp;
}

// TPCObj_PFP_track_passesMIPcut
CC1piPlotVars Var_TPCObj_PFP_track_passesMIPcut(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_passesMIPcut;
  tmp.KeepBelowCut = false;
  tmp.OnlyDaughters = true;
  tmp.TracksNeeded = "exactlytwo";
  tmp.CutValue = 0.5;
  tmp.bins = {2,0,2};
  tmp.histtitle = ";IsMIP;";
  tmp.histname = "isMIP";
  return tmp;
}

// TPCObj_PFP_track_residual_mean_low
CC1piPlotVars Var_TPCObj_PFP_track_residual_mean_low(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_residual_mean;
  tmp.KeepBelowCut = false;
  tmp.OnlyDaughters = true;
  tmp.TracksNeeded = "atleasttwo";
  tmp.CutValue = -0.7;
  tmp.bins = {25,-2.8,0};
  tmp.histtitle = ";<r_{i}>;";
  tmp.histname = "residual_mean_down";
  return tmp;
}

// TPCObj_PFP_track_residual_mean_high
CC1piPlotVars Var_TPCObj_PFP_track_residual_mean_high(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_residual_mean;
  tmp.KeepBelowCut = true;
  tmp.OnlyDaughters = true;
  tmp.TracksNeeded = "atleasttwo";
  tmp.CutValue = 0.7;
  tmp.bins = {25,0,2.8};
  tmp.histtitle = ";<r_{i}>;";
  tmp.histname = "residual_mean_up";
  return tmp;
}

// TPCObj_PFP_track_residual_std_low
CC1piPlotVars Var_TPCObj_PFP_track_residual_std_low(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_residual_std;
  tmp.KeepBelowCut = false;
  tmp.OnlyDaughters = true;
  tmp.TracksNeeded = "atleasttwo";
  tmp.CutValue = 0.;
  tmp.bins = {25,0,2};
  tmp.histtitle = ";#sigma_{r_{i}};";
  tmp.histname = "residual_std_down";
  return tmp;
}

// TPCObj_PFP_track_residual_std_high
CC1piPlotVars Var_TPCObj_PFP_track_residual_std_high(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_residual_std;
  tmp.KeepBelowCut = true;
  tmp.OnlyDaughters = true;
  tmp.TracksNeeded = "atleasttwo";
  tmp.CutValue = 2.5;
  tmp.bins = {25,0,4};
  tmp.histtitle = ";#sigma_{r_{i}};";
  tmp.histname = "residual_std_up";
  return tmp;
}

// TPCObj_PFP_isContained_double
CC1piPlotVars Var_TPCObj_PFP_isContained(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_isContained_double;
  tmp.KeepBelowCut = false;
  tmp.OnlyDaughters = true;
  tmp.TracksNeeded = "all";
  tmp.CutValue = 0.5;
  tmp.bins = {2,0,2};
  tmp.histtitle = ";isContained;";
  tmp.histname = "isContained";
  return tmp;
}

// TPCObj_PFP_lnLmipoverp
CC1piPlotVars Var_TPCObj_PFP_lnLmipoverp(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_lnLmipoverp;
  tmp.KeepBelowCut = false;
  tmp.CutValue = 0.25;
  tmp.bins = {60,-10,10};
  tmp.histtitle = ";ln(L_{MIP})/(L_{p});";
  tmp.histname = "lnLmipoverp";
  return tmp;
}

// TPCObj_PFP_Lmumipovermumipp
CC1piPlotVars Var_TPCObj_PFP_Lmumipovermumipp(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_Lmumipovermumipp;
  tmp.KeepBelowCut = false;
  tmp.CutValue = 0.66;
  tmp.bins = {25,0.5,0.9};
  tmp.histtitle = ";(L_{#mu}+L_{MIP})/(L_{#mu}+L_{MIP}+L_{p});";
  tmp.histname = "Lmumipovermumipp";
  return tmp;
}

// TPCObj_DaughterTracks_Order_dEdxtr
CC1piPlotVars Var_TPCObj_DaughterTracks_Order_dEdxtr(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_DaughterTracks_Order_dEdxtr;
  tmp.bins = {10,0,10};
  tmp.histtitle = ";Daughter PFPs (ordered lowest to highest <dE/dx>_{tr});";
  tmp.histname = "daughterPFPs_bydEdxstart";
  return tmp;
}

// TPCObj_DaughterTracks_Order_dEdxtr_selMIPs
CC1piPlotVars Var_TPCObj_DaughterTracks_Order_dEdxtr_selMIPs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_DaughterTracks_Order_dEdxtr_selMIPs;
  tmp.bins = {10,0,10};
  tmp.histtitle = ";Daughter PFPs classed as MIPs (ordered lowest to highest <dE/dx>_{tr});";
  tmp.histname = "daughterPFPs_bydEdxstart_MIPs";
  return tmp;
}

// TPCObj_DaughterTracks_Order_trklen_selMIPs
CC1piPlotVars Var_TPCObj_DaughterTracks_Order_trklen_selMIPs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_DaughterTracks_Order_trklen_selMIPs;
  tmp.bins = {10,0,10};
  tmp.histtitle = ";Daughter PFPs classed as MIPs (ordered lowest to highest track length);";
  tmp.histname = "daughterPFPs_bytrklen_MIPs";
  return tmp;
}


#endif
