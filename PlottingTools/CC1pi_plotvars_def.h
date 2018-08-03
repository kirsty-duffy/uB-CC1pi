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
  bool PlotOnlyDaughterMIPs;
  bool PlotOnlyContained;
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
  tmp.PlotOnlyDaughterMIPs = false;
  tmp.PlotOnlyContained = false;
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
  tmp.PlotOnlyDaughterMIPs = false;
  tmp.PlotOnlyContained = false;
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
  tmp.PlotOnlyDaughterMIPs = false;
  tmp.PlotOnlyContained = false;
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
  tmp.PlotOnlyDaughterMIPs = false;
  tmp.PlotOnlyContained = false;
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
  tmp.PlotOnlyDaughterMIPs = false;
  tmp.PlotOnlyContained = false;
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
  tmp.PlotOnlyDaughterMIPs = false;
  tmp.PlotOnlyContained = false;
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
  tmp.PlotOnlyDaughterMIPs = false;
  tmp.PlotOnlyContained = false;
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
  tmp.PlotOnlyDaughterMIPs = false;
  tmp.PlotOnlyContained = false;
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
  tmp.PlotOnlyDaughterMIPs = false;
  tmp.PlotOnlyContained = false;
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
  tmp.PlotOnlyDaughterMIPs = false;
  tmp.PlotOnlyContained = false;
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
  tmp.PlotOnlyDaughterMIPs = false;
  tmp.PlotOnlyContained = false;
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
  tmp.PlotOnlyDaughterMIPs = false;
  tmp.PlotOnlyContained = false;
  return tmp;
}

// TPCObj_DaughterTracks_Order_dEdxtr
CC1piPlotVars Var_TPCObj_DaughterTracks_Order_dEdxtr(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_DaughterTracks_Order_dEdxtr;
  tmp.bins = {10,0,10};
  tmp.histtitle = ";Daughter PFPs (ordered lowest to highest <dE/dx>_{tr});";
  tmp.histname = "daughterPFPs_bydEdxstart";
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = false;
  return tmp;
}

// TPCObj_DaughterTracks_Order_dEdxtr_selMIPs
CC1piPlotVars Var_TPCObj_DaughterTracks_Order_dEdxtr_selMIPs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_DaughterTracks_Order_dEdxtr_selMIPs;
  tmp.bins = {10,0,10};
  tmp.histtitle = ";Daughter PFPs classed as MIPs (ordered lowest to highest <dE/dx>_{tr});";
  tmp.histname = "daughterPFPs_bydEdxstart_MIPs";
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = false;
  return tmp;
}

// TPCObj_DaughterTracks_Order_trklen_selMIPs
CC1piPlotVars Var_TPCObj_DaughterTracks_Order_trklen_selMIPs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_DaughterTracks_Order_trklen_selMIPs;
  tmp.bins = {10,0,10};
  tmp.histtitle = ";Daughter PFPs classed as MIPs (ordered lowest to highest track length);";
  tmp.histname = "daughterPFPs_bytrklen_MIPs";
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = false;
  return tmp;
}

// TPCObj_PFP_MCSLLmuMinusLLp: daughter MIPs only
CC1piPlotVars Var_TPCObj_PFP_MCSLLmuMinusLLp_DaughterMIPs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MCSLLmuMinusLLp;
  tmp.bins = {60,-50,10};
  tmp.histtitle = ";MCS LL_{#mu} - LL_{p} (Daughter PFPs classed as MIPs only);";
  tmp.histname = "MCSLLmuMinusLLp_daughterMIPs";
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = false;
  return tmp;
}

// TPCObj_PFP_track_MCSLLmuMinusLLp: contained daughter MIPs only
CC1piPlotVars Var_TPCObj_PFP_track_MCSLLmuMinusLLp_ContDaughterMIPs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MCSLLmuMinusLLp;
  tmp.bins = {60,-50,10};
  tmp.histtitle = ";MCS LL_{#mu} - LL_{p} (Contained daughter PFPs classed as MIPs only);";
  tmp.histname = "MCSLLmuMinusLLp_contdaughterMIPs";
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_track_MomRangeMinusMCS_p: contained daughter MIPs only
CC1piPlotVars Var_TPCObj_PFP_track_MomRangeMinusMCS_p_ContDaughterMIPs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MomRangeMinusMCS_p;
  tmp.bins = {100,-5,5};
  tmp.histtitle = ";Mom. by range - Mom. by MCS, p assumption, GeV (Cont. daughter MIP-like PFPs only);";
  tmp.histname = "MomRangeMinusMCS_p_contdaughterMIPs";
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_track_MomRangeMinusMCS_mu: contained daughter MIPs only
CC1piPlotVars Var_TPCObj_PFP_track_MomRangeMinusMCS_mu_ContDaughterMIPs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MomRangeMinusMCS_mu;
  tmp.bins = {100,-5,5};
  tmp.histtitle = ";Mom. by range - Mom. by MCS, #mu assumption, GeV (Cont. daughter MIP-like PFPs only);";
  tmp.histname = "MomRangeMinusMCS_mu_contdaughterMIPs";
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_track_SimpleCluster_hitLinearity_minimum
CC1piPlotVars Var_TPCObj_PFP_track_SimpleCluster_hitLinearity_minimum(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_SimpleCluster_hitLinearity_minimum;
  tmp.bins = {10,0,1};
  tmp.histtitle = ";Local Linearity Minimum;";
  tmp.histname = "hitLinearity_minimum";
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_track_SimpleCluster_hitLinearity_mean
CC1piPlotVars Var_TPCObj_PFP_track_SimpleCluster_hitLinearity_mean(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_SimpleCluster_hitLinearity_mean;
  tmp.bins = {10,0,1};
  tmp.histtitle = ";Local Linearity Mean;";
  tmp.histname = "hitLinearity_mean";
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_track_SimpleCluster_hitLinearity_truncated_mean
CC1piPlotVars Var_TPCObj_PFP_track_SimpleCluster_hitLinearity_truncated_mean(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_SimpleCluster_hitLinearity_truncated_mean;
  tmp.bins = {10,0,1};
  tmp.histtitle = ";Local Linearity Truncated Mean;";
  tmp.histname = "hitLinearity_truncated_mean";
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

#endif
