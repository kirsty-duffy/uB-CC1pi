#ifndef __CC1PI_PLOTVARS_DEF_H__
#define __CC1PI_PLOTVARS_DEF_H__

#include "CC1pi_treevars.h"

struct CC1piPlotVars{
  std::vector<double> const *Var;

  // For cut variables
  bool isMIPcut;
  bool KeepBelowCut;
  bool OnlyDaughters;
  std::string TracksNeeded;
  double CutValue;

  // For plotting variables
  std::vector<double> bins;
  std::string histtitle;
  std::string histname;
  bool PlotOnlyDaughters;
  bool PlotOnlyDaughterMIPs;
  bool PlotOnlyNotMIPDaughters;
  bool PlotOnlyContained;
  bool PlotOnlyLeadingDaughterMIP;
  bool PlotOnlySecondDaughterMIP;
  bool PlotOnlyPionCandidate;
  bool PlotNotPionCandidate;
  bool PlotOnlyMuMuPairs;

  // Default initialisation of plotting bools only
  CC1piPlotVars(){
    isMIPcut = false;
    KeepBelowCut = false;
    OnlyDaughters = false;
    TracksNeeded = "NA";
    CutValue = -9999;

    PlotOnlyDaughters = false;
    PlotOnlyDaughterMIPs = false;
    PlotOnlyNotMIPDaughters = false;
    PlotOnlyContained = false;
    PlotOnlyLeadingDaughterMIP = false;
    PlotOnlySecondDaughterMIP = false;
    PlotOnlyPionCandidate = false;
    PlotNotPionCandidate = false;
    PlotOnlyMuMuPairs = false;
  }
};

// Reconstructed neutrino daughters
CC1piPlotVars Var_TPCObj_PFP_isDaughter(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_isDaughter_double;
  tmp.KeepBelowCut = false;
  tmp.CutValue = 0.5;
  tmp.histtitle = ";Is Neutrino Daughter?;";
  tmp.histname = "IsDaughter";
  return tmp;
}

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

// TPCObj_PFP_BrokenTrackAngle
CC1piPlotVars Var_TPCObj_AngleBetweenMIPs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_AngleBetweenMIPs;
  tmp.KeepBelowCut = false;
  tmp.TracksNeeded = "NA";
  tmp.CutValue = 0.1;
  tmp.bins = {25,0,3.14};
  tmp.histtitle = ";Angle between MIP candidates [rad];";
  tmp.histname = "AngleBetweenMIPs";
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = false;
  return tmp;
}

// TPCObj_PFP_BrokenTrackAngle
CC1piPlotVars Var_TPCObj_AngleBetweenMIPs_mumupairs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_AngleBetweenMIPs;
  tmp.bins = {32,0,3.2};
  tmp.histtitle = ";Angle between MIP candidates (#mu^{-}#mu^{-} pairs only) [rad];";
  tmp.histname = "AngleBetweenMIPs_mumupairs";
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyMuMuPairs = true;
  return tmp;
}

// TPCObj_PFP_BrokenTrackAngle
CC1piPlotVars Var_TPCObj_AngleBetweenMIPs_high(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_AngleBetweenMIPs;
  tmp.KeepBelowCut = true;
  tmp.TracksNeeded = "NA";
  tmp.CutValue = 2.6;
  tmp.bins = {12,2,3.2};
  tmp.histtitle = ";Angle between MIP candidates [rad];";
  tmp.histname = "AngleBetweenMIPs_high";
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = false;
  return tmp;
}

// TPCObj_PFP_track_perc_used_hits
CC1piPlotVars Var_TPCObj_PFP_track_perc_used_hits(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_perc_used_hits;
  tmp.KeepBelowCut = false;
  tmp.isMIPcut = true;
  tmp.CutValue = 0.7;
  tmp.bins = {5,0.5,1};
  tmp.histtitle = ";Fraction of used hits in cluster;";
  tmp.histname = "perc_used_hits";
  return tmp;
}

// TPCObj_PFP_track_perc_used_hits
CC1piPlotVars Var_TPCObj_PFP_track_perc_used_hits_LeadingMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_perc_used_hits;
  tmp.KeepBelowCut = false;
  tmp.CutValue = 0.7;
  tmp.bins = {25,0,1};
  tmp.histtitle = ";Fraction of used hits in cluster (Leading MIP);";
  tmp.histname = "perc_used_hits_LeadingMIP";
  tmp.PlotOnlyLeadingDaughterMIP = true;
  tmp.PlotOnlyContained = false;
  return tmp;
}

// TPCObj_PFP_track_perc_used_hits
CC1piPlotVars Var_TPCObj_PFP_track_perc_used_hits_SecondMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_perc_used_hits;
  tmp.KeepBelowCut = false;
  tmp.OnlyDaughters = true;
  tmp.TracksNeeded = "atleasttwo";
  tmp.CutValue = 0.7;
  tmp.bins = {25,0,1};
  tmp.histtitle = ";Fraction of used hits in cluster (Second MIP);";
  tmp.histname = "perc_used_hits_SecondMIP";
  tmp.PlotOnlySecondDaughterMIP = true;
  tmp.PlotOnlyContained = false;
  return tmp;
}

// TPCObj_PFP_VtxTrackDist
CC1piPlotVars Var_TPCObj_PFP_VtxTrackDist(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_VtxTrackDist;
  tmp.KeepBelowCut = true;
  tmp.isMIPcut = true;
  tmp.CutValue = 5.;
  tmp.bins = {20,0,20};
  tmp.histtitle = ";Track start distance from reconstructed vertex [cm];";
  tmp.histname = "VtxTrackDist";
  tmp.PlotOnlyDaughters = true;
  return tmp;
}

// TPCObj_PFP_VtxTrackDist
CC1piPlotVars Var_TPCObj_PFP_VtxTrackDist_LeadingMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_VtxTrackDist;
  tmp.KeepBelowCut = true;
  tmp.CutValue = 15.;
  tmp.bins = {25,0,20};
  tmp.histtitle = ";Track start distance from reconstructed vertex [cm] (Leading MIP);";
  tmp.histname = "VtxTrackDist_LeadingMIP";
  tmp.PlotOnlyLeadingDaughterMIP = true;
  tmp.PlotOnlyContained = false;
  return tmp;
}

// TPCObj_PFP_VtxTrackDist
CC1piPlotVars Var_TPCObj_PFP_VtxTrackDist_SecondMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_VtxTrackDist;
  tmp.KeepBelowCut = true;
  tmp.OnlyDaughters = false;
  tmp.TracksNeeded = "atleasttwo";
  tmp.CutValue = 15.;
  tmp.bins = {25,0,20};
  tmp.histtitle = ";Track start distance from reconstructed vertex [cm] (Second MIP);";
  tmp.histname = "VtxTrackDist_SecondMIP";
  tmp.PlotOnlySecondDaughterMIP = true;
  tmp.PlotOnlyContained = false;
  return tmp;
}

// TPCObj_PFP_VtxTrackDist
CC1piPlotVars Var_TPCObj_PFP_VtxTrackDist_mumupairs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_VtxTrackDist;
  tmp.KeepBelowCut = true;
  tmp.OnlyDaughters = false;
  tmp.TracksNeeded = "atleasttwo";
  tmp.CutValue = 15.;
  tmp.bins = {25,0,20};
  tmp.histtitle = ";Track start distance from reconstructed vertex [cm] (Leading MIP, #mu^{-}#mu^{-} pairs only);";
  tmp.histname = "VtxTrackDist_LeadingMIP_mumupairs";
  tmp.PlotOnlyLeadingDaughterMIP = true;
  tmp.PlotOnlyContained = false;
  tmp.PlotOnlyMuMuPairs = true;
  return tmp;
}

// TPCObj_PFP_VtxTrackDist
CC1piPlotVars Var_TPCObj_PFP_VtxTrackDist_SecondMIP_mumupairs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_VtxTrackDist;
  tmp.KeepBelowCut = true;
  tmp.OnlyDaughters = false;
  tmp.TracksNeeded = "atleasttwo";
  tmp.CutValue = 15.;
  tmp.bins = {25,0,20};
  tmp.histtitle = ";Track start distance from reconstructed vertex [cm] (Second MIP, #mu^{-}#mu^{-} pairs only);";
  tmp.histname = "VtxTrackDist_SecondMIP_mumupairs";
  tmp.PlotOnlySecondDaughterMIP = true;
  tmp.PlotOnlyContained = false;
  tmp.PlotOnlyMuMuPairs = true;
  return tmp;
}

// TPCObj_PFP_VtxTrackDist
CC1piPlotVars Var_TPCObj_PFP_track_Chi2Proton_SecondMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_Chi2Proton_plane2;
  tmp.bins = {25,0,200};
  tmp.histtitle = ";#chi^{2}_{proton} (Second MIP);";
  tmp.histname = "chi2p_SecondMIP";
  tmp.PlotOnlySecondDaughterMIP = true;
  tmp.PlotOnlyContained = false;
  return tmp;
}

// TPCObj_PFP_VtxTrackDist
CC1piPlotVars Var_TPCObj_PFP_track_Chi2Proton_SecondMIPcont(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_Chi2Proton_plane2;
  tmp.bins = {25,0,200};
  tmp.histtitle = ";#chi^{2}_{proton} (Second MIP, cont.);";
  tmp.histname = "chi2p_SecondMIPcont";
  tmp.PlotOnlySecondDaughterMIP = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_track_dEdx_truncmean_start
CC1piPlotVars Var_TPCObj_PFP_track_dEdx_truncmean_start(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_dEdx_truncmean_start;
  tmp.KeepBelowCut = true;
  tmp.isMIPcut = true;
  tmp.CutValue = 2.2;
  tmp.bins = {25,0,2.2};
  tmp.histtitle = ";Truncated Mean dE/dx at start of track [MeV/cm];";
  tmp.histname = "dEdx_truncmean_atstart";
  tmp.PlotOnlyDaughterMIPs = false;
  tmp.PlotOnlyContained = false;
  tmp.PlotOnlyDaughters = true;
  return tmp;
}

// TPCObj_PFP_track_dEdx_truncmean_start
CC1piPlotVars Var_TPCObj_PFP_track_dEdx_truncmean_start_lowcut(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_dEdx_truncmean_start;
  tmp.KeepBelowCut = false;
  tmp.isMIPcut = true;
  tmp.CutValue = 1.0;
  tmp.bins = {10,0.5,1.5};
  tmp.histtitle = ";Truncated Mean dE/dx at start of track;";
  tmp.histname = "dEdx_truncmean_atstart_lowcut";
  tmp.PlotOnlyDaughters = true;
  return tmp;
}

// TPCObj_PFP_track_dEdx_truncmean_start
CC1piPlotVars Var_TPCObj_PFP_track_dEdx_stddev_SecondMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_dedx_stddev;
  tmp.bins = {50,0,20};
  tmp.histtitle = ";Standard deviation in dE/dx over track (Second MIP only);";
  tmp.histname = "dEdx_stddev_SecondMIP";
  tmp.PlotOnlySecondDaughterMIP = true;
  return tmp;
}

// TPCObj_PFP_track_dEdx_truncmean_start
CC1piPlotVars Var_TPCObj_PFP_track_dEdx_stddev_start_SecondMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_dEdx_stddev_start;
  tmp.bins = {50,0,50};
  tmp.histtitle = ";Standard deviation in dE/dx at start of track (Second MIP only);";
  tmp.histname = "dEdx_stddev_start_SecondMIP";
  tmp.PlotOnlySecondDaughterMIP = true;
  return tmp;
}

// TPCObj_PFP_track_dEdx_truncmean_start
CC1piPlotVars Var_TPCObj_PFP_track_dEdx_mean_start_SecondMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_dEdx_mean_start;
  tmp.bins = {50,0,10};
  tmp.histtitle = ";Mean dE/dx at start of track (Second MIP only);";
  tmp.histname = "dEdx_mean_start_SecondMIP";
  tmp.PlotOnlySecondDaughterMIP = true;
  return tmp;
}

// TPCObj_PFP_track_dEdx_truncmean_start
CC1piPlotVars Var_TPCObj_PFP_track_dEdx_truncmoverm_start_SecondMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_dEdx_truncmoverm_start;
  tmp.bins = {50,0,1.5};
  tmp.histtitle = ";<dE/dx>_{tr}/<dE/dx> at start of track (Second MIP only);";
  tmp.histname = "dEdx_truncmoverm_start_SecondMIP";
  tmp.PlotOnlySecondDaughterMIP = true;
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
  tmp.isMIPcut = true;
  tmp.CutValue = -0.5;
  tmp.bins = {60,-10,10};
  tmp.histtitle = ";ln(L_{MIP})/(L_{p});";
  tmp.histname = "lnLmipoverp";
  tmp.PlotOnlyDaughterMIPs = false;
  tmp.PlotOnlyContained = false;
  tmp.PlotOnlyDaughters = true;
  return tmp;
}

// TPCObj_PFP_lnLmipoverp
CC1piPlotVars Var_TPCObj_PFP_lnLmipoverp_cont(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_lnLmipoverp;
  tmp.KeepBelowCut = false;
  tmp.CutValue = 0.25;
  tmp.bins = {60,-10,10};
  tmp.histtitle = ";ln(L_{MIP})/(L_{p}) (Cont. MIPs only);";
  tmp.histname = "lnLmipoverp_cont";
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_lnLmipoverp
CC1piPlotVars Var_TPCObj_PFP_lnLmipoverp_SecondMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_lnLmipoverp;
  tmp.bins = {60,-10,10};
  tmp.histtitle = ";ln(L_{MIP})/(L_{p}) (Second MIP only);";
  tmp.histname = "lnLmipoverp_secondMIPonly";
  tmp.PlotOnlySecondDaughterMIP = true;
  tmp.PlotOnlyContained = false;
  return tmp;
}

// TPCObj_PFP_lnLmipoverp
CC1piPlotVars Var_TPCObj_PFP_lnLmipoverp_SecondMIPcont(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_lnLmipoverp;
  tmp.bins = {60,-10,10};
  tmp.histtitle = ";ln(L_{MIP})/(L_{p}) (Second MIP, cont. only);";
  tmp.histname = "lnLmipoverp_secondMIPonly_cont";
  tmp.PlotOnlySecondDaughterMIP = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_Lmumipovermumipp
CC1piPlotVars Var_TPCObj_PFP_Lmumipovermumipp(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_Lmumipovermumipp;
  tmp.KeepBelowCut = false;
  tmp.CutValue = 0.66;
  tmp.bins = {50,0.,1.};
  tmp.histtitle = ";(L_{#mu}+L_{MIP})/(L_{#mu}+L_{MIP}+L_{p}) (MIPs only);";
  tmp.histname = "Lmumipovermumipp";
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = false;
  return tmp;
}

// TPCObj_PFP_Lmumipovermumipp
CC1piPlotVars Var_TPCObj_PFP_Lmumipovermumipp_cont(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_Lmumipovermumipp;
  tmp.KeepBelowCut = false;
  tmp.CutValue = 0.66;
  tmp.bins = {50,0.,1.};
  tmp.histtitle = ";(L_{#mu}+L_{MIP})/(L_{#mu}+L_{MIP}+L_{p}) (Cont. MIPs only);";
  tmp.histname = "Lmumipovermumipp_cont";
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = true;
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

// TPCObj_PFP_MCSLLmuMinusLLp: leading MIP only
CC1piPlotVars Var_TPCObj_PFP_MCSLLpiMinusLLp_LeadingMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MCSLLpiMinusLLp;
  tmp.bins = {88,-2,0.2};
  tmp.histtitle = ";(MCS LL_{#pi} - LL_{p})/MCS LL_{#pi} (Leading MIP only);";
  tmp.histname = "MCSLLmuMinusLLp_leadingMIP";
  tmp.PlotOnlyLeadingDaughterMIP = true;
  tmp.PlotOnlyContained = false;
  return tmp;
}

// TPCObj_PFP_MCSLLmuMinusLLp: second MIP only
CC1piPlotVars Var_TPCObj_PFP_MCSLLpiMinusLLp_SecondMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MCSLLpiMinusLLp;
  tmp.bins = {88,-2,0.2};
  tmp.histtitle = ";(MCS LL_{#pi} - LL_{p})/MCS LL_{#pi} (Second MIP only);";
  tmp.histname = "MCSLLmuMinusLLp_secondMIP";
  tmp.PlotOnlySecondDaughterMIP = true;
  tmp.PlotOnlyContained = false;
  return tmp;
}

// TPCObj_PFP_MCSLLmuMinusLLp
CC1piPlotVars Var_TPCObj_PFP_MCSLLpiMinusLLp_NotPion(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MCSLLpiMinusLLp;
  tmp.bins = {88,-2,0.2};
  tmp.histtitle = ";(MCS LL_{#pi} - LL_{p})/MCS LL_{#pi} (MIP, not #pi candidate);";
  tmp.histname = "MCSLLmuMinusLLp_notPion";
  tmp.PlotNotPionCandidate = true;
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = false;
  return tmp;
}

// TPCObj_PFP_MCSLLmuMinusLLp: leading MIP, contained only
CC1piPlotVars Var_TPCObj_PFP_MCSLLpiMinusLLp_LeadingMIPcont(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MCSLLpiMinusLLp;
  tmp.bins = {88,-2,0.2};
  tmp.histtitle = ";(MCS LL_{#pi} - LL_{p})/MCS LL_{#pi} (Leading MIP, cont. only);";
  tmp.histname = "MCSLLmuMinusLLp_leadingMIPcont";
  tmp.PlotOnlyLeadingDaughterMIP = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_MCSLLmuMinusLLp: second MIP, contained only
CC1piPlotVars Var_TPCObj_PFP_MCSLLpiMinusLLp_SecondMIPcont(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MCSLLpiMinusLLp;
  tmp.bins = {88,-2,0.2};
  tmp.histtitle = ";(MCS LL_{#pi} - LL_{p})/MCS LL_{#pi} (Second MIP, cont. only);";
  tmp.histname = "MCSLLmuMinusLLp_secondMIPcont";
  tmp.PlotOnlySecondDaughterMIP = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_MCSLLmuMinusLLp: second MIP, contained only
CC1piPlotVars Var_TPCObj_PFP_MCSLLpiMinusLLp_NotPioncont(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MCSLLpiMinusLLp;
  tmp.bins = {88,-2,0.2};
  tmp.histtitle = ";(MCS LL_{#pi} - LL_{p})/MCS LL_{#pi} (MIP, not #pi candidate);";
  tmp.histname = "MCSLLmuMinusLLp_notPioncont";
  tmp.PlotNotPionCandidate = true;
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_track_MomRangeMinusMCS_p: leading daughter MIPs only
CC1piPlotVars Var_TPCObj_PFP_track_MomRangeMinusMCS_p_LeadingMIPcont(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MomRangeMinusMCS_p;
  tmp.bins = {60,-2,1};
  tmp.histtitle = ";(Mom. by range - Mom. by MCS)/(Mom. by range), p assumption (Leading MIP, cont. only);";
  tmp.histname = "MomRangeMinusMCS_p_LeadingMIPcont";
  tmp.PlotOnlyLeadingDaughterMIP = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_track_MomRangeMinusMCS_p: second daughter MIPs only
CC1piPlotVars Var_TPCObj_PFP_track_MomRangeMinusMCS_p_SecondMIPcont(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MomRangeMinusMCS_p;
  tmp.bins = {60,-2,2};
  tmp.histtitle = ";(Mom. by range - Mom. by MCS)/(Mom. by range), p assumption (Second MIP, cont. only);";
  tmp.histname = "MomRangeMinusMCS_p_SecondMIPcont";
  tmp.PlotOnlySecondDaughterMIP = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_track_MomRangeMinusMCS_p
CC1piPlotVars Var_TPCObj_PFP_track_MomRangeMinusMCS_p_NotPioncont(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MomRangeMinusMCS_p;
  tmp.bins = {60,-2,2};
  tmp.histtitle = ";(Mom. by range - Mom. by MCS)/(Mom. by range), p assumption (MIP, not #pi candidate, cont. only);";
  tmp.histname = "MomRangeMinusMCS_p_notPioncont";
  tmp.PlotNotPionCandidate = true;
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_track_MomRangeMinusMCS_mu: leading daughter MIP only
CC1piPlotVars Var_TPCObj_PFP_track_MomRangeMinusMCS_mu_LeadingMIPcont(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MomRangeMinusMCS_mu;
  tmp.bins = {60,-1.5,1.5};
  tmp.histtitle = ";(Mom. by range - Mom. by MCS)/(Mom. by range), #mu assumption (Leading MIP, cont. only);";
  tmp.histname = "MomRangeMinusMCS_mu_LeadingMIPcont";
  tmp.PlotOnlyLeadingDaughterMIP= true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_track_MomRangeMinusMCS_mu: second daughter MIP only
CC1piPlotVars Var_TPCObj_PFP_track_MomRangeMinusMCS_mu_SecondMIPcont(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MomRangeMinusMCS_mu;
  tmp.bins = {50,-5,1.5};
  tmp.histtitle = ";(Mom. by range - Mom. by MCS)/(Mom. by range), #mu assumption (Second MIP, cont. only);";
  tmp.histname = "MomRangeMinusMCS_mu_SecondMIPcont";
  tmp.PlotOnlySecondDaughterMIP= true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_track_MomRangeMinusMCS_mu: second daughter MIP only
CC1piPlotVars Var_TPCObj_PFP_track_MomRangeMinusMCS_mu_SecondMIPcont_mumupairs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MomRangeMinusMCS_mu;
  tmp.bins = {50,-5,1.5};
  tmp.histtitle = ";(Mom. by range - Mom. by MCS)/(Mom. by range), #mu assumption (Second MIP, cont., #mu^{-}#mu^{-} pairs only);";
  tmp.histname = "MomRangeMinusMCS_mu_SecondMIPcont_mumupairs";
  tmp.PlotOnlySecondDaughterMIP= true;
  tmp.PlotOnlyContained = true;
  tmp.PlotOnlyMuMuPairs = true;
  return tmp;
}

// TPCObj_PFP_track_MomRangeMinusMCS_mu
CC1piPlotVars Var_TPCObj_PFP_track_MomRangeMinusMCS_mu_NotPioncont(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MomRangeMinusMCS_mu;
  tmp.bins = {50,-5,1.5};
  tmp.histtitle = ";(Mom. by range - Mom. by MCS)/(Mom. by range), #mu assumption (MIP, not #pi candidate, cont. only);";
  tmp.histname = "MomRangeMinusMCS_mu_NotPioncont";
  tmp.PlotNotPionCandidate = true;
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_track_SimpleCluster_hitLinearity_minimum
CC1piPlotVars Var_TPCObj_PFP_track_SimpleCluster_hitLinearity_minimum_LeadingMIPcont(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_SimpleCluster_hitLinearity_minimum;
  tmp.bins = {50,0,1};
  tmp.histtitle = ";Local Linearity Minimum (Leading MIP, cont. only);";
  tmp.histname = "hitLinearity_minimum_LeadingMIPcont";
  tmp.PlotOnlyLeadingDaughterMIP = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_track_SimpleCluster_hitLinearity_minimum
CC1piPlotVars Var_TPCObj_PFP_track_SimpleCluster_hitLinearity_minimum_SecondMIPcont(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_SimpleCluster_hitLinearity_minimum;
  tmp.bins = {50,0,1};
  tmp.histtitle = ";Local Linearity Minimum (Second MIP, cont. only);";
  tmp.histname = "hitLinearity_minimum_SecondMIPcont";
  tmp.PlotOnlySecondDaughterMIP = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_track_SimpleCluster_hitLinearity_minimum
CC1piPlotVars Var_TPCObj_PFP_track_SimpleCluster_hitLinearity_minimum_SecondMIPcont_mumupairs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_SimpleCluster_hitLinearity_minimum;
  tmp.bins = {50,0,1};
  tmp.histtitle = ";Local Linearity Minimum (Second MIP, cont. #mu^{-}#u^{-} pairs only);";
  tmp.histname = "hitLinearity_minimum_SecondMIPcont_mumupairs";
  tmp.PlotOnlySecondDaughterMIP = true;
  tmp.PlotOnlyContained = true;
  tmp.PlotOnlyMuMuPairs = true;
  return tmp;
}

// TPCObj_PFP_track_SimpleCluster_hitLinearity_minimum
CC1piPlotVars Var_TPCObj_PFP_track_SimpleCluster_hitLinearity_minimum_NotPioncont(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_SimpleCluster_hitLinearity_minimum;
  tmp.bins = {50,0,1};
  tmp.histtitle = ";Local Linearity Minimum (MIP, not #pi candidate);";
  tmp.histname = "hitLinearity_minimum_contnotPion";
  tmp.PlotNotPionCandidate = true;
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_track_SimpleCluster_hitLinearity_mean
CC1piPlotVars Var_TPCObj_PFP_track_SimpleCluster_hitLinearity_mean_LeadingMIPcont(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_SimpleCluster_hitLinearity_mean;
  tmp.bins = {50,0,1};
  tmp.histtitle = ";Local Linearity Mean (Leading MIP, cont. only);";
  tmp.histname = "hitLinearity_mean_LeadingMIPcont";
  tmp.PlotOnlyLeadingDaughterMIP = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_track_SimpleCluster_hitLinearity_mean
CC1piPlotVars Var_TPCObj_PFP_track_SimpleCluster_hitLinearity_mean_SecondMIPcont(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_SimpleCluster_hitLinearity_mean;
  tmp.bins = {50,0,1};
  tmp.histtitle = ";Local Linearity Mean (Second MIP, cont. only);";
  tmp.histname = "hitLinearity_mean_SecondMIPcont";
  tmp.PlotOnlySecondDaughterMIP = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_track_SimpleCluster_hitLinearity_mean
CC1piPlotVars Var_TPCObj_PFP_track_SimpleCluster_hitLinearity_mean_SecondMIPcont_mumupairs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_SimpleCluster_hitLinearity_mean;
  tmp.bins = {50,0,1};
  tmp.histtitle = ";Local Linearity Mean (Second MIP, cont. #mu^{-}#mu^{-} pairs only);";
  tmp.histname = "hitLinearity_mean_SecondMIPcont_mumupairs";
  tmp.PlotOnlySecondDaughterMIP = true;
  tmp.PlotOnlyContained = true;
  tmp.PlotOnlyMuMuPairs = true;
  return tmp;
}

// TPCObj_PFP_track_SimpleCluster_hitLinearity_mean
CC1piPlotVars Var_TPCObj_PFP_track_SimpleCluster_hitLinearity_mean_NotPioncont(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_SimpleCluster_hitLinearity_mean;
  tmp.bins = {50,0,1};
  tmp.histtitle = ";Local Linearity Mean (MIP, not #pi candidate);";
  tmp.histname = "hitLinearity_mean_contnotPion";
  tmp.PlotNotPionCandidate = true;
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_track_SimpleCluster_hitLinearity_truncated_mean
CC1piPlotVars Var_TPCObj_PFP_track_SimpleCluster_hitLinearity_truncated_mean_LeadingMIPcont(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_SimpleCluster_hitLinearity_truncated_mean;
  tmp.bins = {50,0,1};
  tmp.histtitle = ";Local Linearity Truncated Mean (Leading MIP, cont. only);";
  tmp.histname = "hitLinearity_truncated_mean_LeadingMIPcont";
  tmp.PlotOnlyLeadingDaughterMIP = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_track_SimpleCluster_hitLinearity_truncated_mean
CC1piPlotVars Var_TPCObj_PFP_track_SimpleCluster_hitLinearity_truncated_mean_SecondMIPcont(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_SimpleCluster_hitLinearity_truncated_mean;
  tmp.bins = {50,0,1};
  tmp.histtitle = ";Local Linearity Truncated Mean (Second MIP, cont. only);";
  tmp.histname = "hitLinearity_truncated_mean_SecondMIPcont";
  tmp.PlotOnlySecondDaughterMIP = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_track_SimpleCluster_hitLinearity_truncated_mean
CC1piPlotVars Var_TPCObj_PFP_track_SimpleCluster_hitLinearity_truncated_mean_SecondMIPcont_mumupairs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_SimpleCluster_hitLinearity_truncated_mean;
  tmp.bins = {50,0,1};
  tmp.histtitle = ";Local Linearity Truncated Mean (Second MIP, cont. #mu^{-}#mu{-} pairs only);";
  tmp.histname = "hitLinearity_truncated_mean_SecondMIPcont_mumupairs";
  tmp.PlotOnlySecondDaughterMIP = true;
  tmp.PlotOnlyContained = true;
  tmp.PlotOnlyMuMuPairs = true;
  return tmp;
}

// TPCObj_PFP_track_SimpleCluster_hitLinearity_truncated_mean
CC1piPlotVars Var_TPCObj_PFP_track_SimpleCluster_hitLinearity_truncated_mean_NotPioncont(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_SimpleCluster_hitLinearity_truncated_mean;
  tmp.bins = {50,0,1};
  tmp.histtitle = ";Local Linearity Truncated Mean (MIP, not #pi candidate);";
  tmp.histname = "hitLinearity_truncated_mean_contnotPion";
  tmp.PlotNotPionCandidate = true;
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = true;
  return tmp;
}

// TPCObj_PFP_track_length_LeadingMIP
CC1piPlotVars Var_TPCObj_PFP_track_nhits_LeadingMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_nhits;
  tmp.KeepBelowCut = false;
  tmp.OnlyDaughters = false;
  tmp.TracksNeeded = "all";
  tmp.CutValue = 10;
  tmp.bins = {400,0,1000};
  tmp.histtitle = ";No. hits in track (leading MIP);";
  tmp.histname = "trknhits_leadingMIP";
  tmp.PlotOnlyLeadingDaughterMIP = true;
  return tmp;
}

// TPCObj_PFP_track_length_SecondMIP
CC1piPlotVars Var_TPCObj_PFP_track_nhits_SecondMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_nhits;
  tmp.bins = {400,0,1000};
  tmp.histtitle = ";No. hits in track (second MIP);";
  tmp.histname = "trknhits_secondMIP";
  tmp.PlotOnlySecondDaughterMIP = true;
  return tmp;
}

// TPCObj_PFP_track_length_notPion
CC1piPlotVars Var_TPCObj_PFP_track_nhits_NotPion(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_nhits;
  tmp.bins = {400,0,1000};
  tmp.histtitle = ";No. hits in track (MIP, not #pi candidate);";
  tmp.histname = "trknhits_MIPnotPion";
  tmp.PlotNotPionCandidate = true;
  tmp.PlotOnlyDaughterMIPs = true;
  return tmp;
}

// TPCObj_PFP_ndaughters_SecondMIP
CC1piPlotVars Var_TPCObj_PFP_ndaughters_SecondMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_ndaughters;
  tmp.bins = {10,0,10};
  tmp.histtitle = ";No. daughter PFPs (second MIP);";
  tmp.histname = "trkndaughters_secondMIP";
  tmp.PlotOnlySecondDaughterMIP = true;
  return tmp;
}

// TPCObj_PFP_ndaughters_SecondMIP
CC1piPlotVars Var_TPCObj_PFP_ndaughters_SecondMIP_mumupairs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_ndaughters;
  tmp.bins = {10,0,10};
  tmp.histtitle = ";No. daughter PFPs (second MIP, #mu^{-}#mu^{-} pairs only);";
  tmp.histname = "trkndaughters_secondMIP_mumupairs";
  tmp.PlotOnlySecondDaughterMIP = true;
  tmp.PlotOnlyMuMuPairs = true;
  return tmp;
}

// TPCObj_PFP_ndaughters_SecondMIP
CC1piPlotVars Var_TPCObj_PFP_ndaughters_isPion(treevars *vars){
  CC1piPlotVars tmp;
  tmp.KeepBelowCut = false;
  tmp.CutValue = 1.5;
  tmp.Var = vars->TPCObj_PFP_ndaughters;
  tmp.bins = {10,0,10};
  tmp.histtitle = ";No. daughter PFPs (#pi candidate);";
  tmp.histname = "trkndaughters_pioncand";
  tmp.PlotOnlyPionCandidate = true;
  return tmp;
}

// TPCObj_PFP_ndaughters_SecondMIP
CC1piPlotVars Var_TPCObj_PFP_ndaughters_notPion(treevars *vars){
  CC1piPlotVars tmp;
  tmp.KeepBelowCut = false;
  tmp.CutValue = 1.5;
  tmp.Var = vars->TPCObj_PFP_ndaughters;
  tmp.bins = {10,0,10};
  tmp.histtitle = ";No. daughter PFPs (MIP, not #pi candidate);";
  tmp.histname = "trkndaughters_notpicand";
  tmp.PlotNotPionCandidate = true;
  tmp.PlotOnlyDaughterMIPs = true;
  return tmp;
}

// TPCObj_PFP_track_length
CC1piPlotVars Var_TPCObj_PFP_track_length(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_length;
  tmp.bins = {10,0,140};
  tmp.histtitle = ";Track length [cm];";
  tmp.histname = "trklength";
  return tmp;
}

// TPCObj_PFP_track_length
CC1piPlotVars Var_TPCObj_PFP_track_length_LeadingMIP_mumupairs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_length;
  tmp.bins = {10,0,140};
  tmp.histtitle = ";Track length (Leading MIP, #mu^{-}#mu^{-} pairs only) [cm];";
  tmp.histname = "trklength_LeadingMIP_mumupairs";
  tmp.PlotOnlyMuMuPairs = true;
  tmp.PlotOnlyLeadingDaughterMIP = true;
  return tmp;
}

// TPCObj_PFP_track_length
CC1piPlotVars Var_TPCObj_PFP_track_length_SecondMIP_mumupairs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_length;
  tmp.bins = {10,0,140};
  tmp.histtitle = ";Track length (Second MIP, #mu^{-}#mu^{-} pairs only) [cm];";
  tmp.histname = "trklength_SecondMIP_mumupairs";
  tmp.PlotOnlyMuMuPairs = true;
  tmp.PlotOnlySecondDaughterMIP = true;
  return tmp;
}

// TPCObj_PFP_track_MCSpi_maxScatter
CC1piPlotVars Var_TPCObj_PFP_track_MCSpi_maxScatter(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MCS_pi_maxScatter;
  tmp.bins = {50,0,1000};
  tmp.histtitle = ";Max. MCS scatter angle [mrad.] (#pi assumption);";
  tmp.histname = "trkMCSpi_maxScatter";
  tmp.PlotOnlyDaughterMIPs = true;
  return tmp;
}

// TPCObj_PFP_track_MCSpi_maxScatter
CC1piPlotVars Var_TPCObj_PFP_track_MCSpi_maxScatter_SecondMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MCS_pi_maxScatter;
  tmp.bins = {50,0,1000};
  tmp.histtitle = ";Max. MCS scatter angle [mrad.] (#pi assumption, Second MIP only);";
  tmp.histname = "trkMCSpi_maxScatter_SecondMIP";
  tmp.PlotOnlySecondDaughterMIP = true;
  return tmp;
}

// TPCObj_PFP_track_MCSpi_maxScatter
CC1piPlotVars Var_TPCObj_PFP_track_MCSpi_maxScatter_SecondMIP_mumupairs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MCS_pi_maxScatter;
  tmp.bins = {50,0,1000};
  tmp.histtitle = ";Max. MCS scatter angle [mrad.] (#pi assumption, Second MIP, #mu^{-}#mu^{-} pairs only);";
  tmp.histname = "trkMCSpi_maxScatter_SecondMIP_mumupairs";
  tmp.PlotOnlySecondDaughterMIP = true;
  tmp.PlotOnlyMuMuPairs = true;
  return tmp;
}

// TPCObj_PFP_track_MCSpi_maxScatter
CC1piPlotVars Var_TPCObj_PFP_track_MCSp_maxScatter(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MCS_p_maxScatter;
  tmp.bins = {50,0,1000};
  tmp.histtitle = ";Max. MCS scatter angle [mrad.] (p assumption);";
  tmp.histname = "trkMCSp_maxScatter";
  tmp.PlotOnlyDaughterMIPs = true;
  return tmp;
}

// TPCObj_PFP_track_MCSpi_maxScatter
CC1piPlotVars Var_TPCObj_PFP_track_MCSp_maxScatter_SecondMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MCS_p_maxScatter;
  tmp.bins = {50,0,1000};
  tmp.histtitle = ";Max. MCS scatter angle [mrad.] (p assumption, Second MIP only);";
  tmp.histname = "trkMCSp_maxScatter_SecondMIP";
  tmp.PlotOnlySecondDaughterMIP = true;
  return tmp;
}

// TPCObj_PFP_track_MCSp_maxScatter
CC1piPlotVars Var_TPCObj_PFP_track_MCSp_maxScatter_SecondMIP_mumupairs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MCS_p_maxScatter;
  tmp.bins = {50,0,1000};
  tmp.histtitle = ";Max. MCS scatter angle [mrad.] (p assumption, Second MIP, #mu^{-}#mu^{-} pairs only);";
  tmp.histname = "trkMCSp_maxScatter_SecondMIP_mumupairs";
  tmp.PlotOnlySecondDaughterMIP = true;
  tmp.PlotOnlyMuMuPairs = true;
  return tmp;
}

// TPCObj_PFP_track_MCSpi_meanScatter
CC1piPlotVars Var_TPCObj_PFP_track_MCSpi_meanScatter(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MCS_pi_meanScatter;
  tmp.bins = {50,0,1000};
  tmp.histtitle = ";Mean MCS scatter angle [mrad.] (#pi assumption);";
  tmp.histname = "trkMCSpi_meanScatter";
  tmp.PlotOnlyDaughterMIPs = true;
  return tmp;
}

// TPCObj_PFP_track_MCSpi_meanScatter
CC1piPlotVars Var_TPCObj_PFP_track_MCSpi_meanScatter_SecondMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MCS_pi_meanScatter;
  tmp.bins = {50,0,1000};
  tmp.histtitle = ";Mean MCS scatter angle [mrad.] (#pi assumption, Second MIP only);";
  tmp.histname = "trkMCSpi_meanScatter_SecondMIP";
  tmp.PlotOnlySecondDaughterMIP = true;
  return tmp;
}

// TPCObj_PFP_track_MCSpi_meanScatter
CC1piPlotVars Var_TPCObj_PFP_track_MCSpi_meanScatter_SecondMIP_mumupairs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MCS_pi_meanScatter;
  tmp.bins = {50,0,1000};
  tmp.histtitle = ";Mean MCS scatter angle [mrad.] (#pi assumption, Second MIP, #mu^{-}#mu^{-} pairs only);";
  tmp.histname = "trkMCSpi_meanScatter_SecondMIP_mumupairs";
  tmp.PlotOnlySecondDaughterMIP = true;
  tmp.PlotOnlyMuMuPairs = true;
  return tmp;
}

// TPCObj_PFP_track_MCSpi_meanScatter
CC1piPlotVars Var_TPCObj_PFP_track_MCSp_meanScatter(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MCS_p_meanScatter;
  tmp.bins = {50,0,1000};
  tmp.histtitle = ";Mean MCS scatter angle [mrad.] (p assumption);";
  tmp.histname = "trkMCSp_meanScatter";
  tmp.PlotOnlyDaughterMIPs = true;
  return tmp;
}

// TPCObj_PFP_track_MCSpi_meanScatter
CC1piPlotVars Var_TPCObj_PFP_track_MCSp_meanScatter_SecondMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MCS_p_meanScatter;
  tmp.bins = {50,0,1000};
  tmp.histtitle = ";Mean MCS scatter angle [mrad.] (p assumption, Second MIP only);";
  tmp.histname = "trkMCSp_meanScatter_SecondMIP";
  tmp.PlotOnlySecondDaughterMIP = true;
  return tmp;
}

// TPCObj_PFP_track_MCSp_meanScatter
CC1piPlotVars Var_TPCObj_PFP_track_MCSp_meanScatter_SecondMIP_mumupairs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_MCS_p_meanScatter;
  tmp.bins = {50,0,1000};
  tmp.histtitle = ";Mean MCS scatter angle [mrad.] (p assumption, Second MIP, #mu^{-}#mu^{-} pairs only);";
  tmp.histname = "trkMCSp_meanScatter_SecondMIP_mumupairs";
  tmp.PlotOnlySecondDaughterMIP = true;
  tmp.PlotOnlyMuMuPairs = true;
  return tmp;
}

// TPCObj_PFP_track_MCSpi_meanScatter
CC1piPlotVars Var_TPCObj_PFP_track_dedx_grminhits(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_dedx_grminhits;
  tmp.KeepBelowCut = false;
  tmp.isMIPcut = true;
  tmp.CutValue = 0.5;
  tmp.histtitle = ";Passes min. no. hits threshold for collection plane;";
  tmp.histname = "grminhits";
  tmp.bins = {2,0,2};
  tmp.PlotOnlyDaughters = true;
  return tmp;
}

CC1piPlotVars Var_TPCObj_PFP_track_dEdx_nhits(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_dEdx_nhits;
  tmp.KeepBelowCut = false;
  tmp.isMIPcut = true;
  tmp.CutValue = 33;
  tmp.histtitle = ";Number of collection plane hits;";
  tmp.histname = "nhits";
  tmp.bins = {25,0,300};
  tmp.PlotOnlyDaughters = true;
  return tmp;
}

// Longest MIP: originally classed as track?
CC1piPlotVars Var_TPCObj_LeadingMIP_PandoraClassedAsTrack(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_PandoraClassedAsTrack_double;
  tmp.KeepBelowCut = false;
  tmp.CutValue = 0.5;
  tmp.TracksNeeded = "LeadingMIP";
  tmp.bins = {2,-1,1};
  tmp.histtitle = ";Pandora classed as track? (1=true, 0=false) (Leading MIP only);";
  tmp.histname = "Pandoraclassedastrack_LeadingMIP";
  tmp.PlotOnlyLeadingDaughterMIP= true;
  return tmp;
}

// Longest MIP: originally classed as track?
CC1piPlotVars Var_TPCObj_SecondMIP_isContained(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_isContained_double;
  tmp.KeepBelowCut = false;
  tmp.OnlyDaughters = true;
  tmp.TracksNeeded = "SecondMIP";
  tmp.CutValue = 0.5;
  tmp.histtitle = ";Second MIP contained;";
  return tmp;
}

// Longest MIP: originally classed as track?
CC1piPlotVars Var_TPCObj_FirstMIP_isContained(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_isContained_double;
  tmp.KeepBelowCut = false;
  tmp.OnlyDaughters = true;
  tmp.TracksNeeded = "LeadingMIP";
  tmp.CutValue = 0.5;
  tmp.histtitle = ";First MIP contained;";
  return tmp;
}

// Longest MIP: originally classed as track?
CC1piPlotVars Var_TPCObj_AllDaughters_isContained(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_isContained_double;
  tmp.KeepBelowCut = false;
  tmp.OnlyDaughters = true;
  tmp.TracksNeeded = "all";
  tmp.CutValue = 0.5;
  tmp.histtitle = ";All daughters contained;";
  return tmp;
}

// Longest MIP: originally classed as track?
CC1piPlotVars Var_TPCObj_AllDaughtersExceptLeadingMIP_isContained(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_isContained_double;
  tmp.KeepBelowCut = false;
  tmp.OnlyDaughters = true;
  tmp.TracksNeeded = "allExceptLeadingMIP";
  tmp.CutValue = 0.5;
  tmp.histtitle = ";All daughters except leading MIP contained;";
  return tmp;
}

// Longest MIP: originally classed as track?
CC1piPlotVars Var_TPCObj_AllTracks_isContained(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_isContained_double;
  tmp.KeepBelowCut = false;
  tmp.OnlyDaughters = false;
  tmp.TracksNeeded = "all";
  tmp.CutValue = 0.5;
  tmp.histtitle = ";All tracks contained;";
  return tmp;
}

// Longest MIP: originally classed as track?
CC1piPlotVars Var_TPCObj_NotMIPs_isContained(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_isContained_double;
  tmp.histtitle = ";Track is contained? (Non-MIP daughters);";
  tmp.bins = {2,0,2};
  tmp.histname = "isContained_nonMIPDaughters";
  tmp.PlotOnlyNotMIPDaughters = true;
  return tmp;
}

//
CC1piPlotVars Var_TPCObj_PFP_trueKE_selMIPs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_trueKE;
  tmp.bins = {100,0,5};
  tmp.histtitle = ";True Kinetic Energy (MIP candidates only);";
  tmp.histname = "TrueKE_MIPs";
  tmp.PlotOnlyDaughterMIPs = true;
  return tmp;
}

//
CC1piPlotVars Var_TPCObj_PFP_trueKE_SecondMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_trueKE;
  tmp.bins = {100,0,5};
  tmp.histtitle = ";True Kinetic Energy (Second MIP candidates only);";
  tmp.histname = "TrueKE_SecondMIP";
  tmp.PlotOnlySecondDaughterMIP = true;
  return tmp;
}

//
CC1piPlotVars Var_TPCObj_PFP_trueEndP_selMIPs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_trueEndP;
  tmp.bins = {100,0,0.1};
  tmp.histtitle = ";True momentum at end of track (MIP candidates only);";
  tmp.histname = "TrueEndP_MIPs";
  tmp.PlotOnlyDaughterMIPs = true;
  return tmp;
}

//
CC1piPlotVars Var_TPCObj_PFP_trueEndP_SecondMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_trueEndP;
  tmp.bins = {100,0,0.1};
  tmp.histtitle = ";True momentum at end of track (MIP candidates only);";
  tmp.histname = "TrueEndP_SecondMIP";
  tmp.PlotOnlySecondDaughterMIP = true;
  return tmp;
}

//
CC1piPlotVars Var_TPCObj_PFP_track_theta_SecondMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_theta;
  tmp.bins = {30,0,3.15};
  tmp.histtitle = ";Track theta (Second MIP only);";
  tmp.histname = "Theta_SecondMIP";
  tmp.PlotOnlySecondDaughterMIP = true;
  return tmp;
}

CC1piPlotVars Var_TPCObj_PFP_track_theta(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_theta;
  tmp.bins = {30,0,3.15};
  tmp.histtitle = ";Track theta [rad];";
  tmp.histname = "Theta";
  tmp.PlotOnlyDaughters = true;
  return tmp;
}

//
CC1piPlotVars Var_TPCObj_PFP_track_phi_SecondMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_phi;
  tmp.bins = {60,-3.15,3.15};
  tmp.histtitle = ";Track phi (Second MIP only);";
  tmp.histname = "Phi_SecondMIP";
  tmp.PlotOnlySecondDaughterMIP = true;
  return tmp;
}

//
CC1piPlotVars Var_TPCObj_PFP_MCP_PDG_mTruePDG(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_MCP_PDG_mTruePDG;
  tmp.bins = {100,-1,1};
  tmp.histtitle = ";MCP PDG - 'True PDG';";
  tmp.histname = "MCPPDG_mTruePDG";
  return tmp;
}

//
CC1piPlotVars Var_TPCObj_PFP_MCP_numdaughters_SecondMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_MCP_numdaughters_notphotons;
  tmp.bins = {100,0,100};
  tmp.histtitle = ";MCP: no. daughters (excl. #gamma);";
  tmp.histname = "MCP_numdaughters_notphotons_SecondMIP";
  tmp.PlotOnlySecondDaughterMIP = true;
  return tmp;
}

//
CC1piPlotVars Var_TPCObj_PFP_MCP_motherIDeq0_SecondMIP(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_MCP_motherIDeq0;
  tmp.bins = {2,0,2};
  tmp.histtitle = ";MCP: mother ID=0? (1=yes, 0=no) (Second MIP only);";
  tmp.histname = "MCP_motherIDeq0_SecondMIP";
  tmp.PlotOnlySecondDaughterMIP = true;
  return tmp;
}

// Truncated mean dE/dx difference for MIP candidates
CC1piPlotVars Var_TPCObj_dEdx_truncmean_MIPdiff(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_dEdx_truncmean_MIPdiff;
  tmp.KeepBelowCut = true;
  tmp.TracksNeeded = "NA";
  tmp.CutValue = 0.3;
  tmp.bins = {10,0,1};
  tmp.histtitle = ";Difference between MIP candidates' truncated mean dE/dx;";
  tmp.histname = "dEdx_truncmean_MIPdiff";
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = false;
  return tmp;
}

// Truncated mean dE/dx difference for MIP candidates
CC1piPlotVars Var_TPCObj_dEdx_truncmean_MIPdiff_mupi(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_dEdx_truncmean_MIPdiff_mupi;
  tmp.bins = {10,0,1};
  tmp.histtitle = ";Difference between MIP candidates' truncated mean dE/dx (true #mu-#pi selected);";
  tmp.histname = "dEdx_truncmean_MIPdiff_mupi";
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = false;
  return tmp;
}

// Truncated mean dE/dx difference for MIP candidates
CC1piPlotVars Var_TPCObj_dEdx_truncmean_MIPdiff_muproton(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_dEdx_truncmean_MIPdiff_muproton;
  tmp.bins = {10,0,1};
  tmp.histtitle = ";Difference between MIP candidates' truncated mean dE/dx (true #mu-p selected);";
  tmp.histname = "dEdx_truncmean_MIPdiff_muproton";
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = false;
  return tmp;
}

// Truncated mean dE/dx difference for MIP candidates
CC1piPlotVars Var_TPCObj_dEdx_truncmean_MIPdiff_other(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_dEdx_truncmean_MIPdiff_other;
  tmp.bins = {10,0,1};
  tmp.histtitle = ";Difference between MIP candidates' truncated mean dE/dx (other);";
  tmp.histname = "dEdx_truncmean_MIPdiff_other";
  tmp.PlotOnlyDaughterMIPs = true;
  tmp.PlotOnlyContained = false;
  return tmp;
}
/*
// BDT score
CC1piPlotVars Var_TPCObj_PFP_track_BDTscore(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_BDTscore;
  tmp.KeepBelowCut = false;
//  tmp.OnlyDaughters = true;
//  tmp.TracksNeeded = "exactlytwo";
  tmp.isMIPcut = true;
  tmp.CutValue = 0.5;
  tmp.bins = {25,-1,1};
  tmp.histtitle = ";BDT score;";
  tmp.histname = "BDTscore";
  return tmp;
}
*/

CC1piPlotVars Var_TPCObj_PFP_track_BDTscore_combined(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_PFP_track_BDTscore_combined;
  tmp.KeepBelowCut = false;
//  tmp.OnlyDaughters = true;
//  tmp.TracksNeeded = "exactlytwo";
  tmp.isMIPcut = true;
  tmp.CutValue = 0.5;
  tmp.bins = {2,0,2};
  tmp.histtitle = ";BDT score (combined);";
  tmp.histname = "BDTscore_combined";
  return tmp;
}

CC1piPlotVars Var_TPCObj_NDaughterPFPs(treevars *vars){
  CC1piPlotVars tmp;
  tmp.Var = vars->TPCObj_NDaughterPFPs;
//  tmp.OnlyDaughters = true;
//  tmp.TracksNeeded = "exactlytwo";
  tmp.bins = {10,0,10};
  tmp.histtitle = ";Number of #nu daughter PFPs in selected TPCObject;";
  tmp.histname = "NDaughterPFPs";
  return tmp;
}

// Decide whether to fill a plot for a given track. This is useful for e.g. plots with PlotOnlyDaughterMIPs=true
bool FillPlotForTrack(CC1piPlotVars *plotvar, treevars *vars, int i_tr){
  bool DoFillPlot = true;

  if (plotvar->PlotOnlyDaughters && !(vars->TPCObj_PFP_isDaughter->at(i_tr))) DoFillPlot = false;

  if (plotvar->PlotOnlyDaughterMIPs && !(vars->TPCObj_PFP_isDaughter->at(i_tr) && vars->TPCObj_PFP_track_passesMIPcut->at(i_tr))) DoFillPlot = false;

  if (plotvar->PlotOnlyContained && !(vars->TPCObj_PFP_track_isContained->at(i_tr))) DoFillPlot = false;

  if (plotvar->PlotOnlyLeadingDaughterMIP && !(i_tr == vars->TPCObj_LeadingMIPtrackIndex)) DoFillPlot = false;

  if (plotvar->PlotOnlySecondDaughterMIP && !(i_tr == vars->TPCObj_SecondMIPtrackIndex)) DoFillPlot = false;

  if (plotvar->PlotOnlyPionCandidate && !(vars->TPCObj_PFP_track_passesPioncut->at(i_tr)==1)) DoFillPlot = false;

  if (plotvar->PlotNotPionCandidate && (vars->TPCObj_PFP_track_passesPioncut->at(i_tr)==1)) DoFillPlot = false;

  if (plotvar->PlotOnlyNotMIPDaughters && !(vars->TPCObj_PFP_isDaughter->at(i_tr) && vars->TPCObj_PFP_track_passesMIPcut->at(i_tr)==false)) DoFillPlot = false;

  if (plotvar->PlotOnlyMuMuPairs){

    // if (plotvar->histname=="AngleBetweenMIPs_mumupairs" && vars->TPCObj_LeadingMIPtrackIndex>=0){
    //   std::cout << vars->TPCObj_LeadingMIPtrackIndex << ", " << vars->TPCObj_SecondMIPtrackIndex << ", " << vars->TPCObj_AngleBetweenMIPs->at(0) << std::endl;
    // }

    if (!(vars->TPCObj_LeadingMIPtrackIndex>=0 && vars->TPCObj_SecondMIPtrackIndex>=0)) DoFillPlot = false;
    else if (!(vars->TPCObj_PFP_truePDG->at(vars->TPCObj_LeadingMIPtrackIndex)==13 && vars->TPCObj_PFP_truePDG->at(vars->TPCObj_SecondMIPtrackIndex)==13)) DoFillPlot = false;
  }

  return DoFillPlot;
}

#endif
