#ifndef __TREEVARS_H__
#define __TREEVARS_H__

#include "../Algorithms/TopologyEnums.h"

#include "TMath.h"
#include "TVector3.h"
#include "TTree.h"
#include "TMVA/Reader.h"

// This code defines the variables we read in from the tree and new derived variables from there.
// The implementation (i.e. the code that actually calculates the new derived variables) is in treevars_header.h


struct treevars{
   // These are the variables that are filled directly from the tree
   std::vector<std::vector<double>> *TPCObj_PFP_LH_fwd_mu = nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_LH_fwd_p = nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_LH_fwd_pi = nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_LH_bwd_mu = nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_LH_bwd_p = nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_LH_bwd_pi = nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_LH_MIP = nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_track_Chi2Proton = nullptr;
   std::vector<double> *TPCObj_PFP_track_rangeE_p = nullptr;
   std::vector<double> *TPCObj_PFP_track_rangeE_mu = nullptr;
   NuIntTopology Truth_topology = kUnknown;
   bool Marco_selected = false;
   std::vector<bool> *TPCObj_PFP_isDaughter = nullptr;
   std::vector<std::vector<int>> *TPCObj_PFP_daughterids = nullptr;
   std::vector<bool> *TPCObj_PFP_PandoraClassedAsTrack = nullptr;
   std::vector<double> *TPCObj_PFP_track_theta = nullptr;
   std::vector<double> *TPCObj_PFP_track_phi = nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_track_start = nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_track_end = nullptr;
   std::vector<double> *TPCObj_PFP_track_residual_mean = nullptr;
   std::vector<double> *TPCObj_PFP_track_residual_std = nullptr;
   std::vector<double> *TPCObj_PFP_track_perc_used_hits = nullptr;
   std::vector<double> *TPCObj_reco_vtx = nullptr;
   std::vector<bool> *TPCObj_PFP_track_isContained = nullptr;
   std::vector<int> *TPCObj_PFP_truePDG = nullptr;
   std::vector<double> *TPCObj_PFP_trueKE = nullptr;
   std::vector<double> *TPCObj_PFP_trueEndP = nullptr;
   std::vector<double> *TPCObj_PFP_track_length = nullptr;
   std::vector<int> *TPCObj_PFP_MCPid = nullptr;
   std::vector<int> *MCP_PDG = nullptr;
   std::vector<int> *MCP_numdaughters = nullptr;
   std::vector<int> *MCP_MotherID = nullptr;
   std::vector<int> *MCP_ID = nullptr;
   std::vector<std::vector<int>> *MCP_DaughterIDs = nullptr;

   std::vector<std::vector<std::vector<double>>> *TPCObj_PFP_track_dedx_perhit=nullptr;
   std::vector<std::vector<std::vector<double>>> *TPCObj_PFP_track_resrange_perhit=nullptr;
   std::vector<double>* TPCObj_PFP_track_nhits=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSmu_fwdMom=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSmu_bwdMom=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSmu_bestMom=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSmu_fwdMomUncert=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSmu_bwdMomUncert=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSmu_bestMomUncert=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSmu_fwdLL=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSmu_bwdLL=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSmu_bestLL=nullptr;
   std::vector<std::vector<double>>* TPCObj_PFP_track_MCSmu_segmentRadLengths=nullptr;
   std::vector<std::vector<double>>* TPCObj_PFP_track_MCSmu_scatterAngles=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSp_fwdMom=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSp_bwdMom=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSp_bestMom=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSp_fwdMomUncert=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSp_bwdMomUncert=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSp_bestMomUncert=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSp_fwdLL=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSp_bwdLL=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSp_bestLL=nullptr;
   std::vector<std::vector<double>>* TPCObj_PFP_track_MCSp_segmentRadLengths=nullptr;
   std::vector<std::vector<double>>* TPCObj_PFP_track_MCSp_scatterAngles=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSpi_fwdMom=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSpi_bwdMom=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSpi_bestMom=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSpi_fwdMomUncert=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSpi_bwdMomUncert=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSpi_bestMomUncert=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSpi_fwdLL=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSpi_bwdLL=nullptr;
   std::vector<double>* TPCObj_PFP_track_MCSpi_bestLL=nullptr;
   std::vector<std::vector<double>>* TPCObj_PFP_track_MCSpi_segmentRadLengths=nullptr;
   std::vector<std::vector<double>>* TPCObj_PFP_track_MCSpi_scatterAngles=nullptr;

   std::vector<std::vector<std::vector<double>>> *TPCObj_PFP_track_SpacepointsXYZ=nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_track_SpacepointsQPlane2=nullptr;

   unsigned int run_num;
   unsigned int subrun_num;
   unsigned int event_num;


   // These are derived quantities - derived from the values above in Calcvars
   std::vector<double> *TPCObj_PFP_isDaughter_double;
   std::vector<double> *TPCObj_PFP_PandoraClassedAsTrack_double;
   std::vector<double> *TPCObj_PFP_LH_p;
   std::vector<double> *TPCObj_PFP_LH_mu;
   std::vector<double> *TPCObj_PFP_LH_pi;
   std::vector<double> *TPCObj_PFP_LH_mip;
   std::vector<double> *TPCObj_PFP_Lmipoverp;
   std::vector<double> *TPCObj_PFP_lnLmipoverp;
   std::vector<double> *TPCObj_PFP_Lmumipovermumipp;
   std::vector<double> *TPCObj_PFP_track_Chi2Proton_plane2;
   std::vector<double> *TPCObj_PFP_BrokenTrackAngle;
   std::vector<double> *TPCObj_PFP_VtxTrackDist;
   std::vector<double> *TPCObj_PFP_isContained_double;
   std::vector<double> *TPCObj_PFP_ndaughters;
   std::vector<double> *TPCObj_PFP_track_dEdx_truncmean_start;
   std::vector<double> *TPCObj_PFP_track_dEdx_mean_start;
   std::vector<double> *TPCObj_PFP_track_dEdx_truncmoverm_start;
   std::vector<double> *TPCObj_PFP_track_dEdx_stddev_start;
   std::vector<double> *TPCObj_PFP_track_dedx_grminhits;
   std::vector<double> *TPCObj_PFP_track_dEdx_nhits;
   std::vector<double> *TPCObj_PFP_track_dedx_stddev;
   std::vector<double> *TPCObj_DaughterTracks_Order_dEdxtr;
   std::vector<double> *TPCObj_DaughterTracks_Order_dEdxtr_selMIPs;
   std::vector<double> *TPCObj_DaughterTracks_Order_trklen_selMIPs;
   std::vector<double> *TPCObj_PFP_track_MCSLLpiMinusLLp;
   std::vector<double> *TPCObj_PFP_track_MomRangeMinusMCS_p;
   std::vector<double> *TPCObj_PFP_track_MomRangeMinusMCS_mu;
   std::vector<double> *TPCObj_PFP_track_MCS_pi_maxScatter;
   std::vector<double> *TPCObj_PFP_track_MCS_p_maxScatter;
   std::vector<double> *TPCObj_PFP_track_MCS_pi_meanScatter;
   std::vector<double> *TPCObj_PFP_track_MCS_p_meanScatter;
   std::vector<double> *TPCObj_PFP_MCP_PDG;
   std::vector<double> *TPCObj_PFP_MCP_numdaughters;
   std::vector<double> *TPCObj_PFP_MCP_numdaughters_notphotons;
   std::vector<double> *TPCObj_PFP_MCP_motherIDeq0;
   std::vector<double> *TPCObj_PFP_MCP_PDG_mTruePDG;
   std::vector<double> *TPCObj_PFP_track_BDTscore;
   std::vector<double> *TPCObj_NDaughterPFPs;
   std::vector<double> *TPCObj_PFP_MCP_trueOrigPDG;

   std::vector<std::vector<std::vector<double>>> *TPCObj_PFP_track_SpacepointsXYZ_Ordered;
   std::vector<std::vector<double>> *TPCObj_PFP_track_SpacepointsQPlane2_Ordered;
   std::vector<std::vector<double>> *TPCObj_PFP_track_SpacepointsdQdsPlane2_Ordered;
   std::vector<std::vector<double>> *TPCObj_PFP_track_Spacepoints_LocalLin;
   std::vector<std::vector<double>> *TPCObj_PFP_track_Spacepoints_Dirz;
   std::vector<std::vector<double>> *TPCObj_PFP_track_Spacepoints_DirCuSum;
   std::vector<std::vector<int>> *TPCObj_PFP_track_Spacepoints_kinkidxs;
   std::vector<double> *TPCObj_PFP_track_Spacepoints_ttest_max;


   // Is the track a MIP? Evaluate MIP cuts and then put this as an input into selection
   // Coded as a double (to fit in with the other code) but should be evaluated as a bool
   std::vector<double> *TPCObj_PFP_track_passesMIPcut;
   std::vector<double> *TPCObj_PFP_track_passesPioncut;

   // Not really a vector, only going to have one entry
   // But easier to store this way so it's compatible with the plotting code
   std::vector<double> *TPCObj_AngleBetweenMIPs;
   std::vector<double> *TPCObj_dEdx_truncmean_MIPdiff;
   std::vector<double> *TPCObj_dEdx_truncmean_MIPdiff_mupi;
   std::vector<double> *TPCObj_dEdx_truncmean_MIPdiff_muproton;
   std::vector<double> *TPCObj_dEdx_truncmean_MIPdiff_other;
   std::vector<double> *TPCObj_mupiContained;
   std::vector<double> *TPCObj_muonLonger;
   std::vector<double> *containment1_pass;
   std::vector<double> *containment2_pass;

   int containment0_muoncandidatePDG;
   int containment0_pioncandidatePDG;
   int containment1_muoncandidatePDG;
   int containment1_pioncandidatePDG;
   int containment2_muoncandidatePDG;
   int containment2_pioncandidatePDG;

   int TPCObj_LeadingMIPtrackIndex;
   int TPCObj_SecondMIPtrackIndex;

   float float_dEdx_truncmean_start;
   float float_VtxTrackDist;
   float float_nhits;
   float float_lnLmipoverp;

};


struct MVAvars {
   double dEdx_truncmean_start;
   double VtxTrackDist;
   double nhits;
   double lnLmipoverp;
};

#endif
