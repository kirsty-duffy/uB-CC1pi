#ifndef __TREEVARS_H__
#define __TREEVARS_H__

#include "../Algorithms/TopologyEnums.h"

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
   std::vector<double> *TPCObj_PFP_track_rangeE_p = nullptr;
   std::vector<double> *TPCObj_PFP_track_rangeE_mu = nullptr;
   NuIntTopology Truth_topology = kUnknown;
   bool Marco_selected = false;
   std::vector<bool> *TPCObj_PFP_isDaughter = nullptr;
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
   std::vector<double> *TPCObj_PFP_track_length = nullptr;

   std::vector<std::vector<std::vector<double>>> *TPCObj_PFP_track_trajPoint_Position=nullptr;
   std::vector<std::vector<std::vector<double>>> *TPCObj_PFP_track_trajPoint_Direction=nullptr;
   std::vector<std::vector<std::vector<double>>> *TPCObj_PFP_track_dedx_perhit=nullptr;
   std::vector<std::vector<std::vector<double>>> *TPCObj_PFP_track_resrange_perhit=nullptr;
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

   std::vector<std::vector<double>>* TPCObj_PFP_track_SimpleCluster_hitTime=nullptr;
   std::vector<std::vector<int>>* TPCObj_PFP_track_SimpleCluster_hitWire=nullptr;
   std::vector<std::vector<int>>* TPCObj_PFP_track_SimpleCluster_hitPlane=nullptr;
   std::vector<std::vector<double>>* TPCObj_PFP_track_SimpleCluster_hitIntegral=nullptr;
   std::vector<std::vector<double>>* TPCObj_PFP_track_SimpleCluster_hitTimeTicks=nullptr;
   std::vector<std::vector<int>>* TPCObj_PFP_track_SimpleCluster_hitWireNo=nullptr;
   std::vector<int>* TPCObj_PFP_track_SimpleCluster_StartIndex=nullptr;
   std::vector<std::vector<double>>* TPCObj_PFP_track_SimpleCluster_hitdQds=nullptr;
   std::vector<std::vector<double>>* TPCObj_PFP_track_SimpleCluster_hitds=nullptr;
   std::vector<std::vector<double>>* TPCObj_PFP_track_SimpleCluster_hitdQdsSlider=nullptr;
   std::vector<std::vector<double>>* TPCObj_PFP_track_SimpleCluster_hitLinearity=nullptr;
   std::vector<bool>* TPCObj_PFP_track_ct_passed_basic=nullptr;
   std::vector<bool>* TPCObj_PFP_track_ct_result_bragg=nullptr;
   std::vector<bool>* TPCObj_PFP_track_ct_result_michel=nullptr;


   // These are derived quantities - derived from the values above in Calcvars
   std::vector<double> *TPCObj_PFP_LH_p;
   std::vector<double> *TPCObj_PFP_LH_mu;
   std::vector<double> *TPCObj_PFP_LH_pi;
   std::vector<double> *TPCObj_PFP_LH_mip;
   std::vector<double> *TPCObj_PFP_Lmipoverp;
   std::vector<double> *TPCObj_PFP_lnLmipoverp;
   std::vector<double> *TPCObj_PFP_Lmumipovermumipp;
   std::vector<double> *TPCObj_PFP_BrokenTrackAngle;
   std::vector<double> *TPCObj_PFP_VtxTrackDist;
   std::vector<double> *TPCObj_PFP_isContained_double;
   std::vector<double> *TPCObj_PFP_track_dEdx_truncmean_start;
   std::vector<double> *TPCObj_DaughterTracks_Order_dEdxtr;
   std::vector<double> *TPCObj_DaughterTracks_Order_dEdxtr_selMIPs;
   std::vector<double> *TPCObj_DaughterTracks_Order_trklen_selMIPs;
   std::vector<double> *TPCObj_PFP_track_MCSLLmuMinusLLp;
   std::vector<double> *TPCObj_PFP_track_MomRangeMinusMCS_p;
   std::vector<double> *TPCObj_PFP_track_MomRangeMinusMCS_mu;
   std::vector<double> *TPCObj_PFP_track_SimpleCluster_hitLinearity_minimum;
   std::vector<double> *TPCObj_PFP_track_SimpleCluster_hitLinearity_mean;
   std::vector<double> *TPCObj_PFP_track_SimpleCluster_hitLinearity_truncated_mean;

   // Is the track a MIP? Evaluate MIP cuts and then put this as an input into selection
   // Coded as a double (to fit in with the other code) but should be evaluated as a bool
   std::vector<double> *TPCObj_PFP_track_passesMIPcut;
};

#endif
