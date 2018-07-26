#ifndef __TREEVARS_H__
#define __TREEVARS_H__


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

   // Is the track a MIP? Evaluate MIP cuts and then put this as an input into selection
   // Coded as a double (to fit in with the other code) but should be evaluated as a bool
   std::vector<double> *TPCObj_PFP_track_passesMIPcut;
};

#endif
