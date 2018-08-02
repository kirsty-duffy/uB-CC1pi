#ifndef __TREEVARS_CXX__
#define __TREEVARS_CXX__


// This code defines the variables we read in from the tree and calculates new derived variables from there.

#include "../../../../larana/larana/TruncatedMean/Algorithm/TruncMean.cxx"
#include "CC1pi_treevars.h"
#include "CC1pi_cuts.h"
#include "CC1pi_MIPcut.cxx"


void settreevars(TTree *intree, treevars *varstoset){
   intree->SetBranchStatus("*",0);
   intree->SetBranchStatus("TPCObj_PFP_LH_fwd_mu",1);
   intree->SetBranchAddress("TPCObj_PFP_LH_fwd_mu", &(varstoset->TPCObj_PFP_LH_fwd_mu));
   intree->SetBranchStatus("TPCObj_PFP_LH_fwd_p",1);
   intree->SetBranchAddress("TPCObj_PFP_LH_fwd_p", &(varstoset->TPCObj_PFP_LH_fwd_p));
   intree->SetBranchStatus("TPCObj_PFP_LH_fwd_pi",1);
   intree->SetBranchAddress("TPCObj_PFP_LH_fwd_pi", &(varstoset->TPCObj_PFP_LH_fwd_pi));
   intree->SetBranchStatus("TPCObj_PFP_LH_bwd_mu",1);
   intree->SetBranchAddress("TPCObj_PFP_LH_bwd_mu", &(varstoset->TPCObj_PFP_LH_bwd_mu));
   intree->SetBranchStatus("TPCObj_PFP_LH_bwd_p",1);
   intree->SetBranchAddress("TPCObj_PFP_LH_bwd_p", &(varstoset->TPCObj_PFP_LH_bwd_p));
   intree->SetBranchStatus("TPCObj_PFP_LH_bwd_pi",1);
   intree->SetBranchAddress("TPCObj_PFP_LH_bwd_pi", &(varstoset->TPCObj_PFP_LH_bwd_pi));
   intree->SetBranchStatus("TPCObj_PFP_LH_MIP",1);
   intree->SetBranchAddress("TPCObj_PFP_LH_MIP", &(varstoset->TPCObj_PFP_LH_MIP));
   intree->SetBranchStatus("TPCObj_PFP_track_rangeE_p",1);
   intree->SetBranchAddress("TPCObj_PFP_track_rangeE_p", &(varstoset->TPCObj_PFP_track_rangeE_p));
   intree->SetBranchStatus("TPCObj_PFP_track_rangeE_mu",1);
   intree->SetBranchAddress("TPCObj_PFP_track_rangeE_mu", &(varstoset->TPCObj_PFP_track_rangeE_mu));
   intree->SetBranchStatus("Truth_topology",1);
   intree->SetBranchAddress("Truth_topology", &(varstoset->Truth_topology));
   intree->SetBranchStatus("Marco_selected",1);
   intree->SetBranchAddress("Marco_selected", &(varstoset->Marco_selected));
   intree->SetBranchStatus("TPCObj_PFP_isDaughter",1);
   intree->SetBranchAddress("TPCObj_PFP_isDaughter", &(varstoset->TPCObj_PFP_isDaughter));
   intree->SetBranchStatus("TPCObj_PFP_track_theta",1);
   intree->SetBranchAddress("TPCObj_PFP_track_theta", &(varstoset->TPCObj_PFP_track_theta));
   intree->SetBranchStatus("TPCObj_PFP_track_phi",1);
   intree->SetBranchAddress("TPCObj_PFP_track_phi", &(varstoset->TPCObj_PFP_track_phi));
   intree->SetBranchStatus("TPCObj_PFP_track_length",1);
   intree->SetBranchAddress("TPCObj_PFP_track_length", &(varstoset->TPCObj_PFP_track_length));
   intree->SetBranchStatus("TPCObj_PFP_track_start",1);
   intree->SetBranchAddress("TPCObj_PFP_track_start", &(varstoset->TPCObj_PFP_track_start));
   intree->SetBranchStatus("TPCObj_PFP_track_end",1);
   intree->SetBranchAddress("TPCObj_PFP_track_end", &(varstoset->TPCObj_PFP_track_end));
   intree->SetBranchStatus("TPCObj_PFP_track_residual_mean",1);
   intree->SetBranchAddress("TPCObj_PFP_track_residual_mean", &(varstoset->TPCObj_PFP_track_residual_mean));
   intree->SetBranchStatus("TPCObj_PFP_track_residual_std",1);
   intree->SetBranchAddress("TPCObj_PFP_track_residual_std", &(varstoset->TPCObj_PFP_track_residual_std));
   intree->SetBranchStatus("TPCObj_PFP_track_perc_used_hits",1);
   intree->SetBranchAddress("TPCObj_PFP_track_perc_used_hits", &(varstoset->TPCObj_PFP_track_perc_used_hits));
   intree->SetBranchStatus("TPCObj_reco_vtx",1);
   intree->SetBranchAddress("TPCObj_reco_vtx", &(varstoset->TPCObj_reco_vtx));
   intree->SetBranchStatus("TPCObj_PFP_track_isContained",1);
   intree->SetBranchAddress("TPCObj_PFP_track_isContained", &(varstoset->TPCObj_PFP_track_isContained));
   intree->SetBranchStatus("TPCObj_PFP_truePDG",1);
   intree->SetBranchAddress("TPCObj_PFP_truePDG", &(varstoset->TPCObj_PFP_truePDG));
   intree->SetBranchStatus("TPCObj_PFP_track_trajPoint_Position",1);
   intree->SetBranchAddress("TPCObj_PFP_track_trajPoint_Position", &(varstoset->TPCObj_PFP_track_trajPoint_Position));
   intree->SetBranchStatus("TPCObj_PFP_track_trajPoint_Direction",1);
   intree->SetBranchAddress("TPCObj_PFP_track_trajPoint_Direction", &(varstoset->TPCObj_PFP_track_trajPoint_Direction));
   intree->SetBranchStatus("TPCObj_PFP_track_dedx_perhit",1);
   intree->SetBranchAddress("TPCObj_PFP_track_dedx_perhit", &(varstoset->TPCObj_PFP_track_dedx_perhit));
   intree->SetBranchStatus("TPCObj_PFP_track_resrange_perhit",1);
   intree->SetBranchAddress("TPCObj_PFP_track_resrange_perhit", &(varstoset->TPCObj_PFP_track_resrange_perhit));

   intree->SetBranchStatus("TPCObj_PFP_track_MCSmu_fwdMom",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSmu_fwdMom", &(varstoset->TPCObj_PFP_track_MCSmu_fwdMom));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSmu_bwdMom",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSmu_bwdMom", &(varstoset->TPCObj_PFP_track_MCSmu_bwdMom));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSmu_bestMom",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSmu_bestMom", &(varstoset->TPCObj_PFP_track_MCSmu_bestMom));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSmu_fwdMomUncert",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSmu_fwdMomUncert", &(varstoset->TPCObj_PFP_track_MCSmu_fwdMomUncert));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSmu_bwdMomUncert",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSmu_bwdMomUncert", &(varstoset->TPCObj_PFP_track_MCSmu_bwdMomUncert));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSmu_bestMomUncert",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSmu_bestMomUncert", &(varstoset->TPCObj_PFP_track_MCSmu_bestMomUncert));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSmu_fwdLL",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSmu_fwdLL", &(varstoset->TPCObj_PFP_track_MCSmu_fwdLL));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSmu_bwdLL",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSmu_bwdLL", &(varstoset->TPCObj_PFP_track_MCSmu_bwdLL));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSmu_bestLL",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSmu_bestLL", &(varstoset->TPCObj_PFP_track_MCSmu_bestLL));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSmu_segmentRadLengths",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSmu_segmentRadLengths", &(varstoset->TPCObj_PFP_track_MCSmu_segmentRadLengths));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSmu_scatterAngles",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSmu_scatterAngles", &(varstoset->TPCObj_PFP_track_MCSmu_scatterAngles));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSp_fwdMom",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSp_fwdMom", &(varstoset->TPCObj_PFP_track_MCSp_fwdMom));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSp_bwdMom",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSp_bwdMom", &(varstoset->TPCObj_PFP_track_MCSp_bwdMom));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSp_bestMom",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSp_bestMom", &(varstoset->TPCObj_PFP_track_MCSp_bestMom));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSp_fwdMomUncert",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSp_fwdMomUncert", &(varstoset->TPCObj_PFP_track_MCSp_fwdMomUncert));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSp_bwdMomUncert",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSp_bwdMomUncert", &(varstoset->TPCObj_PFP_track_MCSp_bwdMomUncert));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSp_bestMomUncert",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSp_bestMomUncert", &(varstoset->TPCObj_PFP_track_MCSp_bestMomUncert));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSp_fwdLL",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSp_fwdLL", &(varstoset->TPCObj_PFP_track_MCSp_fwdLL));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSp_bwdLL",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSp_bwdLL", &(varstoset->TPCObj_PFP_track_MCSp_bwdLL));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSp_bestLL",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSp_bestLL", &(varstoset->TPCObj_PFP_track_MCSp_bestLL));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSp_segmentRadLengths",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSp_segmentRadLengths", &(varstoset->TPCObj_PFP_track_MCSp_segmentRadLengths));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSp_scatterAngles",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSp_scatterAngles", &(varstoset->TPCObj_PFP_track_MCSp_scatterAngles));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSpi_fwdMom",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSpi_fwdMom", &(varstoset->TPCObj_PFP_track_MCSpi_fwdMom));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSpi_bwdMom",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSpi_bwdMom", &(varstoset->TPCObj_PFP_track_MCSpi_bwdMom));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSpi_bestMom",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSpi_bestMom", &(varstoset->TPCObj_PFP_track_MCSpi_bestMom));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSpi_fwdMomUncert",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSpi_fwdMomUncert", &(varstoset->TPCObj_PFP_track_MCSpi_fwdMomUncert));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSpi_bwdMomUncert",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSpi_bwdMomUncert", &(varstoset->TPCObj_PFP_track_MCSpi_bwdMomUncert));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSpi_bestMomUncert",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSpi_bestMomUncert", &(varstoset->TPCObj_PFP_track_MCSpi_bestMomUncert));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSpi_fwdLL",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSpi_fwdLL", &(varstoset->TPCObj_PFP_track_MCSpi_fwdLL));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSpi_bwdLL",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSpi_bwdLL", &(varstoset->TPCObj_PFP_track_MCSpi_bwdLL));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSpi_bestLL",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSpi_bestLL", &(varstoset->TPCObj_PFP_track_MCSpi_bestLL));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSpi_segmentRadLengths",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSpi_segmentRadLengths", &(varstoset->TPCObj_PFP_track_MCSpi_segmentRadLengths));
   intree->SetBranchStatus("TPCObj_PFP_track_MCSpi_scatterAngles",1);
   intree->SetBranchAddress("TPCObj_PFP_track_MCSpi_scatterAngles", &(varstoset->TPCObj_PFP_track_MCSpi_scatterAngles));

   intree->SetBranchStatus("TPCObj_PFP_track_SimpleCluster_hitTime",1);
   intree->SetBranchAddress("TPCObj_PFP_track_SimpleCluster_hitTime", &(varstoset->TPCObj_PFP_track_SimpleCluster_hitTime));
   intree->SetBranchStatus("TPCObj_PFP_track_SimpleCluster_hitWire",1);
   intree->SetBranchAddress("TPCObj_PFP_track_SimpleCluster_hitWire", &(varstoset->TPCObj_PFP_track_SimpleCluster_hitWire));
   intree->SetBranchStatus("TPCObj_PFP_track_SimpleCluster_hitPlane",1);
   intree->SetBranchAddress("TPCObj_PFP_track_SimpleCluster_hitPlane", &(varstoset->TPCObj_PFP_track_SimpleCluster_hitPlane));
   intree->SetBranchStatus("TPCObj_PFP_track_SimpleCluster_hitIntegral",1);
   intree->SetBranchAddress("TPCObj_PFP_track_SimpleCluster_hitIntegral", &(varstoset->TPCObj_PFP_track_SimpleCluster_hitIntegral));
   intree->SetBranchStatus("TPCObj_PFP_track_SimpleCluster_hitTimeTicks",1);
   intree->SetBranchAddress("TPCObj_PFP_track_SimpleCluster_hitTimeTicks", &(varstoset->TPCObj_PFP_track_SimpleCluster_hitTimeTicks));
   intree->SetBranchStatus("TPCObj_PFP_track_SimpleCluster_hitWireNo",1);
   intree->SetBranchAddress("TPCObj_PFP_track_SimpleCluster_hitWireNo", &(varstoset->TPCObj_PFP_track_SimpleCluster_hitWireNo));
   intree->SetBranchStatus("TPCObj_PFP_track_SimpleCluster_StartIndex",1);
   intree->SetBranchAddress("TPCObj_PFP_track_SimpleCluster_StartIndex", &(varstoset->TPCObj_PFP_track_SimpleCluster_StartIndex));
   intree->SetBranchStatus("TPCObj_PFP_track_SimpleCluster_hitdQds",1);
   intree->SetBranchAddress("TPCObj_PFP_track_SimpleCluster_hitdQds", &(varstoset->TPCObj_PFP_track_SimpleCluster_hitdQds));
   intree->SetBranchStatus("TPCObj_PFP_track_SimpleCluster_hitds",1);
   intree->SetBranchAddress("TPCObj_PFP_track_SimpleCluster_hitds", &(varstoset->TPCObj_PFP_track_SimpleCluster_hitds));
   intree->SetBranchStatus("TPCObj_PFP_track_SimpleCluster_hitdQdsSlider",1);
   intree->SetBranchAddress("TPCObj_PFP_track_SimpleCluster_hitdQdsSlider", &(varstoset->TPCObj_PFP_track_SimpleCluster_hitdQdsSlider));
   intree->SetBranchStatus("TPCObj_PFP_track_SimpleCluster_hitLinearity",1);
   intree->SetBranchAddress("TPCObj_PFP_track_SimpleCluster_hitLinearity", &(varstoset->TPCObj_PFP_track_SimpleCluster_hitLinearity));
   intree->SetBranchStatus("TPCObj_PFP_track_ct_passed_basic",1);
   intree->SetBranchAddress("TPCObj_PFP_track_ct_passed_basic", &(varstoset->TPCObj_PFP_track_ct_passed_basic));
   intree->SetBranchStatus("TPCObj_PFP_track_ct_result_bragg",1);
   intree->SetBranchAddress("TPCObj_PFP_track_ct_result_bragg", &(varstoset->TPCObj_PFP_track_ct_result_bragg));
   intree->SetBranchStatus("TPCObj_PFP_track_ct_result_michel",1);
   intree->SetBranchAddress("TPCObj_PFP_track_ct_result_michel", &(varstoset->TPCObj_PFP_track_ct_result_michel));

}

void Calcvars(treevars *vars){

   // Initialise vectors that we are going to fill with calculated values
   int vecsize = vars->TPCObj_PFP_LH_fwd_p->size();

   vars->TPCObj_PFP_LH_p = new std::vector<double>(vecsize);
   vars->TPCObj_PFP_LH_mu = new std::vector<double>(vecsize);
   vars->TPCObj_PFP_LH_pi = new std::vector<double>(vecsize);
   vars->TPCObj_PFP_LH_mip = new std::vector<double>(vecsize);
   vars->TPCObj_PFP_Lmipoverp = new std::vector<double>(vecsize);
   vars->TPCObj_PFP_lnLmipoverp = new std::vector<double>(vecsize);
   vars->TPCObj_PFP_Lmumipovermumipp = new std::vector<double>(vecsize);
   vars->TPCObj_PFP_BrokenTrackAngle = new std::vector<double>(vecsize);
   vars->TPCObj_PFP_VtxTrackDist = new std::vector<double>(vecsize);
   vars->TPCObj_PFP_isContained_double = new std::vector<double>(vecsize);
   vars->TPCObj_PFP_track_dEdx_truncmean_start = new std::vector<double>(vecsize);
   vars->TPCObj_DaughterTracks_Order_dEdxtr = new std::vector<double>(vecsize,-9999);
   vars->TPCObj_DaughterTracks_Order_dEdxtr_selMIPs = new std::vector<double>(vecsize,-9999);
   vars->TPCObj_DaughterTracks_Order_trklen_selMIPs = new std::vector<double>(vecsize,-9999);
   vars->TPCObj_PFP_track_MCSLLmuMinusLLp = new std::vector<double>(vecsize,-9999);
   vars->TPCObj_PFP_track_MomRangeMinusMCS_p = new std::vector<double>(vecsize,-9999);
   vars->TPCObj_PFP_track_MomRangeMinusMCS_mu = new std::vector<double>(vecsize,-9999);

   vars->TPCObj_PFP_track_passesMIPcut = new std::vector<double>(vecsize,-9999.);

   // Just use collection plane for now
   int i_pl = 2;

   // For ordering tracks by various things
   std::vector<std::pair<double, int>> *pair_dEdx_truncmean_index = new std::vector<std::pair<double, int>>(vecsize);
   std::vector<std::pair<double, int>> *pair_trklen_index = new std::vector<std::pair<double, int>>(vecsize);

   // Now calculate the values for all variables
   for (int i_track=0; i_track < vecsize; i_track++){

      // The neutrino and any other PFPs that didn't get reco'd as tracks should get bogus values
      if(vars->TPCObj_PFP_track_theta->at(i_track) == -9999) {

         vars->TPCObj_PFP_LH_p->at(i_track) = -9999;
         vars->TPCObj_PFP_LH_mu->at(i_track) = -9999;
         vars->TPCObj_PFP_LH_pi->at(i_track) = -9999;
         vars->TPCObj_PFP_LH_mip->at(i_track) = -9999;
         vars->TPCObj_PFP_Lmipoverp->at(i_track) = -9999;
         vars->TPCObj_PFP_lnLmipoverp->at(i_track) = -9999;
         vars->TPCObj_PFP_Lmumipovermumipp->at(i_track) = -9999;
         vars->TPCObj_PFP_VtxTrackDist->at(i_track) = -9999;
         vars->TPCObj_PFP_BrokenTrackAngle->at(i_track) = -9999;
         vars->TPCObj_PFP_isContained_double->at(i_track) = -9999;
         vars->TPCObj_PFP_track_dEdx_truncmean_start->at(i_track) = -9999;
         vars->TPCObj_PFP_track_dEdx_truncmean_start->at(i_track) = -9999;
         pair_dEdx_truncmean_index->at(i_track)=std::make_pair(-9999.,i_track);
         pair_trklen_index->at(i_track)=std::make_pair(-9999.,i_track);
         vars->TPCObj_PFP_track_MCSLLmuMinusLLp->at(i_track) = -9999.;
         vars->TPCObj_PFP_track_MomRangeMinusMCS_p->at(i_track) = -9999.;
         vars->TPCObj_PFP_track_MomRangeMinusMCS_mu->at(i_track) = -9999.;
         vars->TPCObj_PFP_track_passesMIPcut->at(i_track) = -9999.;

         continue;
      }

      vars->TPCObj_PFP_LH_p->at(i_track) = std::max(vars->TPCObj_PFP_LH_fwd_p->at(i_track).at(i_pl), vars->TPCObj_PFP_LH_bwd_p->at(i_track).at(i_pl));
      vars->TPCObj_PFP_LH_mu->at(i_track) = std::max(vars->TPCObj_PFP_LH_fwd_mu->at(i_track).at(i_pl), vars->TPCObj_PFP_LH_bwd_mu->at(i_track).at(i_pl));
      vars->TPCObj_PFP_LH_pi->at(i_track) = std::max(vars->TPCObj_PFP_LH_fwd_pi->at(i_track).at(i_pl), vars->TPCObj_PFP_LH_bwd_pi->at(i_track).at(i_pl));
      vars->TPCObj_PFP_LH_mip->at(i_track) = vars->TPCObj_PFP_LH_MIP->at(i_track).at(i_pl);

      // Divisions are only valid if none of the values are <0 (i.e. -999 or -9999: only set for invalid tracks)
      // Otherwise set to -9999
      vars->TPCObj_PFP_Lmipoverp->at(i_track)=-9999;
      vars->TPCObj_PFP_lnLmipoverp->at(i_track)=-9999;
      vars->TPCObj_PFP_Lmumipovermumipp->at(i_track)=-9999;
      if (vars->TPCObj_PFP_LH_mip->at(i_track)>=0 && vars->TPCObj_PFP_LH_p->at(i_track)>=0){
         vars->TPCObj_PFP_Lmipoverp->at(i_track) = vars->TPCObj_PFP_LH_mip->at(i_track) / vars->TPCObj_PFP_LH_p->at(i_track);
         vars->TPCObj_PFP_lnLmipoverp->at(i_track) = TMath::Log(vars->TPCObj_PFP_Lmipoverp->at(i_track));
         if (vars->TPCObj_PFP_LH_mu->at(i_track)>=0){
            vars->TPCObj_PFP_Lmumipovermumipp->at(i_track) = (vars->TPCObj_PFP_LH_mu->at(i_track)+vars->TPCObj_PFP_LH_mip->at(i_track))/(vars->TPCObj_PFP_LH_mu->at(i_track)+vars->TPCObj_PFP_LH_mip->at(i_track)+vars->TPCObj_PFP_LH_p->at(i_track));
         }
      }
      vars->TPCObj_PFP_VtxTrackDist->at(i_track) = std::hypot(std::hypot(vars->TPCObj_reco_vtx->at(0) - vars->TPCObj_PFP_track_start->at(i_track).at(0), vars->TPCObj_reco_vtx->at(1) - vars->TPCObj_PFP_track_start->at(i_track).at(1)), vars->TPCObj_reco_vtx->at(2) - vars->TPCObj_PFP_track_start->at(i_track).at(2));
      vars->TPCObj_PFP_isContained_double->at(i_track) = (double)(vars->TPCObj_PFP_track_isContained->at(i_track));

      // Angle between tracks (using theta and phi)...
      int track1 = i_track;

      // We want non-tracks to get bogus values (-9999)
      // But we want tracks that just don't have other tracks close by to get 0
      double maxangle = -1;

      for (int track2 = track1+1; track2 < vecsize; track2++) {
         if(vars->TPCObj_PFP_track_theta->at(track1) == -9999) continue;
         if(vars->TPCObj_PFP_track_theta->at(track2) == -9999) continue;

         TVector3 v1(1,1,1);
         TVector3 v2(1,1,1);
         v1.SetMag(1);
         v2.SetMag(1);
         v1.SetTheta(vars->TPCObj_PFP_track_theta->at(track1));
         v2.SetTheta(vars->TPCObj_PFP_track_theta->at(track2));
         v1.SetPhi(vars->TPCObj_PFP_track_phi->at(track1));
         v2.SetPhi(vars->TPCObj_PFP_track_phi->at(track2));

         TVector3 start1(vars->TPCObj_PFP_track_start->at(track1).at(0),vars->TPCObj_PFP_track_start->at(track1).at(1),vars->TPCObj_PFP_track_start->at(track1).at(2));
         TVector3 end1(vars->TPCObj_PFP_track_end->at(track1).at(0),vars->TPCObj_PFP_track_end->at(track1).at(1),vars->TPCObj_PFP_track_end->at(track1).at(2));
         TVector3 start2(vars->TPCObj_PFP_track_start->at(track2).at(0),vars->TPCObj_PFP_track_start->at(track2).at(1),vars->TPCObj_PFP_track_start->at(track2).at(2));
         TVector3 end2(vars->TPCObj_PFP_track_end->at(track2).at(0),vars->TPCObj_PFP_track_end->at(track2).at(1),vars->TPCObj_PFP_track_end->at(track2).at(2));
         std::vector<double> distances = { (start1-start2).Mag(),(start1-end2).Mag(), (end1-start2).Mag(), (end1-end2).Mag() };
         double mindist = *std::min_element(distances.begin(),distances.end());

         // Flip tracks if the closest point isn't start-start
         // Not needed: cos is symmetric
         // if(mindist==distances.at(1) || mindist==distances.at(3)) {
         //    v2=-1*v2;
         //    TVector3 oldstart2=start2;
         //    start2=end2;
         //    end2=oldstart2;
         // }
         // if(mindist==distances.at(2) || mindist==distances.at(3)) {
         //    v1=-1*v1;
         //    TVector3 oldstart1=start1;
         //    start1=end1;
         //    end1=oldstart1;
         // }

         double angle = TMath::ACos(v1.Dot(v2));

         if(mindist < 3 && angle > maxangle) maxangle = angle;
      }

      vars->TPCObj_PFP_BrokenTrackAngle->at(i_track) = maxangle;

      // Calculate truncated mean dE/dx at the start of the track
      int nhits_start = 15;
      std::vector<float> dEdx_float;
      //std::cout << "---" << std::endl;
      // Vector order is track/plane/hit
      for (int i=1; i<nhits_start+1; i++){
          // Skip first hit (start from i=1) and last hit
          size_t perhit_size = vars->TPCObj_PFP_track_dedx_perhit->at(i_track).at(i_pl).size();
          // std::cout << "perhit_size = " << perhit_size << std::endl;
          if (i>=perhit_size) continue;
          int index = i;
          if (vars->TPCObj_PFP_track_resrange_perhit->at(i_track).at(i_pl).at(0)<vars->TPCObj_PFP_track_resrange_perhit->at(i_track).at(i_pl).at(perhit_size-1)){ // start of vector is end of track
           index = perhit_size-i;
           // std::cout << "Pushing back residual range " << vars->TPCObj_PFP_track_resrange_perhit->at(i_track).at(i_pl).at(index) << " instead of " << vars->TPCObj_PFP_track_resrange_perhit->at(i_track).at(i_pl).at(i) << std::endl;
          }
          else{
           // std::cout << "Pushing back residual range " << vars->TPCObj_PFP_track_resrange_perhit->at(i_track).at(i_pl).at(index) << " instead of " << vars->TPCObj_PFP_track_resrange_perhit->at(i_track).at(i_pl).at(perhit_size-i) << std::endl;
          }
          dEdx_float.push_back((float)(vars->TPCObj_PFP_track_dedx_perhit->at(i_track).at(i_pl).at(index)));
      }
      TruncMean trm;
      double tmpval = -9999;
      if (dEdx_float.size()>0){
         tmpval = (double)trm.CalcIterativeTruncMean(dEdx_float, 1, 1, 0, 1, 0.1, 1.0);
      } vars->TPCObj_PFP_track_dEdx_truncmean_start->at(i_track) = tmpval;
      pair_dEdx_truncmean_index->at(i_track) = std::make_pair(tmpval,i_track);

      // Make pair for ordering tracks by Length
      pair_trklen_index->at(i_track) = std::make_pair(vars->TPCObj_PFP_track_length->at(i_track),i_track);


      // MCS-related variables: MCS LLmu-LLp, momentum by range - momentum by MCS for a muon, momentum by range - momentum by MCS for a proton
      if (vars->TPCObj_PFP_track_MCSmu_fwdLL->at(i_track)==-9999 || vars->TPCObj_PFP_track_MCSmu_fwdLL->at(i_track)==-999 || vars->TPCObj_PFP_track_MCSp_fwdLL->at(i_track)==-9999 || vars->TPCObj_PFP_track_MCSp_fwdLL->at(i_track)==-999){
         vars->TPCObj_PFP_track_MCSLLmuMinusLLp->at(i_track) = -9999.;
      }
      else{
         vars->TPCObj_PFP_track_MCSLLmuMinusLLp->at(i_track) = vars->TPCObj_PFP_track_MCSmu_fwdLL->at(i_track)-vars->TPCObj_PFP_track_MCSp_fwdLL->at(i_track);
      }

      if (vars->TPCObj_PFP_track_rangeE_p->at(i_track)==-9999 || vars->TPCObj_PFP_track_rangeE_p->at(i_track)==-999 || vars->TPCObj_PFP_track_MCSp_fwdMom->at(i_track)==-9999 || vars->TPCObj_PFP_track_MCSp_fwdMom->at(i_track)==-999){
         vars->TPCObj_PFP_track_MomRangeMinusMCS_p->at(i_track) = -9999;
      }
      else{
         double KE_p = vars->TPCObj_PFP_track_rangeE_p->at(i_track);
         double M_p = 938.272;
         double MomRange_p = TMath::Sqrt((KE_p*KE_p)+(2*M_p*KE_p));
         vars->TPCObj_PFP_track_MomRangeMinusMCS_p->at(i_track) = MomRange_p/1000.-vars->TPCObj_PFP_track_MCSp_fwdMom->at(i_track);
      }

      if (vars->TPCObj_PFP_track_rangeE_mu->at(i_track)==-9999 || vars->TPCObj_PFP_track_rangeE_mu->at(i_track)==-999 || vars->TPCObj_PFP_track_MCSmu_fwdMom->at(i_track)==-9999 || vars->TPCObj_PFP_track_MCSmu_fwdMom->at(i_track)==-999){
         vars->TPCObj_PFP_track_MomRangeMinusMCS_mu->at(i_track) = -9999;
      }
      else{
         double KE_mu = vars->TPCObj_PFP_track_rangeE_mu->at(i_track);
         double M_mu = 105.7;
         double MomRange_mu =  TMath::Sqrt((KE_mu*KE_mu)+(2*M_mu*KE_mu));
         vars->TPCObj_PFP_track_MomRangeMinusMCS_mu->at(i_track) = MomRange_mu/1000.-vars->TPCObj_PFP_track_MCSmu_fwdMom->at(i_track);
      }


      // Evaluate MIP cut (i.e. whether we want to class this track as a MIP). Cut algorithm defined in CutValues_header.h
      vars->TPCObj_PFP_track_passesMIPcut->at(i_track) = (double)EvalMIPCut(vars,i_track);

   } // end loop over tracks in TPCObj (i_track)

   // Now fill TPCObj_DaughterTracks_Order_dEdxtr
   std::sort(pair_dEdx_truncmean_index->begin(),pair_dEdx_truncmean_index->end());
   int i_order=0;
   int i_order_mip=0;
   for (int i_vec=0; i_vec<vecsize; i_vec++){
      int i_track = pair_dEdx_truncmean_index->at(i_vec).second;
      if (vars->TPCObj_PFP_isDaughter->at(i_track) && pair_dEdx_truncmean_index->at(i_vec).first >= 0){
         vars->TPCObj_DaughterTracks_Order_dEdxtr->at(i_track) = i_order;
         if (bool(vars->TPCObj_PFP_track_passesMIPcut->at(i_track))){
            vars->TPCObj_DaughterTracks_Order_dEdxtr_selMIPs->at(i_track) = i_order_mip;
            i_order_mip++;
            //std::cout << " -- Track " << i_track << " is MIP: " << pair_dEdx_truncmean_index->at(i_vec).first << std::endl;
         }
         i_order++;
      }
   }

   // Now fill TPCObj_DaughterTracks_Order_trklen_selMIPs
   std::sort(pair_trklen_index->begin(),pair_trklen_index->end());
   i_order_mip=0;
   for (int i_vec=0; i_vec<vecsize; i_vec++){
      int i_track = pair_trklen_index->at(i_vec).second;
      if (bool(vars->TPCObj_PFP_track_passesMIPcut->at(i_track)) && vars->TPCObj_PFP_isDaughter->at(i_track) && pair_trklen_index->at(i_vec).first >= 0){
         vars->TPCObj_DaughterTracks_Order_trklen_selMIPs->at(i_track) = i_order_mip;
         i_order_mip++;
      }
   }
}


// --------------------------------------------------- //
// Function to clear memory from calculated variables


void Clearvars(treevars *vars){

   delete vars->TPCObj_PFP_LH_p;
   delete vars->TPCObj_PFP_LH_mu;
   delete vars->TPCObj_PFP_LH_pi;
   delete vars->TPCObj_PFP_LH_mip;
   delete vars->TPCObj_PFP_Lmipoverp;
   delete vars->TPCObj_PFP_lnLmipoverp;
   delete vars->TPCObj_PFP_Lmumipovermumipp;
   delete vars->TPCObj_PFP_BrokenTrackAngle;
   delete vars->TPCObj_PFP_VtxTrackDist;
   delete vars->TPCObj_PFP_isContained_double;
   delete vars->TPCObj_PFP_track_dEdx_truncmean_start;
   delete vars->TPCObj_DaughterTracks_Order_dEdxtr;
   delete vars->TPCObj_DaughterTracks_Order_dEdxtr_selMIPs;

}

#endif
