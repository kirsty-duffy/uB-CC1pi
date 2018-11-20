// Struct to hold tree variables for CC1pi analysis

// Note: if you want to add a variable to the tree, you have to do three things:
//   1) add the data object in cc1pianavars
//   2) add a default value for the data object in Clear() (inside cc1pianavars)
//   3) add a line to set the branch in MakeAnaBranches
//   4) add a line to fill the branch in SetReco2Vars *if* it's a variable that's pretty much copied out the existing tree (i.e. we don't have to calculate it)

// Only define this stuff once
#ifndef FILLTREE_H
#define FILLTREE_H

// uboone includes
#include "uboone/AnalysisTree/MCTruth/AssociationsTruth_tool.h"
#include "uboone/AnalysisTree/MCTruth/IMCTruthMatching.h"

// root includes
#include "TFile.h"
#include "TTree.h"
#include "TDirectory.h"

// Marco's code includes
#include "uboone/UBXSec/DataTypes/SelectionResult.h"
#include "uboone/UBXSec/DataTypes/TPCObject.h"
#include "uboone/UBXSec/DataTypes/MCGhost.h"

// larsoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/AnalysisBase/ParticleID.h"
#include "lardataobj/RecoBase/MCSFitResult.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "larreco/RecoAlg/TrackMomentumCalculator.h"

// art includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

// CC1pi method includes
#include "uboone/CC1pi/Algorithms/MIPConsistencyCheck_Marco.h"
#include "uboone/CC1pi/Algorithms/GetTopology.h"
#include "uboone/CC1pi/Algorithms/TopologyEnums.h"
#include "uboone/CC1pi/Algorithms/FVCheck.h"
#include "uboone/CC1pi/Algorithms/PDGEnums.h"
#include "uboone/CC1pi/Algorithms/InputTags.h"
#include "uboone/CC1pi/Algorithms/ShowerRejection.h"
#include "uboone/CC1pi/Algorithms/MapBuilderUtility.h"

// Particle ID includes
#include "uboone/ParticleID/Algorithms/Bragg_Likelihood_Estimator.h"
#include "uboone/ParticleID/Algorithms/uB_PlaneIDBitsetHelperFunctions.h"

struct cc1pianavars{

  // Event variables
  bool isData;
  art::RunNumber_t run_num;
  art::SubRunNumber_t subrun_num;
  art::EventNumber_t event_num;

  NuIntTopology Truth_topology;

  // For likelihood-based PID (specifically for Lmip per hit)
  particleid::Bragg_Likelihood_Estimator braggcalc;

  // Selection result variables
  std::map<std::string,bool> Marco_cutflow;
  bool Marco_selected;
  std::map<std::string,bool> CC1picutflow;

  // Track variables
  std::vector<double> TPCObj_PFP_track_length;
  std::vector<std::vector<double>> TPCObj_PFP_track_start;
  std::vector<std::vector<double>> TPCObj_PFP_track_end;
  std::vector<double> TPCObj_PFP_track_theta;
  std::vector<double> TPCObj_PFP_track_phi;
  std::vector<double> TPCObj_PFP_track_mom;
  std::vector<std::vector<double>> TPCObj_PFP_track_dedx_truncmean;
  std::vector<std::vector<std::vector<double>>> TPCObj_PFP_track_dedx_perhit;
  std::vector<std::vector<std::vector<double>>> TPCObj_PFP_track_resrange_perhit;
  std::vector<double> TPCObj_PFP_track_nhits;
  std::vector<bool> TPCObj_PFP_isMIP;
  std::vector<bool> TPCObj_PFP_track_isContained;
  std::vector<std::vector<double>> TPCObj_PFP_track_AngleBetweenTracks;
  std::vector<double> TPCObj_PFP_track_residual_mean;
  std::vector<double> TPCObj_PFP_track_residual_std;
  std::vector<double> TPCObj_PFP_track_perc_used_hits;
  std::vector<double> TPCObj_PFP_track_MCSmu_fwdMom;
  std::vector<double> TPCObj_PFP_track_MCSmu_bwdMom;
  std::vector<double> TPCObj_PFP_track_MCSmu_bestMom;
  std::vector<double> TPCObj_PFP_track_MCSmu_fwdMomUncert;
  std::vector<double> TPCObj_PFP_track_MCSmu_bwdMomUncert;
  std::vector<double> TPCObj_PFP_track_MCSmu_bestMomUncert;
  std::vector<double> TPCObj_PFP_track_MCSmu_fwdLL;
  std::vector<double> TPCObj_PFP_track_MCSmu_bwdLL;
  std::vector<double> TPCObj_PFP_track_MCSmu_bestLL;
  std::vector<std::vector<double>> TPCObj_PFP_track_MCSmu_segmentRadLengths;
  std::vector<std::vector<double>> TPCObj_PFP_track_MCSmu_scatterAngles;
  std::vector<double> TPCObj_PFP_track_MCSp_fwdMom;
  std::vector<double> TPCObj_PFP_track_MCSp_bwdMom;
  std::vector<double> TPCObj_PFP_track_MCSp_bestMom;
  std::vector<double> TPCObj_PFP_track_MCSp_fwdMomUncert;
  std::vector<double> TPCObj_PFP_track_MCSp_bwdMomUncert;
  std::vector<double> TPCObj_PFP_track_MCSp_bestMomUncert;
  std::vector<double> TPCObj_PFP_track_MCSp_fwdLL;
  std::vector<double> TPCObj_PFP_track_MCSp_bwdLL;
  std::vector<double> TPCObj_PFP_track_MCSp_bestLL;
  std::vector<std::vector<double>> TPCObj_PFP_track_MCSp_segmentRadLengths;
  std::vector<std::vector<double>> TPCObj_PFP_track_MCSp_scatterAngles;
  std::vector<double> TPCObj_PFP_track_MCSpi_fwdMom;
  std::vector<double> TPCObj_PFP_track_MCSpi_bwdMom;
  std::vector<double> TPCObj_PFP_track_MCSpi_bestMom;
  std::vector<double> TPCObj_PFP_track_MCSpi_fwdMomUncert;
  std::vector<double> TPCObj_PFP_track_MCSpi_bwdMomUncert;
  std::vector<double> TPCObj_PFP_track_MCSpi_bestMomUncert;
  std::vector<double> TPCObj_PFP_track_MCSpi_fwdLL;
  std::vector<double> TPCObj_PFP_track_MCSpi_bwdLL;
  std::vector<double> TPCObj_PFP_track_MCSpi_bestLL;
  std::vector<std::vector<double>> TPCObj_PFP_track_MCSpi_segmentRadLengths;
  std::vector<std::vector<double>> TPCObj_PFP_track_MCSpi_scatterAngles;
  std::vector<std::vector<std::vector<double>>> TPCObj_PFP_track_SpacepointsXYZ;
  std::vector<std::vector<double>> TPCObj_PFP_track_SpacepointsQPlane2;

  // Shower variables
  std::vector<double> TPCObj_PFP_shower_length;
  std::vector<std::vector<double>> TPCObj_PFP_shower_start;

  // PFP reco variables
  int TPCObj_NPFPs;
  int TPCObj_NTracks;
  int TPCObj_NShowers;
  std::vector<bool> TPCObj_PFP_PandoraClassedAsTrack;
  std::vector<bool> TPCObj_PFP_PandoraClassedAsShower;
  std::vector<bool> TPCObj_PFP_isDaughter;
  std::vector<std::vector<int>> TPCObj_PFP_daughterids;
  std::vector<int> TPCObj_PFP_id;
  std::vector<double> TPCObj_reco_vtx;

  // PFP true variables
  std::vector<int> TPCObj_PFP_MCPid;
  std::vector<int> TPCObj_PFP_truePDG;
  std::vector<double> TPCObj_PFP_trueE;
  std::vector<double> TPCObj_PFP_trueKE;
  std::vector<double> TPCObj_PFP_trueEndP;
  int TPCObj_origin;
  int TPCObj_origin_extra;

  // PID variables
  std::vector<std::vector<double>> TPCObj_PFP_LH_fwd_mu;
  std::vector<std::vector<double>> TPCObj_PFP_LH_fwd_p;
  std::vector<std::vector<double>> TPCObj_PFP_LH_fwd_pi;
  std::vector<std::vector<double>> TPCObj_PFP_LH_bwd_mu;
  std::vector<std::vector<double>> TPCObj_PFP_LH_bwd_p;
  std::vector<std::vector<double>> TPCObj_PFP_LH_bwd_pi;
  std::vector<std::vector<double>> TPCObj_PFP_LH_MIP;
  std::vector<std::vector<double>> TPCObj_PFP_PIDA;
  std::vector<std::vector<double>> TPCObj_PFP_track_depE;
  std::vector<std::vector<double>> TPCObj_PFP_track_Chi2Proton;
  std::vector<std::vector<double>> TPCObj_PFP_track_Chi2Muon;
  std::vector<std::vector<double>> TPCObj_PFP_track_Chi2Pion;
  std::vector<double> TPCObj_PFP_track_rangeE_mu;
  std::vector<double> TPCObj_PFP_track_rangeE_p;
  std::vector<std::vector<std::vector<double>>> TPCObj_PFP_track_Lmip_perhit;

  // MCParticle variables
  std::vector<int> MCP_PDG;
  std::vector<double> MCP_length;
  std::vector<std::string> MCP_process;
  std::vector<std::string> MCP_endprocess;
  std::vector<int> MCP_numdaughters;
  std::vector<double> MCP_P;
  std::vector<double> MCP_Px;
  std::vector<double> MCP_Py;
  std::vector<double> MCP_Pz;
  std::vector<double> MCP_E;
  std::vector<double> MCP_KE;
  std::vector<bool> MCP_isContained;
  std::vector<int> MCP_ID;
  std::vector<std::vector<int>> MCP_DaughterIDs;
  std::vector<int> MCP_StatusCode;
  std::vector<int> MCP_MotherID;
  std::vector<bool> MCP_isNuDaughter;
  std::vector<std::vector<double>> MCP_StartPosition;
  std::vector<std::vector<double>> MCP_EndPosition;

  // True neutrino variables
  std::vector<double> nu_vtx;
  std::vector<double> nu_vtx_spacecharge;
  bool nu_isCC;
  int nu_PDG;
  double nu_E;
  int nu_MCPID;


  fhicl::ParameterSet pset;

  InputTags *CC1piInputTags;

  // Constructor (pass the fhicl parameters here)
  cc1pianavars(fhicl::ParameterSet const &p);

  void Clear();

  void SetReco2Vars(art::Event &evt);

};

void MakeAnaBranches(TTree *t, cc1pianavars *vars);

#endif
