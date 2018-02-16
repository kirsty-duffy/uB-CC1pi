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
#include "lardataobj/RecoBase/PFParticle.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

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

struct cc1pianavars{

  bool isData;
  art::RunNumber_t run_num;
  art::SubRunNumber_t subrun_num;
  art::EventNumber_t event_num;
  
  std::map<std::string,bool> Marco_cutflow;
  bool Marco_selected;

  NuIntTopology TPCObj_beamnu_topology;
  std::vector<double> TPCObj_PFP_track_length;
  std::vector<std::vector<double>> TPCObj_PFP_track_start;
  std::vector<std::vector<double>> TPCObj_PFP_track_end;
  std::vector<double> TPCObj_PFP_track_dqdx_truncmean;
  std::vector<bool> TPCObj_PFP_isMIP;
  std::vector<double> TPCObj_PFP_shower_length;
  std::vector<std::vector<double>> TPCObj_PFP_shower_start;
  int TPCObj_NPFPs;
  int TPCObj_NTracks;
  int TPCObj_NShowers;
  std::vector<bool> TPCObj_PFP_isTrack;
  std::vector<bool> TPCObj_PFP_isShower;
  std::vector<bool> TPCObj_PFP_isDaughter;
  std::vector<int> TPCObj_PFP_id;
  std::vector<int> TPCObj_PFP_MCPid;
  std::vector<int> TPCObj_PFP_truePDG;
  std::vector<double> TPCObj_PFP_trueE;
  std::vector<double> TPCObj_PFP_trueKE;
  int TPCObj_origin;
  int TPCObj_origin_extra;
  std::vector<double> TPCObj_reco_vtx;

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

  std::vector<double> nu_vtxx;
  std::vector<double> nu_vtxy;
  std::vector<double> nu_vtxz;
  std::vector<bool> nu_isCC;
  std::vector<int> nu_PDG;
  std::vector<double> nu_E;

  std::map<std::string,bool> CC1picutflow;

  fhicl::ParameterSet pset;

  // Constructor (pass the fhicl parameters here)
  cc1pianavars(fhicl::ParameterSet const &p);
  
  void Clear();

  void SetReco2Vars(art::Event &evt);

};

void MakeAnaBranches(TTree *t, cc1pianavars *vars);

#endif
