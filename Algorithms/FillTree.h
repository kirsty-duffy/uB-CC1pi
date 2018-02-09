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

struct cc1pianavars{

  bool isData;
  art::RunNumber_t run_num;
  art::SubRunNumber_t subrun_num;
  art::EventNumber_t event_num;
  
  NuIntTopology topology;
  
  std::map<std::string,bool> cutflow;
  bool isSelected;
  std::vector<double> track_length;
  std::vector<std::vector<double>> track_start;
  std::vector<std::vector<double>> track_end;
  std::vector<double> shower_length;
  int NPFPs;
  int NTracks;
  int NShowers;
  std::vector<bool> Sel_PFP_isTrack;
  std::vector<bool> Sel_PFP_isShower;
  std::vector<bool> Sel_PFP_isDaughter;
  std::vector<bool> Sel_PFP_isMIP;
  std::vector<int> Sel_PFP_ID;
  std::vector<int> Sel_MCP_ID;
  std::vector<int> Sel_MCP_PDG;
  std::vector<double> Sel_MCP_E;
  int tpcobj_origin;
  int tpcobj_origin_extra;

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
  std::vector<bool> MCP_isContained;

  std::vector<double> nu_vtxx;
  std::vector<double> nu_vtxy;
  std::vector<double> nu_vtxz;
  std::vector<bool> nu_isCC;
  std::vector<int> nu_PDG;
  std::vector<double> nu_E;
  
  std::vector<double> dqdx_trunc_uncalib;

  std::map<std::string,bool> CC1picutflow;
  bool PassesCC1piSelec;
  std::string CC1piSelecFailureReason;

  fhicl::ParameterSet pset;

  // Constructor (pass the fhicl parameters here)
  cc1pianavars(fhicl::ParameterSet const &p);
  
  void Clear();

  void SetReco2Vars(art::Event &evt);

};

bool inFV(double x, double y, double z);

void MakeAnaBranches(TTree *t, cc1pianavars *vars);

#endif
