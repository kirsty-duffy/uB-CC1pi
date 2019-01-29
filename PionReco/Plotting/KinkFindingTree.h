#include <map>
#include <vector>
#include <iostream>
#include <algorithm>

#include "TTree.h"

#include "CalcLocalLinearity.h"
#include "MakeLocalLinearityPlots.C"

// Function to make list of primary PFPs
// This is a nasty workaround because it turns out primary doesn't mean what I thought (it means "neutrino", whereas I wanted the first non-neutrino PFPs)
std::vector<int> GetAllPrimaryPFPIDs(std::vector<int> *PFP_primaryPFPid, std::vector<bool> *PFP_isInNuSlice, std::map<int,int> PFP_ID_to_vectorpos){
    std::vector<int> all_primary_PFP_IDs;

    for (size_t i=0; i<PFP_primaryPFPid->size(); i++){
        int primaryPFPID = PFP_primaryPFPid->at(i);
        // Check if we've seen this one before
        // If not, add it to the vector
        if (std::find(all_primary_PFP_IDs.begin(), all_primary_PFP_IDs.end(), primaryPFPID) != all_primary_PFP_IDs.end()){
            continue;
        }

        // Check that primary comes from neutrino slice, if it exists (i.e. if we're looking at Pandora consolidated output)
        // If not, don't keep this primary
        int primary_idx = PFP_ID_to_vectorpos.find(primaryPFPID)->second;
        if (PFP_isInNuSlice && (!PFP_isInNuSlice->at(primary_idx))) continue;

        all_primary_PFP_IDs.push_back(primaryPFPID);
    }

    return all_primary_PFP_IDs;
}

struct KinkFindingTree{
  // Variables to fill from tree
  int nMCPs=0;
  int nu_MCPID=0;
  std::vector<int> *MCP_ID=nullptr;
  std::vector<int> *MCP_MotherID=nullptr;
  std::vector<std::vector<int>> *MCP_DaughterIDs=nullptr;
  std::vector<int> *MCP_PDGcode=nullptr;
  std::vector<std::string> *MCP_EndProcess=nullptr;
  std::vector<int> *MCP_nMatchedPFPs=nullptr;
  std::vector<double> *MCP_totaldepE=nullptr;
  std::vector<std::vector<int>> *MCP_matchedPFP_ID=nullptr;
  std::vector<std::vector<double>> *MCP_matchedPFP_matchedE=nullptr;
  std::vector<double> *MCP_trueStartE=nullptr;
  std::vector<double> *MCP_trueStartP=nullptr;
  std::vector<std::vector<double>> *MCP_Px_eachpoint=nullptr;
  std::vector<std::vector<double>> *MCP_Py_eachpoint=nullptr;
  std::vector<std::vector<double>> *MCP_Pz_eachpoint=nullptr;
  std::vector<std::vector<double>> *MCP_x_eachpoint=nullptr;
  std::vector<std::vector<double>> *MCP_y_eachpoint=nullptr;
  std::vector<std::vector<double>> *MCP_z_eachpoint=nullptr;
  std::vector<std::vector<double>> *MCP_StartXYZ=nullptr;
  std::vector<std::vector<double>> *MCP_EndXYZ=nullptr;

  int n_PFPs=0;
  std::vector<int> *PFP_ID=nullptr;
  std::vector<int> *PFP_TrackShowerPdg=nullptr;
  std::vector<bool> *PFP_IsPrimary=nullptr;
  std::vector<std::pair<int,int>> *PFP_ID_to_bestmatchMCPid=nullptr;
  std::vector<double> *PFP_totaldepE=nullptr;
  std::vector<int> *PFP_primaryPFPid=nullptr;
  std::vector<bool> *PFP_isInNuSlice=nullptr;
  std::vector<std::vector<std::vector<double>>> *PFP_spacepoints_XYZ=nullptr;
  std::vector<std::vector<std::vector<double>>> *PFP_trajpoints_XYZ=nullptr;
  std::vector<double> *PFP_track_length=nullptr;
  std::vector<std::vector<double>> *PFP_track_start=nullptr;
  std::vector<std::vector<double>> *PFP_track_end=nullptr;
  std::vector<double> *PFP_track_theta=nullptr;
  std::vector<double> *PFP_track_phi=nullptr;
  std::vector<std::vector<std::vector<double>>> *PFP_ordered_spacepoints=nullptr;
  std::vector<std::vector<int>> *PFP_ordered_spacepoints_origidx=nullptr;
  std::vector<std::vector<int>> *PFP_goodreco_MCPvectorposes=nullptr;
  std::vector<std::vector<double>> *PFP_trueMCPend=nullptr;

  bool verbose;

  // Things we calculate
  std::map<int,int> MCP_ID_to_vectorpos;
  std::map<int,int> PFP_ID_to_vectorpos;
  double pidepE;
  std::vector<int> all_primary_PFP_IDs;
  std::vector<int> bestmatch_PFP_hierarchy_idxs;
  double bestmatch_matchedE;
  std::vector<int> primary_PFP_idxs;
  int n_p;
  int n_mu;
  int n_visp;
  int primarypfp_pdgcode = -999;
  std::vector<std::vector<int>> PFP_kinkidxs;

  // TH1D *h_mincosth;

  // constructor
  KinkFindingTree();

  // Get maps and do truth selection on depE
  bool Setup();

  // Get variables relating to best-match PFP hierarchy
  void GetPFPHierarchy();

  // Get MCPs associated to a PFP if it is "well reconstructed"
  std::vector<int> GetGoodRecoMCPsfromPFP(int pfpid);

  // Work out if track has a true kink
  double IsTrueKink(std::vector<int> MCP_vectorposes, double threshold);

  // Find reco kinks
  std::vector<int> FindRecoKinks(int pfp_idx, double mincosth, int evtnum);
  std::vector<int> FindRecoKinks_orderedsp(int pfp_idx, double mincosth, int evtnum, std::vector<std::vector<double>>);

  std::vector<int> FindRecoKinks_byangle(int pfp_idx, double mincosth, int evtnum, double kinkcosth, int _slider_window);
  std::vector<int> FindRecoKinks_byangle_orderedsp(int pfp_idx, double mincosth, int evtnum, double kinkcosth, int _slider_window, std::vector<std::vector<double>>);
};

KinkFindingTree::KinkFindingTree(){
  PFP_isInNuSlice = nullptr;
  verbose = false;
  pidepE = 0;
  // h_mincosth = new TH1D("h_mincosth",";Min(cos(#theta))",100,-1,1);
};

void SetKinkFindingTreeVariables(KinkFindingTree *vars, TTree *tr){
  // Set branch addresses
  tr->SetBranchAddress("nMCPs",&(vars->nMCPs));
  tr->SetBranchAddress("nu_MCPID",&(vars->nu_MCPID));
  tr->SetBranchAddress("MCP_ID",&(vars->MCP_ID));
  tr->SetBranchAddress("MCP_MotherID",&(vars->MCP_MotherID));
  tr->SetBranchAddress("MCP_DaughterIDs",&(vars->MCP_DaughterIDs));
  tr->SetBranchAddress("MCP_PDGcode",&(vars->MCP_PDGcode));
  tr->SetBranchAddress("MCP_EndProcess",&(vars->MCP_EndProcess));
  tr->SetBranchAddress("MCP_nMatchedPFPs",&(vars->MCP_nMatchedPFPs));
  tr->SetBranchAddress("MCP_totaldepE",&(vars->MCP_totaldepE));
  tr->SetBranchAddress("MCP_matchedPFP_ID",&(vars->MCP_matchedPFP_ID));
  tr->SetBranchAddress("MCP_matchedPFP_matchedE",&(vars->MCP_matchedPFP_matchedE));
  tr->SetBranchAddress("MCP_trueStartE",&(vars->MCP_trueStartE));
  tr->SetBranchAddress("MCP_trueStartP",&(vars->MCP_trueStartP));
  tr->SetBranchAddress("MCP_Px_eachpoint",&(vars->MCP_Px_eachpoint));
  tr->SetBranchAddress("MCP_Py_eachpoint",&(vars->MCP_Py_eachpoint));
  tr->SetBranchAddress("MCP_Pz_eachpoint",&(vars->MCP_Pz_eachpoint));
  tr->SetBranchAddress("MCP_x_eachpoint",&(vars->MCP_x_eachpoint));
  tr->SetBranchAddress("MCP_y_eachpoint",&(vars->MCP_y_eachpoint));
  tr->SetBranchAddress("MCP_z_eachpoint",&(vars->MCP_z_eachpoint));
  tr->SetBranchAddress("MCP_StartXYZ",&(vars->MCP_StartXYZ));
  tr->SetBranchAddress("MCP_EndXYZ",&(vars->MCP_EndXYZ));
  tr->SetBranchAddress("n_PFPs",&(vars->n_PFPs));
  tr->SetBranchAddress("PFP_ID",&(vars->PFP_ID));
  tr->SetBranchAddress("PFP_TrackShowerPdg",&(vars->PFP_TrackShowerPdg));
  tr->SetBranchAddress("PFP_IsPrimary",&(vars->PFP_IsPrimary));
  tr->SetBranchAddress("PFP_ID_to_bestmatchMCPid",&(vars->PFP_ID_to_bestmatchMCPid));
  tr->SetBranchAddress("PFP_totaldepE",&(vars->PFP_totaldepE));
  tr->SetBranchAddress("PFP_primaryPFPid",&(vars->PFP_primaryPFPid));
  tr->SetBranchAddress("PFP_isInNuSlice",&(vars->PFP_isInNuSlice));
  tr->SetBranchAddress("PFP_spacepoints_XYZ",&(vars->PFP_spacepoints_XYZ));
  tr->SetBranchAddress("PFP_trajpoints_XYZ",&(vars->PFP_trajpoints_XYZ));
  tr->SetBranchAddress("PFP_track_length",&(vars->PFP_track_length));
  tr->SetBranchAddress("PFP_track_start",&(vars->PFP_track_start));
  tr->SetBranchAddress("PFP_track_end",&(vars->PFP_track_end));
  tr->SetBranchAddress("PFP_track_theta",&(vars->PFP_track_theta));
  tr->SetBranchAddress("PFP_track_phi",&(vars->PFP_track_phi));
  tr->SetBranchAddress("PFP_ordered_spacepoints",&(vars->PFP_ordered_spacepoints));
  tr->SetBranchAddress("PFP_ordered_spacepoints_origidx",&(vars->PFP_ordered_spacepoints_origidx));
  tr->SetBranchAddress("PFP_trueMCPend",&(vars->PFP_trueMCPend));
  tr->SetBranchAddress("PFP_goodreco_MCPvectorposes",&(vars->PFP_goodreco_MCPvectorposes));
};

void SetKinkFindingTreeBranches(KinkFindingTree *vars, TTree *tr){
  // Set branch addresses
  tr->Branch("nMCPs",&(vars->nMCPs));
  tr->Branch("MCP_ID",&(vars->MCP_ID));
  tr->Branch("MCP_MotherID",&(vars->MCP_MotherID));
  tr->Branch("MCP_DaughterIDs",&(vars->MCP_DaughterIDs));
  tr->Branch("MCP_PDGcode",&(vars->MCP_PDGcode));
  tr->Branch("MCP_EndProcess",&(vars->MCP_EndProcess));
  tr->Branch("MCP_nMatchedPFPs",&(vars->MCP_nMatchedPFPs));
  tr->Branch("MCP_totaldepE",&(vars->MCP_totaldepE));
  tr->Branch("MCP_matchedPFP_ID",&(vars->MCP_matchedPFP_ID));
  tr->Branch("MCP_matchedPFP_matchedE",&(vars->MCP_matchedPFP_matchedE));
  tr->Branch("MCP_trueStartE",&(vars->MCP_trueStartE));
  tr->Branch("MCP_trueStartP",&(vars->MCP_trueStartP));
  tr->Branch("MCP_Px_eachpoint",&(vars->MCP_Px_eachpoint));
  tr->Branch("MCP_Py_eachpoint",&(vars->MCP_Py_eachpoint));
  tr->Branch("MCP_Pz_eachpoint",&(vars->MCP_Pz_eachpoint));
  tr->Branch("MCP_StartXYZ",&(vars->MCP_StartXYZ));
  tr->Branch("MCP_EndXYZ",&(vars->MCP_EndXYZ));
  tr->Branch("n_PFPs",&(vars->n_PFPs));
  tr->Branch("PFP_ID",&(vars->PFP_ID));
  tr->Branch("PFP_TrackShowerPdg",&(vars->PFP_TrackShowerPdg));
  tr->Branch("PFP_IsPrimary",&(vars->PFP_IsPrimary));
  tr->Branch("PFP_ID_to_bestmatchMCPid",&(vars->PFP_ID_to_bestmatchMCPid));
  tr->Branch("PFP_totaldepE",&(vars->PFP_totaldepE));
  tr->Branch("PFP_primaryPFPid",&(vars->PFP_primaryPFPid));
  tr->Branch("PFP_isInNuSlice",&(vars->PFP_isInNuSlice));
  tr->Branch("PFP_spacepoints_XYZ",&(vars->PFP_spacepoints_XYZ));
  tr->Branch("PFP_trajpoints_XYZ",&(vars->PFP_trajpoints_XYZ));
  tr->Branch("PFP_track_length",&(vars->PFP_track_length));
  tr->Branch("PFP_track_start",&(vars->PFP_track_start));
  tr->Branch("PFP_track_end",&(vars->PFP_track_end));
  tr->Branch("PFP_track_theta",&(vars->PFP_track_theta));
  tr->Branch("PFP_track_phi",&(vars->PFP_track_phi));
};

void SetupKinkFindingTreeForReading(KinkFindingTree *tr){
  tr->nMCPs=0;
  tr->nu_MCPID=0;
  tr->MCP_ID=nullptr;
  tr->MCP_MotherID=nullptr;
  tr->MCP_DaughterIDs=nullptr;
  tr->MCP_PDGcode=nullptr;
  tr->MCP_EndProcess=nullptr;
  tr->MCP_nMatchedPFPs=nullptr;
  tr->MCP_totaldepE=nullptr;
  tr->MCP_matchedPFP_ID=nullptr;
  tr->MCP_matchedPFP_matchedE=nullptr;
  tr->MCP_trueStartE=nullptr;
  tr->MCP_trueStartP=nullptr;
  tr->n_PFPs=0;
  tr->MCP_Px_eachpoint=nullptr;
  tr->MCP_Py_eachpoint=nullptr;
  tr->MCP_Pz_eachpoint=nullptr;
  tr->MCP_x_eachpoint=nullptr;
  tr->MCP_y_eachpoint=nullptr;
  tr->MCP_z_eachpoint=nullptr;
  tr->MCP_StartXYZ=nullptr;
  tr->MCP_EndXYZ=nullptr;
  tr->PFP_ID=nullptr;
  tr->PFP_TrackShowerPdg=nullptr;
  tr->PFP_IsPrimary=nullptr;
  tr->PFP_ID_to_bestmatchMCPid=nullptr;
  tr->PFP_totaldepE=nullptr;
  tr->PFP_primaryPFPid=nullptr;
  tr->PFP_isInNuSlice=nullptr;
  tr->PFP_spacepoints_XYZ=nullptr;
  tr->PFP_trajpoints_XYZ=nullptr;
  tr->PFP_track_length=nullptr;
  tr->PFP_track_start=nullptr;
  tr->PFP_track_end=nullptr;
  tr->PFP_track_theta=nullptr;
  tr->PFP_track_phi=nullptr;
  tr->PFP_goodreco_MCPvectorposes=nullptr;
  tr->PFP_trueMCPend=nullptr;
};

void ClearKinkFindingTree(KinkFindingTree *tr){
  tr->nMCPs=0;
  tr->nu_MCPID=0;
  tr->MCP_ID->clear();
  tr->MCP_MotherID->clear();
  tr->MCP_DaughterIDs->clear();
  tr->MCP_PDGcode->clear();
  tr->MCP_EndProcess->clear();
  tr->MCP_nMatchedPFPs->clear();
  tr->MCP_totaldepE->clear();
  tr->MCP_matchedPFP_ID->clear();
  tr->MCP_matchedPFP_matchedE->clear();
  tr->MCP_trueStartE->clear();
  tr->MCP_trueStartP->clear();
  tr->n_PFPs=0;
  tr->MCP_Px_eachpoint->clear();
  tr->MCP_Py_eachpoint->clear();
  tr->MCP_Pz_eachpoint->clear();
  tr->MCP_x_eachpoint->clear();
  tr->MCP_y_eachpoint->clear();
  tr->MCP_z_eachpoint->clear();
  tr->MCP_StartXYZ->clear();
  tr->MCP_EndXYZ->clear();
  tr->PFP_ID->clear();
  tr->PFP_TrackShowerPdg->clear();
  tr->PFP_IsPrimary->clear();
  tr->PFP_ID_to_bestmatchMCPid->clear();
  tr->PFP_totaldepE->clear();
  tr->PFP_primaryPFPid->clear();
  tr->PFP_isInNuSlice->clear();
  tr->PFP_spacepoints_XYZ->clear();
  tr->PFP_trajpoints_XYZ->clear();
  tr->PFP_track_length->clear();
  tr->PFP_track_start->clear();
  tr->PFP_track_end->clear();
  tr->PFP_track_theta->clear();
  tr->PFP_track_phi->clear();
  tr->PFP_trueMCPend->clear();
  tr->PFP_goodreco_MCPvectorposes->clear();
};

void CopyMCPs(KinkFindingTree *originaltree, KinkFindingTree *newtree){
  newtree->nMCPs = originaltree->nMCPs;
  *(newtree->MCP_ID) = *(originaltree->MCP_ID);
  *(newtree->MCP_MotherID) = *(originaltree->MCP_MotherID);
  *(newtree->MCP_DaughterIDs) = *(originaltree->MCP_DaughterIDs);
  *(newtree->MCP_PDGcode) = *(originaltree->MCP_PDGcode);
  *(newtree->MCP_EndProcess) = *(originaltree->MCP_EndProcess);
  *(newtree->MCP_nMatchedPFPs) = *(originaltree->MCP_nMatchedPFPs);
  *(newtree->MCP_totaldepE) = *(originaltree->MCP_totaldepE);
  *(newtree->MCP_matchedPFP_ID) = *(originaltree->MCP_matchedPFP_ID);
  *(newtree->MCP_matchedPFP_matchedE) = *(originaltree->MCP_matchedPFP_matchedE);
  *(newtree->MCP_trueStartE) = *(originaltree->MCP_trueStartE);
  *(newtree->MCP_trueStartP) = *(originaltree->MCP_trueStartP);
  *(newtree->MCP_Px_eachpoint) = *(originaltree->MCP_Px_eachpoint);
  *(newtree->MCP_Py_eachpoint) = *(originaltree->MCP_Py_eachpoint);
  *(newtree->MCP_Pz_eachpoint) = *(originaltree->MCP_Pz_eachpoint);
  *(newtree->MCP_x_eachpoint) = *(originaltree->MCP_x_eachpoint);
  *(newtree->MCP_y_eachpoint) = *(originaltree->MCP_y_eachpoint);
  *(newtree->MCP_z_eachpoint) = *(originaltree->MCP_z_eachpoint);
  *(newtree->MCP_StartXYZ) = *(originaltree->MCP_StartXYZ);
  *(newtree->MCP_EndXYZ) = *(originaltree->MCP_EndXYZ);
};

void CopyPFP(int i_pfp, KinkFindingTree *originaltree, KinkFindingTree *newtree){
  newtree->PFP_ID->push_back(originaltree->PFP_ID->at(i_pfp));
  newtree->PFP_TrackShowerPdg->push_back(originaltree->PFP_TrackShowerPdg->at(i_pfp));
  newtree->PFP_IsPrimary->push_back(originaltree->PFP_IsPrimary->at(i_pfp));
  newtree->PFP_ID_to_bestmatchMCPid->push_back(originaltree->PFP_ID_to_bestmatchMCPid->at(i_pfp));
  newtree->PFP_totaldepE->push_back(originaltree->PFP_totaldepE->at(i_pfp));
  newtree->PFP_primaryPFPid->push_back(originaltree->PFP_primaryPFPid->at(i_pfp));
  newtree->PFP_spacepoints_XYZ->push_back(originaltree->PFP_spacepoints_XYZ->at(i_pfp));
  newtree->PFP_trajpoints_XYZ->push_back(originaltree->PFP_trajpoints_XYZ->at(i_pfp));
  newtree->PFP_track_length->push_back(originaltree->PFP_track_length->at(i_pfp));
  newtree->PFP_track_start->push_back(originaltree->PFP_track_start->at(i_pfp));
  newtree->PFP_track_end->push_back(originaltree->PFP_track_end->at(i_pfp));
  newtree->PFP_track_theta->push_back(originaltree->PFP_track_theta->at(i_pfp));
  newtree->PFP_track_phi->push_back(originaltree->PFP_track_phi->at(i_pfp));
  if (originaltree->PFP_isInNuSlice){
    newtree->PFP_isInNuSlice->push_back(originaltree->PFP_isInNuSlice->at(i_pfp));
  }
};

bool KinkFindingTree::Setup(){
  // Hack because I messed up in the module
  n_PFPs = PFP_ID->size();

  // Make some maps that will help us later
  MCP_ID_to_vectorpos.clear();
  PFP_ID_to_vectorpos.clear();
  for (int i_mcp=0; i_mcp<nMCPs; i_mcp++){
      MCP_ID_to_vectorpos.insert(std::pair<int,int>(MCP_ID->at(i_mcp), i_mcp));
  }
  for (int i_pfp=0; i_pfp<n_PFPs; i_pfp++){
      PFP_ID_to_vectorpos.insert(std::pair<int,int>(PFP_ID->at(i_pfp),i_pfp));
  }
  return true;
};

void KinkFindingTree::GetPFPHierarchy(){
  // --- Calculate variables relating to PFP hierarchy --- //
  all_primary_PFP_IDs = GetAllPrimaryPFPIDs(PFP_primaryPFPid,PFP_isInNuSlice,PFP_ID_to_vectorpos);

  bestmatch_PFP_hierarchy_idxs.clear();
  bestmatch_matchedE = 0;

  // Need to identify a PFP hierarchy
  // Do this by looking for primary PFPs with best-match MCP in the hierarchy
  // and largest matched energy fraction
  primary_PFP_idxs.clear();
  for (int i_p=0; i_p < n_PFPs; i_p++){
      int pfp_id = PFP_ID->at(i_p);
      if (std::find(all_primary_PFP_IDs.begin(), all_primary_PFP_IDs.end(), pfp_id) == all_primary_PFP_IDs.end()) continue;

      int bestmatch_mcpid = PFP_ID_to_bestmatchMCPid->at(i_p).second;
      if (std::find(MCP_ID->begin(), MCP_ID->end(), bestmatch_mcpid) == MCP_ID->end()) continue;

      primary_PFP_idxs.push_back(i_p);
  }

  for (size_t i_primary=0; i_primary < primary_PFP_idxs.size(); i_primary++){
      int matched_mcpid = PFP_ID_to_bestmatchMCPid->at(primary_PFP_idxs.at(i_primary)).second;
      int matched_mcp_idx = MCP_ID_to_vectorpos.find(matched_mcpid)->second;

      // Now look at PFP hierarchy
      // First, gather all PFPs with this PFP as their primary
      std::vector<int> PFP_hierarchy_idxs;
      for (int i_p=0; i_p < n_PFPs; i_p++){
          if (PFP_primaryPFPid->at(i_p) == PFP_ID->at(primary_PFP_idxs.at(i_primary))){
              PFP_hierarchy_idxs.push_back(i_p);
          }
      }

      // Now match energy between MCP hierarchy and PFP hierarchy
      double hierarchy_pfphierarchy_matchedE = 0;
      // Loop through MCP hierarchy and add matched energy if it comes from a PFP in our PFP hierarchy
      for (int mcp_idx=0; mcp_idx<nMCPs; mcp_idx++){
          for (int pfp_matchidx=0; pfp_matchidx<MCP_nMatchedPFPs->at(mcp_idx); pfp_matchidx++){
              // If the matched PFP is not in the PFP hierarchy, continue
              int pfp_match_id = MCP_matchedPFP_ID->at(mcp_idx).at(pfp_matchidx);
              int pfp_match_trueidx = PFP_ID_to_vectorpos.find(pfp_match_id)->second;
              if (std::find(PFP_hierarchy_idxs.begin(), PFP_hierarchy_idxs.end(), pfp_match_trueidx) == PFP_hierarchy_idxs.end())
                  continue;

              // Now sum up matched energy
              hierarchy_pfphierarchy_matchedE += MCP_matchedPFP_matchedE->at(mcp_idx).at(pfp_matchidx);
          }
      }

      if (hierarchy_pfphierarchy_matchedE > bestmatch_matchedE){
          bestmatch_matchedE = hierarchy_pfphierarchy_matchedE;
          bestmatch_PFP_hierarchy_idxs = PFP_hierarchy_idxs;
      }
  }

  if (bestmatch_PFP_hierarchy_idxs.size()>0 && verbose){
      int matched_mcpid = PFP_ID_to_bestmatchMCPid->at(bestmatch_PFP_hierarchy_idxs.at(0)).second;
      int matched_mcp_idx = MCP_ID_to_vectorpos.find(matched_mcpid)->second;
      std::cout << "Primary PFP matched to PDG code " << MCP_PDGcode->at(matched_mcp_idx) << std::endl;
      primarypfp_pdgcode = MCP_PDGcode->at(matched_mcp_idx);
  }
};

std::vector<int> KinkFindingTree::GetGoodRecoMCPsfromPFP(int pfpid){
  // For a given PFP, find all MCPs associated to it
  // Of those, find all MCPs that have > 90% (tuneable) of their energy in this PFP
  // Final requirement for this to be "good" reco: these MCPs have to account for >90% of the PFP energy (also tuneable)

  // Return a vector of good MCP IDs (where "ID" means position in the vector MCP_blah)
  // If the vector is empty, then this PFP did not have "good" reco matches
  std::vector<int> MCP_vectorpos;

  auto search = PFP_ID_to_vectorpos.find(pfpid);
  int pfp_idx = search->second;

  // The tree saves MCP->all matched PFP mapping but not the other way around, so loop through all MCPs and ask if they match this PFP
  double matchedE = 0;

  for (int i_mcp=0; i_mcp<nMCPs; i_mcp++){
    for (int i_matchpfp=0; i_matchpfp<MCP_nMatchedPFPs->at(i_mcp); i_matchpfp++){

      // Are we looking at the PFP we want? If not, move on to the next one
      if (MCP_matchedPFP_ID->at(i_mcp).at(i_matchpfp)!=pfpid) continue;

      // If we have the right PFP, does the matched energy constitute more than 90% of the MCP energy? If not, move on
      if (MCP_matchedPFP_matchedE->at(i_mcp).at(i_matchpfp) < 0.9*MCP_totaldepE->at(i_mcp)) continue;

      // Check that this MCP accounts for at least 5% of the PFP energy
      if (MCP_matchedPFP_matchedE->at(i_mcp).at(i_matchpfp) < 0.1*PFP_totaldepE->at(pfp_idx)) continue;

      // Don't use electrons because they make terrible tracks
      if (TMath::Abs(MCP_PDGcode->at(i_mcp))==11) continue;

      // If those conditions are passed, put this MCP vectorpos in the vector
      MCP_vectorpos.push_back(i_mcp);

      // Also add the matched energy (we'll use this in the final check)
      matchedE += MCP_matchedPFP_matchedE->at(i_mcp).at(i_matchpfp);
    } // end loop over i_matchpfp
  } // end loop over i_mcp

  // Final check for "good" reconstruction: is the total matched energy from these MCPs more than 90% of the PFP energy?
  if (matchedE < 0.9*PFP_totaldepE->at(pfp_idx)){
    //std::cout << "matchedE = " << matchedE << ", 0.9*PFP_E = " << 0.9*PFP_totaldepE->at(pfp_idx) << ". Not well reconstructed" << std::endl;
    MCP_vectorpos.clear();
  }

  return MCP_vectorpos;
}; // end function

double KinkFindingTree::IsTrueKink(std::vector<int> MCP_vectorposes, double threshold){
  if (MCP_vectorposes.size()>1){
    //std::cout << "Assuming true kink because PFP contains " << MCP_vectorposes.size() << " MCPs" << std::endl;

    double mincostheta_twomcps=1.0;

    for (int i_m=0; i_m<MCP_vectorposes.size()-1; i_m++){
      int mcp_1 = MCP_vectorposes.at(i_m);
      int mcp_2 = MCP_vectorposes.at(i_m+1);
      std::vector<double> Px_1 = MCP_Px_eachpoint->at(mcp_1);
      std::vector<double> Py_1 = MCP_Py_eachpoint->at(mcp_1);
      std::vector<double> Pz_1 = MCP_Pz_eachpoint->at(mcp_1);
      std::vector<double> Px_2 = MCP_Px_eachpoint->at(mcp_2);
      std::vector<double> Py_2 = MCP_Py_eachpoint->at(mcp_2);
      std::vector<double> Pz_2 = MCP_Pz_eachpoint->at(mcp_2);

      // Find last non-(0,0,0) momentum in MCP 1
      double vec1_x=0, vec1_y=0, vec1_z=0;
      for (int i_p=0; i_p<Px_1.size(); i_p++){
        if (!(Px_1.at(i_p)==0 && Py_1.at(i_p) == 0 && Pz_1.at(i_p) ==0)){
          vec1_x = Px_1.at(i_p);
          vec1_y = Py_1.at(i_p);
          vec1_z = Pz_1.at(i_p);
        }
      }

      TVector3 vec1(vec1_x,vec1_y,vec1_z);
      TVector3 vec2(Px_2.at(0),Py_2.at(0),Pz_2.at(0));
      vec1.SetMag(1);
      vec2.SetMag(1);

      if (vec1.Dot(vec2)<mincostheta_twomcps) mincostheta_twomcps = vec1.Dot(vec2);
    }

    // return true;
    return -2.0+mincostheta_twomcps;
  } // end if (MCP_vectorposes.size()>1)

  int mcp_idx = MCP_vectorposes.at(0);
  std::vector<double> Px = MCP_Px_eachpoint->at(mcp_idx);
  std::vector<double> Py = MCP_Py_eachpoint->at(mcp_idx);
  std::vector<double> Pz = MCP_Pz_eachpoint->at(mcp_idx);

  double mincosth = 1.0;

  // std::cout << "Px.size() = " << Px.size() << std::endl;
  // std::cout << "Py.size() = " << Py.size() << std::endl;
  // std::cout << "Pz.size() = " << Pz.size() << std::endl;

  for (int i_pt=0; i_pt < Px.size()-1; i_pt++){

    TVector3 vec1(1,1,1);
    vec1.SetXYZ(Px.at(i_pt),Py.at(i_pt),Pz.at(i_pt));
    TVector3 vec2(1,1,1);
    vec2.SetXYZ(Px.at(i_pt+1),Py.at(i_pt+1),Pz.at(i_pt+1));

    if (vec1.Mag()==0 || vec2.Mag()==0) continue;

    vec1.SetMag(1);
    vec2.SetMag(1);

    double costh = vec1.Dot(vec2);

    // std::cout << Px.at(i_pt) << ", " << Py.at(i_pt) << ", " << Pz.at(i_pt) << std::endl;
    // std::cout << Px.at(i_pt+1) << ", " << Py.at(i_pt+1) << ", " << Pz.at(i_pt+1) << std::endl;
    // std::cout << costh << std::endl;

    if (costh<mincosth) mincosth = costh;

  } // end loop over i_pt

  // std::cout << "Min costheta = " << mincosth << std::endl;
  // h_mincosth->Fill(mincosth);
  return mincosth;

  // if (mincosth<threshold) return true;
  // else return false;
};


std::vector<int> KinkFindingTree::FindRecoKinks(int pfp_idx, double mincosth, int evtnum){

  // Order spacepoints (by closest to reco track start)
  std::vector<std::pair<std::vector<double>,int>> ordered_spacepoints_pair_origidx;
  // if (PFP_track_start->at(pfp_idx).at(0)==-999 && PFP_track_start->at(pfp_idx).at(1)==-999 && PFP_track_start->at(pfp_idx).at(2)==-999){
    ordered_spacepoints_pair_origidx =  OrderSpacepoints(PFP_spacepoints_XYZ->at(pfp_idx));
  // }
  // else{
  //   ordered_spacepoints_pair_origidx =  OrderSpacepoints(PFP_spacepoints_XYZ->at(pfp_idx), PFP_track_start->at(pfp_idx));
  // }

  // Get 1D projection onto track direction
  std::vector<double> proj1d = GetSPDirVec(ordered_spacepoints_pair_origidx,PFP_track_theta->at(pfp_idx),PFP_track_phi->at(pfp_idx));

  // Get vector of cumulative summed distance from mean
  std::vector<double> cusum = GetCuSumVec(proj1d);

  // Get vector of Welch's t-test
  std::vector<double> ttest = GetWelchttestVec(proj1d);

  // Find kinks
  std::vector<int> kinkidxs = *(SplitTracks(proj1d,cusum,ttest));

  std::vector<std::vector<double>> ordered_spacepoints;
  for (int i=0; i<ordered_spacepoints_pair_origidx.size();i++){
    ordered_spacepoints.push_back(ordered_spacepoints_pair_origidx.at(i).first);
  }
  PFP_spacepoints_XYZ->at(pfp_idx) = ordered_spacepoints;

  return kinkidxs;
};


std::vector<int> KinkFindingTree::FindRecoKinks_orderedsp(int pfp_idx, double mincosth, int evtnum, std::vector<std::vector<double>> truekinks=std::vector<std::vector<double>>()){

  // Order spacepoints (by closest to reco track start)
  std::vector<std::pair<std::vector<double>,int>> ordered_spacepoints_pair_origidx;
  for (int i=0; i<PFP_ordered_spacepoints->at(pfp_idx).size();i++){
    std::vector<double> spacepoints_i = PFP_ordered_spacepoints->at(pfp_idx).at(i);
    int idx = PFP_ordered_spacepoints_origidx->at(pfp_idx).at(i);
    ordered_spacepoints_pair_origidx.push_back(std::make_pair(spacepoints_i,idx));
  }

  // Get 1D projection onto track direction
  std::vector<double> proj1d = GetSPDirVec(ordered_spacepoints_pair_origidx,PFP_track_theta->at(pfp_idx),PFP_track_phi->at(pfp_idx));

  // Get vector of cumulative summed distance from mean
  std::vector<double> cusum = GetCuSumVec(proj1d);

  // Get vector of Welch's t-test
  std::vector<double> ttest = GetWelchttestVec(proj1d);

  // Find kinks
  std::vector<int> kinkidxs = *(SplitTracks(proj1d,cusum,ttest));

  // if (kinkidxs.size()>0){
    PlotLocalLinearityDetails(PFP_ordered_spacepoints->at(pfp_idx),proj1d,cusum,ttest,kinkidxs,pfp_idx,PFP_ordered_spacepoints->at(pfp_idx).at(0),truekinks,*MCP_x_eachpoint,*MCP_y_eachpoint,*MCP_z_eachpoint,mincosth,evtnum);
  // }

  return kinkidxs;
};


std::vector<int> KinkFindingTree::FindRecoKinks_byangle(int pfp_idx, double mincosth, int evtnum, double kinkcosth, int _slider_window=20){

  // Order spacepoints (by closest to reco track start)
  std::vector<std::pair<std::vector<double>,int>> ordered_spacepoints_pair_origidx;
  // if (PFP_track_start->at(pfp_idx).at(0)==-999 && PFP_track_start->at(pfp_idx).at(1)==-999 && PFP_track_start->at(pfp_idx).at(2)==-999){
    ordered_spacepoints_pair_origidx =  OrderSpacepoints(PFP_spacepoints_XYZ->at(pfp_idx));
  // }
  // else{
  //   ordered_spacepoints_pair_origidx =  OrderSpacepoints(PFP_spacepoints_XYZ->at(pfp_idx), PFP_track_start->at(pfp_idx));
  // }

  // Get vector of local linearity based on angle between consecutive points
  std::vector<double> locallin_vec = GetLocalLinearityVec(ordered_spacepoints_pair_origidx, _slider_window);

  // Find kinks
  std::vector<int> kinkidxs;
  for (int i=0; i<locallin_vec.size(); i++){
    if (locallin_vec.at(i)<kinkcosth){
      kinkidxs.push_back(i);
    }
  }

  // std::vector<std::vector<double>> ordered_spacepoints;
  // for (int i=0; i<ordered_spacepoints_pair_origidx.size();i++){
  //   ordered_spacepoints.push_back(ordered_spacepoints_pair_origidx.at(i).first);
  // }

  // if (mincosth>0.9 && kinkidxs.size()>0){
    // PlotLocalLinearityDetails(ordered_spacepoints,proj1d,cusum,ttest,kinkidxs,pfp_idx,PFP_track_start->at(pfp_idx),mincosth,evtnum);
  // }

  return kinkidxs;
};


std::vector<int> KinkFindingTree::FindRecoKinks_byangle_orderedsp(int pfp_idx, double mincosth, int evtnum, double kinkcosth, int _slider_window=20,std::vector<std::vector<double>> truekinks=std::vector<std::vector<double>>()){

  // Order spacepoints (by closest to reco track start)
  std::vector<std::pair<std::vector<double>,int>> ordered_spacepoints_pair_origidx;
  for (int i=0; i<PFP_ordered_spacepoints->at(pfp_idx).size();i++){
    std::vector<double> spacepoints_i = PFP_ordered_spacepoints->at(pfp_idx).at(i);
    int idx = PFP_ordered_spacepoints_origidx->at(pfp_idx).at(i);
    ordered_spacepoints_pair_origidx.push_back(std::make_pair(spacepoints_i,idx));
  }

  // Get vector of local linearity based on angle between consecutive points
  std::vector<double> locallin_vec = GetLocalLinearityVec(ordered_spacepoints_pair_origidx, _slider_window);

  // Find kinks
  std::vector<int> kinkidxs;
  for (int i=0; i<locallin_vec.size(); i++){
    if (locallin_vec.at(i)<kinkcosth){
      kinkidxs.push_back(i);
    }
  }

  // if (kinkidxs.size()>0){
    PlotLocalLinearityDetails_byangle(PFP_ordered_spacepoints->at(pfp_idx),locallin_vec,kinkidxs,pfp_idx,PFP_ordered_spacepoints->at(pfp_idx).at(0),truekinks,mincosth,evtnum);
  // }

  return kinkidxs;
};
