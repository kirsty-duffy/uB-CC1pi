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

struct PiRecoTree{
  // Variables to fill from tree
  int PiPlusHierarchy_nMCPs=0;
  std::vector<int> *PiPlusHierarchy_MCP_ID=nullptr;
  std::vector<int> *PiPlusHierarchy_MCP_MotherID=nullptr;
  std::vector<std::vector<int>> *PiPlusHierarchy_MCP_DaughterIDs=nullptr;
  std::vector<int> *PiPlusHierarchy_MCP_PDGcode=nullptr;
  std::vector<std::string> *PiPlusHierarchy_MCP_EndProcess=nullptr;
  std::vector<int> *PiPlusHierarchy_MCP_nMatchedPFPs=nullptr;
  std::vector<double> *PiPlusHierarchy_MCP_totaldepE=nullptr;
  std::vector<std::vector<int>> *PiPlusHierarchy_MCP_matchedPFP_ID=nullptr;
  std::vector<std::vector<double>> *PiPlusHierarchy_MCP_matchedPFP_matchedE=nullptr;
  std::vector<double> *PiPlusHierarchy_MCP_trueStartE=nullptr;
  std::vector<double> *PiPlusHierarchy_MCP_trueStartP=nullptr;
  std::vector<int> *PiPlusHierarchy_LogicalPion_MCPids=nullptr;
  int n_PFPs=0;
  std::vector<std::vector<double>> *PiPlusHierarchy_MCP_Px_eachpoint=nullptr;
  std::vector<std::vector<double>> *PiPlusHierarchy_MCP_Py_eachpoint=nullptr;
  std::vector<std::vector<double>> *PiPlusHierarchy_MCP_Pz_eachpoint=nullptr;
  std::vector<std::vector<double>> *PiPlusHierarchy_MCP_StartXYZ=nullptr;
  std::vector<std::vector<double>> *PiPlusHierarchy_MCP_EndXYZ=nullptr;
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
  PiRecoTree();

  // Get maps and do truth selection on depE
  bool Setup();

  // Get variables relating to logical pion
  void CalcLogicalPionVars();

  // Get variables relating to best-match PFP hierarchy
  void GetPFPHierarchy();

  // Get MCPs associated to a PFP if it is "well reconstructed"
  std::vector<int> GetGoodRecoMCPsfromPFP(int pfpid);

  // Work out if track has a true kink
  double IsTrueKink(std::vector<int> MCP_vectorposes, double threshold);

  // Find reco kinks
  std::vector<int> FindRecoKinks(int pfp_idx, double mincosth, int evtnum);

  std::vector<int> FindRecoKinks_byangle(int pfp_idx, double mincosth, int evtnum, double kinkcosth, int _slider_window);
};

PiRecoTree::PiRecoTree(){
  PFP_isInNuSlice = nullptr;
  verbose = false;
  pidepE = 0;
  // h_mincosth = new TH1D("h_mincosth",";Min(cos(#theta))",100,-1,1);
};

void SetPiRecoTreeVariables(PiRecoTree *vars, TTree *tr){
  // Set branch addresses
  tr->SetBranchAddress("PiPlusHierarchy_nMCPs",&(vars->PiPlusHierarchy_nMCPs));
  tr->SetBranchAddress("PiPlusHierarchy_MCP_ID",&(vars->PiPlusHierarchy_MCP_ID));
  tr->SetBranchAddress("PiPlusHierarchy_MCP_MotherID",&(vars->PiPlusHierarchy_MCP_MotherID));
  tr->SetBranchAddress("PiPlusHierarchy_MCP_DaughterIDs",&(vars->PiPlusHierarchy_MCP_DaughterIDs));
  tr->SetBranchAddress("PiPlusHierarchy_MCP_PDGcode",&(vars->PiPlusHierarchy_MCP_PDGcode));
  tr->SetBranchAddress("PiPlusHierarchy_MCP_EndProcess",&(vars->PiPlusHierarchy_MCP_EndProcess));
  tr->SetBranchAddress("PiPlusHierarchy_MCP_nMatchedPFPs",&(vars->PiPlusHierarchy_MCP_nMatchedPFPs));
  tr->SetBranchAddress("PiPlusHierarchy_MCP_totaldepE",&(vars->PiPlusHierarchy_MCP_totaldepE));
  tr->SetBranchAddress("PiPlusHierarchy_MCP_matchedPFP_ID",&(vars->PiPlusHierarchy_MCP_matchedPFP_ID));
  tr->SetBranchAddress("PiPlusHierarchy_MCP_matchedPFP_matchedE",&(vars->PiPlusHierarchy_MCP_matchedPFP_matchedE));
  tr->SetBranchAddress("PiPlusHierarchy_MCP_trueStartE",&(vars->PiPlusHierarchy_MCP_trueStartE));
  tr->SetBranchAddress("PiPlusHierarchy_MCP_trueStartP",&(vars->PiPlusHierarchy_MCP_trueStartP));
  tr->SetBranchAddress("PiPlusHierarchy_LogicalPion_MCPids",&(vars->PiPlusHierarchy_LogicalPion_MCPids));
  tr->SetBranchAddress("PiPlusHierarchy_MCP_Px_eachpoint",&(vars->PiPlusHierarchy_MCP_Px_eachpoint));
  tr->SetBranchAddress("PiPlusHierarchy_MCP_Py_eachpoint",&(vars->PiPlusHierarchy_MCP_Py_eachpoint));
  tr->SetBranchAddress("PiPlusHierarchy_MCP_Pz_eachpoint",&(vars->PiPlusHierarchy_MCP_Pz_eachpoint));
  tr->SetBranchAddress("PiPlusHierarchy_MCP_StartXYZ",&(vars->PiPlusHierarchy_MCP_StartXYZ));
  tr->SetBranchAddress("PiPlusHierarchy_MCP_EndXYZ",&(vars->PiPlusHierarchy_MCP_EndXYZ));
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
};

bool PiRecoTree::Setup(){
  // Make some maps that will help us later
  MCP_ID_to_vectorpos.clear();
  PFP_ID_to_vectorpos.clear();
  for (int i_mcp=0; i_mcp<PiPlusHierarchy_nMCPs; i_mcp++){
      MCP_ID_to_vectorpos.insert(std::pair<int,int>(PiPlusHierarchy_MCP_ID->at(i_mcp), i_mcp));
  }
  for (int i_pfp=0; i_pfp<n_PFPs; i_pfp++){
      PFP_ID_to_vectorpos.insert(std::pair<int,int>(PFP_ID->at(i_pfp),i_pfp));
  }

  // Truth selection: calculate total dep E from logical pion
  // If it's 0 (actually <1e-12), move on to the next pi+
  pidepE = 0;
  for (size_t i_pi=0; i_pi<PiPlusHierarchy_LogicalPion_MCPids->size(); i_pi++){
      int pi_mcp = MCP_ID_to_vectorpos.find(PiPlusHierarchy_LogicalPion_MCPids->at(i_pi))->second;
      pidepE += PiPlusHierarchy_MCP_totaldepE->at(pi_mcp);
  }
  if (pidepE < 1e-12){
      return false;
  }
  return true;
};

void PiRecoTree::CalcLogicalPionVars(){
  // --- Calculate variables relating to logical pion --- //

  // No. true pi+ in logical pion
  if (verbose) std::cout << PiPlusHierarchy_LogicalPion_MCPids->size() << " pi+ in 'logical pion'" << std::endl;

  // No mu+ or p daughters of logical pion
  int finalPi_MCPID = PiPlusHierarchy_LogicalPion_MCPids->at(PiPlusHierarchy_LogicalPion_MCPids->size()-1);
  int lp_endmcp = MCP_ID_to_vectorpos.find(finalPi_MCPID)->second;
  n_p = 0;
  n_mu = 0;
  n_visp = 0;
  for (size_t i_d=0; i_d < PiPlusHierarchy_MCP_DaughterIDs->at(lp_endmcp).size(); i_d++){
    int daughter_mcp = MCP_ID_to_vectorpos.find(PiPlusHierarchy_MCP_DaughterIDs->at(lp_endmcp).at(i_d))->second;
    if (PiPlusHierarchy_MCP_PDGcode->at(daughter_mcp) == 2212){
        n_p++;
        if (PiPlusHierarchy_MCP_trueStartP->at(daughter_mcp) > 0.3) n_visp++;
    }
    else if (PiPlusHierarchy_MCP_PDGcode->at(daughter_mcp) == -13){
        n_mu++;
    }
  }
  if (verbose){
    std::cout << " --- " << n_p << " proton daughters (" << n_visp << " with p>300 MeV)" << std::endl;
    std::cout << " --- " << n_mu << " mu+ daughters" << std::endl;
  }

};

void PiRecoTree::GetPFPHierarchy(){
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
      if (std::find(PiPlusHierarchy_MCP_ID->begin(), PiPlusHierarchy_MCP_ID->end(), bestmatch_mcpid) == PiPlusHierarchy_MCP_ID->end()) continue;

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
      for (int mcp_idx=0; mcp_idx<PiPlusHierarchy_nMCPs; mcp_idx++){
          for (int pfp_matchidx=0; pfp_matchidx<PiPlusHierarchy_MCP_nMatchedPFPs->at(mcp_idx); pfp_matchidx++){
              // If the matched PFP is not in the PFP hierarchy, continue
              int pfp_match_id = PiPlusHierarchy_MCP_matchedPFP_ID->at(mcp_idx).at(pfp_matchidx);
              int pfp_match_trueidx = PFP_ID_to_vectorpos.find(pfp_match_id)->second;
              if (std::find(PFP_hierarchy_idxs.begin(), PFP_hierarchy_idxs.end(), pfp_match_trueidx) == PFP_hierarchy_idxs.end())
                  continue;

              // Now sum up matched energy
              hierarchy_pfphierarchy_matchedE += PiPlusHierarchy_MCP_matchedPFP_matchedE->at(mcp_idx).at(pfp_matchidx);
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
      std::cout << "Primary PFP matched to PDG code " << PiPlusHierarchy_MCP_PDGcode->at(matched_mcp_idx) << std::endl;
      primarypfp_pdgcode = PiPlusHierarchy_MCP_PDGcode->at(matched_mcp_idx);
  }
};

std::vector<int> PiRecoTree::GetGoodRecoMCPsfromPFP(int pfpid){
  // For a given PFP, find all MCPs associated to it
  // Of those, find all MCPs that have > 90% (tuneable) of their energy in this PFP
  // Final requirement for this to be "good" reco: these MCPs have to account for >90% of the PFP energy (also tuneable)

  // Return a vector of good MCP IDs (where "ID" means position in the vector PiPlusHierarchy_MCP_blah)
  // If the vector is empty, then this PFP did not have "good" reco matches
  std::vector<int> MCP_vectorpos;

  auto search = PFP_ID_to_vectorpos.find(pfpid);
  int pfp_idx = search->second;

  // The tree saves MCP->all matched PFP mapping but not the other way around, so loop through all MCPs and ask if they match this PFP
  double matchedE = 0;

  for (int i_mcp=0; i_mcp<PiPlusHierarchy_nMCPs; i_mcp++){
    for (int i_matchpfp=0; i_matchpfp<PiPlusHierarchy_MCP_nMatchedPFPs->at(i_mcp); i_matchpfp++){

      // Are we looking at the PFP we want? If not, move on to the next one
      if (PiPlusHierarchy_MCP_matchedPFP_ID->at(i_mcp).at(i_matchpfp)!=pfpid) continue;

      // If we have the right PFP, does the matched energy constitute more than 90% of the MCP energy? If not, move on
      if (PiPlusHierarchy_MCP_matchedPFP_matchedE->at(i_mcp).at(i_matchpfp) < 0.9*PiPlusHierarchy_MCP_totaldepE->at(i_mcp)) continue;

      // Check that this MCP accounts for at least 5% of the PFP energy
      if (PiPlusHierarchy_MCP_matchedPFP_matchedE->at(i_mcp).at(i_matchpfp) < 0.1*PFP_totaldepE->at(pfp_idx)) continue;

      // Don't use electrons because they make terrible tracks
      if (TMath::Abs(PiPlusHierarchy_MCP_PDGcode->at(i_mcp))==11) continue;

      // If those conditions are passed, put this MCP vectorpos in the vector
      MCP_vectorpos.push_back(i_mcp);

      // Also add the matched energy (we'll use this in the final check)
      matchedE += PiPlusHierarchy_MCP_matchedPFP_matchedE->at(i_mcp).at(i_matchpfp);
    } // end loop over i_matchpfp
  } // end loop over i_mcp

  // Final check for "good" reconstruction: is the total matched energy from these MCPs more than 90% of the PFP energy?
  if (matchedE < 0.9*PFP_totaldepE->at(pfp_idx)){
    //std::cout << "matchedE = " << matchedE << ", 0.9*PFP_E = " << 0.9*PFP_totaldepE->at(pfp_idx) << ". Not well reconstructed" << std::endl;
    MCP_vectorpos.clear();
  }

  return MCP_vectorpos;
}; // end function

double PiRecoTree::IsTrueKink(std::vector<int> MCP_vectorposes, double threshold){
  if (MCP_vectorposes.size()>1){
    std::cout << "Assuming true kink because PFP contains " << MCP_vectorposes.size() << " MCPs" << std::endl;

    double mincostheta_twomcps=1.0;

    for (int i_m=0; i_m<MCP_vectorposes.size()-1; i_m++){
      int mcp_1 = MCP_vectorposes.at(i_m);
      int mcp_2 = MCP_vectorposes.at(i_m+1);
      std::vector<double> Px_1 = PiPlusHierarchy_MCP_Px_eachpoint->at(mcp_1);
      std::vector<double> Py_1 = PiPlusHierarchy_MCP_Py_eachpoint->at(mcp_1);
      std::vector<double> Pz_1 = PiPlusHierarchy_MCP_Pz_eachpoint->at(mcp_1);
      std::vector<double> Px_2 = PiPlusHierarchy_MCP_Px_eachpoint->at(mcp_2);
      std::vector<double> Py_2 = PiPlusHierarchy_MCP_Py_eachpoint->at(mcp_2);
      std::vector<double> Pz_2 = PiPlusHierarchy_MCP_Pz_eachpoint->at(mcp_2);

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
  std::vector<double> Px = PiPlusHierarchy_MCP_Px_eachpoint->at(mcp_idx);
  std::vector<double> Py = PiPlusHierarchy_MCP_Py_eachpoint->at(mcp_idx);
  std::vector<double> Pz = PiPlusHierarchy_MCP_Pz_eachpoint->at(mcp_idx);

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


std::vector<int> PiRecoTree::FindRecoKinks(int pfp_idx, double mincosth, int evtnum){

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

  // if (mincosth>0.9 && kinkidxs.size()>0){
    // PlotLocalLinearityDetails(ordered_spacepoints,proj1d,cusum,ttest,kinkidxs,pfp_idx,PFP_track_start->at(pfp_idx),mincosth,evtnum);
  // }

  return kinkidxs;
};


std::vector<int> PiRecoTree::FindRecoKinks_byangle(int pfp_idx, double mincosth, int evtnum, double kinkcosth, int _slider_window=20){

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
