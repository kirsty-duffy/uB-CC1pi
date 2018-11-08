#ifndef GETDAUGHTERS_CXX
#define GETDAUGHTERS_CXX

#include "GetDaughters.h"

// Get primary MCParticle (i.e. the MCParticle that is a direct daughter of the neutrino) for any given MCParticle
art::Ptr<simb::MCParticle> MCP_GetPrimary(art::Ptr<simb::MCParticle> test_mcp, std::vector<art::Ptr<simb::MCParticle>> all_mcps, int nu_MCPID){
  // std::cout << "Calling MCP_GetPrimary" << std::endl;
  // Check if this MCParticle's mother is the neutrino. If so, return it
  if (test_mcp->Mother() == nu_MCPID) return test_mcp;

  // If not, get the particle's mother and call this function again
  art::Ptr<simb::MCParticle> mother_mcp;
  for (auto mcp : all_mcps){
    if (mcp->TrackId() != test_mcp->Mother()) continue;
    mother_mcp = mcp;
  }
  if (mother_mcp.isNull()){
    // If it couldn't find the mother, return a null value
    return mother_mcp;
  }

  return MCP_GetPrimary(mother_mcp, all_mcps, nu_MCPID);
}

// Calculate std::map linking MCParticles to the primary MCParticle in that hierarchy
void MCP_CalcFullHierarchyMap(std::vector<art::Ptr<simb::MCParticle>> all_mcps, int nu_MCPID, std::map<art::Ptr<simb::MCParticle>,art::Ptr<simb::MCParticle>> &mcps_to_primaries){
  // std::cout << "Calling MCP_CalcFullHierarchyMap" << std::endl;
  for (auto mcp : all_mcps){
    art::Ptr<simb::MCParticle> primary = MCP_GetPrimary(mcp, all_mcps, nu_MCPID);
    mcps_to_primaries.insert(std::pair<art::Ptr<simb::MCParticle>,art::Ptr<simb::MCParticle>>(mcp,primary));
  }
}

std::vector<art::Ptr<simb::MCParticle>> MCP_GetFullHierarchyFromPrimary(std::vector<art::Ptr<simb::MCParticle>> all_mcps, art::Ptr<simb::MCParticle> primary_mcp, int nu_MCPID){
  // std::cout << "Calling MCP_GetFullHierarchyFromPrimary" << std::endl;
  std::map<art::Ptr<simb::MCParticle>,art::Ptr<simb::MCParticle>> mcps_to_primaries;
  MCP_CalcFullHierarchyMap(all_mcps, nu_MCPID, mcps_to_primaries);

  std::vector<art::Ptr<simb::MCParticle>> hierarchy;
  for (std::map<art::Ptr<simb::MCParticle>,art::Ptr<simb::MCParticle>>::iterator itr = mcps_to_primaries.begin(); itr != mcps_to_primaries.end(); ++itr){
    if (itr->second == primary_mcp){
      hierarchy.push_back(itr->first);
    }
  }

  return hierarchy;
}

// Construct "logical pion": a vector of MCParticles that contains all pions that are true daughters of the preceeding pion (combine together pions from elastic scatters that GEANT splits into two different particles)
std::vector<art::Ptr<simb::MCParticle>> MCP_GetLogicalPion(art::Ptr<simb::MCParticle> test_mcp, std::vector<art::Ptr<simb::MCParticle>> all_mcps, bool fUseMuonsInstead){
  std::vector<art::Ptr<simb::MCParticle>> logicalpion;
  bool keepgoing=true;

  while (keepgoing==true){
    // Start by filling this MCP into the vector
    logicalpion.push_back(test_mcp);

    // Check daughters of test_mcp. If exactly one is a pi+, fill that into the vector and look at its daughters
    art::Ptr<simb::MCParticle> piplus_daughter;
    int n_piplus_daughters = 0;

    for (int i_daughter=0; i_daughter<test_mcp->NumberDaughters(); i_daughter++){
      int daughter_TrackId = test_mcp->Daughter(i_daughter);
      art::Ptr<simb::MCParticle> daughter_mcp;
      for (auto mcp : all_mcps){
        if (mcp->TrackId() == daughter_TrackId){
          daughter_mcp = mcp;
          break;
        }
      }
      if (daughter_mcp.isNull()){
        //std::cout << "Could not find daughter MCP " << daughter_TrackId << " in MCParticle list" << std::endl;
        continue;
      }
      // Now we have the daughter, check it's PDG code
      if (fUseMuonsInstead){
        if (daughter_mcp->PdgCode() == 13){
          n_piplus_daughters++;
          piplus_daughter = daughter_mcp;
        }
      }
      else{
        if (daughter_mcp->PdgCode() == 211){
          n_piplus_daughters++;
          piplus_daughter = daughter_mcp;
        }
      }

    } // end loop over daughters

    // Now we have been through all the daughters, how many were pi+? If exactly one, make this the next test_mcp and keep going. If not exactly one, set keepgoing to false and break out of the loop
    if (n_piplus_daughters == 1){
      test_mcp = piplus_daughter;
    }
    else{
      keepgoing = false;
    }
  }

  return logicalpion;
}

#endif
