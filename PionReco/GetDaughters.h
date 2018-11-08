#ifndef GETDAUGHTERS_H
#define GETDAUGHTERS_H

#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"

#include "nusimdata/SimulationBase/MCParticle.h"

// Get primary MCParticle (i.e. the MCParticle that is a direct daughter of the neutrino) for any given MCParticle
art::Ptr<simb::MCParticle> MCP_GetPrimary(art::Ptr<simb::MCParticle> test_mcp, std::vector<art::Ptr<simb::MCParticle>> all_mcps, int nu_MCPID);

// Calculate std::map linking MCParticles to the primary MCParticle in that hierarchy
void MCP_CalcFullHierarchyMap(std::vector<art::Ptr<simb::MCParticle>> all_mcps, int nu_MCPID, std::map<art::Ptr<simb::MCParticle>,art::Ptr<simb::MCParticle>> &mcps_to_primaries);

std::vector<art::Ptr<simb::MCParticle>> MCP_GetFullHierarchyFromPrimary(std::vector<art::Ptr<simb::MCParticle>> all_mcps, art::Ptr<simb::MCParticle> primary_mcp, int nu_MCPID);

std::vector<art::Ptr<simb::MCParticle>> MCP_GetLogicalPion(art::Ptr<simb::MCParticle> test_mcp, std::vector<art::Ptr<simb::MCParticle>> all_mcps, bool fUseMuonsInstead);

#endif
