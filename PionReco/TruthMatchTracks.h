/**
 * Get true PDG from associations.
 * We do this by using hit <-> MCParticle associations, looping over the
 * hits and finding the MCParticle which contributed the most charge
 * to each hit.
 *
 * .key() is used to get the index in the original collection
 */

 #ifndef TRUTHMATCH_H
 #define TRUTHMATCH_H

 // art includes
 #include "art/Framework/Principal/Event.h"
 #include "art/Framework/Principal/Handle.h"
 #include "canvas/Utilities/InputTag.h"
 #include "canvas/Persistency/Common/FindManyP.h"
 #include "canvas/Persistency/Common/FindMany.h"

 #include "nusimdata/SimulationBase/MCParticle.h"
 #include "lardataobj/RecoBase/Track.h"
 #include "lardataobj/RecoBase/Hit.h"
 #include "lardataobj/RecoBase/PFParticle.h"
 #include "lardataobj/RecoBase/Cluster.h"
 #include "lardataobj/AnalysisBase/BackTrackerMatchingData.h"

 struct MatchingData{
   int recoObject_ID;
   int MCParticle_ID;
   int MCParticle_PDG;
   double recoObject_totalEnergy;
   double MatchedEnergy;
};

void Do_Track_MCP_Matching_BestMatch(const art::Event &e, art::Handle<std::vector<recob::Track>> trackHandle, art::Handle<std::vector<recob::Hit>> hitHandle, std::string fHitTrackAssns, std::string fHitTruthAssns, std::map<int, MatchingData> &trackID_to_matchingData, std::map<int, MatchingData> &MCPID_to_matchingData);

void Do_PFP_MCP_Matching_BestMatch(const art::Event &e, art::Handle<std::vector<recob::PFParticle>> pfpHandle, art::Handle<std::vector<recob::Hit>> hitHandle, std::string fPFPproducerlabel, std::string fHitTruthAssns, std::map<int, MatchingData> &pfpID_to_matchingData, std::map<int, MatchingData> &MCPID_to_matchingData);

 void Do_PFP_MCP_Matching_AllMatches(const art::Event &e, art::Handle<std::vector<simb::MCParticle>> mcpHandle, art::Handle<std::vector<recob::PFParticle>> pfpHandle, art::Handle<std::vector<recob::Hit>> hitHandle, std::string fPFPproducerlabel, std::string fHitTruthAssns, std::map<int, std::vector<MatchingData>> &MCPID_to_many_matchingData);

 #endif
