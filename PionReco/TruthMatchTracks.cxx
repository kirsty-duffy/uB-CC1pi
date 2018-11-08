/**
 * Get true PDG from associations.
 * We do this by using hit <-> MCParticle associations, looping over the
 * hits and finding the MCParticle which contributed the most charge
 * to each hit.
 *
 * .key() is used to get the index in the original collection
 */

 #ifndef TRUTHMATCH_CXX
 #define TRUTHMATCH_CXX

 #include "TruthMatchTracks.h"

 void Do_Track_MCP_Matching_BestMatch(const art::Event &e, art::Handle<std::vector<recob::Track>> trackHandle, art::Handle<std::vector<recob::Hit>> hitHandle, std::string fHitTrackAssns, std::string fHitTruthAssns, std::map<int, MatchingData> &trackID_to_matchingData, std::map<int, MatchingData> &MCPID_to_matchingData){

   // First: get all tracks in a ptr vector
   std::vector<art::Ptr<recob::Track>> trackCollection;
   art::fill_ptr_vector(trackCollection, trackHandle);

   // Get track<->hit associations
   art::FindManyP<recob::Hit>
    hits_from_tracks(trackHandle, e, fHitTrackAssns);

   // Get hit<->MCP associations
   art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitHandle,e,fHitTruthAssns);

   // Loop through all tracks and do matching
   for (auto& track : trackCollection){

     MatchingData matchingData;

     std::unordered_map<int,double> trkide;
     double maxe=-1, tote=0;
     simb::MCParticle const* matched_mcparticle = NULL;

     std::vector<simb::MCParticle const*> particle_vec;
     std::vector<anab::BackTrackerHitMatchingData const*> match_vec;

     std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(track->ID());

     for(size_t i_h=0; i_h<hits_from_track.size(); i_h++){
       particle_vec.clear(); match_vec.clear();
       particles_per_hit.get(hits_from_track[i_h].key(),particle_vec,match_vec);

       for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
         trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy;
         tote += match_vec[i_p]->energy;
         if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){
           maxe = trkide[ particle_vec[i_p]->TrackId() ];
           matched_mcparticle = particle_vec[i_p];
         }
       }//end loop over particles per hit

     } // end loop over hits in track

     matchingData.recoObject_totalEnergy = tote;
     matchingData.MatchedEnergy = maxe;
     matchingData.MCParticle_ID = matched_mcparticle->TrackId();
     matchingData.MCParticle_PDG = matched_mcparticle->PdgCode();
     matchingData.recoObject_ID = track->ID();

     trackID_to_matchingData[track->ID()] = matchingData;
     MCPID_to_matchingData[matched_mcparticle->TrackId()] = matchingData;

   } // end loop over tracks in trackCollection
 }

 // ----------------------------------------------------- //

 void Do_PFP_MCP_Matching_BestMatch(const art::Event &e, art::Handle<std::vector<recob::PFParticle>> pfpHandle, art::Handle<std::vector<recob::Hit>> hitHandle, std::string fPFPproducerlabel, std::string fHitTruthAssns, std::map<int, MatchingData> &pfpID_to_matchingData, std::map<int, MatchingData> &MCPID_to_matchingData){

   // First: get all PFPs in a ptr vector
   std::vector<art::Ptr<recob::PFParticle>> pfpCollection;
   art::fill_ptr_vector(pfpCollection, pfpHandle);

   // Get PFP<->Cluster associations
   art::FindManyP<recob::Cluster> clusters_from_pfps(pfpHandle, e, fPFPproducerlabel);

   // Get Cluster<->hit associations
   art::Handle<std::vector<recob::Cluster>> clusterHandle;
   e.getByLabel(fPFPproducerlabel, clusterHandle);
   art::FindManyP<recob::Hit> hits_from_clusters(clusterHandle, e, fPFPproducerlabel);

   // Get hit<->MCP associations
   art::FindMany<simb::MCParticle,anab::BackTrackerHitMatchingData> particles_per_hit(hitHandle,e,fHitTruthAssns);

   // Loop through all PFPs and do matching
   for (auto& pfp : pfpCollection){
     MatchingData matchingData;

     std::unordered_map<int,double> trkide;
     double maxe=-1, tote=0;
     simb::MCParticle const* matched_mcparticle = NULL;

     std::vector<simb::MCParticle const*> particle_vec;
     std::vector<anab::BackTrackerHitMatchingData const*> match_vec;

     // Get clusters from PFP
     std::vector<art::Ptr<recob::Cluster>> clusters_from_pfp = clusters_from_pfps.at(pfp.key());

     // Loop through clusters. For each cluster, get the hits associated to that cluster, and add them into a vector hits_from_pfp
     std::vector<art::Ptr<recob::Hit>> hits_from_pfp;
     for (auto cluster : clusters_from_pfp){
        std::vector<art::Ptr<recob::Hit>> hits_from_cluster = hits_from_clusters.at(cluster.key());
        hits_from_pfp.insert(hits_from_pfp.end(), hits_from_cluster.begin(), hits_from_cluster.end());
     }

     for(size_t i_h=0; i_h<hits_from_pfp.size(); i_h++){
      particle_vec.clear(); match_vec.clear();
      particles_per_hit.get(hits_from_pfp[i_h].key(),particle_vec,match_vec);

      for(size_t i_p=0; i_p<particle_vec.size(); ++i_p){
         trkide[ particle_vec[i_p]->TrackId() ] += match_vec[i_p]->energy;
         tote += match_vec[i_p]->energy;
         if( trkide[ particle_vec[i_p]->TrackId() ] > maxe ){
           maxe = trkide[ particle_vec[i_p]->TrackId() ];
           matched_mcparticle = particle_vec[i_p];
         }
      }//end loop over particles per hit

     } // end loop over hits in track

     matchingData.recoObject_totalEnergy = tote;
     matchingData.MatchedEnergy = maxe;
     matchingData.recoObject_ID = pfp->Self();
     if (matched_mcparticle != NULL){
        matchingData.MCParticle_ID = matched_mcparticle->TrackId();
        matchingData.MCParticle_PDG = matched_mcparticle->PdgCode();
     }

     pfpID_to_matchingData[pfp->Self()] = matchingData;
     if (matched_mcparticle != NULL)
      MCPID_to_matchingData[matched_mcparticle->TrackId()] = matchingData;

  } // end loop over pfps in pfpCollection
 }

 // ----------------------------------------------------- //

 void Do_PFP_MCP_Matching_AllMatches(const art::Event &e, art::Handle<std::vector<simb::MCParticle>> mcpHandle, art::Handle<std::vector<recob::PFParticle>> pfpHandle, art::Handle<std::vector<recob::Hit>> hitHandle, std::string fPFPproducerlabel, std::string fHitTruthAssns, std::map<int, std::vector<MatchingData>> &MCPID_to_many_matchingData){
   // Get all MCPs in a ptr vector, mcpCollection
   std::vector<art::Ptr<simb::MCParticle>> mcpCollection;
   art::fill_ptr_vector(mcpCollection, mcpHandle);

   // Now get MCP<->hit associations
   art::FindManyP<recob::Hit,anab::BackTrackerHitMatchingData> hits_per_mcp(mcpHandle, e, fHitTruthAssns);

   // Also get hit<->cluster associations
   art::FindManyP<recob::Cluster> clusters_from_hits(hitHandle, e, fPFPproducerlabel);

   // And cluster<->pfp associations
   art::Handle<std::vector<recob::Cluster>> clusterHandle;
   e.getByLabel(fPFPproducerlabel, clusterHandle);
   art::FindManyP<recob::PFParticle> pfps_from_clusters(clusterHandle, e, fPFPproducerlabel);

   // Loop through all MCPs. For each MCP construct a vector of MatchingData, where each entry in the vector corresponds to a cluster/PFParticle that shares the same hits
   for (auto& mcp : mcpCollection){
      std::vector<art::Ptr<recob::Hit>> hit_vec;
      std::vector<anab::BackTrackerHitMatchingData const*> match_vec;

      std::map<int,MatchingData> pfpID_to_matchingData;

      // Loop through hits associated to this MCP
      hits_per_mcp.get(mcp.key(),hit_vec,match_vec);

     for(size_t i_h=0; i_h<hit_vec.size(); i_h++){
        // Find all clusters associated to this hit
        std::vector<art::Ptr<recob::Cluster>> clusters_from_hit = clusters_from_hits.at(hit_vec[i_h].key());

        // Loop through clusters associated to this hit
        for (size_t i_cl=0; i_cl<clusters_from_hit.size(); i_cl++){
         // Find all PFPs associated to this cluster
         std::vector<art::Ptr<recob::PFParticle>> pfps_from_cluster = pfps_from_clusters.at(clusters_from_hit[i_cl].key());

         // Now loop over all PFPs and fill matching data vector
         for (size_t i_pfp=0; i_pfp<pfps_from_cluster.size(); i_pfp++){
            int recoObject_ID = pfps_from_cluster.at(i_pfp).key();
            // int cluster_ID = clusters_from_hit.at(i_cl).key();
            int MCParticle_ID = mcp->TrackId();
            double matchedE = match_vec.at(i_h)->energy;

            // Check if this PFP has been found already. If it has, just add the energy. If it hasn't, make a new matchingData object
            auto search = pfpID_to_matchingData.find(recoObject_ID);
            if (search != pfpID_to_matchingData.end()){
               // Found entry already: just add energy
               search->second.MatchedEnergy += matchedE;
            }
            else{
               // No entry found: make a new one
               MatchingData tmp;
               tmp.recoObject_ID = recoObject_ID;
               tmp.MCParticle_ID = MCParticle_ID;
               tmp.MCParticle_PDG = mcp->PdgCode();
               tmp.recoObject_totalEnergy = -999;
               tmp.MatchedEnergy = matchedE;
               pfpID_to_matchingData.insert(std::pair<int,MatchingData>(recoObject_ID,tmp));
            }

         } // end loop over pfps
        } // end loop over clusters
     } // end loop over hits

     // Now we have found all PFPs associated to this MCP, make a vector of matchingDatas and insert it into our map
     std::vector<MatchingData> matchingData;
     for (std::map<int,MatchingData>::iterator itr = pfpID_to_matchingData.begin(); itr != pfpID_to_matchingData.end(); ++itr){
        matchingData.push_back(itr->second);
     }
     // Finally: insert this into the MCPid<->std::vector<matchingData> map that we want to populate
     MCPID_to_many_matchingData.insert(std::pair<int,std::vector<MatchingData>>(mcp->TrackId(),matchingData));


  } // End loop over MCPs in mcpCollection
 }



 #endif
