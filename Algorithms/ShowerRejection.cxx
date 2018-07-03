#ifndef SHOWERREJECTION_CXX
#define SHOWERREJECTION_CXX

#include "ShowerRejection.h"

std::map<std::string,bool> ShowerRejection(art::Event &evt, InputTags *CC1piInputTags)
{
   std::map<std::string,bool> CC1picutflow;

   art::Handle<std::vector<ubana::SelectionResult>> selection_h;
   evt.getByLabel(CC1piInputTags->fSelectionLabel, selection_h);
   if(!selection_h.isValid()){
      mf::LogError(__PRETTY_FUNCTION__) << "[ShowerRejection] SelectionResult product not found." << std::endl;
      throw std::exception();
   }

   //Get TPCObject
   art::FindManyP<ubana::TPCObject> tpcobject_from_selection(selection_h, evt, CC1piInputTags->fSelectionLabel);
   if(tpcobject_from_selection.at(0).size() == 0) {
      //No TPCObject
      CC1picutflow["residuals_std_up"] = false;
      CC1picutflow["residuals_std_down"] = false;
      CC1picutflow["residuals_mean_up"] = false;
      CC1picutflow["residuals_mean_down"] = false;
      CC1picutflow["perc_used_hits_in_cluster"] = false;
      return CC1picutflow;
   }
   art::Ptr<ubana::TPCObject> tpcobj_candidate = tpcobject_from_selection.at(0).at(0);
   art::Handle<std::vector<ubana::TPCObject>> tpcobj_h;
   evt.getByLabel("TPCObjectMaker", tpcobj_h);

   //Get PFPs (in TPCObject)
   art::FindManyP<recob::PFParticle> pfps_from_tpcobject(tpcobj_h, evt, "TPCObjectMaker");
   std::vector<art::Ptr<recob::PFParticle>> pfps = pfps_from_tpcobject.at(tpcobj_candidate.key());

   //Get PFPs (in event)
   art::Handle<std::vector<recob::PFParticle> > pfp_h;
   evt.getByLabel(CC1piInputTags->fTrackLabel,pfp_h);
   if(!pfp_h.isValid()){
      mf::LogError(__PRETTY_FUNCTION__) << "PFP product not found." << std::endl;
      throw std::exception();
   }

   art::FindManyP<recob::Track> tracks_from_pfps(pfp_h, evt, CC1piInputTags->fTrackLabel);

   // per track cut results
   std::vector<bool> residuals_std_up;
   std::vector<bool> residuals_std_down;
   std::vector<bool> residuals_mean_up;
   std::vector<bool> residuals_mean_down;
   std::vector<bool> perc_used_hits_in_cluster;

   for (auto pfp : pfps) {
      std::vector<art::Ptr<recob::Track>> tracks_pfp = tracks_from_pfps.at(pfp.key());
      if(tracks_pfp.size() > 1) {
         mf::LogError(__PRETTY_FUNCTION__) << "PFP associated to more than one track." << std::endl;
         throw std::exception();
      }
      for (auto track : tracks_pfp){
         std::pair<double, double> residual_mean_std;
         double ratio;
         ShowerCheck(evt, CC1piInputTags, pfp, track, residual_mean_std, ratio);

         double residual_mean = residual_mean_std.first;
         double residual_std  = residual_mean_std.second;

         // Hit Residuals STD Cut
         if (residual_std > CC1piInputTags->fResidualsStdCutUp) {
            residuals_std_up.emplace_back(false);
         } else {
            residuals_std_up.emplace_back(true);
         }
         if (residual_std < CC1piInputTags->fResidualsStdCutDown) {
            residuals_std_down.emplace_back(false);
         } else {
            residuals_std_down.emplace_back(true);
         }

         // Hit Residuals Mean Cut
         if (residual_mean > CC1piInputTags->fResidualsMeanCutUp) {
            residuals_mean_up.emplace_back(false);
         } else {
            residuals_mean_up.emplace_back(true);
         }
         if (residual_mean < CC1piInputTags->fResidualsMeanCutDown) {
            residuals_mean_down.emplace_back(false);
         } else {
            residuals_mean_down.emplace_back(true);
         }

         // Percental of used hits in cluster Cut
         if (ratio < CC1piInputTags->fPercUsedHitsCut) {
            perc_used_hits_in_cluster.emplace_back(false); 
         } else {
            perc_used_hits_in_cluster.emplace_back(true); 
         }

      }
   }

   // Hit Residuals STD Cut
   if (std::count(residuals_std_up.begin(),residuals_std_up.end(),false) > 0) {
      CC1picutflow["residuals_std_up"] = false;
   } else {
      CC1picutflow["residuals_std_up"] = true;
   }
   if (std::count(residuals_std_down.begin(),residuals_std_down.end(),false) > 0) {
      CC1picutflow["residuals_std_down"] = false;
   } else {
      CC1picutflow["residuals_std_down"] = true;
   }

   // Hit Residuals Mean Cut
   if (std::count(residuals_mean_up.begin(),residuals_mean_up.end(),false) > 0) {
      CC1picutflow["residuals_mean_up"] = false;
   } else {
      CC1picutflow["residuals_mean_up"] = true;
   }
   if (std::count(residuals_mean_down.begin(),residuals_mean_down.end(),false) > 0) {
      CC1picutflow["residuals_mean_down"] = false;
   } else {
      CC1picutflow["residuals_mean_down"] = true;
   }

   // Percental of used hits in cluster Cut
   if (std::count(perc_used_hits_in_cluster.begin(),perc_used_hits_in_cluster.end(),false) > 0) {
      CC1picutflow["perc_used_hits_in_cluster"] = false; 
   } else {
      CC1picutflow["perc_used_hits_in_cluster"] = true; 
   }

   return CC1picutflow;
}

void ShowerCheck(art::Event &evt, InputTags *CC1piInputTags, art::Ptr<recob::PFParticle> pfp, art::Ptr<recob::Track> track, std::pair<double, double> &residual_mean_std, double &ratio){

   // Detector info service
   ::detinfo::DetectorProperties const* fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();

   // Geometry service
   ::art::ServiceHandle<geo::Geometry> geo;

   // Get track-hit associations
   lar_pandora::TrackVector allPfParticleTracks;
   lar_pandora::TracksToHits trackToHitsMap;
   lar_pandora::LArPandoraHelper::CollectTracks(evt, (CC1piInputTags->fPFParticleProducer).label(), allPfParticleTracks, trackToHitsMap);

   // Get PFP-spacepoint associations
   lar_pandora::PFParticlesToHits recoParticlesToHits;
   lar_pandora::HitsToPFParticles recoHitsToParticles;
   lar_pandora::LArPandoraHelper::BuildPFParticleHitMaps(evt, (CC1piInputTags->fPFParticleProducer).label(), (CC1piInputTags->fSpacePointProducer).label(), recoParticlesToHits, recoHitsToParticles, lar_pandora::LArPandoraHelper::kUseDaughters, true);

   //
   // Look at residuals
   //

   std::vector<TVector3> hit_v; // a vec of hits from coll plane
   std::vector<TVector3> track_v; // a vec of hits from coll plane

   // Collect hits
   auto iter = trackToHitsMap.find(track);
   if (iter != trackToHitsMap.end()) {
      std::vector<art::Ptr<recob::Hit>> hits = iter->second;
      for (auto hit : hits) {
         if (hit->View() == 2) {
            TVector3 h (hit->WireID().Wire, hit->PeakTime(), 0);
            //std::cout << "emplacing hit with wire " << h.X() << ", and time " << h.Y() << std::endl;
            hit_v.emplace_back(h);
         }
      }
   }
   else {
      mf::LogError(__PRETTY_FUNCTION__) << "Track not found in track to hits map." << std::endl;
      throw std::exception();
   }

   // Collect track points
   for (size_t i = 0; i < track->NumberTrajectoryPoints(); i++) {
      try {
         if (track->HasValidPoint(i)) {
            TVector3 trk_pt = track->LocationAtPoint(i);
            double wire = geo->NearestWire(trk_pt, 2);
            double time = fDetectorProperties->ConvertXToTicks(trk_pt.X(), geo::PlaneID(0,0,2));
            TVector3 p (wire, time, 0.);
            //std::cout << "emplacing track point on wire " << p.X() << ", and time " << p.Y() << std::endl;
            track_v.emplace_back(p);
         }
      } catch (...) {
         continue;
      }
   }

   ubana::TrackQuality _track_quality;
   _track_quality.SetTrackPoints(track_v);
   _track_quality.SetHitCollection(hit_v);
   residual_mean_std = _track_quality.GetResiduals();

   double n_hits_in_cluster = 0;
   auto it = recoParticlesToHits.find(pfp);
   if (it != recoParticlesToHits.end()) {
      for (auto h : it->second) {
         if (h->View() == 2) {
            n_hits_in_cluster++;
         }
      }
   }
   else {
      mf::LogError(__PRETTY_FUNCTION__) << "PFP not found in PFP to spacepoints map." << std::endl;
      throw std::exception();
   }

   ratio = (double)hit_v.size()/n_hits_in_cluster;

}

#endif
