#include "uboone/CC1pi/Algorithms/GetUnusedHitCharges.h"

void GetUnusedHitCharges(art::Event &evt, InputTags *CC1piInputTags, int planenum, TVector3 trackend, /*std::vector<art::Ptr<recob::Hit>> &unused_hits,*/ std::vector<double> &unused_hit_charge, std::vector<int> &unused_hit_wiredist, std::vector<double> &unused_hit_timedist){

  bool debug=false;

  // Get hit collection
  art::Handle<std::vector<recob::Hit>> hit_h; evt.getByLabel(CC1piInputTags->fHitProducer,hit_h);
  std::vector<art::Ptr<recob::Hit>> hits;
  art::fill_ptr_vector(hits,hit_h);

  // Get hits associated to tracks
  art::Handle<std::vector<recob::Track>> track_h; evt.getByLabel(CC1piInputTags->fTrackLabel,track_h);
  std::vector<art::Ptr<recob::Track>> tracks;
  art::fill_ptr_vector(tracks,track_h);
  art::FindManyP<recob::Hit> hits_from_tracks(track_h,evt,CC1piInputTags->fTrackLabel);

  // Make new collection of only hits in the first collection that aren't associated to any tracks
  std::vector<art::Ptr<recob::Hit>> unmatched_hits;
  for (auto hit : hits){
    bool gotonexthit = false;

    // Only consider hits on the correct plane
    if (hit->View()!=planenum) continue;

    for (auto track : tracks){
      std::vector<art::Ptr<recob::Hit>> hits_from_track = hits_from_tracks.at(track.key());

      auto iter = std::find(hits_from_track.begin(),hits_from_track.end(),hit);
      if (iter != hits_from_track.end()) {
        //If hit is found in the track, continue to the next hit
        gotonexthit = true;
        break;
      }
    } // end loop over tracks

    if (gotonexthit) continue;
    // If we get to this point, the hit has not been matched to any of the tracks, so put it into our vector of unmatched hits
    unmatched_hits.push_back(hit);
  } // end loop over hits

  if (debug) std::cout << "Found " << unmatched_hits.size() << " unmatched hits on plane " << planenum << std::endl;

  // Now we have our vector of unmatched hits, calculate distance from the end of the track
  ::art::ServiceHandle<geo::Geometry> geo;
  ::detinfo::DetectorProperties const* fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>(); ;
  //geo::PlaneGeo const& plane = geo->Plane(planenum);
  //geo::WireID trkendwireID;
  /*try {
    trkendwireID = plane.NearestWireID(trackend);
  }
  catch (geo::InvalidWireError const& e) {
    if (!e.hasSuggestedWire()) throw;
    trkendwireID = plane.ClosestWireID(e.suggestedWireID());
  }*/
  unsigned int trkendwireID = geo->NearestWire((TVector3 const)trackend,planenum);

  double trkendtime = fDetectorProperties->ConvertXToTicks(trackend.X(), geo::PlaneID(0,0,planenum));

  if (debug) std::cout << "trkendwireID = " << trkendwireID << ", trkendtime = " << trkendtime << std::endl;

  for (auto hit : unmatched_hits){
    if (debug) std::cout << "- hit integral = " << hit->Integral() << ", wire = " << hit->WireID().Wire << ", time = " << hit->PeakTime() << std::endl;

    unused_hit_charge.push_back(hit->Integral());

    int wiredist = (int)hit->WireID().Wire-(int)trkendwireID;
    unused_hit_wiredist.push_back(wiredist);

    double timedist = hit->PeakTime()-trkendtime;
    unused_hit_timedist.push_back(timedist);

    if (debug) std::cout << "  -- dwire = " << wiredist << ", dt = " << timedist << std::endl;
  } // end loop over unmatched hits
} // end function
