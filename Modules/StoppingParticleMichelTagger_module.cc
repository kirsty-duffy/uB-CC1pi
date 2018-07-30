////////////////////////////////////////////////////////////////////////
// Class:       StoppingParticleMichelTagger
// Plugin Type: producer (art v2_05_00)
// File:        StoppingParticleMichelTagger_module.cc
//
// Generated at Mon July 30 10:18:20 2018 by Kirsty Duffy using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

/**
 * \class StoppingParticleMichelTagger
 *
 * \brief Art producer module that tags stopping and Michel-decaying particles. Modeled on M. Del Tutto's StoppingMuonTagger module, but intended to be applied to tracks produced in neutrino events (as opposed to cosmic muons)
 *
 *
 * \author Kirsty Duffy <kduffy@fnal.gov>
 *
 * \version producer (art v2_05_00)
 *
 * \date 2018/07/30
 *
 * Contact: kduffy@fnal.gov
 *
 * Created on: Mon July 30 10:18:20 2018
 *
 */

#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "art/Framework/Services/Optional/TFileDirectory.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "lardata/Utilities/AssociationUtil.h"

#include <memory>
#include <limits>

// Data products include
#include "lardataobj/MCBase/MCTrack.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Vertex.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "uboone/UBXSec/DataTypes/TPCObject.h"
#include "uboone/UBXSec/DataTypes/SelectionResult.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// LArSoft include
#include "larcore/Geometry/Geometry.h"

// CC1pi method includes
#include "uboone/CC1pi/Algorithms/InputTags.h"

// UBXSec Algorithms include
#include "uboone/LLBasicTool/GeoAlgo/GeoTrajectory.h"
#include "uboone/UBXSec/Algorithms/FiducialVolume.h"
#include "uboone/UBXSec/HitCosmicTag/Base/DataTypes.h"
#include "uboone/UBXSec/HitCosmicTag/Base/CosmicTagManager.h"
#include "uboone/UBXSec/HitCosmicTag/Algorithms/StopMuMichel.h"
#include "uboone/UBXSec/HitCosmicTag/Algorithms/StopMuBragg.h"
#include "uboone/UBXSec/HitCosmicTag/Algorithms/CosmicSimpleMIP.h"

// Root include
#include "TString.h"
#include "TTree.h"

const std::vector<float> endPt1 = {-9999., -9999., -9999.};
const std::vector<float> endPt2 = {-9999., -9999., -9999.};

class StoppingParticleMichelTagger;


class StoppingParticleMichelTagger : public art::EDProducer {
public:
  explicit StoppingParticleMichelTagger(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  StoppingParticleMichelTagger(StoppingParticleMichelTagger const &) = delete;
  StoppingParticleMichelTagger(StoppingParticleMichelTagger &&) = delete;
  StoppingParticleMichelTagger & operator = (StoppingParticleMichelTagger const &) = delete;
  StoppingParticleMichelTagger & operator = (StoppingParticleMichelTagger &&) = delete;

  // Required functions.
  void produce(art::Event & e) override;

private:

  ::cosmictag::CosmicTagManager _ct_manager;

  ::art::ServiceHandle<geo::Geometry> geo;
  ::detinfo::DetectorProperties const* fDetectorProperties;

  // ::ubana::FiducialVolume _fiducial_volume;

  InputTags *CC1piInputTags;

  double _coplanar_cut = 5.; ///< If a track start X and end X difference is below this value, is consiered to be coplanar to a collection plane wire

  bool _debug; ///< Debug flag
  bool _create_tree; ///< If true, creates a tree with info to make plots

  TTree* _tree1;
};


StoppingParticleMichelTagger::StoppingParticleMichelTagger(fhicl::ParameterSet const & p)
  {


  _ct_manager.Configure(p.get<cosmictag::Config_t>("CosmicTagManager"));

  // _fiducial_volume.Configure(p.get<fhicl::ParameterSet>("FiducialVolumeSettings"),
  //                            geo->DetHalfHeight(),
  //                            2.*geo->DetHalfWidth(),
  //                            geo->DetLength());

  // std::cout << "[StoppingParticleMichelTagger] FV: " << std::endl;
  // _fiducial_volume.PrintConfig();

  CC1piInputTags = new InputTags(p);

  // _tpcobject_producer = p.get<std::string>("TPCObjectProducer",  "TPCObjectMaker::CC1pi");
  // _pfp_producer       = p.get<std::string>("PFParticleProducer", "pandoraNu::CC1pi");
  // _cluster_producer   = p.get<std::string>("ClusterProducer",    "pandoraNu::CC1pi");
  // _track_producer     = p.get<std::string>("TrackProducer",      "pandoraNu::CC1pi");

  _coplanar_cut       = p.get<double>("CoplanarCut",   5.);

  _debug = p.get<bool>("DebugMode", true);
  _create_tree = p.get<bool>("CreateTree", true);

  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();

  art::ServiceHandle<art::TFileService> fs;
  if (_create_tree) {
    _tree1 = fs->make<TTree>("StoppingTaggerTree","");
    // _tree1->Branch("run",           &_run,                 "run/I");
    // _tree1->Branch("subrun",        &_subrun,              "subrun/I");
    // _tree1->Branch("event",         &_event,               "event/I");
    // _tree1->Branch("origin",        &_origin,              "origin/I");
    // _tree1->Branch("origin_extra",  &_origin_extra,        "origin_extra/I");
    // _tree1->Branch("fv",            &_fv,                  "fv/O");
    // _tree1->Branch("length",        &_length,              "length/D");
  }

}



void StoppingParticleMichelTagger::produce(art::Event & evt) {

  if (_debug) std::cout <<"[StoppingParticleMichelTagger] Starts." << std::endl;

  if(_create_tree) {
    // _run = evt.id().run();
    // _subrun = evt.id().subRun();
    // _event = evt.id().event();
  }

  // Get result of Marco's selection
  art::Handle<std::vector<ubana::SelectionResult>> selection_h;
  evt.getByLabel(CC1piInputTags->fSelectionLabel, selection_h);
  if(!selection_h.isValid()){
     mf::LogError(__PRETTY_FUNCTION__) << "[FillTree] SelectionResult product not found." << std::endl;
     throw std::exception();
  }
  std::vector<art::Ptr<ubana::SelectionResult>> selection_v;
  art::fill_ptr_vector(selection_v, selection_h);

  // Skip events that are not selected
  // Actually maybe not - will this mess up the alignment of the trees?
  // if (!(selection_v.at(0)->GetSelectionStatus())) continue;

  // Get selected CCnumu candidate TPCObject (there should be only one)
  art::FindManyP<ubana::TPCObject> TPCObject_from_selection(selection_h, evt, CC1piInputTags->fSelectionLabel);
  art::Ptr<ubana::TPCObject> TPCObj_candidate = TPCObject_from_selection.at(0).at(0);
  art::Handle<std::vector<ubana::TPCObject>> TPCObj_h;
  evt.getByLabel(CC1piInputTags->fTPCObjectProducer, TPCObj_h);

  // Get tracks and pfparticles from TPCObject
  art::FindManyP<recob::Track> tracks_from_TPCObject(TPCObj_h, evt, CC1piInputTags->fTPCObjectProducer);
  std::vector<art::Ptr<recob::Track>> tracks = tracks_from_TPCObject.at(TPCObj_candidate.key());

  art::FindManyP<recob::PFParticle> pfps_from_TPCObject(TPCObj_h, evt, CC1piInputTags->fTPCObjectProducer);
  std::vector<art::Ptr<recob::PFParticle>> pfps = pfps_from_TPCObject.at(TPCObj_candidate.key());

  // Get reconstructed neutrino vertex position from TPCObject
  // Actually, no: use the spacepoint closest to the start of the track (if a secondary track is backwards going it might have the wrong end closest to the nu vertex)
  // recob::Vertex TPCObj_nu_vtx = TPCObj_candidate->GetVertex();
  // const double *reco_nu_vtx = TPCObj_nu_vtx->XYZ();
  // if (_debug) std::cout << "[StoppingParticleMichelTagger] Looking at neutrino vertex (" << reco_nu_vtx[0] << "," << reco_nu_vtx[1] << "," << reco_nu_vtx[3] << ")" << std::endl;

  // Map track->pfp, pfp->spacepoint, and pfp->cluster

  // Collect pfparticles and map pfp->spacepoint
  // lar_pandora::PFParticleVector pfp_v;
  // lar_pandora::PFParticlesToSpacePoints pfps_to_spacepoints;
  // lar_pandora::LArPandoraHelper::CollectPFParticles(evt, CC1piInputTags->fPFParticleProducer, pfp_v, pfps_to_spacepoints);

  // Collect pfparticles and map pfp->cluster
  lar_pandora::PFParticleVector pfp_v;
  lar_pandora::PFParticlesToClusters pfps_to_clusters;
  lar_pandora::LArPandoraHelper::CollectPFParticles(evt, CC1piInputTags->fPFParticleProducer.label(), pfp_v, pfps_to_clusters);

  // Collect clusters and map cluster->hit
  lar_pandora::ClusterVector cluster_v;
  lar_pandora::ClustersToHits clusters_to_hits;
  lar_pandora::LArPandoraHelper::CollectClusters(evt, CC1piInputTags->fClusterProducer.label(), cluster_v, clusters_to_hits);

  // Collect tracks and map pfp->track
  // Can't do the same thing for tracks in the CC1pi analysis because we hacked Pandora to make it produce a track for every PFP. The way we made sure this didn't affect Marco's selection was to make sure LArPandoraHelper::CollectTracks ignored our "extra" tracks that would originally have been showers. Therefore, if we want to look at those tracks, we can't use CollectTracks. Instead, I have copied and pasted the contents of LArPandoraHelper::CollectTracks below (without the lines we don't want)
  art::Handle< std::vector<recob::Track> > theTracks;
  evt.getByLabel(CC1piInputTags->fTrackLabel, theTracks);
  if (!theTracks.isValid())
  {
      std::cout << "[StoppingParticleMichelTagger]  Failed to find tracks... " << std::endl;
      return;
  }
  else
  {
      if (_debug) std::cout << "[StoppingParticleMichelTagger]  Found: " << theTracks->size() << " Tracks " << std::endl;
  }
  art::FindOneP<recob::PFParticle> theParticleAssns(theTracks, evt, CC1piInputTags->fTrackLabel);


  // ----- Now the setup is done! ----- //

  // Loop through tracks. For each track we want to run the various HitCosmicTag algorithms and store the output.
  for (auto track_i : tracks) {

  // TODO
  // Store whether the track is deemed to be coplanar (i.e. has start and end x positions within coplanar cut defined above, which would imply that the track is coplanar to a collection plane wire, and we probably don't have good calorimetry). Don't actually do anything with this information here, but store it so we can look at it later.
  // double deltax = track_i->Vertex().Z() - track_i->End().Z();
  // if (_debug) std::cout << "[StoppingParticleMichelTagger] Delta x is " << deltax << std::endl;
  // bool collection_is_coplanar = TMath::Abs(deltax) < _coplanar_cut;

  // Create an approximate start hit on plane 2 from the spacepoint closest to the reconstructed neutrino vertex candidate. Do this by getting the PFP associated with the track (there can only be one), and the clusters/hits/spacepoints associated to that PFP.
  // First, exclude spacepoints outside the TPC. Then get the point closest to the reconstructed neutrino vertex.
  // Turn into a cosmictag::SimpleHit for use with the HitCosmicTag Algorithms
  art::Ptr<recob::PFParticle> pfp_from_track = theParticleAssns.at(track_i.key());

  // std::vector<art::Ptr<recob::SpacePoint>> sp_v;
  // sp_v.clear();

  std::vector<art::Ptr<recob::Hit>> hit_v;
  hit_v.clear();

  // Fill vector of spacepoints
  // auto iter = pfps_to_spacepoints.find(pfp_from_track);
  // if (iter == pfps_to_spacepoints.end()){
  //   if (_debug) std::cout << "[StoppingParticleMichelTagger] Did not find spacepoints to match to PFP for track " << track_i.key() << std::endl;
  //   continue;
  // }
  // sp_v.reserve(sp_v.size() + iter->second.size());
  // sp_v.insert(sp_v.end(), iter->second.begin(), iter->second.end());

  // Find clusters first...
  auto iter2 = pfps_to_clusters.find(pfp_from_track);
  if (iter2 == pfps_to_clusters.end()){
    if (_debug) std::cout << "[StoppingParticleMichelTagger] Did not find clusters to match to PFP for track " << track_i.key() << std::endl;
    continue;
  }

  // ... then fill vector of hits
  for (auto c : iter2->second){
    auto iter3 = clusters_to_hits.find(c);
    if (iter3 == clusters_to_hits.end()){
      if (_debug) std::cout << "[StoppingParticleMichelTagger] Cluster in TPCObject not found by pandora?! " << track_i.key() << std::endl;
      throw std::exception();
    }
    hit_v.reserve(iter3->second.size());
    hit_v.insert(hit_v.begin(), iter3->second.begin(), iter3->second.end());
  }

  // Exclude spacepoints outside the TPC from sp_v
  // std::vector<art::Ptr<recob::SpacePoint>> temp;
  // ::geoalgo::AABox tpcvol(0., (-1.)*geo->DetHalfHeight(), 0., geo->DetHalfWidth()*2, geo->DetHalfHeight(), geo->DetLength());
  //
  // for (auto s : sp_v) {
  //   const double *xyz = s->XYZ();
  //   ::geoalgo::Vector point (xyz[0], xyz[1], xyz[2]);
  //   if (tpcvol.Contain(point)) {
  //     temp.push_back(s);
  //   }
  // }
  // sp_v = temp;

  // Now get the start of the track
  double first_point[3] = {track_i->Start().X(), track_i->Start().Y(), track_i->Start().Z()};

  // Create an approximate start hit on plane 2
  int first_w = geo->NearestWire(first_point, 2);
  double first_t = fDetectorProperties->ConvertXToTicks(first_point[0], geo::PlaneID(0,0,2))/4.;
  if (_debug) std::cout << "[StoppingParticleMichelTagger] First point: wire " << first_w << ", time " << first_t << std::endl;

  // Now make a cosmictag::SimpleHits object for the first point
  cosmictag::SimpleHit start_simple;
  start_simple.time = first_t;
  start_simple.wire = first_w;
  start_simple.plane = 2;

  // Finally, make a vector of cosmictag::SimpleHits for use with the HitCosmicTag algorithms (for hits in collection plane only)
  if (_debug) std::cout << "[StoppingParticleMichelTagger] Now create simple hit vector, size " << hit_v.size() << std::endl;
  std::vector<cosmictag::SimpleHit> simple_hit_v;
  for (auto h : hit_v){
    cosmictag::SimpleHit sh;

    sh.t = fDetectorProperties->ConvertTicksToX(h->PeakTime(), geo::PlaneID(0,0,h->View()));
    sh.w = h->WireID().Wire * geo->WirePitch(geo::PlaneID(0,0,h->View()));

    sh.plane = h->View();
    sh.integral = h->Integral();

    sh.time = h->PeakTime() / 4;
    sh.wire = h->WireID().Wire;

    // Emplace into vector only if collection plane
    if (h->View() == 2) {
      simple_hit_v.emplace_back(sh);
    }
  }

  if (_debug) std::cout << "[StoppingParticleMichelTagger] Simple hit vector size " << simple_hit_v.size() << std::endl;

  // Now get ready to run the HitCosmicTag algorithms.
  // Reset _ct_manager
  _ct_manager.Reset();

  // Give simple hits to the ct manager
  cosmictag::SimpleCluster sc(simple_hit_v);
  _ct_manager.Emplace(std::move(sc));

  // Give start hit to the ct manager
  _ct_manager.SetStartHit(std::move(start_simple));

  // Run the cluster analyser (i.e. the HitCosmicTag algorithms)
  bool passed = _ct_manager.Run();
  if (_debug) _ct_manager.PrintClusterStatus();
  cosmictag::SimpleCluster processed_cluster = _ct_manager.GetCluster();
  if(_debug) std::cout << "[StoppingParticleMichelTagger] Passed basic selection? " << (passed ? "YES" : "NO") << std::endl;

  bool ct_result_michel;
  bool ct_result_bragg;

  // Run the "custom" Michel tagging algorithm
  if (_debug) ((cosmictag::StopMuMichel*)(_ct_manager.GetCustomAlgo("StopMuMichel")))->PrintConfig();
  ct_result_michel = ((cosmictag::StopMuMichel*)(_ct_manager.GetCustomAlgo("StopMuMichel")))->IsStopMuMichel(processed_cluster);
  if(_debug) std::cout << "[StoppingParticleMichelTagger] Is michel? " << (ct_result_michel ? "YES" : "NO") << std::endl;

  // Run the "custom" Bragg peak finding algorithm
  if (_debug) ((cosmictag::StopMuBragg*)(_ct_manager.GetCustomAlgo("StopMuBragg")))->PrintConfig();
  ct_result_bragg = ((cosmictag::StopMuBragg*)(_ct_manager.GetCustomAlgo("StopMuBragg")))->IsStopMuBragg(processed_cluster);
  if(_debug) std::cout << "[StoppingParticleMichelTagger] Found stopping particle (bragg)? " << (ct_result_bragg ? "YES" : "NO") << std::endl;

  // Find some way to save all this information?

  } // End loop over recob::tracks

} // end ::produce



DEFINE_ART_MODULE(StoppingParticleMichelTagger)
