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
#include <vector>

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

  TTree* _tree1;

  int _run;
  int _subrun;
  int _event;

  std::vector<double> TPCObj_PFP_track_length;
  std::vector<double> TPCObj_PFP_track_deltax;
  std::vector<bool> TPCObj_PFP_track_ct_passed_basic;
  std::vector<bool> TPCObj_PFP_track_ct_result_michel;
  std::vector<bool> TPCObj_PFP_track_ct_result_bragg;

  std::vector<std::vector<double>> TPCObj_PFP_track_SimpleCluster_hitTime;
  std::vector<std::vector<int>> TPCObj_PFP_track_SimpleCluster_hitWire;
  std::vector<std::vector<int>> TPCObj_PFP_track_SimpleCluster_hitPlane;
  std::vector<std::vector<double>> TPCObj_PFP_track_SimpleCluster_hitIntegral;
  std::vector<std::vector<double>> TPCObj_PFP_track_SimpleCluster_hitTimeTicks;
  std::vector<std::vector<int>> TPCObj_PFP_track_SimpleCluster_hitWireNo;
  std::vector<int> TPCObj_PFP_track_SimpleCluster_StartIndex;
  std::vector<std::vector<double>> TPCObj_PFP_track_SimpleCluster_hitdQds;
  std::vector<std::vector<double>> TPCObj_PFP_track_SimpleCluster_hitds;
  std::vector<std::vector<double>> TPCObj_PFP_track_SimpleCluster_hitdQdsSlider;
  std::vector<std::vector<double>> TPCObj_PFP_track_SimpleCluster_hitLinearity;
};


StoppingParticleMichelTagger::StoppingParticleMichelTagger(fhicl::ParameterSet const & p)
  {


  _ct_manager.Configure(p.get<cosmictag::Config_t>("CosmicTagManager"));

  CC1piInputTags = new InputTags(p);
  if (_debug) CC1piInputTags->PrintConfig();

  _coplanar_cut       = p.get<double>("CoplanarCut",   5.);

  _debug = p.get<bool>("DebugMode", true);

  fDetectorProperties = lar::providerFrom<detinfo::DetectorPropertiesService>();

  art::ServiceHandle<art::TFileService> fs;
  _tree1 = fs->make<TTree>("StoppingTaggerTree","");
  _tree1->Branch("run",           &_run,                 "run/I");
  _tree1->Branch("subrun",        &_subrun,              "subrun/I");
  _tree1->Branch("event",         &_event,               "event/I");
  _tree1->Branch("TPCObj_PFP_track_length_StoppingModule", &TPCObj_PFP_track_length);
  _tree1->Branch("TPCObj_PFP_track_SimpleCluster_hitTime", &TPCObj_PFP_track_SimpleCluster_hitTime);
  _tree1->Branch("TPCObj_PFP_track_SimpleCluster_hitWire", &TPCObj_PFP_track_SimpleCluster_hitWire);
  _tree1->Branch("TPCObj_PFP_track_SimpleCluster_hitPlane", &TPCObj_PFP_track_SimpleCluster_hitPlane);
  _tree1->Branch("TPCObj_PFP_track_SimpleCluster_hitIntegral", &TPCObj_PFP_track_SimpleCluster_hitIntegral);
  _tree1->Branch("TPCObj_PFP_track_SimpleCluster_hitTimeTicks", &TPCObj_PFP_track_SimpleCluster_hitTimeTicks);
  _tree1->Branch("TPCObj_PFP_track_SimpleCluster_hitWireNo", &TPCObj_PFP_track_SimpleCluster_hitWireNo);
  _tree1->Branch("TPCObj_PFP_track_SimpleCluster_StartIndex", &TPCObj_PFP_track_SimpleCluster_StartIndex);
  _tree1->Branch("TPCObj_PFP_track_SimpleCluster_hitdQds", &TPCObj_PFP_track_SimpleCluster_hitdQds);
  _tree1->Branch("TPCObj_PFP_track_SimpleCluster_hitds", &TPCObj_PFP_track_SimpleCluster_hitds);
  _tree1->Branch("TPCObj_PFP_track_SimpleCluster_hitdQdsSlider", &TPCObj_PFP_track_SimpleCluster_hitdQdsSlider);
  _tree1->Branch("TPCObj_PFP_track_SimpleCluster_hitLinearity", &TPCObj_PFP_track_SimpleCluster_hitLinearity);
  _tree1->Branch("TPCObj_PFP_track_ct_passed_basic", &TPCObj_PFP_track_ct_passed_basic);
  _tree1->Branch("TPCObj_PFP_track_ct_result_michel", &TPCObj_PFP_track_ct_result_michel);
  _tree1->Branch("TPCObj_PFP_track_ct_result_bragg", &TPCObj_PFP_track_ct_result_bragg);

}



void StoppingParticleMichelTagger::produce(art::Event & evt) {

  if (_debug) std::cout <<"[StoppingParticleMichelTagger] Starts." << std::endl;

  _run = evt.id().run();
  _subrun = evt.id().subRun();
  _event = evt.id().event();

  // Make things to put in tree
  TPCObj_PFP_track_length.clear();
  TPCObj_PFP_track_deltax.clear();
  TPCObj_PFP_track_SimpleCluster_hitTime.clear();
  TPCObj_PFP_track_SimpleCluster_hitWire.clear();
  TPCObj_PFP_track_SimpleCluster_hitPlane.clear();
  TPCObj_PFP_track_SimpleCluster_hitIntegral.clear();
  TPCObj_PFP_track_SimpleCluster_hitTimeTicks.clear();
  TPCObj_PFP_track_SimpleCluster_hitWireNo.clear();
  TPCObj_PFP_track_SimpleCluster_StartIndex.clear();
  TPCObj_PFP_track_SimpleCluster_hitdQds.clear();
  TPCObj_PFP_track_SimpleCluster_hitds.clear();
  TPCObj_PFP_track_SimpleCluster_hitdQdsSlider.clear();
  TPCObj_PFP_track_SimpleCluster_hitLinearity.clear();


  // Get result of Marco's selection
  art::Handle<std::vector<ubana::SelectionResult>> selection_h;
  evt.getByLabel(CC1piInputTags->fSelectionLabel, selection_h);
  if(!selection_h.isValid()){
     mf::LogError(__PRETTY_FUNCTION__) << "[FillTree] SelectionResult product not found." << std::endl;
     throw std::exception();
  }
  std::vector<art::Ptr<ubana::SelectionResult>> selection_v;
  art::fill_ptr_vector(selection_v, selection_h);


  // Skip all tracks in events that are not selected
  // This was the best way I could think of to not mess up the alignment of the trees
  if (selection_v.at(0)->GetSelectionStatus()){

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

    // Store whether the track is deemed to be coplanar (i.e. has start and end x positions within coplanar cut defined above, which would imply that the track is coplanar to a collection plane wire, and we probably don't have good calorimetry). Don't actually do anything with this information here, but store it so we can look at it later.
    double deltax = track_i->Vertex().Z() - track_i->End().Z();
    if (_debug) std::cout << "[StoppingParticleMichelTagger] Delta x is " << deltax << std::endl;
    // bool collection_is_coplanar = TMath::Abs(deltax) < _coplanar_cut;

    // Create an approximate start hit on plane 2 from the spacepoint closest to the reconstructed neutrino vertex candidate. Do this by getting the PFP associated with the track (there can only be one), and the clusters/hits/spacepoints associated to that PFP.
    // First, exclude spacepoints outside the TPC. Then get the point closest to the reconstructed neutrino vertex.
    // Turn into a cosmictag::SimpleHit for use with the HitCosmicTag Algorithms
    art::Ptr<recob::PFParticle> pfp_from_track = theParticleAssns.at(track_i.key());

    // std::vector<art::Ptr<recob::SpacePoint>> sp_v;
    // sp_v.clear();

    std::vector<art::Ptr<recob::Hit>> hit_v;
    hit_v.clear();

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
    // Now that all algorithms have been run, get the SimpleCluster and save elements of it
    cosmictag::SimpleCluster sc_tosave = _ct_manager.GetCluster();
    std::vector<cosmictag::SimpleHit> sc_tosave_s_hit_v = sc_tosave._s_hit_v;
    std::vector<double> sc_tosave_s_hit_v_time;
    std::vector<int> sc_tosave_s_hit_v_wire;
    std::vector<int> sc_tosave_s_hit_v_plane;
    std::vector<double> sc_tosave_s_hit_v_integral;
    std::vector<double> sc_tosave_s_hit_v_tticks;
    std::vector<int> sc_tosave_s_hit_v_wno;
    for (size_t i_sh=0; i_sh<sc_tosave_s_hit_v.size(); i_sh++){
      sc_tosave_s_hit_v_time.push_back(sc_tosave_s_hit_v.at(i_sh).time);
      sc_tosave_s_hit_v_wire.push_back(sc_tosave_s_hit_v.at(i_sh).wire);
      sc_tosave_s_hit_v_plane.push_back(sc_tosave_s_hit_v.at(i_sh).plane);
      sc_tosave_s_hit_v_integral.push_back(sc_tosave_s_hit_v.at(i_sh).integral);
      sc_tosave_s_hit_v_tticks.push_back(sc_tosave_s_hit_v.at(i_sh).t);
      sc_tosave_s_hit_v_wno.push_back(sc_tosave_s_hit_v.at(i_sh).w);
    }
    int sc_tosave_start_index = sc_tosave._start_index;
    std::vector<double> sc_tosave_dqds_v = sc_tosave._dqds_v;
    std::vector<double> sc_tosave_ds_v = sc_tosave._ds_v;
    std::vector<double> sc_tosave_dqds_slider = sc_tosave._dqds_slider;
    std::vector<double> sc_tosave_linearity_v = sc_tosave._linearity_v;


    TPCObj_PFP_track_length.emplace_back(track_i->Length());
    TPCObj_PFP_track_deltax.emplace_back(deltax);
    TPCObj_PFP_track_SimpleCluster_hitTime.emplace_back(sc_tosave_s_hit_v_time);
    TPCObj_PFP_track_SimpleCluster_hitWire.emplace_back(sc_tosave_s_hit_v_wire);
    TPCObj_PFP_track_SimpleCluster_hitPlane.emplace_back(sc_tosave_s_hit_v_plane);
    TPCObj_PFP_track_SimpleCluster_hitIntegral.emplace_back(sc_tosave_s_hit_v_integral);
    TPCObj_PFP_track_SimpleCluster_hitTimeTicks.emplace_back(sc_tosave_s_hit_v_tticks);
    TPCObj_PFP_track_SimpleCluster_hitWireNo.emplace_back(sc_tosave_s_hit_v_wno);
    TPCObj_PFP_track_SimpleCluster_StartIndex.emplace_back(sc_tosave_start_index);
    TPCObj_PFP_track_SimpleCluster_hitdQds.emplace_back(sc_tosave_dqds_v);
    TPCObj_PFP_track_SimpleCluster_hitds.emplace_back(sc_tosave_ds_v);
    TPCObj_PFP_track_SimpleCluster_hitdQdsSlider.emplace_back(sc_tosave_dqds_slider);
    TPCObj_PFP_track_SimpleCluster_hitLinearity.emplace_back(sc_tosave_linearity_v);
    TPCObj_PFP_track_ct_passed_basic.emplace_back(passed);
    TPCObj_PFP_track_ct_result_michel.emplace_back(ct_result_michel);
    TPCObj_PFP_track_ct_result_bragg.emplace_back(ct_result_bragg);

    } // End loop over recob::tracks
  } // end if (selection_v.at(0)->GetSelectionStatus())

  _tree1->Fill();

} // end ::produce



DEFINE_ART_MODULE(StoppingParticleMichelTagger)
