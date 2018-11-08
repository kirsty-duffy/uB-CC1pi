////////////////////////////////////////////////////////////////////////
// Class:       ArePionsWellReconstructed
// Plugin Type: analyzer (art v2_05_01)
// File:        ArePionsWellReconstructed_module.cc
//
// Generated at Mon Oct  1 12:16:58 2018 by Kirsty Duffy using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////


#include "art/Framework/Core/EDAnalyzer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/Cluster.h"
#include "lardataobj/RecoBase/SpacePoint.h"

#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"
#include "Pandora/PdgTable.h"

#include "PDGEnums.h"

#include "TruthMatchTracks.h"
#include "TruthMatchTracks.cxx"
#include "GetDaughters.h"
#include "GetDaughters.cxx"

#include "TTree.h"
#include "TVector3.h"

class ArePionsWellReconstructed;


class ArePionsWellReconstructed : public art::EDAnalyzer {
public:
  explicit ArePionsWellReconstructed(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  ArePionsWellReconstructed(ArePionsWellReconstructed const &) = delete;
  ArePionsWellReconstructed(ArePionsWellReconstructed &&) = delete;
  ArePionsWellReconstructed & operator = (ArePionsWellReconstructed const &) = delete;
  ArePionsWellReconstructed & operator = (ArePionsWellReconstructed &&) = delete;

  // Required functions.
  void analyze(art::Event const & e) override;

private:
  // TFile service
  art::ServiceHandle<art::TFileService> tfs;
  TTree *_outtree;

  std::string fTrackLabel;
  std::string fHitLabel;
  std::string fHitTrackAssns;
  std::string fHitTruthAssns;
  std::string fPFPproducerlabel;

  bool fIsVerbose;
  bool fUseMuonsInstead;

  // Variables to fill in tree
  // One entry in tree for each true primary pi+
  int PiPlusHierarchy_nMCPs;
  std::vector<int> PiPlusHierarchy_MCP_ID;
  std::vector<int> PiPlusHierarchy_MCP_MotherID;
  std::vector<std::vector<int>> PiPlusHierarchy_MCP_DaughterIDs;
  std::vector<int> PiPlusHierarchy_MCP_PDGcode;
  std::vector<std::string> PiPlusHierarchy_MCP_EndProcess;
  std::vector<int> PiPlusHierarchy_MCP_nMatchedPFPs;
  std::vector<double> PiPlusHierarchy_MCP_totaldepE;
  std::vector<std::vector<int>> PiPlusHierarchy_MCP_matchedPFP_ID;
  std::vector<std::vector<double>> PiPlusHierarchy_MCP_matchedPFP_matchedE;
  std::vector<double> PiPlusHierarchy_MCP_trueStartE;
  std::vector<double> PiPlusHierarchy_MCP_trueStartP;
  std::vector<int> PiPlusHierarchy_LogicalPion_MCPids;

  int n_PFPs;
  std::vector<int> PFP_ID;
  std::vector<int> PFP_TrackShowerPdg;
  std::vector<bool> PFP_IsPrimary;
  std::vector<std::pair<int,int>> PFP_ID_to_bestmatchMCPid;
  std::vector<double> PFP_totaldepE;
  std::vector<int> PFP_primaryPFPid;
  std::vector<std::vector<std::vector<double>>> PFP_spacepoints_XYZ;
  std::vector<std::vector<std::vector<double>>> PFP_trajpoints_XYZ;
  std::vector<double> PFP_track_length;
  std::vector<std::vector<double>> PFP_track_start;
  std::vector<std::vector<double>> PFP_track_end;

};


ArePionsWellReconstructed::ArePionsWellReconstructed(fhicl::ParameterSet const & p)
  :
  EDAnalyzer(p)  // ,
 // More initializers here.
{
  fHitTrackAssns = p.get<std::string>("HitTrackAssns","pandoraNu::UBXSec");
  fHitTruthAssns = p.get<std::string>("HitTruthAssns","pandoraCosmicHitRemoval::UBXSec");
  fTrackLabel = p.get<std::string>("TrackLabel","pandoraNu::UBXSec");
  fHitLabel = p.get<std::string>("HitLabel","pandoraCosmicHitRemoval::UBXSec");
  fPFPproducerlabel = p.get<std::string>("PFPproducerlabel","pandoraNu::UBXSec");
  fIsVerbose = p.get<bool>("verbose",false);
  fUseMuonsInstead = p.get<bool>("useMuonsInstead",false);

  _outtree = tfs->make<TTree>("PiRecoTree","");
  _outtree->Branch("PiPlusHierarchy_nMCPs",&PiPlusHierarchy_nMCPs);
  _outtree->Branch("PiPlusHierarchy_MCP_ID",&PiPlusHierarchy_MCP_ID);
  _outtree->Branch("PiPlusHierarchy_MCP_MotherID",&PiPlusHierarchy_MCP_MotherID);
  _outtree->Branch("PiPlusHierarchy_MCP_DaughterIDs",&PiPlusHierarchy_MCP_DaughterIDs);
  _outtree->Branch("PiPlusHierarchy_MCP_PDGcode",&PiPlusHierarchy_MCP_PDGcode);
  _outtree->Branch("PiPlusHierarchy_MCP_EndProcess",&PiPlusHierarchy_MCP_EndProcess);
  _outtree->Branch("PiPlusHierarchy_MCP_nMatchedPFPs",&PiPlusHierarchy_MCP_nMatchedPFPs);
  _outtree->Branch("PiPlusHierarchy_MCP_totaldepE",&PiPlusHierarchy_MCP_totaldepE);
  _outtree->Branch("PiPlusHierarchy_MCP_matchedPFP_ID",&PiPlusHierarchy_MCP_matchedPFP_ID);
  _outtree->Branch("PiPlusHierarchy_MCP_matchedPFP_matchedE",&PiPlusHierarchy_MCP_matchedPFP_matchedE);
  _outtree->Branch("PiPlusHierarchy_MCP_trueStartE",&PiPlusHierarchy_MCP_trueStartE);
  _outtree->Branch("PiPlusHierarchy_MCP_trueStartP",&PiPlusHierarchy_MCP_trueStartP);
  _outtree->Branch("PiPlusHierarchy_LogicalPion_MCPids",&PiPlusHierarchy_LogicalPion_MCPids);

  _outtree->Branch("n_PFPs",&n_PFPs);
  _outtree->Branch("PFP_ID",&PFP_ID);
  _outtree->Branch("PFP_TrackShowerPdg",&PFP_TrackShowerPdg);
  _outtree->Branch("PFP_IsPrimary",&PFP_IsPrimary);
  _outtree->Branch("PFP_ID_to_bestmatchMCPid",&PFP_ID_to_bestmatchMCPid);
  _outtree->Branch("PFP_totaldepE",&PFP_totaldepE);
  _outtree->Branch("PFP_primaryPFPid",&PFP_primaryPFPid);
  _outtree->Branch("PFP_spacepoints_XYZ",&PFP_spacepoints_XYZ);
  _outtree->Branch("PFP_trajpoints_XYZ",&PFP_trajpoints_XYZ);
  _outtree->Branch("PFP_track_length",&PFP_track_length);
  _outtree->Branch("PFP_track_start",&PFP_track_start);
  _outtree->Branch("PFP_track_end",&PFP_track_end);
}

void ArePionsWellReconstructed::analyze(art::Event const & evt)
{
  // Get all primary (i.e. direct daughters of the neutrino) MCParticles
  auto const mctruth_h = evt.getValidHandle<std::vector<simb::MCTruth>>("generator"); // Get only GENIE MCtruth
  int nu_MCPID = mctruth_h->at(0).GetNeutrino().Nu().TrackId();

  art::Handle<std::vector<simb::MCParticle>> mcp_h;
  evt.getByLabel("largeant", mcp_h);
  std::vector<art::Ptr<simb::MCParticle>> mcp_v;
  art::fill_ptr_vector(mcp_v, mcp_h);

  std::vector<art::Ptr<simb::MCParticle>> primary_mcps;
  std::vector<art::Ptr<simb::MCParticle>> primary_piplusmcps;
  for (auto mcpar : mcp_v) {
    if (mcpar->Mother() != nu_MCPID) continue;
    primary_mcps.push_back(mcpar);

    if (fUseMuonsInstead && (PDGCode(mcpar->PdgCode()) != kMuMinus)) continue;
    if (!fUseMuonsInstead && (PDGCode(mcpar->PdgCode()) != kPiPlus)) continue;
    primary_piplusmcps.push_back(mcpar);
  }

  if (fIsVerbose) std::cout << "Found " << primary_mcps.size() << " primary neutrino daughters (" << primary_piplusmcps.size() << " pi+)" << std::endl;

  // Now get all PFPs and hits from the event
  art::Handle<std::vector<recob::PFParticle>> pfpHandle;
  evt.getByLabel(fPFPproducerlabel, pfpHandle);
  std::vector<art::Ptr<recob::PFParticle>> pfpCollection;


  // Get mapping from ID to PFParticle
  /*std::map<size_t, art::Ptr<recob::PFParticle> >*/ lar_pandora::PFParticleMap particleMap;
  for (unsigned int i = 0; i < pfpHandle->size(); ++i)
  {
      const art::Ptr<recob::PFParticle> pfParticle(pfpHandle, i);
      if (!particleMap.emplace(pfParticle->Self(), pfParticle).second)
          throw cet::exception("ArePionsWellReconstructed") << "Repeated PFParticles in the input collection!" << std::endl;
  }

  // Ok, here we need a difference between MCC8 and MCC 9
  // In MCC 8, we use pandoraNu, so every PFP in pfpHandle should be put into pfpCollection
  // However, in MCC 9, we use pandora consolidated output, and have to extract only the particles that are daughters of the neutrino PFParticle
  // Tell which case we're considering by looking at fPFPproducerlabel: if it contains "pandora" but not "pandoraNu" or "pandoraCosmic" then assume it's MCC9 (this is not true for all producers, but it is true for the ones we use)
  if ((fPFPproducerlabel.find(std::string("pandora")) != std::string::npos) && (fPFPproducerlabel.find(std::string("pandoraNu")) == std::string::npos) && (fPFPproducerlabel.find(std::string("pandoraCosmic")) == std::string::npos)){ // MCC 9
    if (fIsVerbose) std::cout << "Found PFPproducerlabel " << fPFPproducerlabel << ". Running code for pandora consolidated output." << std::endl;
    // Find the daughters of the neutrino PFParticle
    for (unsigned int i = 0; i < pfpHandle->size(); ++i)
    {
        const art::Ptr<recob::PFParticle> pfParticle(pfpHandle, i);

        // If the particle is primary, it doesn't have a parent
        if (pfParticle->IsPrimary())
            continue;
        // Find the parent particle
        const auto parentIterator = particleMap.find(pfParticle->Parent());
        if (parentIterator == particleMap.end())
            throw cet::exception("ArePionsWellReconstructed") << "PFParticle has parent that's not in the input collection" << std::endl;

        // Check if the parent is primary
        if (!parentIterator->second->IsPrimary())
            continue;

        // Check if the parent has a neutrino PDG code
        const int parentPDG = std::abs(parentIterator->second->PdgCode());
        if (parentPDG != pandora::NU_E && parentPDG != pandora::NU_MU && parentPDG != pandora::NU_TAU)
            continue;
        // If we get to here, then the PFParticle's parent is the reconstructed neutrino! So it's a reconstructed final state
        pfpCollection.push_back(pfParticle);
    }
  }
  else{ // MCC 8
    art::fill_ptr_vector(pfpCollection, pfpHandle);
  }


  art::Handle<std::vector<recob::Hit>> hitHandle;
  evt.getByLabel(fHitLabel, hitHandle);

  // Do truth matching: best match (match 1 MCP to each PFP)
  std::map<int, MatchingData> pfpID_to_matchingData_BestMatch;
  std::map<int, MatchingData> MCPID_to_matchingData_BestMatch;

  Do_PFP_MCP_Matching_BestMatch(evt, pfpHandle, hitHandle, fPFPproducerlabel, fHitTruthAssns, pfpID_to_matchingData_BestMatch, MCPID_to_matchingData_BestMatch);

  // Do truth matching: find all PFP matches for each MCP
  std::map<int, std::vector<MatchingData>> MCPID_to_many_matchingData;
  Do_PFP_MCP_Matching_AllMatches(evt, mcp_h, pfpHandle, hitHandle, fPFPproducerlabel, fHitTruthAssns, MCPID_to_many_matchingData);

  // Get MCP<->hit associations so we can calculate the deposited energy for each MCP
  art::FindManyP<recob::Hit,anab::BackTrackerHitMatchingData> hits_per_mcp(mcp_h, evt, fHitTruthAssns);

  // Get PFP<->spacepoint and PFP<->track associations so we can store the positions of space points and trajectory points
  art::FindManyP<recob::SpacePoint> spacepoints_per_pfp(pfpHandle, evt, fPFPproducerlabel);
  art::FindManyP<recob::Track> tracks_from_pfps(pfpHandle, evt, fPFPproducerlabel);


  // Loop through pi+ primaries and gather all MCParticles in the daughter hierarchy
  for (auto piplusprimary : primary_piplusmcps){

    // Reset tree variables that are filled per pi+ primary particle
    PiPlusHierarchy_nMCPs = 0;
    PiPlusHierarchy_MCP_ID.clear();
    PiPlusHierarchy_MCP_MotherID.clear();
    PiPlusHierarchy_MCP_DaughterIDs.clear();
    PiPlusHierarchy_MCP_PDGcode.clear();
    PiPlusHierarchy_MCP_EndProcess.clear();
    PiPlusHierarchy_MCP_nMatchedPFPs.clear();
    PiPlusHierarchy_MCP_totaldepE.clear();
    PiPlusHierarchy_MCP_matchedPFP_ID.clear();
    PiPlusHierarchy_MCP_matchedPFP_matchedE.clear();
    PiPlusHierarchy_MCP_trueStartE.clear();
    PiPlusHierarchy_MCP_trueStartP.clear();
    PiPlusHierarchy_LogicalPion_MCPids.clear();

    n_PFPs = 0;
    PFP_ID.clear();
    PFP_TrackShowerPdg.clear();
    PFP_IsPrimary.clear();
    PFP_ID_to_bestmatchMCPid.clear();
    PFP_totaldepE.clear();
    PFP_primaryPFPid.clear();
    PFP_spacepoints_XYZ.clear();
    PFP_trajpoints_XYZ.clear();


    // Now get true MCP hierarchy from primary pi+
    std::vector<art::Ptr<simb::MCParticle>> hierarchy = MCP_GetFullHierarchyFromPrimary(mcp_v, piplusprimary, nu_MCPID);

    // Fill variables for all MCPs in hierarchy
    PiPlusHierarchy_nMCPs = hierarchy.size();
    for (auto mcp : hierarchy){
      PiPlusHierarchy_MCP_ID.push_back(mcp->TrackId());
      PiPlusHierarchy_MCP_MotherID.push_back(mcp->Mother());
      PiPlusHierarchy_MCP_PDGcode.push_back(mcp->PdgCode());
      PiPlusHierarchy_MCP_EndProcess.push_back(mcp->EndProcess());
      PiPlusHierarchy_MCP_trueStartE.push_back(mcp->E());
      PiPlusHierarchy_MCP_trueStartP.push_back(mcp->P());

      std::vector<int> daughters;
      for (int i_d=0; i_d<mcp->NumberDaughters(); i_d++){
        daughters.push_back(mcp->Daughter(i_d));
      }
      PiPlusHierarchy_MCP_DaughterIDs.push_back(daughters);

      // Get true deposited energy by MCParticle from truth-matching information
      std::vector<art::Ptr<recob::Hit>> hit_vec;
      std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
      hits_per_mcp.get(mcp.key(),hit_vec,match_vec);
      double depE = 0;
      for (size_t i_matchhit=0; i_matchhit<match_vec.size(); i_matchhit++){
        depE+=match_vec.at(i_matchhit)->energy;
      }
      PiPlusHierarchy_MCP_totaldepE.push_back(depE);

      // Look for matched PFPs
      auto search = MCPID_to_many_matchingData.find(mcp->TrackId());
      if (search != MCPID_to_many_matchingData.end()){
        if (fIsVerbose) std::cout << "particle with energy " << depE << " matched to the following " << search->second.size() <<  " PFPs: " << std::endl;

        PiPlusHierarchy_MCP_nMatchedPFPs.push_back(search->second.size());
      }

      std::vector<int> matchedPFP_ID;
      std::vector<double> matchedPFP_matchedE;
      for (auto mtch : search->second){
        if (fIsVerbose) {
          std::cout << "PFP ID = " << mtch.recoObject_ID << " with matched energy " << mtch.MatchedEnergy << std::endl;
          std::cout << "  (PFP best matched to MCP with PDG " << pfpID_to_matchingData_BestMatch.find(mtch.recoObject_ID)->second.MCParticle_PDG << ", with purity " << pfpID_to_matchingData_BestMatch.find(mtch.recoObject_ID)->second.MatchedEnergy << "/" << pfpID_to_matchingData_BestMatch.find(mtch.recoObject_ID)->second.recoObject_totalEnergy << ")" << std::endl;
        } // end if verbose

        matchedPFP_ID.push_back(mtch.recoObject_ID);
        matchedPFP_matchedE.push_back(mtch.MatchedEnergy);
      }
      PiPlusHierarchy_MCP_matchedPFP_ID.push_back(matchedPFP_ID);
      PiPlusHierarchy_MCP_matchedPFP_matchedE.push_back(matchedPFP_matchedE);

    } // end loop over MCPs in hierarchy

    // Also construct "logical pion": a vector of MCParticles that contains all pions that are true daughters of the preceeding pion (combine together pions from elastic scatters that GEANT splits into two different particles)
    std::vector<art::Ptr<simb::MCParticle>> logicalpion = MCP_GetLogicalPion(piplusprimary,mcp_v,fUseMuonsInstead);

    // Make vector of logical pion MCP IDs so we can navigate it easily later
    for (auto piplus : logicalpion){
      PiPlusHierarchy_LogicalPion_MCPids.push_back(piplus->TrackId());
    }

    // Finally, fill variables related to PFPs in the event
    n_PFPs = pfpCollection.size();
    for (auto pfp : pfpCollection){
      PFP_ID.push_back((int)pfp->Self());
      PFP_TrackShowerPdg.push_back(pfp->PdgCode());

      PFP_IsPrimary.push_back(pfp->IsPrimary());

      // Use LArPandoraHelper functions to get primary PFP
      art::Ptr<recob::PFParticle> primarypfp;
      primarypfp = lar_pandora::LArPandoraHelper::GetFinalStatePFParticle(particleMap, pfp);
      int primarypfpid = -999;
      if (primarypfp) primarypfpid = primarypfp->Self();
      PFP_primaryPFPid.push_back(primarypfpid);

      if (fIsVerbose){
        std::cout << "PFP ID = " << pfp->Self() << ", Primary PFP ID = " << primarypfpid << ". PrimaryPFP->IsPrimary() = " << primarypfp->IsPrimary() << std::endl;
      }

      auto search = pfpID_to_matchingData_BestMatch.find((int)pfp->Self());
      if (search != pfpID_to_matchingData_BestMatch.end()){
        PFP_ID_to_bestmatchMCPid.push_back(std::make_pair(pfp->Self(),search->second.MCParticle_ID));
        PFP_totaldepE.push_back(search->second.recoObject_totalEnergy);
      }
      else{
        PFP_ID_to_bestmatchMCPid.push_back(std::make_pair(pfp->Self(),-999));
        PFP_totaldepE.push_back(-999);
      }


      // Get spacepoints from PFP
      std::vector<art::Ptr<recob::SpacePoint>> sp_vec;
      spacepoints_per_pfp.get(pfp.key(),sp_vec);
      if (fIsVerbose){
        std::cout << "PFP has " << sp_vec.size() << " associated space points" << std::endl;
      }
      std::vector<std::vector<double>> sp_tvector_vec;
      for (size_t i_sp=0; i_sp<sp_vec.size(); i_sp++){
        const double *XYZ = sp_vec.at(i_sp)->XYZ();
        if (XYZ==nullptr) continue;
        std::vector<double> dummytvec = {XYZ[0],XYZ[1],XYZ[2]};
        sp_tvector_vec.push_back(dummytvec);
      } // end loop over spacepoints
      PFP_spacepoints_XYZ.push_back(sp_tvector_vec);

      // Get track from PFP (if it exists)
      std::vector<art::Ptr<recob::Track>> tracks_pfp = tracks_from_pfps.at(pfp.key());
      double track_length = -999;
      std::vector<double> track_start = {-999,-999,-999};
      std::vector<double> track_end = {-999,-999,-999};
      std::vector<std::vector<double>> tr_tvector_vec;
      if(tracks_pfp.size() > 1) {
         mf::LogError(__PRETTY_FUNCTION__) << "PFP associated to more than one track." << std::endl;
         throw std::exception();
      }
      for (auto track : tracks_pfp){
         track_length = track -> Length();
         auto start = track -> Start();
         track_start = {start.X(),start.Y(),start.Z()};
         auto end = track -> End();
         track_end = {end.X(),end.Y(),end.Z()};

         // Now get trajectory points
         size_t ntrajpoints = track -> NumberTrajectoryPoints();
         for (size_t i_trj=0; i_trj<ntrajpoints; i_trj++){
           std::vector<double> dummytvec = {-999,-999,-999};
           dummytvec.at(0) = track->TrajectoryPoint(i_trj).position.X();
           dummytvec.at(1) = track->TrajectoryPoint(i_trj).position.Y();
           dummytvec.at(2) = track->TrajectoryPoint(i_trj).position.Z();
           tr_tvector_vec.push_back(dummytvec);
         } // end loop over traj points
       } // end loop over tracks associated to pfp (there should only be 0 or 1)
      PFP_trajpoints_XYZ.push_back(tr_tvector_vec);
      PFP_track_length.push_back(track_length);
      PFP_track_start.push_back(track_start);
      PFP_track_end.push_back(track_end);


    } // end loop over PFPs

    _outtree->Fill();

  } // end loop over pi+ primaries
}

DEFINE_ART_MODULE(ArePionsWellReconstructed)
