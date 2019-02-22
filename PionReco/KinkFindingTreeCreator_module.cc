////////////////////////////////////////////////////////////////////////
// Class:       KinkFindingTreeCreator
// Plugin Type: analyzer (art v2_05_01)
// File:        KinkFindingTreeCreator_module.cc
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
#include "Plotting/CalcLocalLinearity.h"

#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"

class KinkFindingTreeCreator;


class KinkFindingTreeCreator : public art::EDAnalyzer {
public:
  explicit KinkFindingTreeCreator(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  KinkFindingTreeCreator(KinkFindingTreeCreator const &) = delete;
  KinkFindingTreeCreator(KinkFindingTreeCreator &&) = delete;
  KinkFindingTreeCreator & operator = (KinkFindingTreeCreator const &) = delete;
  KinkFindingTreeCreator & operator = (KinkFindingTreeCreator &&) = delete;

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

  // Variables to fill in tree
  // One entry in tree for each event
  int nMCPs;
  int nu_MCPID;
  std::vector<int> MCP_ID;
  std::vector<int> MCP_MotherID;
  std::vector<std::vector<int>> MCP_DaughterIDs;
  std::vector<int> MCP_PDGcode;
  std::vector<std::string> MCP_EndProcess;
  std::vector<int> MCP_nMatchedPFPs;
  std::vector<double> MCP_totaldepE;
  std::vector<std::vector<int>> MCP_matchedPFP_ID;
  std::vector<std::vector<double>> MCP_matchedPFP_matchedE;
  std::vector<double> MCP_trueStartE;
  std::vector<double> MCP_trueStartP;
  std::vector<std::vector<double>> MCP_Px_eachpoint;
  std::vector<std::vector<double>> MCP_Py_eachpoint;
  std::vector<std::vector<double>> MCP_Pz_eachpoint;
  std::vector<std::vector<double>> MCP_x_eachpoint;
  std::vector<std::vector<double>> MCP_y_eachpoint;
  std::vector<std::vector<double>> MCP_z_eachpoint;
  std::vector<std::vector<double>> MCP_StartXYZ;
  std::vector<std::vector<double>> MCP_EndXYZ;

  int n_PFPs;
  std::vector<int> PFP_ID;
  std::vector<int> PFP_TrackShowerPdg;
  std::vector<bool> PFP_IsPrimary;
  std::vector<std::pair<int,int>> PFP_ID_to_bestmatchMCPid;
  std::vector<double> PFP_totaldepE;
  std::vector<int> PFP_primaryPFPid;
  std::vector<std::vector<std::vector<double>>> PFP_spacepoints_XYZ;
  std::vector<std::vector<std::vector<double>>> PFP_ordered_spacepoints;
  std::vector<std::vector<int>> PFP_ordered_spacepoints_origidx;
  std::vector<std::vector<int>> PFP_goodreco_MCPvectorposes;
  std::vector<std::vector<double>> PFP_trueMCPend;
  std::vector<std::vector<std::vector<double>>> PFP_trajpoints_XYZ;
  std::vector<double> PFP_track_length;
  std::vector<std::vector<double>> PFP_track_start;
  std::vector<std::vector<double>> PFP_track_end;
  std::vector<double> PFP_track_theta;
  std::vector<double> PFP_track_phi;

};


KinkFindingTreeCreator::KinkFindingTreeCreator(fhicl::ParameterSet const & p)
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

  _outtree = tfs->make<TTree>("KinkFindingTree","");
  _outtree->Branch("nMCPs",&nMCPs);
  _outtree->Branch("nu_MCPID",&nu_MCPID);
  _outtree->Branch("MCP_ID",&MCP_ID);
  _outtree->Branch("MCP_MotherID",&MCP_MotherID);
  _outtree->Branch("MCP_DaughterIDs",&MCP_DaughterIDs);
  _outtree->Branch("MCP_PDGcode",&MCP_PDGcode);
  _outtree->Branch("MCP_EndProcess",&MCP_EndProcess);
  _outtree->Branch("MCP_nMatchedPFPs",&MCP_nMatchedPFPs);
  _outtree->Branch("MCP_totaldepE",&MCP_totaldepE);
  _outtree->Branch("MCP_matchedPFP_ID",&MCP_matchedPFP_ID);
  _outtree->Branch("MCP_matchedPFP_matchedE",&MCP_matchedPFP_matchedE);
  _outtree->Branch("MCP_trueStartE",&MCP_trueStartE);
  _outtree->Branch("MCP_trueStartP",&MCP_trueStartP);
  _outtree->Branch("MCP_Px_eachpoint",&MCP_Px_eachpoint);
  _outtree->Branch("MCP_Py_eachpoint",&MCP_Py_eachpoint);
  _outtree->Branch("MCP_Pz_eachpoint",&MCP_Pz_eachpoint);
  _outtree->Branch("MCP_x_eachpoint",&MCP_x_eachpoint);
  _outtree->Branch("MCP_y_eachpoint",&MCP_y_eachpoint);
  _outtree->Branch("MCP_z_eachpoint",&MCP_z_eachpoint);
  _outtree->Branch("MCP_StartXYZ",&MCP_StartXYZ);
  _outtree->Branch("MCP_EndXYZ",&MCP_EndXYZ);

  _outtree->Branch("n_PFPs",&n_PFPs);
  _outtree->Branch("PFP_ID",&PFP_ID);
  _outtree->Branch("PFP_TrackShowerPdg",&PFP_TrackShowerPdg);
  _outtree->Branch("PFP_IsPrimary",&PFP_IsPrimary);
  _outtree->Branch("PFP_ID_to_bestmatchMCPid",&PFP_ID_to_bestmatchMCPid);
  _outtree->Branch("PFP_totaldepE",&PFP_totaldepE);
  _outtree->Branch("PFP_primaryPFPid",&PFP_primaryPFPid);
  _outtree->Branch("PFP_spacepoints_XYZ",&PFP_spacepoints_XYZ);
  _outtree->Branch("PFP_ordered_spacepoints",&PFP_ordered_spacepoints);
  _outtree->Branch("PFP_ordered_spacepoints_origidx",&PFP_ordered_spacepoints_origidx);
  _outtree->Branch("PFP_goodreco_MCPvectorposes",&PFP_goodreco_MCPvectorposes);
  _outtree->Branch("PFP_trueMCPend",&PFP_trueMCPend);
  _outtree->Branch("PFP_trajpoints_XYZ",&PFP_trajpoints_XYZ);
  _outtree->Branch("PFP_track_length",&PFP_track_length);
  _outtree->Branch("PFP_track_start",&PFP_track_start);
  _outtree->Branch("PFP_track_end",&PFP_track_end);
  _outtree->Branch("PFP_track_theta",&PFP_track_theta);
  _outtree->Branch("PFP_track_phi",&PFP_track_phi);
}

void KinkFindingTreeCreator::analyze(art::Event const & evt)
{
  // Get all primary (i.e. direct daughters of the neutrino) MCParticles
  auto const mctruth_h = evt.getValidHandle<std::vector<simb::MCTruth>>("generator"); // Get only GENIE MCtruth
  nu_MCPID = mctruth_h->at(0).GetNeutrino().Nu().TrackId();

  art::Handle<std::vector<simb::MCParticle>> mcp_h;
  evt.getByLabel("largeant", mcp_h);
  std::vector<art::Ptr<simb::MCParticle>> mcp_v;
  art::fill_ptr_vector(mcp_v, mcp_h);

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
          throw cet::exception("KinkFindingTreeCreator") << "Repeated PFParticles in the input collection!" << std::endl;
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
            throw cet::exception("KinkFindingTreeCreator") << "PFParticle has parent that's not in the input collection" << std::endl;

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


  // Reset tree variables that are filled per event
  nMCPs = 0;
  MCP_ID.clear();
  MCP_MotherID.clear();
  MCP_DaughterIDs.clear();
  MCP_PDGcode.clear();
  MCP_EndProcess.clear();
  MCP_nMatchedPFPs.clear();
  MCP_totaldepE.clear();
  MCP_matchedPFP_ID.clear();
  MCP_matchedPFP_matchedE.clear();
  MCP_trueStartE.clear();
  MCP_trueStartP.clear();
  MCP_Px_eachpoint.clear();
  MCP_Py_eachpoint.clear();
  MCP_Pz_eachpoint.clear();
  MCP_x_eachpoint.clear();
  MCP_y_eachpoint.clear();
  MCP_z_eachpoint.clear();
  MCP_StartXYZ.clear();
  MCP_EndXYZ.clear();

  n_PFPs = 0;
  PFP_ID.clear();
  PFP_TrackShowerPdg.clear();
  PFP_IsPrimary.clear();
  PFP_ID_to_bestmatchMCPid.clear();
  PFP_totaldepE.clear();
  PFP_primaryPFPid.clear();
  PFP_spacepoints_XYZ.clear();  PFP_ordered_spacepoints.clear();    PFP_ordered_spacepoints_origidx.clear();
  PFP_goodreco_MCPvectorposes.clear();
  PFP_trueMCPend.clear();
  PFP_trajpoints_XYZ.clear();
  PFP_track_length.clear();
  PFP_track_start.clear();
  PFP_track_end.clear();
  PFP_track_theta.clear();
  PFP_track_phi.clear();



  // Fill variables for all MCPs
  nMCPs = mcp_v.size();
  for (auto mcp : mcp_v){

    MCP_ID.push_back(mcp->TrackId());
    MCP_MotherID.push_back(mcp->Mother());
    MCP_PDGcode.push_back(mcp->PdgCode());
    MCP_EndProcess.push_back(mcp->EndProcess());
    MCP_trueStartE.push_back(mcp->E());
    MCP_trueStartP.push_back(mcp->P());

    std::vector<double> px, py, pz;
    std::vector<double> x, y, z;
    for (size_t i_traj=0; i_traj<mcp->NumberTrajectoryPoints(); i_traj++){
      px.push_back(mcp->Px(i_traj));
      py.push_back(mcp->Py(i_traj));
      pz.push_back(mcp->Pz(i_traj));
      x.push_back(mcp->Position(i_traj).X());
      y.push_back(mcp->Position(i_traj).Y());
      z.push_back(mcp->Position(i_traj).Z());
    }
    MCP_Px_eachpoint.push_back(px);
    MCP_Py_eachpoint.push_back(py);
    MCP_Pz_eachpoint.push_back(pz);
    MCP_x_eachpoint.push_back(x);
    MCP_y_eachpoint.push_back(y);
    MCP_z_eachpoint.push_back(z);

    std::vector<int> daughters;
    for (int i_d=0; i_d<mcp->NumberDaughters(); i_d++){
      daughters.push_back(mcp->Daughter(i_d));
    }
    MCP_DaughterIDs.push_back(daughters);

    TLorentzVector start = mcp->Position();
    std::vector<double> startxyz{start.X(),start.Y(),start.Z()};
    MCP_StartXYZ.push_back(startxyz);

    TLorentzVector end = mcp->EndPosition();
    std::vector<double> endxyz{end.X(),end.Y(),end.Z()};
    MCP_EndXYZ.push_back(endxyz);

    // Get true deposited energy by MCParticle from truth-matching information
    std::vector<art::Ptr<recob::Hit>> hit_vec;
    std::vector<anab::BackTrackerHitMatchingData const*> match_vec;
    hits_per_mcp.get(mcp.key(),hit_vec,match_vec);
    double depE = 0;
    for (size_t i_matchhit=0; i_matchhit<match_vec.size(); i_matchhit++){
      depE+=match_vec.at(i_matchhit)->energy;
    }
    MCP_totaldepE.push_back(depE);

    // Look for matched PFPs
    auto search = MCPID_to_many_matchingData.find(mcp->TrackId());
    if (search != MCPID_to_many_matchingData.end()){
      if (fIsVerbose) std::cout << "particle with energy " << depE << " matched to the following " << search->second.size() <<  " PFPs: " << std::endl;

      MCP_nMatchedPFPs.push_back(search->second.size());
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
    MCP_matchedPFP_ID.push_back(matchedPFP_ID);
    MCP_matchedPFP_matchedE.push_back(matchedPFP_matchedE);

  }  // end loop over MCPs


  // Finally, fill variables related to PFPs in the event
  n_PFPs = pfpCollection.size();
  for (auto pfp : pfpCollection){

    // Is this a primary PFP? By that I mean, is its parent a neutrino PFP? (pfParticle->IsPrimary() means "is it a neutrino?" so that's not what we want)
    // Find the parent particle
    bool isprimarypfp=false;
    const auto parentIterator = particleMap.find(pfp->Parent());
    if (parentIterator == particleMap.end()){
      isprimarypfp = false;
    }
    else{
      // Check if the parent has a neutrino PDG code
      const int parentPDG = std::abs(parentIterator->second->PdgCode());
      if (parentIterator!=particleMap.end() && parentIterator->second->IsPrimary() && (parentPDG==pandora::NU_E || parentPDG==pandora::NU_MU || parentPDG==pandora::NU_TAU)){
        isprimarypfp = true;
      }
      else {
        isprimarypfp = false;
      }
    }


    // Is this PFP a good kink candidate? That means: is it a primary PFP, does it come from a pi+, and is it truth-matched to two particles that are mother-daughter related with 90% completeness?
    // If not, don't use it

    // 1) is it primary?
    if (!isprimarypfp) continue;

    // 2) is it truth-matched to a pi+ from a neutrino?
    int mcpid=-999;
    auto search1 = pfpID_to_matchingData_BestMatch.find((int)pfp->Self());
    if (search1 != pfpID_to_matchingData_BestMatch.end()){
      if (PDGCode(search1->second.MCParticle_PDG) != kPiPlus) continue;

      mcpid=search1->second.MCParticle_ID;
      art::Ptr<simb::MCParticle> matchedmcp;
      for (auto mcp : mcp_v){
        if (mcp->TrackId() == mcpid){
          matchedmcp = mcp;
          break;
        }
      }
      if (matchedmcp->Mother() != nu_MCPID) continue;
    }
    else{
      continue;
    }

    // 3) is it truth-matched to exactly two particles, of which it captures 90% of the energy?
    std::vector<int> MatchedMCP_vectorposes;
    // We have MCP->all matched PFP mapping but not the other way around, so loop through all MCPs and ask if they match this PFP
    for (size_t i_mcp=0; i_mcp < mcp_v.size(); i_mcp++){
      auto mcp_tmp = mcp_v.at(i_mcp);
      auto matchingdata_mcp_tmp = MCPID_to_many_matchingData.find((int)mcp_tmp->TrackId());
      if (matchingdata_mcp_tmp==MCPID_to_many_matchingData.end()) continue;
      for (size_t i_matchpfp=0; i_matchpfp<matchingdata_mcp_tmp->second.size(); i_matchpfp++){

        auto matchingdat = matchingdata_mcp_tmp->second.at(i_matchpfp);

        // Are we looking at the PFP we want? If not, move on to the next one
        if (matchingdat.recoObject_ID != (int)pfp->Self()) continue;

        // std::cout << "PFP " << pfp->Self() << ": MCP with PDG " << MCP_PDGcode.at(i_mcp) << " with matched energy " << matchingdat.MatchedEnergy << ". True daughter? " << bool(MCP_MotherID.at(i_mcp)==mcpid) << " Is primary pi+? " << bool(MCP_ID.at(i_mcp)==mcpid) << std::endl;

        // If we have the right PFP, does the matched energy constitute more than 90% of the MCP energy? If not, move on
        if (matchingdat.MatchedEnergy < 0.9*MCP_totaldepE.at(i_mcp)) continue;

        // Only want MCPs that are true muons/pions/protons(with p>30GeV) and daughters of original pi+ (or the original pi+)
        if (!(MCP_MotherID.at(i_mcp)==mcpid || MCP_ID.at(i_mcp)==mcpid)) continue;
        if (!(PDGCode(MCP_PDGcode.at(i_mcp))==kMuMinus || PDGCode(MCP_PDGcode.at(i_mcp))==kPiPlus || (PDGCode(MCP_PDGcode.at(i_mcp))==kProton && MCP_trueStartP.at(i_mcp)>0.03) )) continue;

        // std::cout << "  -- keep" << std::endl;

        MatchedMCP_vectorposes.push_back(i_mcp);

      } // end loop over i_matchpfp
    } // end loop over i_mcp

    // std::cout << " -- Matched " << MatchedMCP_vectorposes.size() << " MCParticles -- " << std::endl;
    if (MatchedMCP_vectorposes.size() != 2) continue;


    // 4) Are those two particles mother-daughter related?
    int mcp_vecpos_1 = MatchedMCP_vectorposes.at(0);
    int mcp_vecpos_2 = MatchedMCP_vectorposes.at(1);
    std::vector<double> trueendxyz{-999,-999,-999};
    if (mcp_v.at(mcp_vecpos_1)->Mother() == mcp_v.at(mcp_vecpos_2)->TrackId()){
      TLorentzVector end = mcp_v.at(mcp_vecpos_2)->EndPosition();
      trueendxyz = std::vector<double> {end.X(),end.Y(),end.Z()};
      // do nothing
    }
    else if (mcp_v.at(mcp_vecpos_2)->Mother() == mcp_v.at(mcp_vecpos_1)->TrackId()){
      TLorentzVector end = mcp_v.at(mcp_vecpos_1)->EndPosition();
      trueendxyz = std::vector<double> {end.X(),end.Y(),end.Z()};
      // do nothing
    }
    else{
      continue;
    }

    // If it passed all these requirements, save the PFP!
    PFP_IsPrimary.push_back(isprimarypfp);
    PFP_ID.push_back((int)pfp->Self());
    PFP_TrackShowerPdg.push_back(pfp->PdgCode());
    PFP_goodreco_MCPvectorposes.push_back(MatchedMCP_vectorposes);
    PFP_trueMCPend.push_back(trueendxyz);


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

    // Save ordered spacepoints
    std::vector<std::pair<std::vector<double>,int>> ordered_spacepoints_pair_origidx =  OrderSpacepoints(sp_tvector_vec);
    std::vector<std::vector<double>> ordered_spacepoints;
    std::vector<int> ordered_spacepoints_origidx;
    for (size_t i_sp=0; i_sp<ordered_spacepoints_pair_origidx.size(); i_sp++){
      ordered_spacepoints.push_back(ordered_spacepoints_pair_origidx.at(i_sp).first);
      ordered_spacepoints_origidx.push_back(ordered_spacepoints_pair_origidx.at(i_sp).second);
    }
    PFP_ordered_spacepoints.push_back(ordered_spacepoints);
    PFP_ordered_spacepoints_origidx.push_back(ordered_spacepoints_origidx);

    // Get track from PFP (if it exists)
    std::vector<art::Ptr<recob::Track>> tracks_pfp = tracks_from_pfps.at(pfp.key());
    double track_length = -999;
    double track_theta = -999;
    double track_phi = -999;
    std::vector<double> track_start = {-999,-999,-999};
    std::vector<double> track_end = {-999,-999,-999};
    std::vector<std::vector<double>> tr_tvector_vec;
    if(tracks_pfp.size() > 1) {
       mf::LogError(__PRETTY_FUNCTION__) << "PFP associated to more than one track." << std::endl;
       throw std::exception();
    }
    for (auto track : tracks_pfp){
       track_length = track -> Length();
       track_theta = track -> Theta();
       track_phi = track -> Phi();
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
    PFP_track_theta.push_back(track_theta);
    PFP_track_phi.push_back(track_phi);


  } // end loop over PFPs

  if (PFP_IsPrimary.size()>0) _outtree->Fill();

}

DEFINE_ART_MODULE(KinkFindingTreeCreator)
