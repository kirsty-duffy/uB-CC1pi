// Only define this stuff once
#ifndef FILLTREE_CXX
#define FILLTREE_CXX

// Only need to include header file
// All other includes should be in there
#include "FillTree.h"

cc1pianavars::cc1pianavars(fhicl::ParameterSet const &p){
   pset = p;
}


void cc1pianavars::Clear(){

   // Set default values (that will be used if the event isn't filled)

   isData = false;
   run_num = -9999;
   subrun_num = -9999;
   event_num = -9999;
//   pot = 0;

   cutflow.clear();
   isSelected = false;
   track_length.clear();
   shower_length.clear();
   NPFPs = -9999;
   NTracks = -9999;
   NShowers = -9999;
   Sel_PFP_isTrack.clear();
   Sel_PFP_isShower.clear();
   Sel_PFP_ID.clear();
   Sel_MCP_ID.clear();
   Sel_MCP_PDG.clear();
   Sel_MCP_E.clear();
   tpcobj_origin = -9999;
   tpcobj_origin_extra = -9999;

   MCP_PDG.clear();
   MCP_length.clear();
   MCP_process.clear();
   MCP_endprocess.clear();
   MCP_numdaughters.clear();
   MCP_P.clear();
   MCP_Px.clear();
   MCP_Py.clear();
   MCP_Pz.clear();
   MCP_E.clear();

   nu_vtxx.clear();
   nu_vtxy.clear();
   nu_vtxz.clear();
   nu_isCC.clear();
   nu_PDG.clear();
   nu_E.clear();

   MIPConsistency.clear();
   dqdx_trunc_uncalib.clear();

   PassesCC1piSelec = false;
   CC1piSelecFailureReason = "";
}

void cc1pianavars::endSubRun(art::SubRun &sr) {
   // Note: the entire subrun's POT is recorded in the tree for every event.
   // You must only add it once per subrun to get the correct number.

   art::Handle<sumdata::POTSummary> potsum_h;

   // MC
   if (!isData) {
      if(sr.getByLabel("generator", potsum_h)) {
         pot = potsum_h->totpot;
      }
      else {
         pot = 0;
         //Should this raise an error?
      }
   }

   // Data
   else {
      if(sr.getByLabel("beamdata", "bnbETOR860", potsum_h)){
         pot = potsum_h->totpot;
      }
      else {
         pot = 0;
      }
   }

}

void cc1pianavars::SetReco2Vars(art::Event &evt){

   // Set all the values that are in the reco2 file
   // This just cleans the module up - move all these lines to here instead of in the module

   isData = evt.isRealData();
   run_num = evt.run();
   subrun_num = evt.subRun();
   event_num = evt.event();


   art::Handle<std::vector<ubana::SelectionResult>> selection_h;
   evt.getByLabel("UBXSec", selection_h);
   if(!selection_h.isValid()){
      mf::LogError(__PRETTY_FUNCTION__) << "SelectionResult product not found." << std::endl;
      throw std::exception();
   }
   std::vector<art::Ptr<ubana::SelectionResult>> selection_v;
   art::fill_ptr_vector(selection_v, selection_h);

   cutflow = selection_v.at(0)->GetCutFlowStatus();

   if (selection_v.at(0)->GetSelectionStatus()) {

      isSelected = true;

      //Get TPCObject
      art::FindManyP<ubana::TPCObject> tpcobject_from_selection(selection_h, evt, "UBXSec");
      art::Ptr<ubana::TPCObject> tpcobj_candidate = tpcobject_from_selection.at(0).at(0);
      art::Handle<std::vector<ubana::TPCObject>> tpcobj_h;
      evt.getByLabel("TPCObjectMaker", tpcobj_h);

      tpcobj_origin = tpcobj_candidate -> GetOrigin();
      tpcobj_origin_extra = tpcobj_candidate -> GetOriginExtra();

      //Get PFPs (in TPCObject)
      art::FindManyP<recob::PFParticle> pfps_from_tpcobject(tpcobj_h, evt, "TPCObjectMaker");
      std::vector<art::Ptr<recob::PFParticle>> pfps = pfps_from_tpcobject.at(tpcobj_candidate.key());
      NPFPs = pfps.size();

      //Get tracks (in TPCObject)
      art::FindManyP<recob::Track> tracks_from_tpcobject(tpcobj_h, evt, "TPCObjectMaker");
      std::vector<art::Ptr<recob::Track>> tracks = tracks_from_tpcobject.at(tpcobj_candidate.key());
      NTracks = tracks.size();

      // Get calo objects (for dqdx)
      // For now use uncalibrated "calo" variables
      // Calibrated dqdx will be available in MCC 8.6, but currently has a bug
      // Use association: pandoraNucalo, art::Assns<recob::Track,anab::Calorimetry,void>
      // (Couldn't find another way to do this except going back to track handle -.-)
      // also hardcode that this is for pandoraNu only
      art::InputTag _caloTag = "pandoraNucalo";
      art::InputTag _trackTag = "pandoraNu";
      auto const& track_h = evt.getValidHandle<std::vector<recob::Track>>(_trackTag);
      art::FindManyP<anab::Calorimetry> calos_from_tracks(track_h, evt, _caloTag);

      for (auto track : tracks) {
         track_length.emplace_back(track -> Length());

         unsigned int trkid = track->ID();
         std::vector<art::Ptr<anab::Calorimetry>> calos = calos_from_tracks.at(trkid);
         double dqdx_truncmean = GetDqDxTruncatedMean(calos); // this function is in MIPConsistency_Marco
         dqdx_trunc_uncalib.emplace_back(dqdx_truncmean);

         MIPConsistency.emplace_back(IsMIP(track->Length(), dqdx_truncmean));
      }

      //Get showers (in TPCObject)
      art::FindManyP<recob::Shower> showers_from_tpcobject(tpcobj_h, evt, "TPCObjectMaker");
      std::vector<art::Ptr<recob::Shower>> showers = showers_from_tpcobject.at(tpcobj_candidate.key());
      NShowers = showers.size();

      for (auto shower : showers) {
         shower_length.emplace_back(shower -> Length());
      }

      //Get MCGhosts (in event)
      art::Handle<std::vector<ubana::MCGhost> > ghost_h;
      evt.getByLabel("RecoTrueMatcher",ghost_h);
      if(!ghost_h.isValid()){
         mf::LogError(__PRETTY_FUNCTION__) << "MCGhost product not found." << std::endl;
         throw std::exception();
      }

      //Get PFPs (in event)
      art::Handle<std::vector<recob::PFParticle> > pfp_h;
      evt.getByLabel("pandoraNu::UBXSec",pfp_h);
      if(!pfp_h.isValid()){
         mf::LogError(__PRETTY_FUNCTION__) << "PFP product not found." << std::endl;
         throw std::exception();
      }

      //Get associations
      //Could just use "pfps" and not bother also getting all the non-selected PFPs from the event?
      art::FindManyP<ubana::MCGhost>   mcghost_from_pfp (pfp_h,   evt, "RecoTrueMatcher");
      art::FindManyP<simb::MCParticle> mcp_from_mcghost (ghost_h, evt, "RecoTrueMatcher");

      for (auto pfp : pfps) {

         Sel_PFP_isTrack.emplace_back(lar_pandora::LArPandoraHelper::IsTrack(pfp));
         Sel_PFP_isShower.emplace_back(lar_pandora::LArPandoraHelper::IsShower(pfp));
         Sel_PFP_ID.emplace_back(pfp -> Self());

         auto mcghosts = mcghost_from_pfp.at(pfp.key());
         if (mcghosts.size() != 1) {
            Sel_MCP_ID.emplace_back(-9999);
            Sel_MCP_PDG.emplace_back(-9999);
            Sel_MCP_E.emplace_back(-9999);
            continue;
         }
         auto mcghost = mcghosts[0];

         auto mcps = mcp_from_mcghost.at(mcghost.key());
         if(mcps.size() != 1) {
            Sel_MCP_ID.emplace_back(-9999);
            Sel_MCP_PDG.emplace_back(-9999);
            Sel_MCP_E.emplace_back(-9999);
            continue;
         }
         auto mcp = mcps[0];

         Sel_MCP_ID.emplace_back(mcp -> TrackId());
         Sel_MCP_PDG.emplace_back(mcp -> PdgCode());
         Sel_MCP_E.emplace_back(mcp -> E());

      }

   }

   // Get all MCParticles
   art::Handle<std::vector<simb::MCParticle>> mcp_h;
   evt.getByLabel("largeant", mcp_h);
   std::vector<art::Ptr<simb::MCParticle>> mcp_v;
   art::fill_ptr_vector(mcp_v, mcp_h);

   bool recorded = false;

   //backtracker replacement boilerplate
   std::unique_ptr<truth::IMCTruthMatching> fMCTruthMatching;
   const fhicl::ParameterSet& truthParams = pset.get<fhicl::ParameterSet>("MCTruthMatching");
   fMCTruthMatching = std::unique_ptr<truth::IMCTruthMatching>(new truth::AssociationsTruth(truthParams));
   fMCTruthMatching->Rebuild(evt);

   //Loop over MCParticles...
   for (auto mcpar : mcp_v) {

      art::Ptr<simb::MCTruth> mc_truth = fMCTruthMatching -> TrackIDToMCTruth(mcpar -> TrackId());
      if (mc_truth -> Origin() != simb::kBeamNeutrino) continue;
      if (mcpar -> Mother() != 0) continue;
      if (mcpar -> StatusCode() != 1) continue;

      MCP_PDG.emplace_back(mcpar -> PdgCode());

      MCP_P.emplace_back(mcpar -> P());
      MCP_Px.emplace_back(mcpar -> Px());
      MCP_Py.emplace_back(mcpar -> Py());
      MCP_Pz.emplace_back(mcpar -> Pz());
      MCP_E.emplace_back(mcpar -> E());

      simb::MCTrajectory traj = mcpar -> Trajectory();
      MCP_length.emplace_back(traj.TotalLength());

      MCP_process.emplace_back(mcpar -> Process());
      MCP_endprocess.emplace_back(mcpar -> EndProcess());
      MCP_numdaughters.emplace_back(mcpar -> NumberDaughters());

      //Record neutrino info (once)
      if(!recorded) {

         if (!(mc_truth -> NeutrinoSet())) continue;

         recorded = true;

         simb::MCNeutrino nu = mc_truth -> GetNeutrino();
         simb::MCParticle neutrino = nu.Nu();

         if (nu.CCNC() == 0) nu_isCC.emplace_back(true);
         else if (nu.CCNC() == 1) nu_isCC.emplace_back(false);

         nu_PDG.emplace_back(neutrino.PdgCode());

         nu_E.emplace_back(neutrino.E());

         const TLorentzVector& vertex = neutrino.Position(0);
         double MC_vertex[4];
         vertex.GetXYZT(MC_vertex);
         nu_vtxx.emplace_back(MC_vertex[0]);
         nu_vtxy.emplace_back(MC_vertex[1]);
         nu_vtxz.emplace_back(MC_vertex[2]);

      }
   }
}


void MakeAnaBranches(TTree *t, cc1pianavars *vars){

   // Note: it's very important that we use the syntax &(vars->evtnum) to make sure
   // we get the value from reference, otherwise it won't work!

   t -> Branch("isData", &(vars->isData));
   t -> Branch("run_num", &(vars->run_num));
   t -> Branch("subrun_num", &(vars->subrun_num));
   t -> Branch("event_num", &(vars->event_num));
   t -> Branch("pot", &(vars->pot));

   t -> Branch("cutflow", &(vars->cutflow));
   t -> Branch("isSelected", &(vars->isSelected));
   t -> Branch("track_length", &(vars->track_length));
   t -> Branch("shower_length", &(vars->shower_length));
   t -> Branch("NPFPs", &(vars->NPFPs));
   t -> Branch("NTracks", &(vars->NTracks));
   t -> Branch("NShowers", &(vars->NShowers));
   t -> Branch("Sel_PFP_isTrack", &(vars->Sel_PFP_isTrack));
   t -> Branch("Sel_PFP_isShower", &(vars->Sel_PFP_isShower));
   t -> Branch("Sel_PFP_ID", &(vars->Sel_PFP_ID));
   t -> Branch("Sel_MCP_ID", &(vars->Sel_MCP_ID));
   t -> Branch("Sel_MCP_PDG", &(vars->Sel_MCP_PDG));
   t -> Branch("Sel_MCP_E", &(vars->Sel_MCP_E));
   t -> Branch("tpcobj_origin", &(vars->tpcobj_origin));
   t -> Branch("tpcobj_origin_extra", &(vars->tpcobj_origin_extra));

   t -> Branch("MCP_PDG", &(vars->MCP_PDG));
   t -> Branch("MCP_length", &(vars->MCP_length));
   t -> Branch("MCP_process", &(vars->MCP_process));
   t -> Branch("MCP_endprocess", &(vars->MCP_endprocess));
   t -> Branch("MCP_numdaughters", &(vars->MCP_numdaughters));
   t -> Branch("MCP_P", &(vars->MCP_P));
   t -> Branch("MCP_Px", &(vars->MCP_Px));
   t -> Branch("MCP_Py", &(vars->MCP_Py));
   t -> Branch("MCP_Pz", &(vars->MCP_Pz));
   t -> Branch("MCP_E", &(vars->MCP_E));

   t -> Branch("nu_vtxx", &(vars->nu_vtxx));
   t -> Branch("nu_vtxy", &(vars->nu_vtxy));
   t -> Branch("nu_vtxz", &(vars->nu_vtxz));
   t -> Branch("nu_isCC", &(vars->nu_isCC));
   t -> Branch("nu_PDG", &(vars->nu_PDG));
   t -> Branch("nu_E", &(vars->nu_E));

   t -> Branch("MIPConsistency", &(vars->MIPConsistency));
   t -> Branch("dqdx_trunc_uncalib", &(vars->dqdx_trunc_uncalib));

   t -> Branch("PassesCC1piSelec", &(vars->PassesCC1piSelec));
   t -> Branch("CC1piSelecFailureReason", &(vars->CC1piSelecFailureReason));

}

#endif
