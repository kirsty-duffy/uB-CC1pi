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

   Marco_cutflow.clear();
   Marco_selected = false;

   Truth_topology = kUnknown;
   
   TPCObj_PFP_track_length.clear();
   TPCObj_PFP_track_start.clear();
   TPCObj_PFP_track_end.clear();
   TPCObj_PFP_track_dqdx_truncmean.clear();
   TPCObj_PFP_isMIP.clear();
   TPCObj_PFP_shower_length.clear();
   TPCObj_PFP_shower_start.clear();
   TPCObj_NPFPs = -9999;
   TPCObj_NTracks = -9999;
   TPCObj_NShowers = -9999;
   TPCObj_PFP_isTrack.clear();
   TPCObj_PFP_isShower.clear();
   TPCObj_PFP_isDaughter.clear();
   TPCObj_PFP_id.clear();
   TPCObj_PFP_MCPid.clear();
   TPCObj_PFP_truePDG.clear();
   TPCObj_PFP_trueE.clear();
   TPCObj_PFP_trueKE.clear();
   TPCObj_origin = -9999;
   TPCObj_origin_extra = -9999;
   TPCObj_reco_vtx.clear();

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
   MCP_KE.clear();
   MCP_isContained.clear();

   nu_vtx.clear();
   nu_vtx_spacecharge.clear();
   nu_isCC = false;
   nu_PDG = -9999;
   nu_E = -9999;

   CC1picutflow.clear();
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

   Marco_cutflow = selection_v.at(0)->GetCutFlowStatus();

   // Get topology from MC truth (MC only)
   // Do this for all events (not just those selected by Marco)
   // For events not selected by Marco we don't have TPCObj_origin information so we
   // can only set neutrino interaction topology (including kUnknown) or outFV
   // Events that are selected and cosmic/mixed/unknown origin will be overwritten later
   if (!isData){
     
     auto const mctruth_h = evt.getValidHandle<std::vector<simb::MCTruth>>("generator"); // Get only GENIE MCtruth
     simb::MCTruth top_mctruth = mctruth_h->at(0);
     simb::MCNeutrino top_nu = top_mctruth.GetNeutrino();
     simb::MCParticle top_neutrino = top_nu.Nu();

     const TLorentzVector& top_vertex = top_neutrino.Position(0);
     double top_vertexXYZT[4];
     top_vertex.GetXYZT(top_vertexXYZT);


     if(inFV(top_vertexXYZT[0], top_vertexXYZT[1], top_vertexXYZT[2])){ // If true vertex is in the fiducial volume
        Truth_topology = GetTopology(mctruth_h);
     }
     else{ // outFV topology
        Truth_topology = kOutFV; 
     }
   } // if !isData


   // Now fill most variables only if the event is selected by Marco
   if (selection_v.at(0)->GetSelectionStatus()) {

      Marco_selected = true;

      //Get TPCObject
      art::FindManyP<ubana::TPCObject> TPCObject_from_selection(selection_h, evt, "UBXSec");
      art::Ptr<ubana::TPCObject> TPCObj_candidate = TPCObject_from_selection.at(0).at(0);
      art::Handle<std::vector<ubana::TPCObject>> TPCObj_h;
      evt.getByLabel("TPCObjectMaker", TPCObj_h);

      TPCObj_origin = TPCObj_candidate -> GetOrigin();
      TPCObj_origin_extra = TPCObj_candidate -> GetOriginExtra();

      // Get topology from MC truth (MC only)
      // Overwrite what was done above for selected cosmic/mixed/unknown origin events only
      if (!isData){

         // For selected event, check TPC object origin for cosmic/mixed/unknown
         // Only overwrite events that don't have TPC_origin == 0
         // For TPC_origin == 0 events we want to keep the Truth_topology assigned above
         if (TPCObj_origin == 1){
            Truth_topology = kCosmic;
         }
         else if (TPCObj_origin == 2){
            Truth_topology = kMixed;
         }
         else if (TPCObj_origin != 0){ // if selected object has some other non-neutrino origin, set it to unknown
            Truth_topology = kUnknown;
         }
      } // if !isData


      //Get PFPs (in TPCObject)
      art::FindManyP<recob::PFParticle> pfps_from_TPCObject(TPCObj_h, evt, "TPCObjectMaker");
      std::vector<art::Ptr<recob::PFParticle>> pfps = pfps_from_TPCObject.at(TPCObj_candidate.key());
      TPCObj_NPFPs = pfps.size();

      //Get tracks (in TPCObject)
      art::FindManyP<recob::Track> tracks_from_TPCObject(TPCObj_h, evt, "TPCObjectMaker");
      std::vector<art::Ptr<recob::Track>> tracks = tracks_from_TPCObject.at(TPCObj_candidate.key());
      TPCObj_NTracks = tracks.size();

      // Get calo objects (for dqdx)
      // For now use uncalibrated "calo" variables
      // Calibrated dqdx will be available in MCC 8.6, but currently has a bug
      // also hardcode that this is for pandoraNu only
      // Use association: pandoraNucalo, art::Assns<recob::Track,anab::Calorimetry,void>
      // (Couldn't find another way to do this except going back to track handle -.-)
      art::InputTag _caloTag = "pandoraNucalo";
      art::InputTag _trackTag = "pandoraNu";
      auto const& track_h = evt.getValidHandle<std::vector<recob::Track>>(_trackTag);
      art::FindManyP<anab::Calorimetry> calos_from_tracks(track_h, evt, _caloTag);

      //Get showers (in TPCObject)
      art::FindManyP<recob::Shower> showers_from_TPCObject(TPCObj_h, evt, "TPCObjectMaker");
      std::vector<art::Ptr<recob::Shower>> showers = showers_from_TPCObject.at(TPCObj_candidate.key());
      TPCObj_NShowers = showers.size();

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
      art::FindManyP<ubana::MCGhost>   mcghost_from_pfp (pfp_h,   evt, "RecoTrueMatcher");
      art::FindManyP<simb::MCParticle> mcp_from_mcghost (ghost_h, evt, "RecoTrueMatcher");

      //Find neutrino
      int nuID = -1;
      for (auto pfp : pfps) {
         if(lar_pandora::LArPandoraHelper::IsNeutrino(pfp)) {
            nuID = (int)pfp->Self();
            break;
         }
      }

      art::FindManyP<recob::Track> tracks_from_pfps(pfp_h, evt, "pandoraNu::UBXSec");
      art::FindManyP<recob::Shower> showers_from_pfps(pfp_h, evt, "pandoraNu::UBXSec");

      for (auto pfp : pfps) {

         TPCObj_PFP_isTrack.emplace_back(lar_pandora::LArPandoraHelper::IsTrack(pfp));
         TPCObj_PFP_isShower.emplace_back(lar_pandora::LArPandoraHelper::IsShower(pfp));
         TPCObj_PFP_isDaughter.emplace_back(pfp->Parent()==(size_t)nuID);
         TPCObj_PFP_id.emplace_back(pfp -> Self());

         //Set default values for track/shower specific variables
         double track_length = -9999;
         std::vector<double> track_start = {-9999, -9999, -9999};
         std::vector<double> track_end = {-9999, -9999, -9999};
         double track_dqdx_truncmean = -9999;
         bool isMIP = false;
         double shower_length = -9999;
         std::vector<double> shower_start = {-9999, -9999, -9999};

//         if(lar_pandora::LArPandoraHelper::IsTrack(pfp)) {
            std::vector<art::Ptr<recob::Track>> tracks_pfp = tracks_from_pfps.at(pfp.key());
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

               unsigned int trkid = track->ID();
               std::vector<art::Ptr<anab::Calorimetry>> calos = calos_from_tracks.at(trkid);
               track_dqdx_truncmean = GetDqDxTruncatedMean(calos); // this function is in MIPConsistency_Marco

               // Apply David's calibration: constant factor applied to trunc mean dqdx
               // to convert it to electrons/cm, and compare MC and data
               // Multiply MC by 198 and data by 243
               // !!! Note that this should change when real calibration is available !!!
               if (evt.isRealData()){ // Data: multiply by 243
                  track_dqdx_truncmean *= 243.;
               }
               else{ // MC: multiply by 198
                  track_dqdx_truncmean *= 198.;
               }

               isMIP = IsMIP(track_length, track_dqdx_truncmean);
            }

//         } // IsTrack

//         else if(lar_pandora::LArPandoraHelper::IsShower(pfp)) {
            std::vector<art::Ptr<recob::Shower>> showers_pfp = showers_from_pfps.at(pfp.key());
            if(showers_pfp.size() > 1) {
               mf::LogError(__PRETTY_FUNCTION__) << "PFP associated to more than one showers." << std::endl;
               throw std::exception();
            }
            for (auto shower : showers_pfp) {
               shower_length = shower -> Length();
               auto start = shower -> ShowerStart();
               shower_start = {start.X(),start.Y(),start.Z()};
            }

//         } // IsShower

         // Fill track/shower specific variables
         TPCObj_PFP_track_length.emplace_back(track_length);
         TPCObj_PFP_track_start.emplace_back(track_start);
         TPCObj_PFP_track_end.emplace_back(track_end);
         TPCObj_PFP_track_dqdx_truncmean.emplace_back(track_dqdx_truncmean);
         TPCObj_PFP_isMIP.emplace_back(isMIP);
         TPCObj_PFP_shower_length.emplace_back(shower_length);
         TPCObj_PFP_shower_start.emplace_back(shower_start);


         // Do reco-truth matching...

         auto mcghosts = mcghost_from_pfp.at(pfp.key());
         if (mcghosts.size() != 1) {
            TPCObj_PFP_MCPid.emplace_back(-9999);
            TPCObj_PFP_truePDG.emplace_back(-9999);
            TPCObj_PFP_trueE.emplace_back(-9999);
            TPCObj_PFP_trueKE.emplace_back(-9999);
            continue;
         }
         auto mcghost = mcghosts[0];

         auto mcps = mcp_from_mcghost.at(mcghost.key());
         if(mcps.size() != 1) {
            TPCObj_PFP_MCPid.emplace_back(-9999);
            TPCObj_PFP_truePDG.emplace_back(-9999);
            TPCObj_PFP_trueE.emplace_back(-9999);
            TPCObj_PFP_trueKE.emplace_back(-9999);
            continue;
         }
         auto mcp = mcps[0];

         TPCObj_PFP_MCPid.emplace_back(mcp -> TrackId());
         TPCObj_PFP_truePDG.emplace_back(mcp -> PdgCode());
         TPCObj_PFP_trueE.emplace_back(mcp -> E());
         TPCObj_PFP_trueKE.emplace_back(mcp -> E() - mcp -> Mass());

      } // loop over pfps


      // Get neutrino candidate vertex from TPCObject
      recob::Vertex TPCObj_nu_vtx = TPCObj_candidate->GetVertex();
      double reco_nu_vtx[3];
      TPCObj_nu_vtx.XYZ(reco_nu_vtx);
      TPCObj_reco_vtx = {reco_nu_vtx[0], reco_nu_vtx[1], reco_nu_vtx[2]};
      // We eventually need to correct for X position (time offset) and space charge
      // See issues #8 and #10 on github

   } // if selected by Marco

   // Get all MCParticles
   art::Handle<std::vector<simb::MCParticle>> mcp_h;
   evt.getByLabel("largeant", mcp_h);
   std::vector<art::Ptr<simb::MCParticle>> mcp_v;
   art::fill_ptr_vector(mcp_v, mcp_h);

   bool nu_recorded = false;

   // backtracker replacement boilerplate
   std::unique_ptr<truth::IMCTruthMatching> fMCTruthMatching;
   const fhicl::ParameterSet& truthParams = pset.get<fhicl::ParameterSet>("MCTruthMatching");
   fMCTruthMatching = std::unique_ptr<truth::IMCTruthMatching>(new truth::AssociationsTruth(truthParams));
   fMCTruthMatching->Rebuild(evt);

   // Loop over all MCParticles in the event...
   // Note: Only filled for particles from true beam neutrinos!
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
      MCP_KE.emplace_back(mcpar -> E() - mcpar -> Mass());

      simb::MCTrajectory traj = mcpar -> Trajectory();
      MCP_length.emplace_back(traj.TotalLength());
      if(inFV(traj.X(0),traj.Y(0),traj.Z(0)) && inFV(traj.X(traj.size()-1),traj.Y(traj.size()-1),traj.Z(traj.size()-1))) MCP_isContained.emplace_back(true);
      else MCP_isContained.emplace_back(false);

      MCP_process.emplace_back(mcpar -> Process());
      MCP_endprocess.emplace_back(mcpar -> EndProcess());
      MCP_numdaughters.emplace_back(mcpar -> NumberDaughters());

      //Record neutrino info (once)
      //Note: Potential issues for events with multiple neutrinos
      if(!nu_recorded) {

         if (!(mc_truth -> NeutrinoSet())) continue;

         nu_recorded = true;

         simb::MCNeutrino nu = mc_truth -> GetNeutrino();
         simb::MCParticle neutrino = nu.Nu();

         if (nu.CCNC() == 0) nu_isCC = true;
         else if (nu.CCNC() == 1) nu_isCC = false;

         nu_PDG = neutrino.PdgCode();

         nu_E = neutrino.E();

         const TLorentzVector& vertex = neutrino.Position(0);
         double MC_vertex[4];
         vertex.GetXYZT(MC_vertex);
         nu_vtx = {MC_vertex[0], MC_vertex[1], MC_vertex[2]};

         // Space Charge correction
         auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
         std::vector<double> sce_corr = SCE->GetPosOffsets(nu_vtx[0], nu_vtx[1], nu_vtx[2]);
         std::cout << "SCE correction in x, y, z = " << sce_corr.at(0)
            << ", " << sce_corr.at(1) 
            << ", " << sce_corr.at(2) << std::endl;
         nu_vtx_spacecharge = {nu_vtx[0] - sce_corr.at(0), nu_vtx[1] + sce_corr.at(1), nu_vtx[2] + sce_corr.at(2)};

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

   t -> Branch("Marco_cutflow", &(vars->Marco_cutflow));
   t -> Branch("Marco_selected", &(vars->Marco_selected));

   t -> Branch("Truth_topology", &(vars->Truth_topology), "Truth_topology/I");
   
   t -> Branch("TPCObj_PFP_track_length", &(vars->TPCObj_PFP_track_length));
   t -> Branch("TPCObj_PFP_track_start", &(vars->TPCObj_PFP_track_start));
   t -> Branch("TPCObj_PFP_track_end", &(vars->TPCObj_PFP_track_end));
   t -> Branch("TPCObj_PFP_track_dqdx_truncmean", &(vars->TPCObj_PFP_track_dqdx_truncmean));
   t -> Branch("TPCObj_PFP_isMIP", &(vars->TPCObj_PFP_isMIP));
   t -> Branch("TPCObj_PFP_shower_length", &(vars->TPCObj_PFP_shower_length));
   t -> Branch("TPCObj_PFP_shower_start", &(vars->TPCObj_PFP_shower_start));
   t -> Branch("TPCObj_NPFPs", &(vars->TPCObj_NPFPs));
   t -> Branch("TPCObj_NTracks", &(vars->TPCObj_NTracks));
   t -> Branch("TPCObj_NShowers", &(vars->TPCObj_NShowers));
   t -> Branch("TPCObj_PFP_isTrack", &(vars->TPCObj_PFP_isTrack));
   t -> Branch("TPCObj_PFP_isShower", &(vars->TPCObj_PFP_isShower));
   t -> Branch("TPCObj_PFP_isDaughter", &(vars->TPCObj_PFP_isDaughter));
   t -> Branch("TPCObj_PFP_id", &(vars->TPCObj_PFP_id));
   t -> Branch("TPCObj_PFP_MCPid", &(vars->TPCObj_PFP_MCPid));
   t -> Branch("TPCObj_PFP_truePDG", &(vars->TPCObj_PFP_truePDG));
   t -> Branch("TPCObj_PFP_trueE", &(vars->TPCObj_PFP_trueE));
   t -> Branch("TPCObj_PFP_trueKE", &(vars->TPCObj_PFP_trueKE));
   t -> Branch("TPCObj_origin", &(vars->TPCObj_origin));
   t -> Branch("TPCObj_origin_extra", &(vars->TPCObj_origin_extra));
   t -> Branch("TPCObj_reco_vtx", &(vars->TPCObj_reco_vtx));

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
   t -> Branch("MCP_KE", &(vars->MCP_KE));
   t -> Branch("MCP_isContained", &(vars->MCP_isContained));

   t -> Branch("nu_vtx", &(vars->nu_vtx));
   t -> Branch("nu_vtx_spacecharge", &(vars->nu_vtx_spacecharge));
   t -> Branch("nu_isCC", &(vars->nu_isCC));
   t -> Branch("nu_PDG", &(vars->nu_PDG));
   t -> Branch("nu_E", &(vars->nu_E));

   t -> Branch("CC1picutflow", &(vars->CC1picutflow));

}

#endif
