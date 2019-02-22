// Only define this stuff once
#ifndef FILLTREE_CXX
#define FILLTREE_CXX

// Only need to include header file
// All other includes should be in there
#include "FillTree.h"

cc1pianavars::cc1pianavars(fhicl::ParameterSet const &p){
   pset = p;

   CC1piInputTags = new InputTags(p);

  fhicl::ParameterSet p_bragg = p.get<fhicl::ParameterSet>("BraggAlgoConfig");

  braggcalc.configure(p_bragg);
  braggcalc.checkRange=false;
  braggcalc.nHitsToDrop = 0;
  braggcalc.printConfiguration();
}


void cc1pianavars::Clear(){

   // Set default values (that will be used if the event isn't filled)

   // Event variables
   isData = false;
   run_num = -9999;
   subrun_num = -9999;
   event_num = -9999;

   Truth_topology = kUnknown;

   // Selection result variables
   Marco_cutflow.clear();
   Marco_selected = false;
   CC1picutflow.clear();

   // Track variables
   TPCObj_PFP_track_length.clear();
   TPCObj_PFP_track_start.clear();
   TPCObj_PFP_track_end.clear();
   TPCObj_PFP_track_theta.clear();
   TPCObj_PFP_track_phi.clear();
   TPCObj_PFP_track_mom.clear();
   TPCObj_PFP_track_dedx_truncmean.clear();
   TPCObj_PFP_track_dedx_perhit.clear();
   TPCObj_PFP_track_nhits.clear();
   TPCObj_PFP_track_resrange_perhit.clear();
   TPCObj_PFP_isMIP.clear();
   TPCObj_PFP_track_isContained.clear();
   TPCObj_PFP_track_AngleBetweenTracks.clear();
   TPCObj_PFP_track_residual_mean.clear();
   TPCObj_PFP_track_residual_std.clear();
   TPCObj_PFP_track_perc_used_hits.clear();
   TPCObj_PFP_track_MCSmu_fwdMom.clear();
   TPCObj_PFP_track_MCSmu_bwdMom.clear();
   TPCObj_PFP_track_MCSmu_bestMom.clear();
   TPCObj_PFP_track_MCSmu_fwdMomUncert.clear();
   TPCObj_PFP_track_MCSmu_bwdMomUncert.clear();
   TPCObj_PFP_track_MCSmu_bestMomUncert.clear();
   TPCObj_PFP_track_MCSmu_fwdLL.clear();
   TPCObj_PFP_track_MCSmu_bwdLL.clear();
   TPCObj_PFP_track_MCSmu_bestLL.clear();
   TPCObj_PFP_track_MCSmu_segmentRadLengths.clear();
   TPCObj_PFP_track_MCSmu_scatterAngles.clear();
   TPCObj_PFP_track_MCSp_fwdMom.clear();
   TPCObj_PFP_track_MCSp_bwdMom.clear();
   TPCObj_PFP_track_MCSp_bestMom.clear();
   TPCObj_PFP_track_MCSp_fwdMomUncert.clear();
   TPCObj_PFP_track_MCSp_bwdMomUncert.clear();
   TPCObj_PFP_track_MCSp_bestMomUncert.clear();
   TPCObj_PFP_track_MCSp_fwdLL.clear();
   TPCObj_PFP_track_MCSp_bwdLL.clear();
   TPCObj_PFP_track_MCSp_bestLL.clear();
   TPCObj_PFP_track_MCSp_segmentRadLengths.clear();
   TPCObj_PFP_track_MCSp_scatterAngles.clear();
   TPCObj_PFP_track_MCSpi_fwdMom.clear();
   TPCObj_PFP_track_MCSpi_bwdMom.clear();
   TPCObj_PFP_track_MCSpi_bestMom.clear();
   TPCObj_PFP_track_MCSpi_fwdMomUncert.clear();
   TPCObj_PFP_track_MCSpi_bwdMomUncert.clear();
   TPCObj_PFP_track_MCSpi_bestMomUncert.clear();
   TPCObj_PFP_track_MCSpi_fwdLL.clear();
   TPCObj_PFP_track_MCSpi_bwdLL.clear();
   TPCObj_PFP_track_MCSpi_bestLL.clear();
   TPCObj_PFP_track_MCSpi_segmentRadLengths.clear();
   TPCObj_PFP_track_MCSpi_scatterAngles.clear();
   TPCObj_PFP_track_SpacepointsXYZ.clear();
   TPCObj_PFP_track_SpacepointsQPlane2.clear();
   TPCObj_PFP_track_unusedhits_charge_plane0.clear();
   TPCObj_PFP_track_unusedhits_endwiredist_plane0.clear();
   TPCObj_PFP_track_unusedhits_endtimedist_plane0.clear();
   TPCObj_PFP_track_unusedhits_charge_plane1.clear();
   TPCObj_PFP_track_unusedhits_endwiredist_plane1.clear();
   TPCObj_PFP_track_unusedhits_endtimedist_plane1.clear();
   TPCObj_PFP_track_unusedhits_charge_plane2.clear();
   TPCObj_PFP_track_unusedhits_endwiredist_plane2.clear();
   TPCObj_PFP_track_unusedhits_endtimedist_plane2.clear();

   // Shower variables
   TPCObj_PFP_shower_length.clear();
   TPCObj_PFP_shower_start.clear();

   // PFP reco variables
   TPCObj_NPFPs = -9999;
   TPCObj_NTracks = -9999;
   TPCObj_NShowers = -9999;
   TPCObj_PFP_PandoraClassedAsTrack.clear();
   TPCObj_PFP_PandoraClassedAsShower.clear();
   TPCObj_PFP_isDaughter.clear();
   TPCObj_PFP_daughterids.clear();
   TPCObj_PFP_id.clear();
   TPCObj_reco_vtx.clear();
   TPCObj_reco_vtx_nospacecharge.clear();

   // PFP true variables
   TPCObj_PFP_MCPid.clear();
   TPCObj_PFP_truePDG.clear();
   TPCObj_PFP_trueE.clear();
   TPCObj_PFP_trueKE.clear();
   TPCObj_PFP_trueEndP.clear();
   TPCObj_origin = -9999;
   TPCObj_origin_extra = -9999;

   // PID variables
   TPCObj_PFP_LH_fwd_mu.clear();
   TPCObj_PFP_LH_fwd_p.clear();
   TPCObj_PFP_LH_fwd_pi.clear();
   TPCObj_PFP_LH_bwd_mu.clear();
   TPCObj_PFP_LH_bwd_p.clear();
   TPCObj_PFP_LH_bwd_pi.clear();
   TPCObj_PFP_LH_MIP.clear();
   TPCObj_PFP_PIDA.clear();
   TPCObj_PFP_track_depE.clear();
   TPCObj_PFP_track_Chi2Proton.clear();
   TPCObj_PFP_track_Chi2Muon.clear();
   TPCObj_PFP_track_Chi2Pion.clear();
   TPCObj_PFP_track_rangeE_mu.clear();
   TPCObj_PFP_track_rangeE_p.clear();
   TPCObj_PFP_track_Lmip_perhit.clear();

   // MCParticle variables
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
   MCP_ID.clear();
   MCP_DaughterIDs.clear();
   MCP_StatusCode.clear();
   MCP_MotherID.clear();
   MCP_StartPosition.clear();
   MCP_EndPosition.clear();

   // True neutrino variables
   nu_vtx.clear();
   nu_vtx_spacecharge.clear();
   nu_isCC = false;
   nu_PDG = -9999;
   nu_E = -9999;
   nuint_mode = -9999;
   nuint_W = -9999;
   nuint_Qsq = -9999;

}

void cc1pianavars::SetReco2Vars(art::Event &evt){
   // Set all the values that are in the reco2 file
   // This just cleans the module up - move all these lines to here instead of in the module

   isData = evt.isRealData();
   run_num = evt.run();
   subrun_num = evt.subRun();
   event_num = evt.event();

   //std::cout << "[CC1pi::FillTree::SetReco2Vars] Using isData = " << isData << std::endl;

   art::Handle<std::vector<ubana::SelectionResult>> selection_h;
   evt.getByLabel(CC1piInputTags->fSelectionLabel, selection_h);
   if(!selection_h.isValid()){
      mf::LogError(__PRETTY_FUNCTION__) << "[FillTree] SelectionResult product not found." << std::endl;
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

      nu_MCPID = top_neutrino.TrackId();

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
      art::FindManyP<ubana::TPCObject> TPCObject_from_selection(selection_h, evt, CC1piInputTags->fSelectionLabel);
      art::Ptr<ubana::TPCObject> TPCObj_candidate = TPCObject_from_selection.at(0).at(0);
      art::Handle<std::vector<ubana::TPCObject>> TPCObj_h;
      evt.getByLabel(CC1piInputTags->fTPCObjectProducer, TPCObj_h);

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
      art::FindManyP<recob::PFParticle> pfps_from_TPCObject(TPCObj_h, evt, CC1piInputTags->fTPCObjectProducer);
      std::vector<art::Ptr<recob::PFParticle>> pfps = pfps_from_TPCObject.at(TPCObj_candidate.key());
      TPCObj_NPFPs = pfps.size();

      //Get tracks (in TPCObject)
      art::FindManyP<recob::Track> tracks_from_TPCObject(TPCObj_h, evt, CC1piInputTags->fTPCObjectProducer);
      std::vector<art::Ptr<recob::Track>> tracks = tracks_from_TPCObject.at(TPCObj_candidate.key());
      TPCObj_NTracks = tracks.size();

      // Get calo objects (for dedx)
      // Use association: pandoraNucalo, art::Assns<recob::Track,anab::Calorimetry,void>
      // (Couldn't find another way to do this except going back to track handle -.-)
      auto const& track_h = evt.getValidHandle<std::vector<recob::Track>>(CC1piInputTags->fTrackLabel);
      art::FindManyP<anab::Calorimetry> calos_from_tracks(track_h, evt, CC1piInputTags->fCalorimetryLabel);

      // Also get track-PID association objects
      art::FindManyP<anab::ParticleID> trackPIDAssn(track_h, evt, CC1piInputTags->fPIDProducer);
      // std::cout << "[CC1pi] trackPIDAssn.IsValid() = " << trackPIDAssn.isValid() << std::endl;
      // std::cout << "[CC1pi] trackPIDAssn.size() = " << trackPIDAssn.size() << std::endl;
      art::FindManyP<anab::ParticleID> trackPIDAssnforChi2(track_h, evt, CC1piInputTags->fPIDLabelChi2);

      //Get showers (in TPCObject)
      art::FindManyP<recob::Shower> showers_from_TPCObject(TPCObj_h, evt, CC1piInputTags->fTPCObjectProducer);
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
      evt.getByLabel(CC1piInputTags->fTrackLabel,pfp_h);
      if(!pfp_h.isValid()){
         mf::LogError(__PRETTY_FUNCTION__) << "PFP product not found." << std::endl;
         throw std::exception();
      }

      //Get associations
      art::FindManyP<ubana::MCGhost>   mcghost_from_pfp (pfp_h,   evt, "RecoTrueMatcher");
      art::FindManyP<simb::MCParticle> mcp_from_mcghost (ghost_h, evt, "RecoTrueMatcher");

      // Get MCS fit results for all tracks in event (we'll only use the result for the tracks that are in the selected TPC object but there's no way to do this without getting them all)
      art::Handle<std::vector<recob::MCSFitResult> > mcsfitresultHandle_mu;
      evt.getByLabel(CC1piInputTags->fMCSMuProducer,mcsfitresultHandle_mu);
      if(!mcsfitresultHandle_mu.isValid()){
        std::cout << "[UBXSec] MCSFitResult product " << CC1piInputTags->fMCSMuProducer.label() << " not found..." << std::endl;
        //throw std::exception();
      }
      std::vector<art::Ptr<recob::MCSFitResult>> mcsfitresultVector_mu;
      art::fill_ptr_vector(mcsfitresultVector_mu, mcsfitresultHandle_mu);

      art::Handle<std::vector<recob::MCSFitResult> > mcsfitresultHandle_p;
      evt.getByLabel(CC1piInputTags->fMCSPProducer,mcsfitresultHandle_p);
      if(!mcsfitresultHandle_p.isValid()){
        std::cout << "[UBXSec] MCSFitResult product " << CC1piInputTags->fMCSPProducer.label() << " not found..." << std::endl;
        //throw std::exception();
      }
      std::vector<art::Ptr<recob::MCSFitResult>> mcsfitresultVector_p;
      art::fill_ptr_vector(mcsfitresultVector_p, mcsfitresultHandle_p);

      art::Handle<std::vector<recob::MCSFitResult> > mcsfitresultHandle_pi;
      evt.getByLabel(CC1piInputTags->fMCSPiProducer,mcsfitresultHandle_pi);
      if(!mcsfitresultHandle_pi.isValid()){
        std::cout << "[UBXSec] MCSFitResult product " << CC1piInputTags->fMCSPiProducer.label() << " not found..." << std::endl;
        //throw std::exception();
      }
      std::vector<art::Ptr<recob::MCSFitResult>> mcsfitresultVector_pi;
      art::fill_ptr_vector(mcsfitresultVector_pi, mcsfitresultHandle_pi);

      // there is no association between tracks and the MCSResult stored in the file,
      // so going to build a map between each track ID and an art::Ptr to the
      // recob::MCSFitResult for quick finding
      MapBuilderUtil mbutil;
      std::vector<art::Ptr<recob::Track>> trackVector;
      art::fill_ptr_vector(trackVector, track_h);
      std::map<int, art::Ptr<recob::MCSFitResult> > trackIdMcsFitMap_mu = mbutil.buildTrackIdMcsResultMap(trackVector, mcsfitresultVector_mu);
      std::map<int, art::Ptr<recob::MCSFitResult> > trackIdMcsFitMap_p = mbutil.buildTrackIdMcsResultMap(trackVector, mcsfitresultVector_p);
      std::map<int, art::Ptr<recob::MCSFitResult> > trackIdMcsFitMap_pi = mbutil.buildTrackIdMcsResultMap(trackVector, mcsfitresultVector_pi);


      //Find neutrino
      int nuID = -1;
      for (auto pfp : pfps) {
         if(lar_pandora::LArPandoraHelper::IsNeutrino(pfp)) {
            nuID = (int)pfp->Self();
            break;
         }
      }

      art::FindManyP<recob::Track> tracks_from_pfps(pfp_h, evt, CC1piInputTags->fTrackLabel);
      art::FindManyP<recob::Shower> showers_from_pfps(pfp_h, evt, CC1piInputTags->fTrackLabel);
      art::FindManyP<recob::SpacePoint> spacepoints_from_pfps(pfp_h, evt, CC1piInputTags->fTrackLabel);

      art::Handle<std::vector<recob::SpacePoint>> sp_h;
      evt.getByLabel(CC1piInputTags->fTrackLabel, sp_h);
      art::FindManyP<recob::Hit> hits_from_spacepoints(sp_h, evt, CC1piInputTags->fTrackLabel);


      std::vector<TVector3> vtxDirs = {};

      for (auto pfp : pfps) {

         //Set default values for track/shower specific variables
         double track_length = -9999;
         std::vector<double> track_start = {-9999, -9999, -9999};
         std::vector<double> track_end = {-9999, -9999, -9999};
         double track_theta = -9999.;
         double track_phi = -9999.;
         double track_mom = -9999.;
         std::vector<double> track_dedx_truncmean(3,-9999.);
         std::vector<std::vector<double>> track_dedx_perhit(3);
         std::vector<std::vector<double>> track_resrange_perhit(3);
         double nhits = -9999;
         bool isMIP = false;
         std::vector<std::vector<double>> track_trajpoint_position;
         std::vector<std::vector<double>> track_trajpoint_direction;
         std::vector<std::vector<double>> noBragg_MIP_perhit(3);
         std::vector<double> Bragg_fwd_mu(3,-999.);
         std::vector<double> Bragg_fwd_p(3,-999.);
         std::vector<double> Bragg_fwd_pi(3,-999.);
         std::vector<double> Bragg_bwd_mu(3,-999.);
         std::vector<double> Bragg_bwd_p(3,-999.);
         std::vector<double> Bragg_bwd_pi(3,-999.);
         std::vector<double> noBragg_MIP(3,-999.);
         std::vector<double> PIDAval(3,-999.);
         std::vector<double> track_depE(3,-999.);
         std::vector<double> track_Chi2Proton(3,-999.);
         std::vector<double> track_Chi2Muon(3,-999.);
         std::vector<double> track_Chi2Pion(3,-999.);
         double track_rangeE_mu = -9999;
         double track_rangeE_p = -9999;
         bool isContained = false;
         std::pair<double,double> residual_mean_std(-9999,-9999);
         double perc_used_hits = -9999;
         double MCSmu_fwdMom = -9999;
         double MCSmu_bwdMom = -9999;
         double MCSmu_bestMom = -9999;
         double MCSmu_fwdMomUncert = -9999;
         double MCSmu_bwdMomUncert = -9999;
         double MCSmu_bestMomUncert = -9999;
         double MCSmu_fwdLL = -9999;
         double MCSmu_bwdLL = -9999;
         double MCSmu_bestLL = -9999;
         std::vector<double> MCSmu_segmentRadLengths(1,-9999.);
         std::vector<double> MCSmu_scatterAngles(1,-9999.);
         double MCSp_fwdMom = -9999;
         double MCSp_bwdMom = -9999;
         double MCSp_bestMom = -9999;
         double MCSp_fwdMomUncert = -9999;
         double MCSp_bwdMomUncert = -9999;
         double MCSp_bestMomUncert = -9999;
         double MCSp_fwdLL = -9999;
         double MCSp_bwdLL = -9999;
         double MCSp_bestLL = -9999;
         std::vector<double> MCSp_segmentRadLengths(1,-9999.);
         std::vector<double> MCSp_scatterAngles(1,-9999.);
         double MCSpi_fwdMom = -9999;
         double MCSpi_bwdMom = -9999;
         double MCSpi_bestMom = -9999;
         double MCSpi_fwdMomUncert = -9999;
         double MCSpi_bwdMomUncert = -9999;
         double MCSpi_bestMomUncert = -9999;
         double MCSpi_fwdLL = -9999;
         double MCSpi_bwdLL = -9999;
         double MCSpi_bestLL = -9999;
         std::vector<double> MCSpi_segmentRadLengths(1,-9999.);
         std::vector<double> MCSpi_scatterAngles(1,-9999.);
         std::vector<std::vector<double>> spacepoints_XYZ;
         std::vector<double> spacepoints_qplane2;
         std::vector<double> unusedhits_charge_plane0;
         std::vector<int> unusedhits_endwiredist_plane0;
         std::vector<double> unusedhits_endtimedist_plane0;
         std::vector<double> unusedhits_charge_plane1;
         std::vector<int> unusedhits_endwiredist_plane1;
         std::vector<double> unusedhits_endtimedist_plane1;
         std::vector<double> unusedhits_charge_plane2;
         std::vector<int> unusedhits_endwiredist_plane2;
         std::vector<double> unusedhits_endtimedist_plane2;


         double shower_length = -9999;
         std::vector<double> shower_start = {-9999, -9999, -9999};

         std::vector<art::Ptr<recob::Track>> tracks_pfp = tracks_from_pfps.at(pfp.key());
         if(tracks_pfp.size() > 1) {
            mf::LogError(__PRETTY_FUNCTION__) << "PFP associated to more than one track." << std::endl;
            throw std::exception();
         }
         for (auto track : tracks_pfp){
            track_length = track -> Length();
            track_theta = track -> Theta();
            track_phi = track -> Phi();
            //if (track->HasMomentum()) track_mom = track -> VertexMomentum(); // Commented out because I don't think this is right, we need to think more about if we can use momentum

            // space charge correction
            auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();

            auto start = track -> Start();
            std::vector<double> sce_corr_start = SCE->GetPosOffsets(start.X(), start.Y(), start.Z());
            track_start = {start.X() + sce_corr_start.at(0), start.Y() - sce_corr_start.at(1), start.Z() - sce_corr_start.at(2)};

            auto end = track -> End();
            std::vector<double> sce_corr_end = SCE->GetPosOffsets(end.X(), end.Y(), end.Z());
            track_end = {end.X() + sce_corr_end.at(0), end.Y() - sce_corr_end.at(1), end.Z() - sce_corr_end.at(2)};

            vtxDirs.emplace_back(track -> VertexDirection());

            //unsigned int trkid = track->ID();

            // Test containment
            if(inFV(track_start.at(0),track_start.at(1),track_start.at(2)) && inFV(track_end.at(0),track_end.at(1),track_end.at(2))) isContained = true;
            else isContained = false;

            //Number of 3D hits
            nhits = track -> NPoints();

            // Get calorimetry information
            // Only fill when there is valid calorimetry information (of course)
            if (calos_from_tracks.isValid()){
               std::vector< art::Ptr<anab::Calorimetry> > caloFromTrack = calos_from_tracks.at(track->ID());

               // for time being, only use Y plane calorimetry
               art::Ptr< anab:: Calorimetry > calo_plane0;
               art::Ptr< anab:: Calorimetry > calo_plane1;
               art::Ptr< anab:: Calorimetry > calo_plane2;
               for (auto c : caloFromTrack){
                  int planenum = c->PlaneID().Plane;

                  if (planenum < 0 || planenum > 2){
                     std::cout << "[CC1pi::FillTree] No calorimetry information for plane number " << planenum << std::endl;
                     continue;
                  }

                  track_dedx_perhit.at(planenum) = c->dEdx();
                  track_resrange_perhit.at(planenum) = c->ResidualRange();

                  std::vector<double> trkpitchvec = c->TrkPitchVec();

                  double track_depE_tmp = 0;
                  std::vector<double> dEdx_dummy = {0.};
                  std::vector<double> rr_dummy = {0.};

                  for (size_t i_hit=0; i_hit < c->dEdx().size(); i_hit++){
                     track_depE_tmp += c->dEdx().at(i_hit)*trkpitchvec.at(i_hit);

                     double Lmip = -9999;
                     dEdx_dummy.at(0) = c->dEdx().at(i_hit);
                     rr_dummy.at(0) = c->ResidualRange().at(i_hit);
                     Lmip = braggcalc.getLikelihood(dEdx_dummy,rr_dummy,0,1,planenum);
                     noBragg_MIP_perhit.at(planenum).emplace_back(Lmip);
                  }
                  track_depE.at(planenum) = track_depE_tmp;
               }
            }  // end if(caloFromTracks.isValid())



            // Store isMIP
            //isMIP = IsMIP(trackPIDAssn, trkid);
            isMIP = IsMIP(track, evt);

            // Now get PID information from track association
            // Note that this relies on you having the feature branch of lardataobj lardataobj/feature/kduffy_pidrefactor_v1_11_00_04 checked out, otherwise you won't be able to read these products
            if (trackPIDAssn.isValid()){

               std::vector<art::Ptr<anab::ParticleID>> trackPID = trackPIDAssn.at(track->ID());
               if (trackPID.size() != 0){
                  // Only fill variables if the track-PID association exists
                  std::vector<anab::sParticleIDAlgScores> AlgScoresVec = trackPID.at(0)->ParticleIDAlgScores();
                  // Loop through AlgScoresVec and find the variables we want
                  for (size_t i_algscore=0; i_algscore<AlgScoresVec.size(); i_algscore++){

                     anab::sParticleIDAlgScores AlgScore = AlgScoresVec.at(i_algscore);
                     int planeid = UBPID::uB_getSinglePlane(AlgScore.fPlaneID);

                     if (planeid < 0 || planeid > 2){
                       std::cout << "[CC1pi::FillTree] No ParticleID information for planeid " << planeid << std::endl;
                       continue;
                     }


                     if (AlgScore.fAlgName == "BraggPeakLLH" && anab::kVariableType(AlgScore.fVariableType) == anab::kLikelihood){
                       if (AlgScore.fAssumedPdg == 0) noBragg_MIP.at(planeid) = AlgScore.fValue;
                        if (anab::kTrackDir(AlgScore.fTrackDir) == anab::kForward){
                           if (AlgScore.fAssumedPdg == 13)   Bragg_fwd_mu.at(planeid) = AlgScore.fValue;
                           if (AlgScore.fAssumedPdg == 2212) Bragg_fwd_p.at(planeid) =  AlgScore.fValue;
                           if (AlgScore.fAssumedPdg == 211) Bragg_fwd_pi.at(planeid) =  AlgScore.fValue;
                        }// if fTrackDir == anab::kForward
                        else if (anab::kTrackDir(AlgScore.fTrackDir) == anab::kBackward){
                           if (AlgScore.fAssumedPdg == 13)   Bragg_bwd_mu.at(planeid) = AlgScore.fValue;
                           if (AlgScore.fAssumedPdg == 2212) Bragg_bwd_p.at(planeid) =  AlgScore.fValue;
                           if (AlgScore.fAssumedPdg == 211) Bragg_bwd_pi.at(planeid) =  AlgScore.fValue;
                        } // if fTrackDir == anab::kBackward
                     } // if fAlName = BraggPeakLLH && fVariableType == anab::kLikelihood

                     if (AlgScore.fAlgName == "PIDA_median" && anab::kVariableType(AlgScore.fVariableType) == anab::kPIDA){
                        PIDAval.at(planeid) = AlgScore.fValue;
                     }// if AlgName = PIDA && fVariableType == anab::kPIDA


                     if (AlgScore.fAlgName == "TruncatedMean"){
                        if (anab::kVariableType(AlgScore.fVariableType) == anab::kdEdxtruncmean) track_dedx_truncmean.at(planeid) = AlgScore.fValue;
                     }// if AlgName = TruncatedMean

                  } // end loop through AlgScoresVec
               } // end if trackPID.size() != 0

               // else{
               // std::cout << "[CC1pi] Found valid trackPIDAssn but trackPID.size() == " << trackPID.size() << std::endl;
               // }

            } // if (trackPIDAssn.isValid())

            // else{
            //    std::cout << "[CC1pi] Could not find valid trackPIDAssn" << std::endl;
            // }

            // Get Chi2 objects from a different ParticleID data product
            if (!trackPIDAssnforChi2.isValid()){
               std::cout << "[ParticleIDValidation] trackPIDAssnforChi2.isValid() == false. Not filling Chi2 variables." << std::endl;
               track_Chi2Proton = {-9999,-9999,-9999};
               track_Chi2Pion   = {-9999,-9999,-9999};
               track_Chi2Muon   = {-9999,-9999,-9999};
            }
            else{
               std::vector<art::Ptr<anab::ParticleID>> trackPIDforChi2 = trackPIDAssnforChi2.at(track->ID());

               for (size_t i_chi2pid=0; i_chi2pid<trackPIDforChi2.size(); i_chi2pid++){
                  //std::cout << "trackPIDforChi2.at(" << i_chi2pid << ")->PlaneID().Plane = " << trackPIDforChi2.at(i_chi2pid)->PlaneID().Plane << std::endl;
                  int planeid = trackPIDforChi2.at(i_chi2pid)->PlaneID().Plane;
                  if (planeid < 0 || planeid > 2){
                     std::cout << "[CC1pi::FillTree] No Chi2 ParticleID information for planeid " << planeid << std::endl;
                     continue;
                  }

                  track_Chi2Proton.at(planeid) = trackPIDforChi2.at(i_chi2pid)->Chi2Proton();
                  track_Chi2Pion.at(planeid)   = trackPIDforChi2.at(i_chi2pid)->Chi2Pion();
                  track_Chi2Muon.at(planeid)   = trackPIDforChi2.at(i_chi2pid)->Chi2Muon();
               } // end loop over i_plane
            } // end else

            // Get energy estimation by range (code for momentum by range copied from analysistree, then convert momentum to energy)
            // Calculations only exist in TrackMomentumCalculator for muons and protons
            // TrackMomentumCalculator returns GeV, multiply by 1000 to get MeV
            trkf::TrackMomentumCalculator trkm;
            double track_rangeP_mu = trkm.GetTrackMomentum(track->Length(),kMuMinus)*1000.;
            double track_rangeP_p = trkm.GetTrackMomentum(track->Length(),kProton)*1000.;

            // Now convert P->E
            // From TrackMomentumCalculator::GetTrackMomentum: P = TMath::Sqrt((KE*KE)+(2*M*KE))
            // P = TMath::Sqrt((E*E)-(M*M)) and E = KE+M
            // => KE = TMath::Sqrt((P*P)+(M*M))-M
            // TrackMometumCalculator uses Muon_M = 105.7 MeV, Proton_M = 938.272 MeV so use these values here
            track_rangeE_mu = TMath::Sqrt((track_rangeP_mu*track_rangeP_mu)+(105.7*105.7)) - 105.7;
            track_rangeE_p = TMath::Sqrt((track_rangeP_p*track_rangeP_p)+(938.272*938.272)) - 938.272;


            // Shower Rejection
            ShowerCheck(evt, CC1piInputTags, pfp, track, residual_mean_std, perc_used_hits);


            // Now save MCS fit results
            art::Ptr<recob::MCSFitResult> mcsFitResult;

            for (auto const& x : trackIdMcsFitMap_mu){
              if (track->ID() == x.first){
                mcsFitResult = x.second;
                break;
              }
            }
            MCSmu_fwdMom = mcsFitResult->fwdMomentum();
            MCSmu_bwdMom = mcsFitResult->bwdMomentum();
            MCSmu_bestMom = mcsFitResult->bestMomentum();
            MCSmu_fwdMomUncert = mcsFitResult->fwdMomUncertainty();
            MCSmu_bwdMomUncert = mcsFitResult->bwdMomUncertainty();
            MCSmu_bestMomUncert = mcsFitResult->bestMomUncertainty();
            MCSmu_fwdLL = mcsFitResult->fwdLogLikelihood();
            MCSmu_bwdLL = mcsFitResult->fwdLogLikelihood();
            MCSmu_bestLL = mcsFitResult->fwdLogLikelihood();
            MCSmu_segmentRadLengths = mcsFitResult->segmentRadLengths();
            MCSmu_scatterAngles = mcsFitResult->scatterAngles();

            for (auto const& x : trackIdMcsFitMap_p){
              if (track->ID() == x.first){
                mcsFitResult = x.second;
                break;
              }
            }
            MCSp_fwdMom = mcsFitResult->fwdMomentum();
            MCSp_bwdMom = mcsFitResult->bwdMomentum();
            MCSp_bestMom = mcsFitResult->bestMomentum();
            MCSp_fwdMomUncert = mcsFitResult->fwdMomUncertainty();
            MCSp_bwdMomUncert = mcsFitResult->bwdMomUncertainty();
            MCSp_bestMomUncert = mcsFitResult->bestMomUncertainty();
            MCSp_fwdLL = mcsFitResult->fwdLogLikelihood();
            MCSp_bwdLL = mcsFitResult->fwdLogLikelihood();
            MCSp_bestLL = mcsFitResult->fwdLogLikelihood();
            MCSp_segmentRadLengths = mcsFitResult->segmentRadLengths();
            MCSp_scatterAngles = mcsFitResult->scatterAngles();

            for (auto const& x : trackIdMcsFitMap_pi){
              if (track->ID() == x.first){
                mcsFitResult = x.second;
                break;
              }
            }
            MCSpi_fwdMom = mcsFitResult->fwdMomentum();
            MCSpi_bwdMom = mcsFitResult->bwdMomentum();
            MCSpi_bestMom = mcsFitResult->bestMomentum();
            MCSpi_fwdMomUncert = mcsFitResult->fwdMomUncertainty();
            MCSpi_bwdMomUncert = mcsFitResult->bwdMomUncertainty();
            MCSpi_bestMomUncert = mcsFitResult->bestMomUncertainty();
            MCSpi_fwdLL = mcsFitResult->fwdLogLikelihood();
            MCSpi_bwdLL = mcsFitResult->fwdLogLikelihood();
            MCSpi_bestLL = mcsFitResult->fwdLogLikelihood();
            MCSpi_segmentRadLengths = mcsFitResult->segmentRadLengths();
            MCSpi_scatterAngles = mcsFitResult->scatterAngles();

            // Save information about unused hits
            GetUnusedHitCharges(evt, CC1piInputTags, 0, track->End(), unusedhits_charge_plane0, unusedhits_endwiredist_plane0, unusedhits_endtimedist_plane0);
            GetUnusedHitCharges(evt, CC1piInputTags, 1, track->End(), unusedhits_charge_plane1, unusedhits_endwiredist_plane1, unusedhits_endtimedist_plane1);
            GetUnusedHitCharges(evt, CC1piInputTags, 2, track->End(), unusedhits_charge_plane2, unusedhits_endwiredist_plane2, unusedhits_endtimedist_plane2);

         } // end loop over tracks

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

         // Fill track/shower specific variables
         TPCObj_PFP_track_length.emplace_back(track_length);
         TPCObj_PFP_track_start.emplace_back(track_start);
         TPCObj_PFP_track_end.emplace_back(track_end);
         TPCObj_PFP_track_theta.emplace_back(track_theta);
         TPCObj_PFP_track_phi.emplace_back(track_phi);
         TPCObj_PFP_track_mom.emplace_back(track_mom);
         TPCObj_PFP_track_dedx_truncmean.emplace_back(track_dedx_truncmean);
         TPCObj_PFP_track_dedx_perhit.emplace_back(track_dedx_perhit);
         TPCObj_PFP_track_nhits.emplace_back(nhits);
         TPCObj_PFP_track_resrange_perhit.emplace_back(track_resrange_perhit);
         TPCObj_PFP_isMIP.emplace_back(isMIP);
         TPCObj_PFP_track_isContained.emplace_back(isContained);
         TPCObj_PFP_track_residual_mean.emplace_back(residual_mean_std.first);
         TPCObj_PFP_track_residual_std.emplace_back(residual_mean_std.second);
         TPCObj_PFP_track_perc_used_hits.emplace_back(perc_used_hits);

         TPCObj_PFP_shower_length.emplace_back(shower_length);
         TPCObj_PFP_shower_start.emplace_back(shower_start);

         TPCObj_PFP_LH_fwd_mu.emplace_back(Bragg_fwd_mu);
         TPCObj_PFP_LH_fwd_p.emplace_back(Bragg_fwd_p);
         TPCObj_PFP_LH_fwd_pi.emplace_back(Bragg_fwd_pi);
         TPCObj_PFP_LH_bwd_mu.emplace_back(Bragg_bwd_mu);
         TPCObj_PFP_LH_bwd_p.emplace_back(Bragg_bwd_p);
         TPCObj_PFP_LH_bwd_pi.emplace_back(Bragg_bwd_pi);
         TPCObj_PFP_LH_MIP.emplace_back(noBragg_MIP);
         TPCObj_PFP_PIDA.emplace_back(PIDAval);
         TPCObj_PFP_track_depE.emplace_back(track_depE);
         TPCObj_PFP_track_Chi2Proton.emplace_back(track_Chi2Proton);
         TPCObj_PFP_track_Chi2Muon.emplace_back(track_Chi2Muon);
         TPCObj_PFP_track_Chi2Pion.emplace_back(track_Chi2Pion);
         TPCObj_PFP_track_rangeE_mu.emplace_back(track_rangeE_mu);
         TPCObj_PFP_track_rangeE_p.emplace_back(track_rangeE_p);
         TPCObj_PFP_track_Lmip_perhit.emplace_back(noBragg_MIP_perhit);


        TPCObj_PFP_track_MCSmu_fwdMom.emplace_back(MCSmu_fwdMom);
        TPCObj_PFP_track_MCSmu_bwdMom.emplace_back(MCSmu_bwdMom);
        TPCObj_PFP_track_MCSmu_bestMom.emplace_back(MCSmu_bestMom);
        TPCObj_PFP_track_MCSmu_fwdMomUncert.emplace_back(MCSmu_fwdMomUncert);
        TPCObj_PFP_track_MCSmu_bwdMomUncert.emplace_back(MCSmu_bwdMomUncert);
        TPCObj_PFP_track_MCSmu_bestMomUncert.emplace_back(MCSmu_bestMomUncert);
        TPCObj_PFP_track_MCSmu_fwdLL.emplace_back(MCSmu_fwdLL);
        TPCObj_PFP_track_MCSmu_bwdLL.emplace_back(MCSmu_bwdLL);
        TPCObj_PFP_track_MCSmu_bestLL.emplace_back(MCSmu_bestLL);
        TPCObj_PFP_track_MCSmu_segmentRadLengths.emplace_back(MCSmu_segmentRadLengths);
        TPCObj_PFP_track_MCSmu_scatterAngles.emplace_back(MCSmu_scatterAngles);
        TPCObj_PFP_track_MCSp_fwdMom.emplace_back(MCSp_fwdMom);
        TPCObj_PFP_track_MCSp_bwdMom.emplace_back(MCSp_bwdMom);
        TPCObj_PFP_track_MCSp_bestMom.emplace_back(MCSp_bestMom);
        TPCObj_PFP_track_MCSp_fwdMomUncert.emplace_back(MCSp_fwdMomUncert);
        TPCObj_PFP_track_MCSp_bwdMomUncert.emplace_back(MCSp_bwdMomUncert);
        TPCObj_PFP_track_MCSp_bestMomUncert.emplace_back(MCSp_bestMomUncert);
        TPCObj_PFP_track_MCSp_fwdLL.emplace_back(MCSp_fwdLL);
        TPCObj_PFP_track_MCSp_bwdLL.emplace_back(MCSp_bwdLL);
        TPCObj_PFP_track_MCSp_bestLL.emplace_back(MCSp_bestLL);
        TPCObj_PFP_track_MCSp_segmentRadLengths.emplace_back(MCSp_segmentRadLengths);
        TPCObj_PFP_track_MCSp_scatterAngles.emplace_back(MCSp_scatterAngles);
        TPCObj_PFP_track_MCSpi_fwdMom.emplace_back(MCSpi_fwdMom);
        TPCObj_PFP_track_MCSpi_bwdMom.emplace_back(MCSpi_bwdMom);
        TPCObj_PFP_track_MCSpi_bestMom.emplace_back(MCSpi_bestMom);
        TPCObj_PFP_track_MCSpi_fwdMomUncert.emplace_back(MCSpi_fwdMomUncert);
        TPCObj_PFP_track_MCSpi_bwdMomUncert.emplace_back(MCSpi_bwdMomUncert);
        TPCObj_PFP_track_MCSpi_bestMomUncert.emplace_back(MCSpi_bestMomUncert);
        TPCObj_PFP_track_MCSpi_fwdLL.emplace_back(MCSpi_fwdLL);
        TPCObj_PFP_track_MCSpi_bwdLL.emplace_back(MCSpi_bwdLL);
        TPCObj_PFP_track_MCSpi_bestLL.emplace_back(MCSpi_bestLL);
        TPCObj_PFP_track_MCSpi_segmentRadLengths.emplace_back(MCSpi_segmentRadLengths);
        TPCObj_PFP_track_MCSpi_scatterAngles.emplace_back(MCSpi_scatterAngles);

        TPCObj_PFP_track_unusedhits_charge_plane0.emplace_back(unusedhits_charge_plane0);
        TPCObj_PFP_track_unusedhits_endwiredist_plane0.emplace_back(unusedhits_endwiredist_plane0);
        TPCObj_PFP_track_unusedhits_endtimedist_plane0.emplace_back(unusedhits_endtimedist_plane0);
        TPCObj_PFP_track_unusedhits_charge_plane1.emplace_back(unusedhits_charge_plane1);
        TPCObj_PFP_track_unusedhits_endwiredist_plane1.emplace_back(unusedhits_endwiredist_plane1);
        TPCObj_PFP_track_unusedhits_endtimedist_plane1.emplace_back(unusedhits_endtimedist_plane1);
        TPCObj_PFP_track_unusedhits_charge_plane2.emplace_back(unusedhits_charge_plane2);
        TPCObj_PFP_track_unusedhits_endwiredist_plane2.emplace_back(unusedhits_endwiredist_plane2);
        TPCObj_PFP_track_unusedhits_endtimedist_plane2.emplace_back(unusedhits_endtimedist_plane2);


         // Fill non-specific variables
         TPCObj_PFP_PandoraClassedAsTrack.emplace_back(lar_pandora::LArPandoraHelper::IsTrack(pfp));
         TPCObj_PFP_PandoraClassedAsShower.emplace_back(lar_pandora::LArPandoraHelper::IsShower(pfp));
         TPCObj_PFP_isDaughter.emplace_back(pfp->Parent()==(size_t)nuID && track_length != -9999); //If Pandora failed to make a track, don't count the PFP as a daughter (even if the associated shower technically is)
         std::vector<int> daughterids(pfp->Daughters().begin(),pfp->Daughters().end());
         TPCObj_PFP_daughterids.emplace_back(daughterids);
         TPCObj_PFP_id.emplace_back(pfp -> Self());

         // Save spacepoints
         // Get spacepoints from PFP
         std::vector<art::Ptr<recob::SpacePoint>> sp_vec;
         spacepoints_from_pfps.get(pfp.key(),sp_vec);
         for (size_t i_sp=0; i_sp<sp_vec.size(); i_sp++){
           const double *XYZ = sp_vec.at(i_sp)->XYZ();
           if (XYZ==nullptr) continue;
           std::vector<double> dummytvec = {XYZ[0],XYZ[1],XYZ[2]};
           spacepoints_XYZ.push_back(dummytvec);

           // we also want to store the charge per spacepoint (so we can do non-track-dependent dE/dx calculations). Since we only trust the plane 2 calorimetry in MCC8, let's only use the plane 2 integrated hit charge. If no plane 2 hit exists for a given spacepoint, set the charge to -9999.
           double charge=-9999;
           std::vector<art::Ptr<recob::Hit>> hits_from_sp;
           hits_from_spacepoints.get(sp_vec.at(i_sp).key(),hits_from_sp);
           int nhits_plane2=0;
           for (auto hit_sp : hits_from_sp){
             if (hit_sp->View()==1){ // plane 2 hits only
               nhits_plane2++;
               charge = hit_sp->Integral();
             }
           } // end loop over hits from spacepoint
           if (nhits_plane2>1){
             std::cout << "[CC1pi::FillTree] Found " << nhits_plane2 << " plane 2 hits associated with one spacepoint. Not using charge from any of them." << std::endl;
             charge = -9999;
           }
           spacepoints_qplane2.push_back(charge);
         } // end loop over spacepoints
         TPCObj_PFP_track_SpacepointsXYZ.push_back(spacepoints_XYZ);
         TPCObj_PFP_track_SpacepointsQPlane2.push_back(spacepoints_qplane2);

         // Do reco-truth matching...

         auto mcghosts = mcghost_from_pfp.at(pfp.key());
         if (mcghosts.size() != 1) {
            TPCObj_PFP_MCPid.emplace_back(-9999);
            TPCObj_PFP_truePDG.emplace_back(-9999);
            TPCObj_PFP_trueE.emplace_back(-9999);
            TPCObj_PFP_trueKE.emplace_back(-9999);
            TPCObj_PFP_trueEndP.emplace_back(-9999);
            continue;
         }
         auto mcghost = mcghosts[0];

         auto mcps = mcp_from_mcghost.at(mcghost.key());
         if(mcps.size() != 1) {
            TPCObj_PFP_MCPid.emplace_back(-9999);
            TPCObj_PFP_truePDG.emplace_back(-9999);
            TPCObj_PFP_trueE.emplace_back(-9999);
            TPCObj_PFP_trueKE.emplace_back(-9999);
            TPCObj_PFP_trueEndP.emplace_back(-9999);
            continue;
         }
         auto mcp = mcps[0];

         TPCObj_PFP_MCPid.emplace_back(mcp -> TrackId());
         TPCObj_PFP_truePDG.emplace_back(mcp -> PdgCode());
         TPCObj_PFP_trueE.emplace_back(mcp -> E());
         TPCObj_PFP_trueKE.emplace_back(mcp -> E() - mcp -> Mass());
         // TVector3 trueEndP(mcp -> EndPx(), mcp -> EndPy(), mcp -> EndPz());
         // TPCObj_PFP_trueEndP.emplace_back(trueEndP.Mag());
         TPCObj_PFP_trueEndP.emplace_back(mcp->EndE() - mcp->Mass());
         // std::cout << "trueEndP = " << mcp->EndPx() << ", " << mcp->EndPy() << ", " << mcp->EndPz() << std::endl;
         // std::cout << "trueEndP.Mag() = " << trueEndP.Mag() << std::endl;
         // std::cout << "EndE = " << mcp -> EndE() << std::endl;
         // std::cout << "EndE - Mass = " << mcp->EndE() - mcp -> Mass() << std::endl;

      } // loop over pfps

      // Calculate angle between tracks
      for (size_t track1 = 0; track1 < vtxDirs.size(); track1++) {
         std::vector<double> angles = {};
         for (size_t track2 = 0; track2 < vtxDirs.size(); track2++) {
            angles.emplace_back(vtxDirs.at(track1).Angle(vtxDirs.at(track2)));
         }
         TPCObj_PFP_track_AngleBetweenTracks.emplace_back(angles);
      }

      // Get neutrino candidate vertex from TPCObject
      recob::Vertex TPCObj_nu_vtx = TPCObj_candidate->GetVertex();
      double reco_nu_vtx[3];
      TPCObj_nu_vtx.XYZ(reco_nu_vtx);
      TPCObj_reco_vtx_nospacecharge = {reco_nu_vtx[0], reco_nu_vtx[1], reco_nu_vtx[2]};

      // Space Charge correction
      auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
      std::vector<double> sce_corr = SCE->GetPosOffsets(reco_nu_vtx[0], reco_nu_vtx[1], reco_nu_vtx[2]);
      TPCObj_reco_vtx = {reco_nu_vtx[0] + sce_corr.at(0), reco_nu_vtx[1] - sce_corr.at(1), reco_nu_vtx[2] - sce_corr.at(2)};

      // We eventually need to correct for X position (time offset)
      // See issue #8 on github

   } // if selected

   // Get all MCParticles (only for MC)
   if (!isData){
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
//         if (mcpar -> Mother() != 0) continue;
//         if (mcpar -> StatusCode() != 1) continue;

         MCP_PDG.emplace_back(mcpar -> PdgCode());

         MCP_P.emplace_back(mcpar -> P());
         MCP_Px.emplace_back(mcpar -> Px());
         MCP_Py.emplace_back(mcpar -> Py());
         MCP_Pz.emplace_back(mcpar -> Pz());
         MCP_E.emplace_back(mcpar -> E());
         MCP_KE.emplace_back(mcpar -> E() - mcpar -> Mass());

         simb::MCTrajectory traj = mcpar -> Trajectory();
         MCP_length.emplace_back(traj.TotalLength());
         // Note: Currently just checks start and end point are in FV. Should technically check all trajectory points.
         if(inFV(traj.X(0),traj.Y(0),traj.Z(0)) && inFV(traj.X(traj.size()-1),traj.Y(traj.size()-1),traj.Z(traj.size()-1))) MCP_isContained.emplace_back(true);
         else MCP_isContained.emplace_back(false);

         MCP_process.emplace_back(mcpar -> Process());
         MCP_endprocess.emplace_back(mcpar -> EndProcess());
         MCP_numdaughters.emplace_back(mcpar -> NumberDaughters());

         MCP_ID.emplace_back(mcpar -> TrackId());
         MCP_StatusCode.emplace_back(mcpar -> StatusCode());
         MCP_MotherID.emplace_back(mcpar -> Mother());

         if (mcpar->Mother() == nu_MCPID)
          MCP_isNuDaughter.emplace_back(true);
         else
          MCP_isNuDaughter.emplace_back(false);

         auto StartPosition = mcpar -> Position();
         auto EndPosition = mcpar -> EndPosition();
         std::vector<double> StartPosition_vect = {StartPosition.X(),StartPosition.Y(),StartPosition.Z()};
         std::vector<double> EndPosition_vect = {EndPosition.X(),EndPosition.Y(),EndPosition.Z()};
         MCP_StartPosition.emplace_back(StartPosition_vect);
         MCP_EndPosition.emplace_back(EndPosition_vect);

         std::vector<int> DaughterIDs;
         for (int daughter_i = 0; daughter_i < mcpar -> NumberDaughters(); daughter_i++) {
            DaughterIDs.emplace_back(mcpar->Daughter(daughter_i));
         }
         MCP_DaughterIDs.emplace_back(DaughterIDs);

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

            nuint_mode = nu.Mode();
            nuint_W = nu.W();
            nuint_Qsq = nu.QSqr();

            const TLorentzVector& vertex = neutrino.Position(0);
            double MC_vertex[4];
            vertex.GetXYZT(MC_vertex);
            nu_vtx = {MC_vertex[0], MC_vertex[1], MC_vertex[2]};

            // Space Charge correction
            auto const* SCE = lar::providerFrom<spacecharge::SpaceChargeService>();
            std::vector<double> sce_corr = SCE->GetPosOffsets(nu_vtx[0], nu_vtx[1], nu_vtx[2]);
            //std::cout << "SCE correction in x, y, z = " << sce_corr.at(0)
            //   << ", " << sce_corr.at(1)
            //   << ", " << sce_corr.at(2) << std::endl;
            nu_vtx_spacecharge = {nu_vtx[0] - sce_corr.at(0), nu_vtx[1] + sce_corr.at(1), nu_vtx[2] + sce_corr.at(2)};
         }
      } // end loop over MCParticles
   } // end if (!isData)
}


void MakeAnaBranches(TTree *t, cc1pianavars *vars){

   // Note: it's very important that we use the syntax &(vars->evtnum) to make sure
   // we get the value from reference, otherwise it won't work!

   t -> Branch("isData", &(vars->isData));
   t -> Branch("run_num", &(vars->run_num));
   t -> Branch("subrun_num", &(vars->subrun_num));
   t -> Branch("event_num", &(vars->event_num));

   t -> Branch("Truth_topology", &(vars->Truth_topology), "Truth_topology/I");

   t -> Branch("Marco_cutflow", &(vars->Marco_cutflow));
   t -> Branch("Marco_selected", &(vars->Marco_selected));
   t -> Branch("CC1picutflow", &(vars->CC1picutflow));

   t -> Branch("TPCObj_PFP_track_length", &(vars->TPCObj_PFP_track_length));
   t -> Branch("TPCObj_PFP_track_start", &(vars->TPCObj_PFP_track_start));
   t -> Branch("TPCObj_PFP_track_end", &(vars->TPCObj_PFP_track_end));
   t -> Branch("TPCObj_PFP_track_theta", &(vars->TPCObj_PFP_track_theta));
   t -> Branch("TPCObj_PFP_track_phi", &(vars->TPCObj_PFP_track_phi));
   t -> Branch("TPCObj_PFP_track_mom", &(vars->TPCObj_PFP_track_mom));
   t -> Branch("TPCObj_PFP_track_dedx_truncmean", &(vars->TPCObj_PFP_track_dedx_truncmean));
   t -> Branch("TPCObj_PFP_track_dedx_perhit",&(vars->TPCObj_PFP_track_dedx_perhit));
   t -> Branch("TPCObj_PFP_track_nhits", &(vars->TPCObj_PFP_track_nhits));
   t -> Branch("TPCObj_PFP_track_resrange_perhit",&(vars->TPCObj_PFP_track_resrange_perhit));
   t -> Branch("TPCObj_PFP_isMIP", &(vars->TPCObj_PFP_isMIP));
   t -> Branch("TPCObj_PFP_track_isContained", &(vars->TPCObj_PFP_track_isContained));
   t -> Branch("TPCObj_PFP_track_AngleBetweenTracks", &(vars->TPCObj_PFP_track_AngleBetweenTracks));
   t -> Branch("TPCObj_PFP_track_residual_mean", &(vars->TPCObj_PFP_track_residual_mean));
   t -> Branch("TPCObj_PFP_track_residual_std", &(vars->TPCObj_PFP_track_residual_std));
   t -> Branch("TPCObj_PFP_track_perc_used_hits", &(vars->TPCObj_PFP_track_perc_used_hits));

   t -> Branch("TPCObj_PFP_track_MCSmu_fwdMom", &(vars->TPCObj_PFP_track_MCSmu_fwdMom));
   t -> Branch("TPCObj_PFP_track_MCSmu_bwdMom", &(vars->TPCObj_PFP_track_MCSmu_bwdMom));
   t -> Branch("TPCObj_PFP_track_MCSmu_bestMom", &(vars->TPCObj_PFP_track_MCSmu_bestMom));
   t -> Branch("TPCObj_PFP_track_MCSmu_fwdMomUncert", &(vars->TPCObj_PFP_track_MCSmu_fwdMomUncert));
   t -> Branch("TPCObj_PFP_track_MCSmu_bwdMomUncert", &(vars->TPCObj_PFP_track_MCSmu_bwdMomUncert));
   t -> Branch("TPCObj_PFP_track_MCSmu_bestMomUncert", &(vars->TPCObj_PFP_track_MCSmu_bestMomUncert));
   t -> Branch("TPCObj_PFP_track_MCSmu_fwdLL", &(vars->TPCObj_PFP_track_MCSmu_fwdLL));
   t -> Branch("TPCObj_PFP_track_MCSmu_bwdLL", &(vars->TPCObj_PFP_track_MCSmu_bwdLL));
   t -> Branch("TPCObj_PFP_track_MCSmu_bestLL", &(vars->TPCObj_PFP_track_MCSmu_bestLL));
   t -> Branch("TPCObj_PFP_track_MCSmu_segmentRadLengths", &(vars->TPCObj_PFP_track_MCSmu_segmentRadLengths));
   t -> Branch("TPCObj_PFP_track_MCSmu_scatterAngles", &(vars->TPCObj_PFP_track_MCSmu_scatterAngles));
   t -> Branch("TPCObj_PFP_track_MCSp_fwdMom", &(vars->TPCObj_PFP_track_MCSp_fwdMom));
   t -> Branch("TPCObj_PFP_track_MCSp_bwdMom", &(vars->TPCObj_PFP_track_MCSp_bwdMom));
   t -> Branch("TPCObj_PFP_track_MCSp_bestMom", &(vars->TPCObj_PFP_track_MCSp_bestMom));
   t -> Branch("TPCObj_PFP_track_MCSp_fwdMomUncert", &(vars->TPCObj_PFP_track_MCSp_fwdMomUncert));
   t -> Branch("TPCObj_PFP_track_MCSp_bwdMomUncert", &(vars->TPCObj_PFP_track_MCSp_bwdMomUncert));
   t -> Branch("TPCObj_PFP_track_MCSp_bestMomUncert", &(vars->TPCObj_PFP_track_MCSp_bestMomUncert));
   t -> Branch("TPCObj_PFP_track_MCSp_fwdLL", &(vars->TPCObj_PFP_track_MCSp_fwdLL));
   t -> Branch("TPCObj_PFP_track_MCSp_bwdLL", &(vars->TPCObj_PFP_track_MCSp_bwdLL));
   t -> Branch("TPCObj_PFP_track_MCSp_bestLL", &(vars->TPCObj_PFP_track_MCSp_bestLL));
   t -> Branch("TPCObj_PFP_track_MCSp_segmentRadLengths", &(vars->TPCObj_PFP_track_MCSp_segmentRadLengths));
   t -> Branch("TPCObj_PFP_track_MCSp_scatterAngles", &(vars->TPCObj_PFP_track_MCSp_scatterAngles));
   t -> Branch("TPCObj_PFP_track_MCSpi_fwdMom", &(vars->TPCObj_PFP_track_MCSpi_fwdMom));
   t -> Branch("TPCObj_PFP_track_MCSpi_bwdMom", &(vars->TPCObj_PFP_track_MCSpi_bwdMom));
   t -> Branch("TPCObj_PFP_track_MCSpi_bestMom", &(vars->TPCObj_PFP_track_MCSpi_bestMom));
   t -> Branch("TPCObj_PFP_track_MCSpi_fwdMomUncert", &(vars->TPCObj_PFP_track_MCSpi_fwdMomUncert));
   t -> Branch("TPCObj_PFP_track_MCSpi_bwdMomUncert", &(vars->TPCObj_PFP_track_MCSpi_bwdMomUncert));
   t -> Branch("TPCObj_PFP_track_MCSpi_bestMomUncert", &(vars->TPCObj_PFP_track_MCSpi_bestMomUncert));
   t -> Branch("TPCObj_PFP_track_MCSpi_fwdLL", &(vars->TPCObj_PFP_track_MCSpi_fwdLL));
   t -> Branch("TPCObj_PFP_track_MCSpi_bwdLL", &(vars->TPCObj_PFP_track_MCSpi_bwdLL));
   t -> Branch("TPCObj_PFP_track_MCSpi_bestLL", &(vars->TPCObj_PFP_track_MCSpi_bestLL));
   t -> Branch("TPCObj_PFP_track_MCSpi_segmentRadLengths", &(vars->TPCObj_PFP_track_MCSpi_segmentRadLengths));
   t -> Branch("TPCObj_PFP_track_MCSpi_scatterAngles", &(vars->TPCObj_PFP_track_MCSpi_scatterAngles));
   t -> Branch("TPCObj_PFP_track_SpacepointsXYZ", &(vars->TPCObj_PFP_track_SpacepointsXYZ));
   t -> Branch("TPCObj_PFP_track_SpacepointsQPlane2", &(vars->TPCObj_PFP_track_SpacepointsQPlane2));
   t -> Branch("TPCObj_PFP_track_unusedhits_charge_plane0", &(vars->TPCObj_PFP_track_unusedhits_charge_plane0));
   t -> Branch("TPCObj_PFP_track_unusedhits_endwiredist_plane0", &(vars->TPCObj_PFP_track_unusedhits_endwiredist_plane0));
   t -> Branch("TPCObj_PFP_track_unusedhits_endtimedist_plane0", &(vars->TPCObj_PFP_track_unusedhits_endtimedist_plane0));
   t -> Branch("TPCObj_PFP_track_unusedhits_charge_plane1", &(vars->TPCObj_PFP_track_unusedhits_charge_plane1));
   t -> Branch("TPCObj_PFP_track_unusedhits_endwiredist_plane1", &(vars->TPCObj_PFP_track_unusedhits_endwiredist_plane1));
   t -> Branch("TPCObj_PFP_track_unusedhits_endtimedist_plane1", &(vars->TPCObj_PFP_track_unusedhits_endtimedist_plane1));
   t -> Branch("TPCObj_PFP_track_unusedhits_charge_plane2", &(vars->TPCObj_PFP_track_unusedhits_charge_plane2));
   t -> Branch("TPCObj_PFP_track_unusedhits_endwiredist_plane2", &(vars->TPCObj_PFP_track_unusedhits_endwiredist_plane2));
   t -> Branch("TPCObj_PFP_track_unusedhits_endtimedist_plane2", &(vars->TPCObj_PFP_track_unusedhits_endtimedist_plane2));


   t -> Branch("TPCObj_PFP_shower_length", &(vars->TPCObj_PFP_shower_length));
   t -> Branch("TPCObj_PFP_shower_start", &(vars->TPCObj_PFP_shower_start));

   t -> Branch("TPCObj_NPFPs", &(vars->TPCObj_NPFPs));
   t -> Branch("TPCObj_NTracks", &(vars->TPCObj_NTracks));
   t -> Branch("TPCObj_NShowers", &(vars->TPCObj_NShowers));
   t -> Branch("TPCObj_PFP_PandoraClassedAsTrack", &(vars->TPCObj_PFP_PandoraClassedAsTrack));
   t -> Branch("TPCObj_PFP_PandoraClassedAsShower", &(vars->TPCObj_PFP_PandoraClassedAsShower));
   t -> Branch("TPCObj_PFP_isDaughter", &(vars->TPCObj_PFP_isDaughter));
   t -> Branch("TPCObj_PFP_daughterids", &(vars->TPCObj_PFP_daughterids));
   t -> Branch("TPCObj_PFP_id", &(vars->TPCObj_PFP_id));
   t -> Branch("TPCObj_reco_vtx", &(vars->TPCObj_reco_vtx));
   t -> Branch("TPCObj_reco_vtx_nospacecharge", &(vars->TPCObj_reco_vtx_nospacecharge));

   t -> Branch("TPCObj_PFP_MCPid", &(vars->TPCObj_PFP_MCPid));
   t -> Branch("TPCObj_PFP_truePDG", &(vars->TPCObj_PFP_truePDG));
   t -> Branch("TPCObj_PFP_trueE", &(vars->TPCObj_PFP_trueE));
   t -> Branch("TPCObj_PFP_trueKE", &(vars->TPCObj_PFP_trueKE));
   t -> Branch("TPCObj_PFP_trueEndP", &(vars->TPCObj_PFP_trueEndP));
   t -> Branch("TPCObj_origin", &(vars->TPCObj_origin));
   t -> Branch("TPCObj_origin_extra", &(vars->TPCObj_origin_extra));

   t -> Branch("TPCObj_PFP_LH_fwd_mu", &(vars->TPCObj_PFP_LH_fwd_mu));
   t -> Branch("TPCObj_PFP_LH_fwd_p", &(vars->TPCObj_PFP_LH_fwd_p));
   t -> Branch("TPCObj_PFP_LH_fwd_pi", &(vars->TPCObj_PFP_LH_fwd_pi));
   t -> Branch("TPCObj_PFP_LH_bwd_mu", &(vars->TPCObj_PFP_LH_bwd_mu));
   t -> Branch("TPCObj_PFP_LH_bwd_p", &(vars->TPCObj_PFP_LH_bwd_p));
   t -> Branch("TPCObj_PFP_LH_bwd_pi", &(vars->TPCObj_PFP_LH_bwd_pi));
   t -> Branch("TPCObj_PFP_LH_MIP", &(vars->TPCObj_PFP_LH_MIP));
   t -> Branch("TPCObj_PFP_PIDA", &(vars->TPCObj_PFP_PIDA));
   t -> Branch("TPCObj_PFP_track_depE", &(vars->TPCObj_PFP_track_depE));
   t -> Branch("TPCObj_PFP_track_Chi2Proton", &(vars->TPCObj_PFP_track_Chi2Proton));
   t -> Branch("TPCObj_PFP_track_Chi2Muon", &(vars->TPCObj_PFP_track_Chi2Muon));
   t -> Branch("TPCObj_PFP_track_Chi2Pion", &(vars->TPCObj_PFP_track_Chi2Pion));
   t -> Branch("TPCObj_PFP_track_rangeE_mu", &(vars->TPCObj_PFP_track_rangeE_mu));
   t -> Branch("TPCObj_PFP_track_rangeE_p", &(vars->TPCObj_PFP_track_rangeE_p));
   t -> Branch("TPCObj_PFP_track_Lmip_perhit",&(vars->TPCObj_PFP_track_Lmip_perhit));

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
   t -> Branch("MCP_ID", &(vars->MCP_ID));
   t -> Branch("MCP_DaughterIDs", &(vars->MCP_DaughterIDs));
   t -> Branch("MCP_StatusCode", &(vars->MCP_StatusCode));
   t -> Branch("MCP_MotherID", &(vars->MCP_MotherID));
   t -> Branch("MCP_isNuDaughter", &(vars->MCP_isNuDaughter));
   t -> Branch("MCP_StartPosition", &(vars->MCP_StartPosition));
   t -> Branch("MCP_EndPosition", &(vars->MCP_EndPosition));

   t -> Branch("nu_vtx", &(vars->nu_vtx));
   t -> Branch("nu_vtx_spacecharge", &(vars->nu_vtx_spacecharge));
   t -> Branch("nu_isCC", &(vars->nu_isCC));
   t -> Branch("nu_PDG", &(vars->nu_PDG));
   t -> Branch("nu_E", &(vars->nu_E));
   t -> Branch("nu_MCPID", &(vars->nu_MCPID));
   t -> Branch("nuint_mode", &(vars->nuint_mode));
   t -> Branch("nuint_W", &(vars->nuint_W));
   t -> Branch("nuint_Qsq", &(vars->nuint_Qsq));

}

#endif
