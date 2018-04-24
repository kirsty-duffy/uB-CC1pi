#ifndef TWOTRACKCHECK_CXX
#define TWOTRACKCHECK_CXX

#include "TwoTrackCheck.h"

std::map<std::string,bool> TwoTrackCheck(art::Event &evt)
{
   std::map<std::string,bool> CC1picutflow;

   art::Handle<std::vector<ubana::SelectionResult>> selection_h;
   evt.getByLabel("UBXSec", selection_h);
   if(!selection_h.isValid()){
      mf::LogError(__PRETTY_FUNCTION__) << "SelectionResult product not found." << std::endl;
      throw std::exception();
   }

   //Get TPCObject
   art::FindManyP<ubana::TPCObject> tpcobject_from_selection(selection_h, evt, "UBXSec");
   if(tpcobject_from_selection.at(0).size() == 0) {
      //No TPCObject
      CC1picutflow["TwoTrackCut"] = false;
      CC1picutflow["TwoMIPCut"] = false;
      CC1picutflow["ExactlyTwoMIPCut"] = false;
      return CC1picutflow;
   }
   art::Ptr<ubana::TPCObject> tpcobj_candidate = tpcobject_from_selection.at(0).at(0);
   art::Handle<std::vector<ubana::TPCObject>> tpcobj_h;
   evt.getByLabel("TPCObjectMaker", tpcobj_h);

   //Get PFPs (in TPCObject)
   art::FindManyP<recob::PFParticle> pfps_from_tpcobject(tpcobj_h, evt, "TPCObjectMaker");
   std::vector<art::Ptr<recob::PFParticle>> pfps = pfps_from_tpcobject.at(tpcobj_candidate.key());

   //Find neutrino
   int nuID = -1;
   for (auto pfp : pfps) {
      if(lar_pandora::LArPandoraHelper::IsNeutrino(pfp)) {
         nuID = (int)pfp->Self();
         break;
      }
   }

   if(nuID == -1) {
      mf::LogError(__PRETTY_FUNCTION__) << "No neutrino PFP in TPCObject" << std::endl;
      throw std::exception();
   }

   //Find track-like daughters of the neutrino
   std::vector<art::Ptr<recob::PFParticle>> daughter_track_pfps;
   for (auto pfp : pfps) {
      if(pfp->Parent()==(size_t)nuID) {
         daughter_track_pfps.emplace_back(pfp);
      }
   }

   if(daughter_track_pfps.size() < 2) {
      CC1picutflow["TwoTrackCut"] = false;
      CC1picutflow["TwoMIPCut"] = false;
      CC1picutflow["ExactlyTwoMIPCut"] = false;
      return CC1picutflow;
   }
   else CC1picutflow["TwoTrackCut"] = true;

   //Get PFPs (in event)
   art::Handle<std::vector<recob::PFParticle> > pfp_h;
   evt.getByLabel("pandoraNu::UBXSec",pfp_h);
   if(!pfp_h.isValid()){
      mf::LogError(__PRETTY_FUNCTION__) << "PFP product not found." << std::endl;
      throw std::exception();
   }

   art::FindManyP<recob::Track> tracks_from_pfps(pfp_h, evt, "pandoraNu::UBXSec");

   // Also get track-PID association objects
   auto const& track_h = evt.getValidHandle<std::vector<recob::Track>>("pandoraNu::UBXSec");
   art::FindManyP<anab::ParticleID> trackPIDAssn(track_h, evt, "pid");

   // Every PFP could in theory have more than one track
   // At least for pandoraNu, however, it "should" always be one-to-one
   // So raise an error if this isn't the case
   int MIPs = 0;
   for (auto pfp : daughter_track_pfps) {
      std::vector<art::Ptr<recob::Track>> daughter_tracks = tracks_from_pfps.at(pfp.key());
      if(daughter_tracks.size() > 1) {
         mf::LogError(__PRETTY_FUNCTION__) << "PFP associated to more than one track." << std::endl;
         throw std::exception();
      }
      for (auto track : daughter_tracks){
         if(IsMIP(trackPIDAssn,track->ID())) MIPs++;
      }
   }

   if(MIPs < 2) {
      CC1picutflow["TwoMIPCut"] = false;
      CC1picutflow["ExactlyTwoMIPCut"] = false;
      return CC1picutflow;
   }
   else CC1picutflow["TwoMIPCut"] = true;

   if(MIPs != 2) {
      CC1picutflow["ExactlyTwoMIPCut"] = false;
      return CC1picutflow;
   }
   else CC1picutflow["ExactlyTwoMIPCut"] = true;

   return CC1picutflow;
}

#endif
