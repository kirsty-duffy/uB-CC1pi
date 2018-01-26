#include "TwoTrackCheck.h"

bool TwoTrackCheck(art::Event &evt, bool checkMIPs)
{
   art::Handle<std::vector<ubana::SelectionResult>> selection_h;
   evt.getByLabel("UBXSec", selection_h);
   if(!selection_h.isValid()){
      mf::LogError(__PRETTY_FUNCTION__) << "SelectionResult product not found." << std::endl;
      throw std::exception();
   }

   //Get TPCObject
   art::FindManyP<ubana::TPCObject> tpcobject_from_selection(selection_h, evt, "UBXSec");
   if(tpcobject_from_selection.at(0).size() == 0) return false; //No TPCObject
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
     if(lar_pandora::LArPandoraHelper::IsTrack(pfp) && pfp->Parent()==(size_t)nuID) {
         daughter_track_pfps.emplace_back(pfp);
      }
   }

   if(daughter_track_pfps.size() < 2) return false;

   else if(!checkMIPs) return true;

   else {
      //Get PFPs (in event)
      art::Handle<std::vector<recob::PFParticle> > pfp_h;
      evt.getByLabel("pandoraNu::UBXSec",pfp_h);
      if(!pfp_h.isValid()){
         mf::LogError(__PRETTY_FUNCTION__) << "PFP product not found." << std::endl;
         throw std::exception();
      }

      art::FindManyP<recob::Track> tracks_from_pfps(pfp_h, evt, "pandoraNu");

      // This is a little complicated. Every PFP could in theory have more than one track
      // (is that true?)
      // Anyway, right now if any track associated with a daughter PFP is MIP-like, that
      // counts for one MIP-like track
      int MIPs = 0;
      for (auto pfp : daughter_track_pfps) {
	std::vector<art::Ptr<recob::Track>> daughter_tracks = tracks_from_pfps.at(pfp.key());
	for (auto track : daughter_tracks){
	  if(IsMIP(track, evt))
	    {
	      MIPs++;
	      break; // break out of the loop over tracks so the maximum one PFP can contribute is 1
	    }
	}
      }

      if(MIPs >= 2) return true;
      else return false;
   }
}
