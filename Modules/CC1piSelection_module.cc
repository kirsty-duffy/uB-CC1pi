////////////////////////////////////////////////////////////////////////
// Class:       CC1piSelection
// Plugin Type: producer (art v2_05_01)
// File:        CC1piSelection_module.cc
//
// Generated at Fri Jan 19 09:15:19 2018 by Kirsty Duffy using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

// Standard art includes
#include "art/Framework/Core/EDProducer.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "art/Framework/Services/Optional/TFileService.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"

// Larsoft data objects includes
#include "larcoreobj/SummaryData/POTSummary.h"

#include <memory>

// ROOT includes
#include "TTree.h"

// Include CC1pi files from elsewhere
#include "uboone/CC1pi/Algorithms/FillTree.h"
#include "uboone/CC1pi/Algorithms/TwoTrackCheck.h"

class CC1piSelection;


class CC1piSelection : public art::EDProducer {
   public:
      explicit CC1piSelection(fhicl::ParameterSet const & p);
      // The compiler-generated destructor is fine for non-base
      // classes without bare pointers or other resource use.

      // Plugins should not be copied or assigned.
      CC1piSelection(CC1piSelection const &) = delete;
      CC1piSelection(CC1piSelection &&) = delete;
      CC1piSelection & operator = (CC1piSelection const &) = delete;
      CC1piSelection & operator = (CC1piSelection &&) = delete;

      // Required functions.
      void produce(art::Event & e) override;
      void endSubRun(art::SubRun &sr) override;

   private:

      // Declare member data here.

      // TFile service
      art::ServiceHandle<art::TFileService> tfs;

      // Output trees
      TTree *_outtree;
      TTree *_pot_tree;

      // Variables for output tree
      cc1pianavars *anavars;

      // CC1pi selection variables
      bool _PassesCC1piSelec;
      std::string _CC1piSelecFailureReason;

      // Other variables set in module
      bool _isData;
      double _sr_pot;

};


CC1piSelection::CC1piSelection(fhicl::ParameterSet const & p)
   // :
   // Initialize member data here.
{
   // Call appropriate produces<>() functions here.
   // Things we want to produce (i.e. add to the event):
   //  - Flag for whether it passes our selection
   //  - Failure reason if it doesn't
   produces< bool >();
   produces< std::string >();


   // Start out assuming everything passes
   _PassesCC1piSelec = true;
   _CC1piSelecFailureReason = "";

   // Instantiate struct to hold variables (and pass fhicl parameters)
   anavars = new cc1pianavars(p);

   // Instantiate tree and make branches
   _outtree = tfs->make<TTree>("outtree","");
   MakeAnaBranches(_outtree,anavars);

   _pot_tree = tfs->make<TTree>("pottree","");
   _pot_tree->Branch("pot", &_sr_pot, "pot/D");
   _isData = false;
}

void CC1piSelection::produce(art::Event & evt)
{
   // Check if we're looking at data or MC (important for POT counting in endSubRun
   _isData = evt.isRealData();

   // Set anavars values to default
   // (this will prevent bugs if any variables don't get set for a given event)
   anavars->Clear();

   // ----- Cut: must pass Marco's selection ----- //
   art::Handle<std::vector<ubana::SelectionResult>> selection_h;
   evt.getByLabel("UBXSec", selection_h);
   if(!selection_h.isValid()){
      mf::LogError(__PRETTY_FUNCTION__) << "SelectionResult product not found." << std::endl;
      throw std::exception();
   }
   std::vector<art::Ptr<ubana::SelectionResult>> selection_v;
   art::fill_ptr_vector(selection_v, selection_h);

   bool PassesMarcosSelec = selection_v.at(0)->GetSelectionStatus();
   if (!PassesMarcosSelec){
      _PassesCC1piSelec = false;
      _CC1piSelecFailureReason = selection_v.at(0)->GetFailureReason();;
   }

   if (_PassesCC1piSelec) // don't bother evaluating next cut if it's already failed
   {
      // ----- Cut: minimum two tracks (don't have to be MIPs) ----- //
      bool TwoTracks = TwoTrackCheck(evt, false, false);
      if (!TwoTracks){
         _PassesCC1piSelec = false;
         _CC1piSelecFailureReason="TwoTrackCut";
      }
   }

   if (_PassesCC1piSelec) // don't bother evaluating next cut if it's already failed
      // This means that this cut will only cut out events that have >=2 tracks but *not* MIP-like
      // Events that have <2 tracks will already have failed
   {
      // ----- Cut: minimum two MIP-consistent tracks ----- //
      bool TwoMIPTracks = TwoTrackCheck(evt, true, false);
      if (!TwoMIPTracks){
         _PassesCC1piSelec = false;
         _CC1piSelecFailureReason="TwoMIPCut";
      }
   }

   if (_PassesCC1piSelec) // don't bother evaluating next cut if it's already failed
      // This means that this cut will only cut out events that have >2 MIP-like tracks
      // Events that have <2 tracks total or with <2 MIP-like tracks will already have failed
   {
      // ----- Cut: exactly two MIP-consistent tracks ----- //
      bool TwoMIPTracks = TwoTrackCheck(evt, true, true);
      if (!TwoMIPTracks){
         _PassesCC1piSelec = false;
         _CC1piSelecFailureReason="ExactlyTwoMIPCut";
      }
   }


   // ----- Almost at the end: fill tree ------ //

   // Set anavars values that are already in the reco2 file
   anavars->SetReco2Vars(evt);

   // Set anavars for _PassesCC1piSelec and _CC1piSelecFailureReason
   anavars->PassesCC1piSelec = _PassesCC1piSelec;
   anavars->CC1piSelecFailureReason = _CC1piSelecFailureReason;

   _outtree -> Fill();


   // ----- Finally: add things to event ------ //
   std::unique_ptr<bool> PassesCC1piSelec = std::make_unique<bool>(_PassesCC1piSelec);
   std::unique_ptr<std::string> CC1piSelecFailureReason = std::make_unique<std::string>(_CC1piSelecFailureReason); 
   evt.put(std::move(PassesCC1piSelec));
   evt.put(std::move(CC1piSelecFailureReason));

}

void CC1piSelection::endSubRun(art::SubRun &sr) {
   // Note: the entire subrun's POT is recorded in the tree for every event.
   // You must only add it once per subrun to get the correct number.

   std::cout << "Hello I'm ending a subRun and setting POT" << std::endl;
   art::Handle<sumdata::POTSummary> potsum_h;

   // MC
   if (!_isData) {
      if(sr.getByLabel("generator", potsum_h)) {
         _sr_pot = potsum_h->totpot;
      }
      else {
         _sr_pot = 0.;
         //Should this raise an error?
      }
   }

   // Data
   else {
      if(sr.getByLabel("beamdata", "bnbETOR860", potsum_h)){
         _sr_pot = potsum_h->totpot;
      }
      else {
         _sr_pot = 0.;
      }
   }

   _pot_tree->Fill();

}

DEFINE_ART_MODULE(CC1piSelection)
