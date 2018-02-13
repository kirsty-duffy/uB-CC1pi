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
#include "uboone/CC1pi/Algorithms/FVCheck.h"

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
      std::map<std::string,bool> _CC1picutflow;
      bool _PassesCC1piSelec;

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
   //  - Cutflow map of whether it passes each cut
   //  - Flag for whether it passes our selection
   produces< std::map<std::string,bool> >();
   produces< bool >();
   produces< std::string >();


   // Start out assuming everything passes
   _PassesCC1piSelec = true;

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
   _CC1picutflow = selection_v.at(0)->GetCutFlowStatus();

   std::map<std::string,bool> TwoTrackcutflow = TwoTrackCheck(evt);
   _CC1picutflow.insert(TwoTrackcutflow.begin(), TwoTrackcutflow.end());

   // Check if "passing" events still pass the full FV cut
   if (PassesMarcosSelec) {
      if(!FVCheck(evt)) {
         PassesMarcosSelec = false;
         _CC1picutflow["fiducial_volume"] = false;
      }
   }

   if (!PassesMarcosSelec){
      _PassesCC1piSelec = false;
      _CC1picutflow["MarcosSelec"] = false;
   }
   else {
      _CC1picutflow["MarcosSelec"] = true;

      if(_CC1picutflow["TwoTrackCut"] == false || _CC1picutflow["TwoMIPCut"] == false || _CC1picutflow["ExactlyTwoMIPCut"] == false) _PassesCC1piSelec = false;
      else _PassesCC1piSelec = true;
   }

   // ----- Almost at the end: fill tree ------ //

   // Set anavars values that are already in the reco2 file
   anavars->SetReco2Vars(evt);

   // Set anavars for _CC1picutflow and _PassesCC1piSelec
   anavars->CC1picutflow = _CC1picutflow;
   anavars->PassesCC1piSelec = _PassesCC1piSelec;

   _outtree -> Fill();


   // ----- Finally: add things to event ------ //
   std::unique_ptr<std::map<std::string,bool>> CC1picutflow = std::make_unique<std::map<std::string,bool>>(_CC1picutflow);
   std::unique_ptr<bool> PassesCC1piSelec = std::make_unique<bool>(_PassesCC1piSelec);

   evt.put(std::move(CC1picutflow));
   evt.put(std::move(PassesCC1piSelec));

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
