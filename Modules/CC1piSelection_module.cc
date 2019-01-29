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
#include "uboone/CC1pi/Algorithms/InputTags.h"
#include "uboone/CC1pi/Algorithms/ShowerRejection.h"

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

      // Input Tags
      InputTags *CC1piInputTags;

      // CC1pi selection variables
      std::map<std::string,bool> _CC1picutflow;

      // Other variables set in module
      bool _isData;
      double _sr_pot;
      unsigned int _run_num;
      unsigned int _subrun_num;

};


CC1piSelection::CC1piSelection(fhicl::ParameterSet const & p)
   // :
   // Initialize member data here.
{
   // Call appropriate produces<>() functions here.
   // Things we want to produce (i.e. add to the event):
   //  - Cutflow map of whether it passes each cut (possibly also pseudo-cuts such as "passes the whole selection")
   produces< std::map<std::string,bool> >();

   // Instantiate struct to hold variables (and pass fhicl parameters)
   anavars = new cc1pianavars(p);

   // Instantiate struct to hold input tags
   CC1piInputTags = new InputTags(p);

   // Instantiate tree and make branches
   _outtree = tfs->make<TTree>("outtree","");
   MakeAnaBranches(_outtree,anavars);

   _pot_tree = tfs->make<TTree>("pottree","");
   _pot_tree->Branch("pot", &_sr_pot, "pot/D");
   _pot_tree->Branch("run_num", &_run_num);
   _pot_tree->Branch("subrun_num", &_subrun_num);
   _isData = false;
}

void CC1piSelection::produce(art::Event & evt)
{
   // Check if we're looking at data or MC (important for POT counting in endSubRun)
   _isData = evt.isRealData();
   std::cout << "[CC1pi] Using isData = " << _isData << std::endl;

   // Check TPC boundary
   std::cout << "TPC Boundary (according to Geometry): " << std::endl;
   auto const* geo = lar::providerFrom<geo::Geometry>();
   for (size_t i = 0; i < geo -> NTPC(); i++){
      double local[3] = {0.,0.,0.};
      double world[3] = {0.,0.,0.};
      const geo::TPCGeo &tpc = geo->TPC(i);
      tpc.LocalToWorld(local,world);
      std::cout << "minx = " << world[0]-geo->DetHalfWidth(i) << std::endl;
      std::cout << "maxx = " << world[0]+geo->DetHalfWidth(i) << std::endl;
      std::cout << "miny = " << world[1]-geo->DetHalfHeight(i) << std::endl;
      std::cout << "maxy = " << world[1]+geo->DetHalfHeight(i) << std::endl;
      std::cout << "minz = " << world[2]-geo->DetLength(i)/2. << std::endl;
      std::cout << "maxz = " << world[2]+geo->DetLength(i)/2. << std::endl;
   }


   // Set anavars values to default
   // (this will prevent bugs if any variables don't get set for a given event)
   anavars->Clear();

   // ----- Cut: must pass Marco's selection ----- //
   art::Handle<std::vector<ubana::SelectionResult>> selection_h;
   evt.getByLabel(CC1piInputTags->fSelectionLabel, selection_h);
   if(!selection_h.isValid()){
      mf::LogError(__PRETTY_FUNCTION__) << "[CC1piSelection] SelectionResult product not found." << std::endl;
      throw std::exception();
   }
   std::vector<art::Ptr<ubana::SelectionResult>> selection_v;
   art::fill_ptr_vector(selection_v, selection_h);

   bool PassesMarcosSelec = selection_v.at(0)->GetSelectionStatus();
   _CC1picutflow = selection_v.at(0)->GetCutFlowStatus();

   std::map<std::string,bool> TwoTrack_cutflow = TwoTrackCheck(evt, CC1piInputTags);
   _CC1picutflow.insert(TwoTrack_cutflow.begin(), TwoTrack_cutflow.end());

   std::map<std::string,bool> ShowerRejection_cutflow = ShowerRejection(evt, CC1piInputTags);
   _CC1picutflow.insert(ShowerRejection_cutflow.begin(), ShowerRejection_cutflow.end());

   // Check if "passing" events still pass the full FV cut
   if (PassesMarcosSelec) {
      if(!FVCheck(evt, CC1piInputTags)) {
         PassesMarcosSelec = false;
         _CC1picutflow["fiducial_volume"] = false;
      }
   }

   if (PassesMarcosSelec) _CC1picutflow["MarcosSelec"] = true;
   else _CC1picutflow["MarcosSelec"] = false;

   _CC1picutflow["CC1piSelec"] = true;
   for (std::map<std::string,bool>::iterator iter=_CC1picutflow.begin(); iter!=_CC1picutflow.end(); ++iter) {
      if(iter->second == false) _CC1picutflow["CC1piSelec"] = false;
   }

   // ----- Almost at the end: fill tree ------ //

   // Set anavars values that are already in the reco2 file
   anavars->SetReco2Vars(evt);

   // Set anavars for _CC1picutflow
   anavars->CC1picutflow = _CC1picutflow;

   _outtree -> Fill();


   // ----- Finally: add things to event ------ //
   std::unique_ptr<std::map<std::string,bool>> CC1picutflow = std::make_unique<std::map<std::string,bool>>(_CC1picutflow);

   evt.put(std::move(CC1picutflow));

}

void CC1piSelection::endSubRun(art::SubRun &sr) {
   // Note: the entire subrun's POT is recorded in the tree for every event.
   // You must only add it once per subrun to get the correct number.

   _run_num = sr.run();
   _subrun_num = sr.subRun();

   std::cout << "Hello I'm ending a subRun and setting POT" << std::endl;
   art::Handle<sumdata::POTSummary> potsum_h;

   // MC
   if (!_isData) {
      std::cout << "!_isData" << std::endl;
      if(sr.getByLabel("generator", potsum_h)) {
         std::cout << "POT found: " << potsum_h->totpot << std::endl;
         _sr_pot = potsum_h->totpot;
      }
      else {
         std::cout << "POT not found. Setting to 0." << std::endl;
         _sr_pot = 0.;
         //Should this raise an error?
      }
   }

   // Data
   else {
      std::cout << "_isData" << std::endl;
      if(sr.getByLabel("beamdata", "bnbETOR860", potsum_h)){
         std::cout << "POT found: " << potsum_h->totpot << std::endl;
         _sr_pot = potsum_h->totpot;
      }
      else {
         std::cout << "POT not found. Setting to 0." << std::endl;
         _sr_pot = 0.;
      }
   }

   std::cout << "Filling POT tree" << std::endl;
   _pot_tree->Fill();
}

DEFINE_ART_MODULE(CC1piSelection)
