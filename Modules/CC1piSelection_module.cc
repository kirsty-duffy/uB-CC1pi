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

#include <memory>

// ROOT includes
#include "TTree.h"

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

private:

  // Declare member data here.

  // TFile service
  art::ServiceHandle<art::TFileService> tfs;
  
  // Output tree
  TTree *_outtree;

  // Variables for output tree
  cc1pianavars *anavars = new cc1pianavars();

};


CC1piSelection::CC1piSelection(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.

  // Things we want to produce (i.e. add to the event):
  //  - Flag for whether it passes our selection
  //  - Failure reason if it doesn't

  // Instantiate tree and make branches
  _outtree = tfs->make<TTree>("outtree","");
  MakeAnaBranches(_outtree,anavars);
}

void CC1piSelection::produce(art::Event & e)
{
  // Implementation of required member function here.
  // This is where we'll call our functions
  // MIP consistency (actually probably call this in SelectTwoMIPTracks)
  // SelectTwoMIPTracks

  // At this point we'll want to define our new data products (flag for whether
  // this event passes our selection and failure reason. Does Marco already have a
  // data product for failure reason that we can use/write to?)
  // These will have to sync up with the things we have in our produces<>() functions
  // in the constructor CC1piSelection::CC1piSelection(fhicl::ParameterSet const & p)
  
  // FillTree
  // Make sure the new flag and failure reasons get added to the tree too
  // Do we also want to add things like PIDA or MIP consistency for each track?
  // If so we'll have to think about how to fill those.

  // Set anavars values to default
  // (this will prevent bugs if any variables don't get set for a given event)
  anavars->Clear();

  // Set anavars values that are already in the reco2 file
  anavars->SetReco2Vars(e);
  
  // This won't compile but is an example of giving it a value we calculated (not copied from artroot file)
  anavars->MIPConsistency = MIPConsistency;
}

DEFINE_ART_MODULE(CC1piSelection)
