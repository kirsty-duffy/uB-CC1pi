////////////////////////////////////////////////////////////////////////
// Class:       MinTwoTrackFilter
// Plugin Type: filter (art v2_05_01)
// File:        MinTwoTrackFilter_module.cc
//
// Generated at Wed Feb 21 10:04:34 2018 by Kirsty Duffy using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

// Standard art includes
#include "art/Framework/Core/EDFilter.h"
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

// Include CC1pi files from elsewhere
#include "uboone/CC1pi/Algorithms/TwoTrackCheck.h"
#include "uboone/CC1pi/Algorithms/FVCheck.h"

class MinTwoTrackFilter;


class MinTwoTrackFilter : public art::EDFilter {
public:
  explicit MinTwoTrackFilter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.

  // Plugins should not be copied or assigned.
  MinTwoTrackFilter(MinTwoTrackFilter const &) = delete;
  MinTwoTrackFilter(MinTwoTrackFilter &&) = delete;
  MinTwoTrackFilter & operator = (MinTwoTrackFilter const &) = delete;
  MinTwoTrackFilter & operator = (MinTwoTrackFilter &&) = delete;

  // Required functions.
  bool filter(art::Event & evt) override;

private:

  // Declare member data here.

};


MinTwoTrackFilter::MinTwoTrackFilter(fhicl::ParameterSet const & p)
// :
// Initialize member data here.
{
  // Call appropriate produces<>() functions here.
}

bool MinTwoTrackFilter::filter(art::Event & evt)
{

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

   // Check if "passing" events still pass Marco's full FV cut
   if (PassesMarcosSelec) {
      if(!FVCheck(evt)) {
         PassesMarcosSelec = false;
      }
   }
   

   // ----- Cut: must have >= 2 tracks ----- //
   std::map<std::string,bool> _TwoTrackCheck = TwoTrackCheck(evt);

   // Now decide whether to keep the event
   // Keep only events which pass Marco's selection and the TwoTrackCut check
   // Note that "TwoTrackCut" is true for events with >= 2 tracks that are reconstructed
   // as daughters of the neutrino in Marco's TPCObject
   if (_TwoTrackCheck["TwoTrackCut"] == true && PassesMarcosSelec){
     return true;
   }
   else{
     return false;
   }
}

DEFINE_ART_MODULE(MinTwoTrackFilter)
