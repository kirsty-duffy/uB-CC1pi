////////////////////////////////////////////////////////////////////////
// Class:       TrueCC1piFilter
// Plugin Type: filter (art v2_05_01)
// File:        TrueCC1piFilter_module.cc
//
// Generated at Thu Feb  1 11:31:03 2018 by Kirsty Duffy using cetskelgen
// from cetlib version v1_21_00.
////////////////////////////////////////////////////////////////////////

// standard module/art includes
#include "art/Framework/Core/EDFilter.h"
#include "art/Framework/Core/ModuleMacros.h"
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

// standard C++ includes
#include <iostream>
#include <string>
#include <memory>

// Larsoft object includes
#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"

class TrueCC1piFilter;


class TrueCC1piFilter : public art::EDFilter {
public:
  explicit TrueCC1piFilter(fhicl::ParameterSet const & p);
  // The compiler-generated destructor is fine for non-base
  // classes without bare pointers or other resource use.
  
  // Plugins should not be copied or assigned.
  TrueCC1piFilter(TrueCC1piFilter const &) = delete;
  TrueCC1piFilter(TrueCC1piFilter &&) = delete;
  TrueCC1piFilter & operator = (TrueCC1piFilter const &) = delete;
  TrueCC1piFilter & operator = (TrueCC1piFilter &&) = delete;
  
  // Required functions.
  bool filter(art::Event & e) override;
  
private:
  
  // Declare member data here.
  bool _isCC1pi;
};


TrueCC1piFilter::TrueCC1piFilter(fhicl::ParameterSet const & p) :
  _isCC1pi(false)
  // Initialize member data here.
{
// Call appropriate produces<>() functions here.
}

bool TrueCC1piFilter::filter(art::Event & evt)
{
  // Check if we're looking at data or MC
  // Give an error if data - this filter won't work!
  if (evt.isRealData()){
    std::cout << "[ERROR : TrueCC1piFilter] Attempting to run on real data! Returning false for all events" << std::endl;
    return false;
  }
  
  // Get all MCParticles and associated MCtruths
  auto const& mcp_handle = evt.getValidHandle<std::vector<simb::MCParticle>>("largeant");
  art::FindManyP<simb::MCTruth> MCtruth_from_mcp(mcp_handle, evt, "largeant");
  art::Ptr<simb::MCTruth> muon_MCTruth;
  
  // Loop through MCParticles and count number of primary muons and pions
  int primary_muons = 0;
  int primary_pions = 0;
  for (unsigned int i_mcp=0; i_mcp < mcp_handle->size(); i_mcp++){
    auto mcp = mcp_handle->at(i_mcp);
    auto mctruths = MCtruth_from_mcp.at(i_mcp);
    
    
    // Check every variable is doing what I think it does
    if (mctruths.size() != 1) std::cout << "[ERROR : TrueCC1piFilter] mctruths.size = " << mctruths.size() << " but code assumes it will be 1" << std::endl;
    
    if (mcp.Process() != "primary" && mcp.Mother() == 0) std::cout << "[ERROR : TrueCC1piFilter] mcp.Process = " << mcp.Process() << " and mother = neutrino" << std::endl;
    
    if (mcp.Process() == "primary" && mcp.Mother() != 0) std::cout << "[ERROR : TrueCC1piFilter] mcp.Process = " << mcp.Process() << " and mother = " << mcp.Mother() << std::endl;
    
    
    // Only look at primary particles
    if (mcp.Process() != "primary"){
      //std::cout << "Not primary, skipping event: mcp.Process = " << mcp.Process() << std::endl;
      continue;
    }

    // Only look at particles coming from the neutrino
    if (mcp.Mother() != 0){
      continue;
    }

    // Only look at true beam neutrino interactions
    if (mctruths.at(0)->Origin() != simb::kBeamNeutrino){
      //std::cout << "[ERROR : TrueCC1piFilter] Not neutrino event. Origin: " << mctruths.at(0)->Origin() << std::endl;
      continue;
    }

    // Only look at particles with status code = 1
    if (mcp.StatusCode() != 1){
      std::cout << "MCP status != 1" << std::endl << mcp << std::endl;
      continue;
    }

    
    // Now: count the mu-s and pi+s    
    if (mcp.PdgCode() == 13){ // true mu-
      primary_muons++;
      muon_MCTruth = mctruths.at(0);
    }

    if (mcp.PdgCode() == 211){ // true pi+
      primary_pions++;
    }
  } // end loop over MCParticles

  // Is true CC1pi+?
  if (primary_muons == 1 && primary_pions == 1){
    _isCC1pi = true;
    //std::cout << "Selected event mode: " << muon_MCTruth->GetNeutrino().Mode() << std::endl;
  }
  else{
    _isCC1pi = false;
  }
  
  return _isCC1pi;
}
  
DEFINE_ART_MODULE(TrueCC1piFilter)
