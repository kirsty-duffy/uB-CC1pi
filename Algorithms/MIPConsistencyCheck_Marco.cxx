// MIP Consistency Check copied from Marco's code on 25th Jan 2018
// This is not intended to be the final MIP consistency cut, just a placeholder
// until we can put in something more concrete.

// It uses a hardcoded array for the cut in length-dqdx space. I have copied
// this from Marco's github (fcl file UBXSec/job/muoncandidatefinder.fcl) but
// put the hardcoded array directly in the header file

// -- Kirsty Duffy, Fermilab, 25/01/2018

// Check if we've already defined a MIP consistency function (don't redefine!)
#ifndef MIPCONSISTENCY_CXX
#define MIPCONSISTENCY_CXX

// Don't need any include statements other than to include MIPConsistencyCheck_PIDA.h
// All the other necessary includes should be in there
#include "MIPConsistencyCheck_Marco.h"


bool IsMIP(art::Ptr<recob::Track> track, art::Event &evt)
{
  // Method relies on track length and truncated mean of dqdx
  double length, dqdx_truncmean;

  // Get length directly from track object
  length = track->Length();

  // Now get dqdx truncated mean
  // For now use uncalibrated "calo" variables
  // Calibrated dqdx will be available in MCC 8.6, but currently has a bug
  // use association: pandoraNucalo, art::Assns<recob::Track,anab::Calorimetry,void>
  // also hardcode that this is for pandoraNu only
  art::InputTag _trackTag = "pandoraNu";
  art::InputTag _caloTag = "pandoraNucalo";

  auto const& track_h = evt.getValidHandle<std::vector<recob::Track>>(_trackTag);
  //std::cout << "tracks.size() = " << track_h->size() << std::endl;
  
  art::FindManyP<anab::Calorimetry> calos_from_tracks(track_h, evt, _caloTag);
  //std::cout << "calos_from_tracks.size() = " << calos_from_tracks.size() << std::endl;
  unsigned int trkid = track->ID();
  std::vector<art::Ptr<anab::Calorimetry>> calos = calos_from_tracks.at(trkid);

  // Got calorimetry object, now work out truncated mean (plane 2 only)
  dqdx_truncmean = GetDqDxTruncatedMean(calos);
  
  // Now evaluate whether it passes the cut!

  // First: all long tracks are MIPs
  if (length > 1000)
    return true;
  
  // If length is negative something has gone very wrong 
  if (length < 0) {
    std::cout << "[MIPConsistencyCheck_Marco] Track length < 0?!" << std::endl;
    return false;
  }

  // Now evaluate Marco's cut
  int l = std::round(length);
  double dqdx_cut = _dqdx_cutvals.at(l);
  
  std::cout << "[MIPConsistencyCheck_Marco] Track length is " << length << ", dqdx_cut is " << dqdx_cut << std::endl;
  std::cout << "[MIPConsistencyCheck_Marco] \t Truncated mean dQ/dx for candidate is: " << dqdx_truncmean << std::endl;
  std::cout << "[MIPConsistencyCheck_Marco] \t MIP consistent ? : " << (dqdx_truncmean <= dqdx_cut ? "YES" : "NO") << std::endl;
  
  if (dqdx_truncmean <= dqdx_cut){
    return true;
  }
  
  return false;
  
}

// ------------------------------------------------------------------------------- //

double GetDqDxTruncatedMean(std::vector<art::Ptr<anab::Calorimetry>> calos)
{
  for (auto c : calos) {
    if (!c) continue; // avoid art errors if c doesn't exist
    if (!c->PlaneID().isValid) continue; // check has valid plane
    int planenum = c->PlaneID().Plane;
    if (planenum != 2) continue; // only use calorimetry from plane 2
    
    std::vector<double> dqdx_v = c->dQdx(); 
    
    return GetDqDxTruncatedMean(dqdx_v);
  }
  
  return -9999;
}

// ------------------------------------------------------------------------------- //

double GetDqDxTruncatedMean(std::vector<double> dqdx_v)
{
  double result = -9999;
  double n = 1.;
  
  if (dqdx_v.size() == 0)
    return result;
  
  //for (auto q : dqdx_v) {
  //std::cout << "dqdx before trim: " << q * 198 << std::endl;
  //}
  
  double median = GetMedian_dqdx(dqdx_v);
  double std    = GetSTD_dqdx(dqdx_v);
  
  // std::cout << "median dqdx = " << median << std::endl;
  //std::cout << "std dqdx =    " << std << std::endl;
  
  std::vector<double> dqdx_v_trimmed;
  dqdx_v_trimmed.clear();
  
  for (auto q : dqdx_v) {
    if (q > median - n * std && 
	q < median + n * std) {
      dqdx_v_trimmed.emplace_back(q);
      //std::cout << "dqdx after trim: " << q * 198 << std::endl;
    }
  }
  
  result = GetMean_dqdx(dqdx_v_trimmed);
  
  return result;
}

// ------------------------------------------------------------------------------- //

double GetMedian_dqdx(std::vector<double> dqdx_v)
{
  double median = -9999;
  size_t size = dqdx_v.size();
  
  if (size == 0)
    return median;
  
  std::sort(dqdx_v.begin(), dqdx_v.end());
  if (size % 2 == 0){
    median = (dqdx_v[size/2 - 1] + dqdx_v[size/2]) / 2;
  }
  else{
    median = dqdx_v[size/2];
  }
  
  return median;
}

// ------------------------------------------------------------------------------- //

double GetMean_dqdx(std::vector<double> dqdx_v)
{
  double mean = -9999;
  size_t size = dqdx_v.size();

  if (size == 0)
    return mean;

  double sum = 0;

  for (auto v : dqdx_v) {
    sum += v;
  }

  mean = sum/(double)size;

  return mean;
}

// ------------------------------------------------------------------------------- //

double GetSTD_dqdx(std::vector<double> dqdx_v)
{

  double variance = -1;

  double sum = 0;
  double sum2 = 0;
  size_t size = dqdx_v.size();

  if (size == 0)
    return -9999;

  for (auto value : dqdx_v) {

    sum  += value;
    sum2 += value*value;

  }  

  variance = sum2/(double)size - (sum/(double)size)*(sum/(double)size);

  if (variance > 0){
    return std::sqrt(variance);
  }
  else {
    return -9999;
  }
}

#endif
