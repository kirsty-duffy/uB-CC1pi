#ifndef SHOWERREJECTION_H
#define SHOWERREJECTION_H

// Marco's code includes
#include "uboone/UBXSec/DataTypes/TPCObject.h"
#include "uboone/UBXSec/DataTypes/SelectionResult.h"
#include "uboone/UBXSec/Algorithms/TrackQuality.h"

// larsoft includes
#include "larcore/Geometry/Geometry.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardataobj/RecoBase/PFParticle.h"
#include "lardataobj/RecoBase/SpacePoint.h"
#include "lardataobj/RecoBase/Track.h"
#include "larpandora/LArPandoraInterface/LArPandoraHelper.h"

// art includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

// CC1pi method includes
#include "uboone/CC1pi/Algorithms/InputTags.h"

std::map<std::string,bool> ShowerRejection(art::Event &evt, InputTags *CC1piInputTags);

void ShowerCheck(art::Event &evt, InputTags *CC1piInputTags, art::Ptr<recob::PFParticle> pfp, art::Ptr<recob::Track> track, std::pair<double, double> &residual_mean_std, double &ratio);

#endif
