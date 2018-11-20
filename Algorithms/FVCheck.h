#ifndef FVCHECK_H
#define FVCHECK_H

// Marco's code includes
#include "uboone/UBXSec/DataTypes/TPCObject.h"
#include "uboone/UBXSec/DataTypes/SelectionResult.h"

// larsoft includes
#include "lardataobj/RecoBase/Vertex.h"
#include "larevt/SpaceChargeServices/SpaceChargeService.h"

// art includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

// CC1pi method includes
#include "uboone/CC1pi/Algorithms/InputTags.h"

bool FVCheck(art::Event &evt, InputTags *CC1piInputTags);

bool inFV(double x, double y, double z);

#endif
