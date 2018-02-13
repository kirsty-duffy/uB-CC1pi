// Marco's code includes
#include "uboone/UBXSec/DataTypes/TPCObject.h"
#include "uboone/UBXSec/DataTypes/SelectionResult.h"

// larsoft includes
#include "lardataobj/RecoBase/Vertex.h"

// art includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "messagefacility/MessageLogger/MessageLogger.h"
#include "canvas/Persistency/Common/FindManyP.h"

bool FVCheck(art::Event &evt);

bool inFV(double x, double y, double z);
