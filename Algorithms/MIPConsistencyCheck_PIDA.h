// Check if we've already defined a MIP consistency header file (don't redefine)
#ifndef MIPCONSISTENCY_H
#define MIPCONSISTENCY_H

// art includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"

// larsoft includes
#include "lardataobj/RecoBase/Track.h"


// Now declare the function
bool IsMIP(art::Ptr<recob::Track> track, art::Event &evt);

#endif
