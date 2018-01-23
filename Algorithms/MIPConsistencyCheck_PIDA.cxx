// Check if we've already defined a MIP consistency function (don't redefine!)
#ifndef MIPCONSISTENCY_CXX
#define MIPCONSISTENCY_CXX

// Don't need any include statements other than to include MIPConsistencyCheck_PIDA.h
// All the other necessary includes should be in there
#include "MIPConsistencyCheck_PIDA.h"

bool IsMIP(art::Ptr<recob::Track> track, art::Event &evt)
{
  // do stuff and return whether the particle is a MIP
  return true;
}

#endif
