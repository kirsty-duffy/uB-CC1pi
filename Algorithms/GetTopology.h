#ifndef GETTOPOLOGY_H
#define GETTOPOLOGY_H

#include "uboone/CC1pi/Algorithms/TopologyEnums.h"

#include "nusimdata/SimulationBase/MCParticle.h"
#include "nusimdata/SimulationBase/MCTruth.h"


#include "art/Framework/Principal/Handle.h"

#include <vector>
#include <iostream>

NuIntTopology GetTopology(art::ValidHandle<std::vector<simb::MCTruth>> mc_truths);

#endif
