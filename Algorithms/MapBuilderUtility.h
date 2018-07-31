//Copied from A. Lister: https://github.com/admlw/NuMuSelection/blob/master/Algos/mapBuilderUtility.h

#ifndef MAPBUILDERUTILITY_H
#define MAPBUILDERUTILITY_H

// art includes
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"

// larsoft includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/MCSFitResult.h"

// cpp includes
#include <map>
#include <utility>

class MapBuilderUtil{

  public:

    /**
     * Builds a map between a track ID and a art::Ptr to a MCSFitResult
     */
    std::map<int, art::Ptr<recob::MCSFitResult>> buildTrackIdMcsResultMap(std::vector< art::Ptr<recob::Track> > trackPtrVec, std::vector< art::Ptr<recob::MCSFitResult> > mcsPtrVec);

};

#endif
