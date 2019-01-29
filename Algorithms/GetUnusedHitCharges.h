#ifndef GETUNUSEDHITCHARGES_H
#define GETUNUSEDHITCHARGES_H

// root includes
#include "TVector.h"

// art includes
#include "art/Framework/Principal/Event.h"
#include "canvas/Persistency/Common/Ptr.h"
#include "canvas/Persistency/Common/PtrVector.h"
#include "canvas/Persistency/Common/FindManyP.h"

// larsoft includes
#include "larcore/Geometry/Geometry.h"
#include "larcore/Geometry/GeometryCore.h"
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/RecoBase/Hit.h"
#include "lardata/DetectorInfoServices/DetectorPropertiesService.h"

// CC1pi method includes
#include "uboone/CC1pi/Algorithms/InputTags.h"

void GetUnusedHitCharges(art::Event &evt, InputTags *CC1piInputTags, int planenum, TVector3 trackend, /*std::vector<art::Ptr<recob::Hit>> &unused_hits,*/ std::vector<double> &unused_hit_charge, std::vector<int> &unused_hit_wiredist, std::vector<double> &unused_hit_timedist);

#endif
