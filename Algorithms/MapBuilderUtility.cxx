#include "MapBuilderUtility.h"

std::map< int, art::Ptr<recob::MCSFitResult> >MapBuilderUtil::buildTrackIdMcsResultMap(std::vector<art::Ptr<recob::Track> > trackPtrVec, std::vector< art::Ptr<recob::MCSFitResult> > mcsPtrVec){

  // if the track ptr vector and the mcs ptr vector aren't the same size that's a problem
  size_t nTracks = trackPtrVec.size();
  size_t nMcs    = mcsPtrVec.size();

  if (nTracks != nMcs) throw std::exception();

  // create map
  std::map< int, art::Ptr<recob::MCSFitResult> > thisMap;

  for (size_t i = 0; i < nTracks; i++){

    art::Ptr<recob::Track> thisTrack = trackPtrVec.at(i);
    art::Ptr<recob::MCSFitResult> thisMcs = mcsPtrVec.at(i);

    // the casting is necessary... not sure why
    thisMap.insert( std::make_pair< int, art::Ptr<recob::MCSFitResult> >(thisTrack->ID(), (art::Ptr<recob::MCSFitResult>)thisMcs));

  }

return thisMap;

}
