#ifndef INPUTTAGS_H
#define INPUTTAGS_H

// art includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

struct InputTags{
   art::InputTag fTrackLabel;
   art::InputTag fCalorimetryLabel;
   art::InputTag fSelectionLabel;

   // Constructor (pass the fhicl parameters here)
   InputTags(fhicl::ParameterSet const &p);
};

#endif
