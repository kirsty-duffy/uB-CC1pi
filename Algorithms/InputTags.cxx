#include "InputTags.h"

InputTags::InputTags(fhicl::ParameterSet const &p){
   fTrackLabel       = p.get<std::string>("TrackLabel");
   fCalorimetryLabel = p.get<std::string>("CalorimetryLabel");
   fSelectionLabel   = p.get<std::string>("SelectionLabel");
   fPIDLabelChi2     = p.get<std::string>("ParticleIdChi2Label");
}
