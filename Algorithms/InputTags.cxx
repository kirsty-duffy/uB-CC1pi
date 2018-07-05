#include "InputTags.h"

InputTags::InputTags(fhicl::ParameterSet const &p){
   fTrackLabel           = p.get<std::string>("TrackLabel");
   fCalorimetryLabel     = p.get<std::string>("CalorimetryLabel");
   fSelectionLabel       = p.get<std::string>("SelectionLabel");
   fPIDLabelChi2         = p.get<std::string>("ParticleIdChi2Label");
   fPFParticleProducer   = p.get<std::string>("PFParticleProducer");
   fSpacePointProducer   = p.get<std::string>("SpacePointProducer");
   fResidualsStdCutUp    = p.get<double>("ResidualsStdCutUp");
   fResidualsStdCutDown  = p.get<double>("ResidualsStdCutDown");
   fResidualsMeanCutDown = p.get<double>("ResidualsMeanCutDown");
   fResidualsMeanCutUp   = p.get<double>("ResidualsMeanCutUp");
   fPercUsedHitsCut      = p.get<double>("PercUsedHitsCut");
}
