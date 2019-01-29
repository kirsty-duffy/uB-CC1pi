#ifndef INPUTTAGS_H
#define INPUTTAGS_H

// art includes
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

struct InputTags{
   art::InputTag fTrackLabel;
   art::InputTag fCalorimetryLabel;
   art::InputTag fSelectionLabel;
   art::InputTag fPIDLabelChi2;
   art::InputTag fPIDProducer;
   art::InputTag fPFParticleProducer;
   art::InputTag fClusterProducer;
   art::InputTag fHitProducer;
   art::InputTag fSpacePointProducer;
   art::InputTag fTPCObjectProducer;
   art::InputTag fMCSMuProducer;
   art::InputTag fMCSPProducer;
   art::InputTag fMCSPiProducer;
   double fResidualsStdCutUp;
   double fResidualsStdCutDown;
   double fResidualsMeanCutDown;
   double fResidualsMeanCutUp;
   double fPercUsedHitsCut;

   // Constructor (pass the fhicl parameters here)
   InputTags(fhicl::ParameterSet const &p);
   void PrintConfig();
};

#endif
