// MIP Consistency Check using fit of Bragg peak to theoretical prediction for various particle types
//
// Reads in anab::ParticleID object that is produced by ParticleID/Modules/ParticleId_module.cc. Note that to run this you will need to check out the ParticleID package from github () and install it in uboonecode/uboone/ParticleID
//
// To read the data object produced by the ParticleID class you will need to download the feature branch of lardataobj: lardataobj/feature/kduffy_pidrefactor_v1_11_00_04. You can do this by:
  // - cd $MRB_SOURCE
  // - mrb g lardataobj
  // - cd lardataobj
  // - git checkout -b feature/kduffy_pidrefactor_v1_11_00_04 origin/feature/kduffy_pidrefactor_v1_11_00_04
// Both the ParticleID and lardataobj packages should be built on uboonecode v06_26_01_12 (this is true on 18th April 2018)

// -- Kirsty Duffy, Fermilab, 18/04/2018

// Check if we've already defined a MIP consistency header file (don't redefine)
#ifndef MIPCONSISTENCY_H
#define MIPCONSISTENCY_H

// art includes
#include "art/Framework/Principal/Event.h"
#include "art/Framework/Principal/Handle.h"
#include "art/Framework/Principal/Run.h"
#include "art/Framework/Principal/SubRun.h"
#include "canvas/Persistency/Common/FindManyP.h"
#include "canvas/Utilities/InputTag.h"
#include "fhiclcpp/ParameterSet.h"

// larsoft includes
#include "lardataobj/RecoBase/Track.h"
#include "lardataobj/AnalysisBase/Calorimetry.h"
#include "lardataobj/AnalysisBase/ParticleID.h"

// standard C++ includes
#include <iostream>

// Now declare the functions
bool IsMIP(art::FindManyP<anab::ParticleID> trackPIDAssn, int TrackID);

#endif
