#include "InputTags.h"

InputTags::InputTags(fhicl::ParameterSet const &p){
   fTrackLabel           = p.get<std::string>("TrackLabel");
   fCalorimetryLabel     = p.get<std::string>("CalorimetryLabel");
   fSelectionLabel       = p.get<std::string>("SelectionLabel");
   fPIDLabelChi2         = p.get<std::string>("ParticleIdChi2Label");
   fPIDProducer          = p.get<std::string>("ParticleIdProducer");
   fPFParticleProducer   = p.get<std::string>("PFParticleProducer");
   fClusterProducer   = p.get<std::string>("ClusterProducer");
   fSpacePointProducer   = p.get<std::string>("SpacePointProducer");
   fTPCObjectProducer    = p.get<std::string>("TPCObjectProducer");
   fResidualsStdCutUp    = p.get<double>("ResidualsStdCutUp");
   fResidualsStdCutDown  = p.get<double>("ResidualsStdCutDown");
   fResidualsMeanCutDown = p.get<double>("ResidualsMeanCutDown");
   fResidualsMeanCutUp   = p.get<double>("ResidualsMeanCutUp");
   fPercUsedHitsCut      = p.get<double>("PercUsedHitsCut");

   std::cout << "[CC1pi::InputTags] Using configuration: " << std::endl
      << "\t <<<< fTrackLabel = " << fTrackLabel << std::endl
      << "\t <<<< fCalorimetryLabel = " << fCalorimetryLabel << std::endl
      << "\t <<<< fSelectionLabel = " << fSelectionLabel << std::endl
      << "\t <<<< fPIDLabelChi2 = " << fPIDLabelChi2 << std::endl
      << "\t <<<< fPIDProducer = " << fPIDProducer << std::endl
      << "\t <<<< fPFParticleProducer = " << fPFParticleProducer  << std::endl
      << "\t <<<< fClusterProducer = " << fClusterProducer  << std::endl
      << "\t <<<< fSpacePointProducer = " << fSpacePointProducer << std::endl
      << "\t <<<< fTPCObjectProducer = " << fTPCObjectProducer << std::endl
      << "\t <<<< fResidualsStdCutUp = " << fResidualsStdCutUp << std::endl
      << "\t <<<< fResidualsStdCutDown = " << fResidualsStdCutDown << std::endl
      << "\t <<<< fResidualsMeanCutUp = " << fResidualsMeanCutUp << std::endl
      << "\t <<<< fResidualsMeanCutDown = " << fResidualsMeanCutDown << std::endl
      << "\t <<<< fPercUsedHitsCut = " << fPercUsedHitsCut << std::endl << std::endl;
}
