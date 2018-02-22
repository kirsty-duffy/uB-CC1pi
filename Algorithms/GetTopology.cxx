#ifndef GETTOPOLOGY_CXX
#define GETTOPOLOGY_CXX

#include "GetTopology.h"


NuIntTopology GetTopology(art::ValidHandle<std::vector<simb::MCTruth>> mc_truths)
{
  // Check whether MCtruth exists
  if (mc_truths->size() == 0){
    std::cout << "[ERROR: GetTopology] mc_truths->size() = " << mc_truths->size() << std::endl;
    return kUnknown;
  }
  
  // Check whether there is more than one MCtruth - if there is this will be wrong so issue an error!
  if (mc_truths->size() > 1){
    std::cout << "[ERROR: GetTopology] mc_truths->size() = " << mc_truths->size() << std::endl;
    return kUnknown;
  }

  // Otherwise, get first MCtruth object
  simb::MCTruth mctruth = mc_truths->at(0);

  // If mctruth origin is not beam neutrino, return unknown
  // (I don't think this should happen if we're looking at genie mctruth, but just in case)
  if (mctruth.Origin() != simb::kBeamNeutrino){
    std::cout << "[ERROR: GetTopology] mctruth.Origin() = " << mctruth.Origin() << std::endl;
    return kUnknown;
  }

  // Loop through particles in MCtruth and count pions and protons
  // I have been reliably informed (by Andy Furmanski) that the particles defined here should
  // only be true neutrino daughters
  int npiplus = 0;
  int npiminus = 0;
  int npizero = 0;
  int nprotons = 0;
  
  int npart = mctruth.NParticles();
  for (int ipart=0; ipart < npart; ipart++){
    simb::MCParticle particle = mctruth.GetParticle(ipart);
    if (particle.StatusCode() != 1){continue;} // only look at detectable particles (that make it out of the nucleus)
    
    if (particle.PdgCode()      ==  211){npiplus++;}
    else if (particle.PdgCode() == -211){npiminus++;}
    else if (particle.PdgCode() ==  111){npizero++;}
    else if (particle.PdgCode() == 2212){nprotons++;}
  }// loop over particles in MCtruth (ipart)

  // Now work out if we're looking at CC or NC and get topology
  NuIntTopology inttype = kUnknown;
  
  simb::MCNeutrino neutrino = mctruth.GetNeutrino();
  if (neutrino.CCNC() == simb::kCC && abs(neutrino.Nu().PdgCode())==14){ // If numu CC
    if (npiplus+npiminus+npizero==0){ // numu CC0pi
      switch (nprotons){
      case 0:
	inttype = kCC0pi0p;
	break;
      case 1:
	inttype = kCC0pi1p;
	break;
      default:
	inttype = kCC0piNp;
      }
    }
    else if (npiplus+npiminus+npizero==1){ // numu CC1pi
      if (npiplus == 1){ // numu CCpiplus
	switch (nprotons){
	case 0:
	  inttype = kCC1piplus0p;
	  break;
	case 1:
	  inttype = kCC1piplus1p;
	  break;
	default:
	  inttype = kCC1piplusNp;
	}
      }
      else if (npiminus == 1){ // numu CCpiminus
	switch (nprotons){
	case 0:
	  inttype = kCC1piminus0p;
	  break;
	case 1:
	  inttype = kCC1piminus1p;
	  break;
	default:
	  inttype = kCC1piminusNp;
	}
      }
      else if (npizero == 1){ // CCpizero
	switch (nprotons){
	case 0:
	  inttype = kCC1pizero0p;
	  break;
	case 1:
	  inttype = kCC1pizero1p;
	  break;
	default:
	  inttype = kCC1pizeroNp;
	}
      }
    }
    else if (npiplus+npiminus+npizero > 1){
      switch (nprotons){
      case 0:
	inttype = kCCmultipi0p;
	break;
      case 1:
	inttype = kCCmultipi1p;
	break;
      default:
	inttype = kCCmultipiNp;
      }
    }
    else{
      inttype = kCCother;
    }
  } // end if numu CC
  else if (neutrino.CCNC() == simb::kCC && abs(neutrino.Nu().PdgCode())==12){ // if nue CC
    inttype = kCCNue;
  } // end if nue CC
  else if (neutrino.CCNC() == simb::kNC){
    inttype = kNC;
  } // end if NC
  else{
     inttype = kUnknown;
  }

  // Check it's doing sensible things
  /*if (neutrino.CCNC() == simb::kCC){std::cout << "CC interaction. ";}
  if (neutrino.CCNC() == simb::kNC){std::cout << "NC interaction. ";}
  std::cout << "Mode is " << neutrino.Mode() << ", topology is " << topologyenum2str(inttype) << std::endl;*/
  
  return inttype;
}

#endif
