#ifndef TOPOLOGYENUMS
#define TOPOLOGYENUMS

#include <string>
#include <iostream>

// Note: kCosmic, kMixed, and kOutFV will never be set by GetTopology.cxx
// They will be kUnknown and you must check for them yourself
enum NuIntTopology{
  kCC0pi0p,
  kCC0pi1p,
  kCC0piNp,
  kCC1piplus0p,
  kCC1piplus1p,
  kCC1piplusNp,
  kCC1piminus0p,
  kCC1piminus1p,
  kCC1piminusNp,
  kCC1pizero0p,
  kCC1pizero1p,
  kCC1pizeroNp,
  kCCmultipi0p,
  kCCmultipi1p,
  kCCmultipiNp,
  kCCother,
  kNC,
  kCosmic,
  kMixed,
  kOutFV,
  kUnknown
};

inline std::string topologyenum2str(NuIntTopology topology)
{
  std::string returnString = "";
  
  switch(topology){
  case kCC0pi0p:
    returnString = "CC 0#pi 0p";
    break;
  case kCC0pi1p:
    returnString = "CC 0#pi 1p";
    break;
  case kCC0piNp:
    returnString = "CC 0#pi Np (N>1)";
    break;
  case kCC1piplus0p:
    returnString = "CC 1#pi^{+} 0p";
    break;
  case kCC1piplus1p:
    returnString = "CC 1#pi^{+} 1p";
    break;
  case kCC1piplusNp:
    returnString = "CC 1#pi^{+} Np (N>1)";
    break;
  case kCC1piminus0p:
    returnString = "CC 1#pi^{-} 0p";
    break;
  case kCC1piminus1p:
    returnString = "CC 1#pi^{-} 1p";
    break;
  case kCC1piminusNp:
    returnString = "CC 1#pi^{-} Np (N>1)";
    break;
  case kCC1pizero0p:
    returnString = "CC 1pi^{0} 0p";
    break;
  case kCC1pizero1p:
    returnString = "CC 1pi^{0} 1p";
    break;
  case kCC1pizeroNp:
    returnString = "CC 1pi^{0} Np (N>1)";
    break;
  case kCCmultipi0p:
    returnString = "CC multi#pi 0p";
    break;
  case kCCmultipi1p:
    returnString = "CC multi#pi 1p";
    break;
  case kCCmultipiNp:
    returnString = "CC multi#pi Np (N>1)";
    break;
  case kCCother:
    returnString = "CC other";
    break;
  case kNC:
    returnString = "NC";
    break;
  case kCosmic:
    returnString = "Cosmic";
    break;
  case kMixed:
    returnString = "Mixed";
    break;
  case kOutFV:
    returnString = "Out of FV";
    break;
  default:
    std::cout << "[ERROR: TopologyEnums.h] Could not find string conversion for " << topology << std::endl;
    returnString = "Unknown";
    break;
  }

  return returnString;
}


inline std::string topologyenum2str_coarse(NuIntTopology topology)
{
  std::string returnString = "";
  
  switch(topology){
  case kCC0pi0p:
  case kCC0pi1p:
  case kCC0piNp:
    returnString = "CC 0#pi";
    break;
  case kCC1piplus0p:
  case kCC1piplus1p:
  case kCC1piplusNp:
    returnString = "CC 1#pi^{+}";
    break;
  case kCC1piminus0p:
  case kCC1piminus1p:
  case kCC1piminusNp:
    returnString = "CC 1#pi^{-}";
    break;
  case kCC1pizero0p:
  case kCC1pizero1p:
  case kCC1pizeroNp:
    returnString = "CC 1pi^{0}";
    break;
  case kCCmultipi0p:
  case kCCmultipi1p:
  case kCCmultipiNp:
    returnString = "CC multi#pi";
    break;
  case kCCother:
    returnString = "CC other";
    break;
  case kNC:
    returnString = "NC";
    break;
  case kCosmic:
    returnString = "Cosmic";
    break;
  case kMixed:
    returnString = "Mixed";
    break;
  case kOutFV:
    returnString = "Out of FV";
    break;
  default:
    std::cout << "[ERROR: TopologyEnums.h] Could not find string conversion for " << topology << std::endl;
    returnString = "Unknown";
    break;
  }

  return returnString;
}

#endif
