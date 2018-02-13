#ifndef PDGENUMS
#define PDGENUMS

enum PDGCode{
    kNuMu        = 14,
    kNuMuBar     = -14,
    kNuE         = 12,
    kNuEBar      = -12,
    kNuTau       = 16,
    kNuTauBar    = -16,
    kMuMinus     = 13,
    kMuPlus      = -13,
    kPiMinus     = -211,
    kPiPlus      = 211,
    kPiZero      = 111,
    kElectron    = 11,
    kPositron    = -11,
    kTauMinus    = 15,
    kTauPlus     = -15,
    kPhoton      = 22,
    kProton      = 2212,
    kNeutron     = 2112,
    kArgon       = 1000180400,
    kBindino     = 2000000101,
    kPDGUnknown  = -9999
};

inline std::string PDGenum2str(PDGCode type)
{
  std::string returnString = "";

  switch(type){
    case kNuMu:
      returnString = "#nu_{#mu}";
      break;
    case kNuMuBar:
      returnString = "#bar{#nu}_{#mu}";
      break;
    case kNuE:
      returnString = "#nu_{e}";
      break;
    case kNuEBar:
      returnString = "#bar{#nu}_{e}";
      break;
    case kNuTau:
      returnString = "#nu_{#tau}";
      break;
    case kNuTauBar:
      returnString = "#bar{#nu}_{#tau}";
      break;
    case kMuMinus:
      returnString = "#mu^{-}";
      break;
    case kMuPlus:
      returnString = "#mu^{+}";
      break;
    case kPiMinus:
      returnString = "#pi^{-}";
      break;
    case kPiPlus:
      returnString = "#pi^{+}";
      break;
    case kPiZero:
      returnString = "#pi^{0}";
      break;
    case kElectron:
      returnString = "e^{-}";
      break;
    case kPositron:
      returnString = "e^{+}";
      break;
    case kTauMinus:
      returnString = "#tau^{-}";
      break;
    case kTauPlus:
      returnString = "#tau^{+}";
      break;
    case kPhoton:
      returnString = "#gamma";
      break;
    case kProton:
      returnString = "p";
      break;
    case kNeutron:
      returnString = "n";
      break;
    case kArgon:
      returnString = "Ar-40";
      break;
  case kBindino:
      returnString = "Bindino";
      break;
  default:
    std::cout << "[ERROR: PDGEnums.h] Could not find string conversion for " << type << std::endl;
    returnString = "Unknown";
    break;
  }

  return returnString;
}


#endif


