#ifndef __CC1PI_MIPCUTS_CXX__
#define __CC1PI_MIPCUTS_CXX__

// #include <vector>

#include "CC1pi_plotvars_def.h"
#include "CC1pi_cuts.h"

// -------------------------------------------------- //

// Definition of MIP cut
double EvalMIPCut(treevars *vars, int i_tr){

   double isMIP=1;

   // Get MIP cuts
   std::vector<CC1piPlotVars> MIPCutVars = GetMIPCutVars(vars);

   // Loop over MIP cuts and apply each one to the track
   for (size_t i_mipcut=0; i_mipcut < MIPCutVars.size(); i_mipcut++){
      double value = MIPCutVars.at(i_mipcut).Var->at(i_tr);
      double cutval = MIPCutVars.at(i_mipcut).CutValue;
      bool MIPlow = MIPCutVars.at(i_mipcut).KeepBelowCut;

      if (value == -9999 || value == -999) {
         isMIP=-9999;
         break;
      }
      else if (MIPlow && (value > cutval)){
         isMIP=0;
         break;
      }
      else if (!MIPlow && (value < cutval)){
         isMIP=0;
         break;
      }
   }
   return isMIP;
};

#endif
