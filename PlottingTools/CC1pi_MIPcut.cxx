#ifndef __CC1PI_MIPCUTS_CXX__
#define __CC1PI_MIPCUTS_CXX__

// #include <vector>

#include "CC1pi_plotvars_def.h"
#include "CC1pi_cuts.h"

// -------------------------------------------------- //

// Definition of MIP cut
double doEvalMIPCut(std::vector<CC1piPlotVars> MIPCutVars, treevars *vars, int i_tr){

   double passed = 1;

   // Loop over MIP cuts and apply each one to the track
   for (size_t i_mipcut=0; i_mipcut < MIPCutVars.size(); i_mipcut++){

      double value = MIPCutVars.at(i_mipcut).Var->at(i_tr);
      bool MIPlow = MIPCutVars.at(i_mipcut).KeepBelowCut;
      double cutval = MIPCutVars.at(i_mipcut).CutValue;

      // std::cout << " MIP cut " << i_mipcut << " cutval " << cutval << "(KeepBelowCut = " << MIPlow << ")" ", value = " << value << std::endl;

      if (value == -9999 || value == -999) {
         return -9999;
      }
      else if (!std::isnormal(value) && value!=0){
         // std::cout << "[WARNING :: CC1pi_MIPcut.cxx] not-normal value for " << MIPCutVars.at(i_mipcut).histtitle << std::endl;
         return -9999;
      }
      else if (MIPlow && (value > cutval)){
         // std::cout << "value > cutval so fails" << std::endl;
         passed = 0;
      }
      else if (!MIPlow && (value < cutval)){
         // std::cout << "value < cutval so fails" << std::endl;
         passed = 0;
      }
   }

   // If it makes it all the way to the end of the loop, it must have passed every cut. Return 1.
   return passed;
};

// -------------------------------------------------- //

// Definition of MIP cut
double EvalMIPCut(treevars *vars, int i_tr, std::vector<CC1piPlotVars> *MIPCutVars){
   std::vector<CC1piPlotVars> MIPCutVars_tmp;
   // Get MIP cuts
   if (MIPCutVars == nullptr){
      // std::cout << "Getting MIP Cut vars from scratch" << std::endl;
      MIPCutVars_tmp = GetMIPCutVars(vars);
   }
   else{
      MIPCutVars_tmp = (*MIPCutVars);
   }
   // std::cout << "Got MIP Cut vars" << std::endl;

   return doEvalMIPCut(MIPCutVars_tmp, vars, i_tr);
};

#endif
