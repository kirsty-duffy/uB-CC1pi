#ifndef __CC1PI_PIONCUT_CXX__
#define __CC1PI_PIONCUT_CXX__

// #include <vector>

#include "CC1pi_plotvars_def.h"
#include "CC1pi_cuts.h"

// -------------------------------------------------- //

// Definition of MIP cut
double EvalPionCut(treevars *vars, int i_tr){

   double isPion=1;

   // Get Pion cuts
   std::vector<CC1piPlotVars> PionCutVars = GetPionCutVars(vars);

   // Loop over Pion cuts and apply each one to the track
   for (size_t i_cut=0; i_cut < PionCutVars.size(); i_cut++){
      double value = PionCutVars.at(i_cut).Var->at(i_tr);
      double cutval = PionCutVars.at(i_cut).CutValue;
      bool Pionlow = PionCutVars.at(i_cut).KeepBelowCut;

      if (value == -9999 || value == -999) {
         isPion=-9999;
         break;
      }
      else if (Pionlow && (value > cutval)){
         isPion=0;
         break;
      }
      else if (!Pionlow && (value < cutval)){
         isPion=0;
         break;
      }
   }
   return isPion;
};

#endif
