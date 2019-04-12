#ifndef __CC1PI_PIONCUT_CXX__
#define __CC1PI_PIONCUT_CXX__

// #include <vector>

#include "CC1pi_plotvars_def.h"
#include "CC1pi_cuts.h"

// -------------------------------------------------- //

// Definition of MIP cut
double EvalPionCut(treevars *vars, int i_tr){

   double isPion=1;

   // std::cout << "isdaughter = " << vars->TPCObj_PFP_isDaughter->at(i_tr) << std::endl;
   // std::cout << "ismip = " << vars->TPCObj_PFP_track_passesMIPcut->at(i_tr) << std::endl;
   // std::cout << "iscontained = " << vars->TPCObj_PFP_track_isContained->at(i_tr) << std::endl;

   // Must be a reconstructed daughter of the Neutrino
   double isdaughter = vars->TPCObj_PFP_isDaughter->at(i_tr);
   if (isdaughter == -9999 || isdaughter == -999){
     return -9999;
   }
   else if (isdaughter==0){
      return 0;
   }

   // Must be selected as a MIP
   double ismip = vars->TPCObj_PFP_track_passesMIPcut->at(i_tr);
   if (ismip == -9999 || ismip == -999){
     return -9999;
   }
   else if (ismip==0){
     return 0;
   }

   // A muon candidate must have been selected
   if (vars->TPCObj_MuCandtrackIndex < 0){
      return 0;
   }

   // Must be pion candidate
   if (vars->TPCObj_PiCandtrackIndex != i_tr){
      return 0;
   }

   // Mu-pi BDT cut
   if (vars->TPCObj_PFP_track_isContained->at(i_tr)){
      // Contained BDT cut

      // double mupibdt = vars->TPCObj_PFP_track_mupiBDTscore_cont->at(i_tr);
      //
      // // std::cout << "mupibdt = " << mupibdt << std::endl;
      //
      // if (mupibdt == -9999 || mupibdt == -999){
      //   return -9999;
      // }
      // else if (mupibdt <= 0.2){
      //    return 0;
      // }
   }
   else{
      // Exiting BDT cut
      double mupibdt = vars->TPCObj_PFP_track_mupiBDTscore_exit->at(i_tr);

      // std::cout << "mupibdt = " << mupibdt << std::endl;

      if (mupibdt == -9999 || mupibdt == -999){
        return -9999;
      }
      else if (mupibdt <= 0.8){
         return 0;
      }
   }

   return isPion;
};

#endif
