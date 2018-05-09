// Class to hold Theoretical predictions for dEdx vs residual range
// For protons, kaons, pions, and muons
// Taken from Bruce Baller's spread sheet: MicroBooNE DocDB-6572

// Converted to TGraphs, where x is Res Range(cm), and y is dE/dx(MeV/cm)

// To use: include this header file, and then do something like this
// particleid::Theory_dEdx_resrange blah;
// blah.g_ThdEdxRR_Proton->Draw("ALP");

// To get the predicted dEdx for any given length, use
// blah.g_ThdEdxRR_Proton->Eval(length,0,"S");
// (this will create an internal TSpline3 object to properly interpolate between the
// nearest two points on the TGraph and give you the most accurate result)

#ifndef THEORY_DEDX_RESRANGE_H
#define THEORY_DEDX_RESRANGE_H

#include "TGraph.h"

namespace particleid{

  class Theory_dEdx_resrange{

  public:
    Theory_dEdx_resrange();
    int GetNPoints(){return npoints;}
    
    TGraph *g_ThdEdxRR_Proton;
    TGraph *g_ThdEdxRR_Kaon;
    TGraph *g_ThdEdxRR_Pion;
    TGraph *g_ThdEdxRR_Muon;
    TGraph *g_ThdEdxRR_MuonNoBragg;

  protected:
    const int npoints = 107;
    
  };

}


#endif
