#ifndef PROTONDEDXRRPLOTS_H
#define PROTONDEDXRRPLOTS_H

#include "CC1pi_treevars.h"
#include "CC1pi_treevars.cxx"

void MakeProtondEdxrrPlot_SingleTrack(TCanvas *c1, treevars *vars, int trackIndex){

  // Get the things you want out of vars
  std::vector<double> dEdx = vars->TPCObj_PFP_track_dedx_perhit->at(trackIndex).at(2); // collection plane only
  std::vector<double> resrg = vars->TPCObj_PFP_track_resrange_perhit->at(trackIndex).at(2); // collection plane only

  // Now make the TGraphs we want to fill...
  int npoints = dEdx.size();
  if (npoints==0) return;

  // ... TGraph for hit time vs hit wire
  TGraph *g = new TGraph(npoints);
  g->SetTitle(std::string(std::string("<dE/dx>tr at start = ")+std::to_string(vars->TPCObj_PFP_track_dEdx_truncmean_start->at(trackIndex))+std::string(";Residual range [cm];dE/dx")).c_str());

  for (size_t i_hit=0; i_hit<npoints; i_hit++){
    g->SetPoint(i_hit,resrg.at(i_hit),dEdx.at(i_hit));
  }

  // Draw plots and style on TCanvas
  g->SetMarkerStyle(kPlus);
  g->Draw("AP");
}

#endif
