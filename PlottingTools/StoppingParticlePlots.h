#ifndef STOPPINGPARTICLEPLOTS_H
#define STOPPINGPARTICLEPLOTS_H

#include "CC1pi_treevars.h"
#include "CC1pi_treevars.cxx"
#include "../Algorithms/PDGEnums.h"

void MakeStoppingParticlePlots_SingleTrack(TCanvas *c1, treevars *vars, int trackIndex){

  // Get the things you want out of vars
  std::vector<double> track_SimpleCluster_hitTime = vars->TPCObj_PFP_track_SimpleCluster_hitTime->at(trackIndex);
  std::vector<int> track_SimpleCluster_hitWire = vars->TPCObj_PFP_track_SimpleCluster_hitWire->at(trackIndex);
  // std::vector<double> track_SimpleCluster_hitTimeTicks = vars->TPCObj_PFP_track_SimpleCluster_hitTimeTicks->at(trackIndex);
  // std::vector<int> track_SimpleCluster_hitWireNo = vars->TPCObj_PFP_track_SimpleCluster_hitWireNo->at(trackIndex);
  int track_SimpleCluster_StartIndex = vars->TPCObj_PFP_track_SimpleCluster_StartIndex->at(trackIndex);
  std::vector<double> track_SimpleCluster_hitdQds = vars->TPCObj_PFP_track_SimpleCluster_hitdQds->at(trackIndex);
  std::vector<double> track_SimpleCluster_hitdQdsSlider = vars->TPCObj_PFP_track_SimpleCluster_hitdQdsSlider->at(trackIndex);
  std::vector<double> track_SimpleCluster_hitLinearity = vars->TPCObj_PFP_track_SimpleCluster_hitLinearity->at(trackIndex);
  bool track_ct_passed_basic = vars->TPCObj_PFP_track_ct_passed_basic->at(trackIndex);
  bool track_ct_result_bragg = vars->TPCObj_PFP_track_ct_result_bragg->at(trackIndex);
  bool track_ct_result_michel = vars->TPCObj_PFP_track_ct_result_michel->at(trackIndex);

  // Now make the TGraphs we want to fill...
  int npoints = track_SimpleCluster_hitTime.size();
  // std::cout << "npoints = " << npoints << ", StartIndex = " << track_SimpleCluster_StartIndex << std::endl;

  // ... TGraph for hit time vs hit wire
  TGraph *g_timevswire_sm = new TGraph(npoints);
  g_timevswire_sm->SetTitle(";hit Wire (cm);hit Time (cm)");
  g_timevswire_sm->SetLineWidth(2);

  // ... TGraph for hit time vs hit wire for first hit
  TGraph *g_timevswire_Start = new TGraph(1);
  g_timevswire_Start->SetPoint(0,track_SimpleCluster_hitWire.at(0),track_SimpleCluster_hitTime.at(0));
  g_timevswire_Start->SetMarkerColor(kRed);
  g_timevswire_Start->SetMarkerStyle(29);

  // ... TGraph for hit dQds vs hit number in ordered vector
  TGraph *g_dQdsvshitno_raw = new TGraph(npoints);
  g_dQdsvshitno_raw->SetMarkerColor(kGray+2);
  g_dQdsvshitno_raw->SetLineColor(kGray);
  g_dQdsvshitno_raw->SetLineWidth(2);
  g_dQdsvshitno_raw->SetTitle(";hit no.;hit dQ/ds (e^{-}/cm)");

  // ... TGraph for smoothed hit dQds vs hit number in ordered vector (to overlay and just check the smoothing)
  TGraph *g_dQdsvshitno_sm = new TGraph(npoints);
  g_dQdsvshitno_sm->SetTitle(";hit no.;hit dQ/ds (e^{-}/cm)");

  // ... TGraph for hit linearity vs hit number in ordered vector
  TGraph *g_linearityvshitno = new TGraph(npoints);
  g_linearityvshitno->SetTitle(";hit no.;Local linearity at hit");

  for (size_t i_hit=0; i_hit<npoints; i_hit++){
    g_timevswire_sm->SetPoint(i_hit,track_SimpleCluster_hitWire.at(i_hit),track_SimpleCluster_hitTime.at(i_hit));
    g_dQdsvshitno_raw->SetPoint(i_hit,i_hit,track_SimpleCluster_hitdQds.at(i_hit));
    g_dQdsvshitno_sm->SetPoint(i_hit,i_hit,track_SimpleCluster_hitdQdsSlider.at(i_hit));
    g_linearityvshitno->SetPoint(i_hit,i_hit,track_SimpleCluster_hitLinearity.at(i_hit));
  }

  // Draw plots and style on TCanvas
  c1->Divide(1,4);

  c1->cd(1);
  TPaveText *pt = new TPaveText(0.05,0.51,0.95,0.95,"NB");
  pt->SetLineColor(kWhite);
  pt->SetFillStyle(0);
  pt->SetTextFont(132);
  pt->AddText(std::string(std::string("True PDG: ")+PDGenum2str((PDGCode)vars->TPCObj_PFP_truePDG->at(trackIndex))).c_str());
  pt->AddText(TString::Format("Classed as MIP? %d",(int)vars->TPCObj_PFP_track_passesMIPcut->at(trackIndex)).Data());
  pt->AddText(TString::Format("Passes basic cosmic tagging? %d",track_ct_passed_basic).Data());
  pt->AddText(TString::Format("Identified Michel? %d",track_ct_result_michel).Data());
  pt->AddText(TString::Format("Identified Bragg (no Michel)? %d",track_ct_result_bragg).Data());
  pt->Draw();

  TLegend *l = new TLegend(0.05,0.05,0.95,0.49);
  l->SetNColumns(3);
  l->SetLineColor(kWhite);
  l->SetFillStyle(0);
  l->AddEntry(g_timevswire_Start,"Start hit identified by algorithm","p");
  l->AddEntry(g_timevswire_sm,"Smoothed result","lp");
  l->AddEntry(g_dQdsvshitno_raw,"Raw result (dQ/ds vs hit no. only)","lp");
  l->Draw();

  c1->cd(2);
  g_timevswire_sm->Draw("AP");
  g_timevswire_Start->Draw("same P");

  c1->cd(3);
  g_dQdsvshitno_raw->Draw("AP");
  g_dQdsvshitno_sm->Draw("same P");

  c1->cd(4);
  g_linearityvshitno->Draw("AP");
}

// --------------------------------------------------- //

void MakeStoppingParticlePlots_AllTracks(std::string mcfile){

   gStyle->SetTitleX(0.5);
   gStyle->SetTitleAlign(23);
   gStyle->SetOptTitle(1);
   gStyle->SetTitleBorderSize(0.);

   TFile *f_bnbcos = new TFile(mcfile.c_str(), "read");
   TTree *t_bnbcos = (TTree*)f_bnbcos->Get("cc1piselec/outtree");
   TTree *t_bnbcos_friend = (TTree*)f_bnbcos->Get("StoppingParticleTagger/StoppingTaggerTree");

   if (t_bnbcos->GetEntries() != t_bnbcos_friend->GetEntries()){
     std::cout << "ERROR: cc1piselec/outtree has " << t_bnbcos->GetEntries() << " entries, and StoppingParticleTagger/StoppingTaggerTree has " << t_bnbcos_friend->GetEntries() << " entries. Cannot make friend tree, exiting." << std::endl;
     return;
   }
   t_bnbcos->AddFriend(t_bnbcos_friend);

   treevars mc_vars;
   settreevars(t_bnbcos,&mc_vars);

   // Loop through the tree and make plots. Each entry in the tree represents a different TPCObject, so loop through these and make one plot per track, so use this with care!
   for (size_t i=0; i<t_bnbcos->GetEntries(); i++){
     t_bnbcos->GetEntry(i);
     for (size_t i_tr=0; i_tr < mc_vars.TPCObj_PFP_track_SimpleCluster_StartIndex->size(); i_tr++){
       if (mc_vars.TPCObj_PFP_track_length->at(i_tr)<0){
         std::cout << "Skipping track " << i_tr << " in event " << i << std::endl;
         continue;
       }
       TCanvas *c1 = new TCanvas("","",400,500);

       MakeStoppingParticlePlots_SingleTrack(c1, &mc_vars, (int)i_tr);

       c1->Print(std::string(std::string("StoppingParticlePlots_TPCObj")+i+std::string("_track")+i_tr+std::string(".pdf")).c_str());
     } // end loop over tracks in TPCObject
   } // end loop over TPCObjects (entries in tree)
}

#endif
