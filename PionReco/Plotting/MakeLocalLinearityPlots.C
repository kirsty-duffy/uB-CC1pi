#include "TCanvas.h"
#include "TTree.h"
#include "TFile.h"
#include "TGraph2D.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TString.h"
#include "TLegend.h"
#include "TPaveText.h"
#include "TGraph.h"


void PlotLocalLinearityDetails(std::vector<std::vector<double>> spacepoints_vec, std::vector<double> proj1d_vec, std::vector<double> cusum_vec, std::vector<double> ttest_vec, std::vector<int> kinkidxs_vec, int pfp_idx, std::vector<double> startpoint_vec, std::vector<std::vector<double>> truekink_vec=std::vector<double>(), double mincosth=-999, int evtnum=-999){

  // Make TGraphs to fill
  TGraph *spacepoints_yz = new TGraph();
  spacepoints_yz->SetMarkerColor(kViolet);
  spacepoints_yz->SetMarkerStyle(7);
  TGraph *spacepoints_xz = new TGraph();
  spacepoints_xz->SetMarkerColor(kViolet);
  spacepoints_xz->SetMarkerStyle(7);

  TGraph *startpoint_yz = new TGraph(1);
  startpoint_yz->SetMarkerColor(kRed);
  startpoint_yz->SetMarkerStyle(29);
  TGraph *startpoint_xz = new TGraph(1);
  startpoint_xz->SetMarkerColor(kRed);
  startpoint_xz->SetMarkerStyle(29);

  TGraph *ttest = new TGraph();

  TGraph *breakpoints_yz = new TGraph();
  breakpoints_yz->SetMarkerStyle(7);
  breakpoints_yz->SetMarkerColor(kOrange-3);
  TGraph *breakpoints_xz = new TGraph();
  breakpoints_xz->SetMarkerStyle(7);
  breakpoints_xz->SetMarkerColor(kOrange-3);

  TGraph *breakpoints_idx = new TGraph();
  breakpoints_idx->SetMarkerStyle(7);
  breakpoints_idx->SetMarkerColor(kOrange-3);

  TGraph *truekink_yz = new TGraph();
  truekink_yz->SetMarkerColor(kBlue);
  truekink_yz->SetMarkerStyle(29);
  TGraph *truekink_xz = new TGraph();
  truekink_xz->SetMarkerColor(kBlue);
  truekink_xz->SetMarkerStyle(29);

  TGraph *dirs = new TGraph();

  TGraph *cusum = new TGraph();

  // Reset canvas
  TCanvas *c1 = new TCanvas("c1","c1",500,700);
  c1->Clear();
  c1->Divide(1,5);

  // Fill spacepoints plots
  double minx = 9999;
  double maxx = -9999;
  double miny = 9999;
  double maxy = -9999;
  double minz = 9999;
  double maxz = -9999;

  // Draw spacepoints
  for (size_t i_sp=0; i_sp<spacepoints_vec.size(); i_sp++){
    double x = spacepoints_vec.at(i_sp).at(0);
    double y = spacepoints_vec.at(i_sp).at(1);
    double z = spacepoints_vec.at(i_sp).at(2);

    spacepoints_yz->SetPoint(i_sp,z,y);
    spacepoints_xz->SetPoint(i_sp,z,x);

     if (x<minx) minx = x;
     if (x>maxx) maxx = x;
     if (y<miny) miny = y;
     if (y>maxy) maxy = y;
     if (z<minz) minz = z;
     if (z>maxz) maxz = z;
  } // end loop over spacepoints

  startpoint_yz->SetPoint(0,startpoint_vec.at(2),startpoint_vec.at(1));
  startpoint_xz->SetPoint(0,startpoint_vec.at(2),startpoint_vec.at(0));

  // Fill local linearity plots
  for (size_t i_pt=0; i_pt < ttest_vec.size(); i_pt++){

    cusum->SetPoint(i_pt,i_pt,cusum_vec.at(i_pt));

    ttest->SetPoint(i_pt,i_pt,ttest_vec.at(i_pt));

    dirs->SetPoint(i_pt,i_pt,proj1d_vec.at(i_pt));
  }

  for (int i_kink=0; i_kink<kinkidxs_vec.size(); i_kink++){
    int kink_idx = kinkidxs_vec.at(i_kink);
    double x = spacepoints_vec.at(kink_idx).at(0);
    double y = spacepoints_vec.at(kink_idx).at(1);
    double z = spacepoints_vec.at(kink_idx).at(2);

    breakpoints_yz->SetPoint(i_kink,z,y);
    breakpoints_xz->SetPoint(i_kink,z,x);
    breakpoints_idx->SetPoint(i_kink,kink_idx,0);
  }

  for (int i_true=0; i_true<truekink_vec.size(); i_true++){
    double x = truekink_vec.at(i_true).at(0);
    double y = truekink_vec.at(i_true).at(1);
    double z = truekink_vec.at(i_true).at(2);
    truekink_yz->SetPoint(i_kink,z,y);
    truekink_yx->SetPoint(i_kink,z,x);
  }




  // Draw plot
  // Have to draw first TGraph, then set axis ranges, then draw again because it's the only way ROOT will do what I want (the TAxis objects don't exist until you draw the TGraph so you can't set ranges)
  c1->cd(2);
  gPad->Divide(2);
  if (spacepoints_yz->GetN()>0){
    c1->cd(2); gPad->cd(1);
    spacepoints_yz->Draw("Ap");
    spacepoints_yz->GetYaxis()->SetTitle("Y (cm)");
    spacepoints_yz->GetXaxis()->SetTitle("Z (cm)");
    if (truekink_yz->GetN()>0) truekink_yz->Draw("same p");
    spacepoints_yz->Draw("same p"); // draw again to make sure spacepoints are on the top
    startpoint_yz->Draw("same p");
    if (breakpoints_yz->GetN()>0) breakpoints_yz->Draw("same p");
    c1->cd(2); gPad->cd(2);
    spacepoints_xz->GetYaxis()->SetTitle("X (cm)");
    spacepoints_xz->GetXaxis()->SetTitle("Z (cm)");
    spacepoints_xz->Draw("Ap");
    if (truekink_yz->GetN()>0) truekink_xz->Draw("same p");
    spacepoints_xz->Draw("same p"); // draw again to make sure spacepoints are on the top
    startpoint_xz->Draw("same p");
    if (breakpoints_xz->GetN()>0) breakpoints_xz->Draw("same p");
  }
  else{
    //return false;
  }

  c1->cd(3);
  dirs->GetXaxis()->SetTitle("Spacepoint index in PFP (coloured PFP only)");
  dirs->GetYaxis()->SetTitle("Direction relative to track start");
  dirs->Draw("APL");
  if (breakpoints_idx->GetN()>0) breakpoints_idx->Draw("P same");

  // Draw local linearity plot
  c1->cd(4);
  ttest->GetXaxis()->SetTitle("Spacepoint index in PFP (coloured PFP only)");
  ttest->GetYaxis()->SetTitle("Welch's t-test statistic");
  // ttest->GetYaxis()->SetRangeUser(0,1.05);
  ttest->Draw("APL");
  if (breakpoints_idx->GetN()>0) breakpoints_idx->Draw("P same");

  c1->cd(5);
  cusum->GetXaxis()->SetTitle("Spacepoint index in PFP (coloured PFP only)");
  cusum->GetYaxis()->SetTitle("Cumulative sum of differences between values and average");
  cusum->Draw("APL");
  if (breakpoints_idx->GetN()>0) breakpoints_idx->Draw("P same");

  // Make legend
  c1->cd(1);
  TPaveText *pt = new TPaveText(0.05,0.41,0.95,0.95,"NB");
  pt->SetLineColor(kWhite);
  pt->SetFillStyle(0);
  pt->SetTextFont(132);
  pt->AddText(TString::Format("PFP %d",pfp_idx).Data());
  pt->AddText(TString::Format("%d reconstructed kinks found",(int)kinkidxs_vec.size()).Data());
  pt->AddText(TString::Format("Smallest true cos(theta): %.2f",mincosth).Data());
  pt->Draw();

  TLegend *l = new TLegend(0.05,0.05,0.95,0.39);
  l->SetNColumns(3);
  l->SetLineColor(kWhite);
  l->SetFillStyle(0);
  l->SetTextFont(132);
  // l->AddEntry(nuvtx_yz,"Reco #nu vertex","p");
  l->AddEntry(startpoint_yz,"Linearity start point","p");
  l->AddEntry(spacepoints_yz,"Studied PFP","p");
  l->AddEntry(breakpoints_yz,"Reconstructed kink(s)","p");
  l->AddEntry(truekink_yz,"True kink(s)","p");
  l->Draw();

  c1->Update();
  c1->Draw();

  c1->Print(TString::Format("mincosth_%.2f_%dkinks_evt_%d_pfp_%d.pdf",mincosth,(int)kinkidxs_vec.size(),evtnum,pfp_idx).Data());

  c1->Delete();

  // return true;

} // end function: PlotLocalLinearityDetails


// void MakeLocalLinearityPlots(std::string inputfile="/uboone/app/users/kduffy/CC1pi/CC1pi_uboonecode/srcs/uboonecode/uboone/CC1pi/ignore/2018-11-08/CC1pi_out.root", std::string outputdir="./", int nplots=50, int _slider_window=5){
//   gStyle->SetOptStat(0);
//   gStyle->SetOptTitle(0);
//
//   TFile *fin = new TFile(inputfile.c_str(),"READ");
//   TTree *tr = (TTree*)fin->Get("cc1piselec/outtree");
//
//   treevars mc_vars;
//   settreevars(tr,&mc_vars);
//
//   // Dummy, just to stop segfaults
//
//   std::string containedBookMVAType = "BDTG";
//   std::string containedBookMVALoc = "/uboone/app/users/ddevitt/LArSoft_v06_26_01_14_uboonecode_v06_26_01_22/srcs/uboonecode/uboone/CC1pi/MVA/dataset_contained/weights/TMVAClassification_BDTG.weights.xml";
//   std::string uncontainedBookMVAType = "BDTG";
//   std::string uncontainedBookMVALoc = "/uboone/app/users/ddevitt/LArSoft_v06_26_01_14_uboonecode_v06_26_01_22/srcs/uboonecode/uboone/CC1pi/MVA/dataset_uncontained/weights/TMVAClassification_BDTG.weights.xml";
//   TMVA::Reader fReader_contained("");
//   fReader_contained.AddVariable("dEdx_truncmean_start", &(mc_vars.float_dEdx_truncmean_start));
//   fReader_contained.AddVariable("VtxTrackDist", &(mc_vars.float_VtxTrackDist));
//   fReader_contained.AddVariable("nhits", &(mc_vars.float_nhits));
//   fReader_contained.AddVariable("lnLmipoverp", &(mc_vars.float_lnLmipoverp));
//   fReader_contained.BookMVA(containedBookMVAType.c_str(), containedBookMVALoc.c_str());
//   TMVA::Reader fReader_uncontained("");
//   fReader_uncontained.AddVariable("dEdx_truncmean_start", &(mc_vars.float_dEdx_truncmean_start));
//   fReader_uncontained.AddVariable("VtxTrackDist", &(mc_vars.float_VtxTrackDist));
//   fReader_uncontained.AddVariable("nhits", &(mc_vars.float_nhits));
//   fReader_uncontained.AddVariable("lnLmipoverp", &(mc_vars.float_lnLmipoverp));
//   fReader_uncontained.BookMVA(uncontainedBookMVAType.c_str(), uncontainedBookMVALoc.c_str());
//
//   TCanvas *c1 = new TCanvas("c1","",400,500);
//
//   int plots_made = 0;
//   for (int i_evt=0; i_evt<tr->GetEntries(); i_evt++){
//     tr->GetEntry(i_evt);
//
//     Calcvars(&mc_vars, &fReader_contained, &fReader_uncontained);
//
//     for (int i_track=0; i_track<mc_vars.TPCObj_PFP_isDaughter->size(); i_track++){
//
//       // Only make plots for muons and pions, and only if there are reconstructed spacepoints
//       if (mc_vars.TPCObj_PFP_track_SpacepointsXYZ_Ordered->at(i_track).size()==0) continue;
//       if (!(mc_vars.TPCObj_PFP_truePDG->at(i_track)==211)) continue;
//
//       bool makeplot = PlotLocalLinearityDetails(_slider_window, &mc_vars, c1, i_track, false);
//
//       if (!makeplot) continue;
//
//       TString savename = TString::Format("%s/linearity_%d_evt%d_pfp%d.pdf",outputdir.c_str(),plots_made,i_evt,i_track);
//       c1->Print(savename.Data());
//
//       plots_made++;
//       if (plots_made>=nplots) return;
//
//     } // end loop over tracks in event (i_track)
//   } // end loop over events in tree (i_evt)
// } // end function: MakeLocalLinearityPlots
