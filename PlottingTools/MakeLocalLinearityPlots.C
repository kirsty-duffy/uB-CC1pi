#include "CC1pi_treevars.h"
#include "CC1pi_treevars.cxx"
#include "CalcLocalLinearity.h"
#include "../Algorithms/PDGEnums.h"

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


bool PlotLocalLinearityDetails(int _slider_window, treevars *vars, TCanvas *c1, int i_track, bool isSelected){

  // Make TGraphs to fill
  TGraph *spacepoints_yz = new TGraph();
  spacepoints_yz->SetMarkerColor(kViolet);
  spacepoints_yz->SetMarkerStyle(7);
  TGraph *spacepoints_xz = new TGraph();
  spacepoints_xz->SetMarkerColor(kViolet);
  spacepoints_xz->SetMarkerStyle(7);
  TGraph *spacepoints_yz_bkg = new TGraph();
  spacepoints_yz_bkg->SetMarkerColor(kGray);
  spacepoints_yz_bkg->SetMarkerStyle(6);
  TGraph *spacepoints_xz_bkg = new TGraph();
  spacepoints_xz_bkg->SetMarkerColor(kGray);
  spacepoints_xz_bkg->SetMarkerStyle(6);

  TGraph *startpoint_yz = new TGraph(1);
  startpoint_yz->SetMarkerColor(kRed);
  startpoint_yz->SetMarkerStyle(29);
  TGraph *startpoint_xz = new TGraph(1);
  startpoint_xz->SetMarkerColor(kRed);
  startpoint_xz->SetMarkerStyle(29);

  TGraph *nuvtx_yz = new TGraph(1);
  nuvtx_yz->SetMarkerColor(kGreen+3);
  nuvtx_yz->SetMarkerStyle(29);
  TGraph *nuvtx_xz = new TGraph(1);
  nuvtx_xz->SetMarkerColor(kGreen+3);
  nuvtx_xz->SetMarkerStyle(29);

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

  // TGraph *hitcharge = new TGraph();
  TGraph *dirs = new TGraph();

  TGraph *cusum = new TGraph();

  // Reset canvas
  c1->Clear();
  c1->Divide(1,5);

  // Fill spacepoints plots
  double minx = 9999;
  double maxx = -9999;
  double miny = 9999;
  double maxy = -9999;
  double minz = 9999;
  double maxz = -9999;

  std::vector<double> LocalLin;
  std::vector<int> PFPboundaries;

  // Draw spacepoints from all PFPs in grey, except from the one we're currently looking at (draw that in a brighter colour)
   int i_pt = 0;
   int i_pt_bkg = 0;
   for (size_t i_pfp=0; i_pfp<vars->TPCObj_PFP_track_SpacepointsXYZ_Ordered->size(); i_pfp++){
     //std::cout << i_pfph << std::endl;

    for (size_t i_sp=0; i_sp<vars->TPCObj_PFP_track_SpacepointsXYZ_Ordered->at(i_pfp).size(); i_sp++){
      double x = vars->TPCObj_PFP_track_SpacepointsXYZ_Ordered->at(i_pfp).at(i_sp).at(0);
      double y = vars->TPCObj_PFP_track_SpacepointsXYZ_Ordered->at(i_pfp).at(i_sp).at(1);
      double z = vars->TPCObj_PFP_track_SpacepointsXYZ_Ordered->at(i_pfp).at(i_sp).at(2);

      spacepoints_yz_bkg->SetPoint(i_pt_bkg,z,y);
      spacepoints_xz_bkg->SetPoint(i_pt_bkg,z,x);
      i_pt_bkg++;

      if (i_pfp == i_track){
        spacepoints_yz->SetPoint(i_pt,z,y);
        spacepoints_xz->SetPoint(i_pt,z,x);
        i_pt++;
      }

       if (x<minx) minx = x;
       if (x>maxx) maxx = x;
       if (y<miny) miny = y;
       if (y>maxy) maxy = y;
       if (z<minz) minz = z;
       if (z>maxz) maxz = z;

       if (i_pfp==i_track && i_sp==0){ // set start point
         startpoint_yz->SetPoint(0,z,y);
         startpoint_xz->SetPoint(0,z,x);
       }
    } // end loop over spacepoints
  } // end loop over other pfps

  // Set neutrino vertex point
  double x = vars->TPCObj_reco_vtx->at(0);
  double y = vars->TPCObj_reco_vtx->at(1);
  double z = vars->TPCObj_reco_vtx->at(2);
  nuvtx_yz->SetPoint(0,z,y);
  nuvtx_xz->SetPoint(0,z,x);

  if (x<minx) minx = x;
  if (x>maxx) maxx = x;
  if (y<miny) miny = y;
  if (y>maxy) maxy = y;
  if (z<minz) minz = z;
  if (z>maxz) maxz = z;

  // Fill local linearity plot and hit charge plot
  int i_pt_charge = 0;
  for (size_t i_pt=0; i_pt < vars->TPCObj_PFP_track_Spacepoints_LocalLin->at(i_track).size(); i_pt++){

    cusum->SetPoint(i_pt,i_pt,vars->TPCObj_PFP_track_Spacepoints_DirCuSum->at(i_track).at(i_pt));

    ttest->SetPoint(i_pt,i_pt,vars->TPCObj_PFP_track_Spacepoints_LocalLin->at(i_track).at(i_pt));

    dirs->SetPoint(i_pt,i_pt,vars->TPCObj_PFP_track_Spacepoints_Dirz->at(i_track).at(i_pt));

    // double charge = vars->TPCObj_PFP_track_SpacepointsdQdsPlane2_Ordered->at(i_track).at(i_pt);
    // if (charge != -9999){
    //   hitcharge->SetPoint(i_pt_charge,i_pt,charge);
    //   i_pt_charge++;
    // }
  }

  for (int i_kink=0; i_kink<vars->TPCObj_PFP_track_Spacepoints_kinkidxs->at(i_track).size(); i_kink++){
    int kink_idx = vars->TPCObj_PFP_track_Spacepoints_kinkidxs->at(i_track).at(i_kink);
    double x = vars->TPCObj_PFP_track_SpacepointsXYZ_Ordered->at(i_track).at(kink_idx).at(0);
    double y = vars->TPCObj_PFP_track_SpacepointsXYZ_Ordered->at(i_track).at(kink_idx).at(1);
    double z = vars->TPCObj_PFP_track_SpacepointsXYZ_Ordered->at(i_track).at(kink_idx).at(2);

    breakpoints_yz->SetPoint(i_kink,z,y);
    breakpoints_xz->SetPoint(i_kink,z,x);
    breakpoints_idx->SetPoint(i_kink,kink_idx,0);
  }




  // Draw plot
  // Have to draw first TGraph, then set axis ranges, then draw again because it's the only way ROOT will do what I want (the TAxis objects don't exist until you draw the TGraph so you can't set ranges)
  c1->cd(2);
  gPad->Divide(2);
  if (spacepoints_yz_bkg->GetN()>0){
    c1->cd(2); gPad->cd(1);
    // spacepoints_yz_bkg->Draw("Ap");
    // spacepoints_yz_bkg->GetYaxis()->SetRangeUser(miny-5,maxy+5);
    // spacepoints_yz_bkg->GetXaxis()->SetRangeUser(minz-5,maxz+5);
    spacepoints_yz_bkg->GetYaxis()->SetTitle("Y (cm)");
    spacepoints_yz_bkg->GetXaxis()->SetTitle("Z (cm)");
    spacepoints_yz_bkg->Draw("Ap");
    // c1->Update();
    spacepoints_yz->Draw("same p");
    nuvtx_yz->Draw("same p");
    startpoint_yz->Draw("same p");
    breakpoints_yz->Draw("same p");
    c1->cd(2); gPad->cd(2);
    // spacepoints_xz_bkg->Draw("Ap");
    // spacepoints_xz_bkg->GetYaxis()->SetRangeUser(minx-5,maxx+5);
    // spacepoints_xz_bkg->GetXaxis()->SetRangeUser(minz-5,maxz+5);
    spacepoints_xz_bkg->GetYaxis()->SetTitle("X (cm)");
    spacepoints_xz_bkg->GetXaxis()->SetTitle("Z (cm)");
    spacepoints_xz_bkg->Draw("Ap");
    // c1->Update();
    spacepoints_xz->Draw("same p");
    nuvtx_xz->Draw("same p");
    startpoint_xz->Draw("same p");
    breakpoints_xz->Draw("same p");
  }
  else{
    return false;
  }

  c1->cd(3);
  dirs->GetXaxis()->SetTitle("Spacepoint index in PFP (coloured PFP only)");
  dirs->GetYaxis()->SetTitle("Direction relative to track start");
  dirs->Draw("APL");
  breakpoints_idx->Draw("P same");

  // Draw local linearity plot
  c1->cd(4);
  ttest->GetXaxis()->SetTitle("Spacepoint index in PFP (coloured PFP only)");
  ttest->GetYaxis()->SetTitle("Welch's t-test statistic");
  // ttest->GetYaxis()->SetRangeUser(0,1.05);
  ttest->Draw("APL");
  breakpoints_idx->Draw("P same");

  // c1->cd(4);
  // hitcharge->GetXaxis()->SetTitle("Spacepoint index in PFP (coloured PFP only)");
  // hitcharge->GetYaxis()->SetTitle("Plane 2 dQ/ds");
  // hitcharge->Draw("APL");

  c1->cd(5);
  cusum->GetXaxis()->SetTitle("Spacepoint index in PFP (coloured PFP only)");
  cusum->GetYaxis()->SetTitle("Cumulative sum of differences between values and average");
  cusum->Draw("APL");
  breakpoints_idx->Draw("P same");

  // Make legend
  c1->cd(1);
  TPaveText *pt = new TPaveText(0.05,0.41,0.95,0.95,"NB");
  pt->SetLineColor(kWhite);
  pt->SetFillStyle(0);
  pt->SetTextFont(132);
  pt->AddText(TString::Format("PFP matched to true %s",PDGenum2str((PDGCode)vars->TPCObj_PFP_truePDG->at(i_track)).c_str()).Data());
  pt->AddText(TString::Format("%.0f reconstructed daughters",vars->TPCObj_PFP_ndaughters->at(i_track)).Data());
  pt->AddText(TString::Format("True interaction type: %s",topologyenum2str((NuIntTopology)vars->Truth_topology).c_str()).Data());
  if (isSelected==true){
    pt->AddText("Event passes selection cuts");
  }
  else if (isSelected==false){
    pt->AddText("Event is not selected");
  }
  pt->Draw();

  TLegend *l = new TLegend(0.05,0.05,0.95,0.39);
  l->SetNColumns(3);
  l->SetLineColor(kWhite);
  l->SetFillStyle(0);
  l->SetTextFont(132);
  l->AddEntry(nuvtx_yz,"Reco #nu vertex","p");
  l->AddEntry(startpoint_yz,"Linearity start point","p");
  l->AddEntry(spacepoints_yz,"Studied PFP","p");
  l->AddEntry(spacepoints_yz_bkg,"Other PFPs","p");
  l->Draw();

  c1->Update();
  c1->Draw();

  return true;

} // end function: PlotLocalLinearityDetails


void MakeLocalLinearityPlots(std::string inputfile="/uboone/app/users/kduffy/CC1pi/CC1pi_uboonecode/srcs/uboonecode/uboone/CC1pi/ignore/2018-11-08/CC1pi_out.root", std::string outputdir="./", int nplots=50, int _slider_window=5){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile *fin = new TFile(inputfile.c_str(),"READ");
  TTree *tr = (TTree*)fin->Get("cc1piselec/outtree");

  treevars mc_vars;
  settreevars(tr,&mc_vars);

  // Dummy, just to stop segfaults

  std::string BookMVAType = "BDTG";
  std::string BookMVALoc = "/uboone/app/users/ddevitt/LArSoft_v06_26_01_14_uboonecode_v06_26_01_22/srcs/uboonecode/uboone/CC1pi/MVA/dataset_NeutrinoOnly/weights/TMVAClassification_BDTG.weights.xml";
  TMVA::Reader fReader_contained("");
  TMVA::Reader fReader("");
  fReader.AddVariable("dEdx_truncmean_start", &(mc_vars.float_dEdx_truncmean_start));
  fReader.AddVariable("VtxTrackDist", &(mc_vars.float_VtxTrackDist));
  fReader.AddVariable("nhits", &(mc_vars.float_nhits));
  fReader.AddVariable("lnLmipoverp", &(mc_vars.float_lnLmipoverp));
  fReader.BookMVA(BookMVAType.c_str(), BookMVALoc.c_str());

  TCanvas *c1 = new TCanvas("c1","",400,500);

  int plots_made = 0;
  for (int i_evt=0; i_evt<tr->GetEntries(); i_evt++){
    tr->GetEntry(i_evt);

    Calcvars(&mc_vars, &fReader);

    for (int i_track=0; i_track<mc_vars.TPCObj_PFP_isDaughter->size(); i_track++){

      // Only make plots for muons and pions, and only if there are reconstructed spacepoints
      if (mc_vars.TPCObj_PFP_track_SpacepointsXYZ_Ordered->at(i_track).size()==0) continue;
      if (!(mc_vars.TPCObj_PFP_truePDG->at(i_track)==211)) continue;

      bool makeplot = PlotLocalLinearityDetails(_slider_window, &mc_vars, c1, i_track, false);

      if (!makeplot) continue;

      TString savename = TString::Format("%s/linearity_%d_evt%d_pfp%d.pdf",outputdir.c_str(),plots_made,i_evt,i_track);
      c1->Print(savename.Data());

      plots_made++;
      if (plots_made>=nplots) return;

    } // end loop over tracks in event (i_track)
  } // end loop over events in tree (i_evt)
} // end function: MakeLocalLinearityPlots
