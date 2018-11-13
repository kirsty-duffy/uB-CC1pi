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

  TGraph *linearity = new TGraph();

  TGraph *hitcharge = new TGraph();

  // Reset canvas
  c1->Clear();
  c1->Divide(1,4);

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
   for (size_t i_pfp=0; i_pfp<vars->TPCObj_PFP_track_SpacepointsXYZ_Ordered->size(); i_pfp++){
     //std::cout << i_pfph << std::endl;

    int i_pt = 0;
    int i_pt_bkg = 0;
    for (size_t i_sp=0; i_sp<vars->TPCObj_PFP_track_SpacepointsXYZ_Ordered->at(i_pfp).size(); i_sp++){
      double x = vars->TPCObj_PFP_track_SpacepointsXYZ_Ordered->at(i_pfp).at(i_sp).at(0);
      double y = vars->TPCObj_PFP_track_SpacepointsXYZ_Ordered->at(i_pfp).at(i_sp).at(1);
      double z = vars->TPCObj_PFP_track_SpacepointsXYZ_Ordered->at(i_pfp).at(i_sp).at(2);

      if (i_pfp == i_track){
        spacepoints_yz->SetPoint(i_pt+1,z,y);
        spacepoints_xz->SetPoint(i_pt+1,z,x);
        i_pt++;
      }
      spacepoints_yz_bkg->SetPoint(i_pt_bkg+1,z,y);
      spacepoints_xz_bkg->SetPoint(i_pt_bkg+1,z,x);
      i_pt_bkg++;

       if (x<minx) minx = x;
       if (x>maxx) maxx = x;
       if (y<miny) miny = y;
       if (y>maxy) maxy = y;
       if (z<minz) minz = z;
       if (z>maxz) maxz = z;

       if (i_pfp==i_track && i_sp==0){ // set start point
         startpoint_yz->SetPoint(1,z,y);
         startpoint_xz->SetPoint(1,z,x);
       }
    } // end loop over spacepoints
  } // end loop over other pfps

  // Set neutrino vertex point
  double x = vars->TPCObj_reco_vtx->at(0);
  double y = vars->TPCObj_reco_vtx->at(1);
  double z = vars->TPCObj_reco_vtx->at(2);
  nuvtx_yz->SetPoint(1,z,y);
  nuvtx_xz->SetPoint(1,z,x);

  if (x<minx) minx = x;
  if (x>maxx) maxx = x;
  if (y<miny) miny = y;
  if (y>maxy) maxy = y;
  if (z<minz) minz = z;
  if (z>maxz) maxz = z;

  // Fill local linearity plot and hit charge plot
  int i_pt_charge = 0;
  for (size_t i_pt=0; i_pt < vars->TPCObj_PFP_track_Spacepoints_LocalLin->at(i_track).size(); i_pt++){
    linearity->SetPoint(i_pt+1,i_pt,vars->TPCObj_PFP_track_Spacepoints_LocalLin->at(i_track).at(i_pt));
    double charge = vars->TPCObj_PFP_track_SpacepointsdQdsPlane2_Ordered->at(i_track).at(i_pt);
    if (charge != -9999){
      hitcharge->SetPoint(i_pt_charge+1,i_pt,charge);
      i_pt_charge++;
    }
  }


  // Draw plot
  // Have to draw first TGraph, then set axis ranges, then draw again because it's the only way ROOT will do what I want (the TAxis objects don't exist until you draw the TGraph so you can't set ranges)
  c1->cd(2);
  gPad->Divide(2);
  if (spacepoints_yz_bkg->GetN()>0){
    c1->cd(2); gPad->cd(1);
    spacepoints_yz_bkg->Draw("Ap");
    spacepoints_yz_bkg->GetYaxis()->SetRangeUser(miny-5,maxy+5);
    spacepoints_yz_bkg->GetXaxis()->SetRangeUser(minz-5,maxz+5);
    spacepoints_yz_bkg->GetYaxis()->SetTitle("Y (cm)");
    spacepoints_yz_bkg->GetXaxis()->SetTitle("Z (cm)");
    spacepoints_yz_bkg->Draw("Ap");
    c1->Update();
    spacepoints_yz->Draw("same p");
    nuvtx_yz->Draw("same p");
    startpoint_yz->Draw("same p");
    c1->cd(2); gPad->cd(2);
    spacepoints_xz_bkg->Draw("Ap");
    spacepoints_xz_bkg->GetYaxis()->SetRangeUser(minx-5,maxx+5);
    spacepoints_xz_bkg->GetXaxis()->SetRangeUser(minz-5,maxz+5);
    spacepoints_xz_bkg->GetYaxis()->SetTitle("X (cm)");
    spacepoints_xz_bkg->GetXaxis()->SetTitle("Z (cm)");
    spacepoints_xz_bkg->Draw("Ap");
    c1->Update();
    spacepoints_xz->Draw("same p");
    nuvtx_xz->Draw("same p");
    startpoint_xz->Draw("same p");
  }
  else{
    return false;
  }

  // Draw local linearity plot
  c1->cd(3);
  linearity->GetXaxis()->SetTitle("Spacepoint index in PFP (coloured PFP only)");
  linearity->GetYaxis()->SetTitle("Local linearity");
  linearity->GetYaxis()->SetRangeUser(0,1.05);
  linearity->Draw("APL");

  c1->cd(4);
  hitcharge->GetXaxis()->SetTitle("Spacepoint index in PFP (coloured PFP only)");
  hitcharge->GetYaxis()->SetTitle("Plane 2 dQ/ds");
  hitcharge->Draw("APL");

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

  std::string containedBookMVAType = "BDTG";
  std::string containedBookMVALoc = "/uboone/app/users/ddevitt/LArSoft_v06_26_01_10/srcs/uboonecode/uboone/CC1pi/MVA/dataset_contained/weights/TMVAClassification_BDTG.weights.xml";
  std::string uncontainedBookMVAType = "BDTG";
  std::string uncontainedBookMVALoc = "/uboone/app/users/ddevitt/LArSoft_v06_26_01_10/srcs/uboonecode/uboone/CC1pi/MVA/dataset_uncontained/weights/TMVAClassification_BDTG.weights.xml";
  TMVA::Reader fReader_contained("");
  fReader_contained.AddVariable("dEdx_truncmean_start", &(mc_vars.float_dEdx_truncmean_start));
  fReader_contained.AddVariable("VtxTrackDist", &(mc_vars.float_VtxTrackDist));
  fReader_contained.AddVariable("nhits", &(mc_vars.float_nhits));
  fReader_contained.AddVariable("lnLmipoverp", &(mc_vars.float_lnLmipoverp));
  fReader_contained.BookMVA(containedBookMVAType.c_str(), containedBookMVALoc.c_str());
  TMVA::Reader fReader_uncontained("");
  fReader_uncontained.AddVariable("dEdx_truncmean_start", &(mc_vars.float_dEdx_truncmean_start));
  fReader_uncontained.AddVariable("VtxTrackDist", &(mc_vars.float_VtxTrackDist));
  fReader_uncontained.AddVariable("nhits", &(mc_vars.float_nhits));
  fReader_uncontained.AddVariable("lnLmipoverp", &(mc_vars.float_lnLmipoverp));
  fReader_uncontained.BookMVA(uncontainedBookMVAType.c_str(), uncontainedBookMVALoc.c_str());

  TCanvas *c1 = new TCanvas("c1","",400,500);

  int plots_made = 0;
  for (int i_evt=0; i_evt<tr->GetEntries(); i_evt++){
    tr->GetEntry(i_evt);
    Calcvars(&mc_vars, &fReader_contained, &fReader_uncontained);

    for (int i_track=0; i_track<mc_vars.TPCObj_PFP_isDaughter->size(); i_track++){

      // Only make plots for muons and pions, and only if there are reconstructed spacepoints
      if (mc_vars.TPCObj_PFP_track_SpacepointsXYZ_Ordered->at(i_track).size()==0) continue;
      if (!(mc_vars.TPCObj_PFP_truePDG->at(i_track)==13 ||mc_vars.TPCObj_PFP_truePDG->at(i_track)==211)) continue;

      bool makeplot = PlotLocalLinearityDetails(_slider_window, &mc_vars, c1, i_track, false);

      if (!makeplot) continue;

      TString savename = TString::Format("%s/linearity_%d_evt%d_pfp%d.pdf",outputdir.c_str(),plots_made,i_evt,i_track);
      c1->Print(savename.Data());

      plots_made++;
      if (plots_made>=nplots) break;

    } // end loop over tracks in event (i_track)
  } // end loop over events in tree (i_evt)
} // end function: MakeLocalLinearityPlots





/*void MakeLocalLinearitySummaryPlots(std::string inputfile="/uboone/data/users/kduffy/CC1pi/book/v06_26_01_13_PionReco_v1/pionrecovalid_merged.root", std::string outputdir="./", int _slider_window=20){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile *fin = new TFile(inputfile.c_str(),"OPEN");
  TTree *tr = (TTree*)fin->Get("pirecovalid/PiRecoTree");

  PiRecoTree PiRecoTree;
  SetPiRecoTreeVariables(&PiRecoTree,tr);


  TCanvas *c1 = new TCanvas();
  TH1D *mean_lin = new TH1D("mean_lin","True #pi^{+}: first PFP;Mean local linearity;Arb. units",50,0,1);
  TH1D *mean_lin_expectkink = new TH1D("mean_lin_expectkink","True #pi^{+}: first PFP;Mean local linearity;Arb. units",50,0,1);
  TH1D *mean_lin_other = new TH1D("mean_lin_other","True #pi^{+}: first PFP;Mean local linearity;Arb. units",50,0,1);
  TH1D *min_lin = new TH1D("min_lin","True #pi^{+}: first PFP;Minimum local linearity;Arb. units",50,0,1);
  TH1D *min_lin_expectkink = new TH1D("min_lin_expectkink","True #pi^{+}: first PFP;Minimum local linearity;Arb. units",50,0,1);
  TH1D *min_lin_other = new TH1D("min_lin_other","True #pi^{+}: first PFP;Minimum local linearity;Arb. units",50,0,1);

  TH2D *lin2d_expectkink = new TH2D("lin2d_expectkink","True #pi^{+}: first PFP (single PFP, multiple MCPs);Minimum local linearity;Mean local linearity",50,0,1,50,0,1);
  TH2D *lin2d_other = new TH2D("lin2d_expectkink","True #pi^{+}: first PFP (other tracks);Minimum local linearity;Mean local linearity",50,0,1,50,0,1);


  mean_lin_expectkink->SetLineColor(kRed);
  min_lin_expectkink->SetLineColor(kRed);

  for (int i_evt=0; i_evt<tr->GetEntries(); i_evt++){
    tr->GetEntry(i_evt);

    bool goodevt = PiRecoTree.Setup();
    if (!goodevt) continue;

    PiRecoTree.CalcLogicalPionVars();
    PiRecoTree.GetPFPHierarchy();

    std::vector<double> LocalLin;

     std::cout << "Reconstructed " << PiRecoTree.bestmatch_PFP_hierarchy_idxs.size() << " PFPs in hierarchy" << std::endl;

     for (size_t i_pfph=0; i_pfph<PiRecoTree.bestmatch_PFP_hierarchy_idxs.size(); i_pfph++){

       if (i_pfph>0) continue;

       int pfp_idx=PiRecoTree.bestmatch_PFP_hierarchy_idxs.at(i_pfph);

         // Order spacepoints
         std::vector<std::vector<double>> ordered_spacepoints = OrderSpacepoints(PiRecoTree.PFP_spacepoints_XYZ->at(pfp_idx),PiRecoTree.PFP_track_start->at(pfp_idx));

        // Get local linearity
        std::vector<double> tmp = GetLocalLinearityVec(ordered_spacepoints, _slider_window);
        LocalLin.insert(LocalLin.end(), tmp.begin(), tmp.end());
      }

    if (LocalLin.size()==0) continue;
    // Fill plots
    mean_lin->Fill(accumulate(LocalLin.begin(),LocalLin.end(),0.0)/LocalLin.size());
    min_lin->Fill(*std::min_element(LocalLin.begin(), LocalLin.end()));

    if (PiRecoTree.bestmatch_PFP_hierarchy_idxs.size()==1 && (PiRecoTree.n_mu+PiRecoTree.n_visp+PiRecoTree.PiPlusHierarchy_LogicalPion_MCPids->size())>1){
      mean_lin_expectkink->Fill(accumulate(LocalLin.begin(),LocalLin.end(),0.0)/LocalLin.size());
      min_lin_expectkink->Fill(*std::min_element(LocalLin.begin(), LocalLin.end()));
      lin2d_expectkink->Fill(*std::min_element(LocalLin.begin(), LocalLin.end()), accumulate(LocalLin.begin(),LocalLin.end(),0.0)/LocalLin.size());
    }
    else{
      mean_lin_other->Fill(accumulate(LocalLin.begin(),LocalLin.end(),0.0)/LocalLin.size());
      min_lin_other->Fill(*std::min_element(LocalLin.begin(), LocalLin.end()));
      lin2d_other->Fill(*std::min_element(LocalLin.begin(), LocalLin.end()), accumulate(LocalLin.begin(),LocalLin.end(),0.0)/LocalLin.size());
    }

  } // end loop over events in tree (i_evt)

  // Draw plots
  TLegend *l = new TLegend(0.13,0.7,0.5,0.89);
  l->SetLineColor(kWhite);
  l->SetFillStyle(0);
  l->SetTextFont(132);
  // l->AddEntry(mean_lin,"All tracks","l");
  l->AddEntry(mean_lin_expectkink,"Single PFP and multiple MCParticles","l");
  l->AddEntry(mean_lin,"Other pions","l");

  // mean_lin->Draw();
  mean_lin_other->DrawNormalized("");
  mean_lin_expectkink->DrawNormalized("same");
  l->Draw();
  TString savename = TString::Format("%s/mean_summary.pdf",outputdir.c_str());
  c1->Print(savename.Data());

  // min_lin->Draw();
  min_lin_expectkink->DrawNormalized("");
  min_lin_other->DrawNormalized("same");
  l->Draw();
  savename = TString::Format("%s/min_summary.pdf",outputdir.c_str());
  c1->Print(savename.Data());

  lin2d_expectkink->GetYaxis()->SetRangeUser(0.7,1);
  lin2d_other->GetYaxis()->SetRangeUser(0.7,1);

  lin2d_expectkink->Draw("colz");
  savename = TString::Format("%s/expectkink_2d_summary.pdf",outputdir.c_str());
  c1->Print(savename.Data());

  lin2d_other->Draw("colz");
  savename = TString::Format("%s/other_2d_summary.pdf",outputdir.c_str());
  c1->Print(savename.Data());

  lin2d_expectkink->SetLineColor(kRed);
  lin2d_expectkink->Draw("box");
  lin2d_other->Draw("box same");
  savename = TString::Format("%s/boxcomparison_2d_summary.pdf",outputdir.c_str());
  c1->Print(savename.Data());

}






void CheckSpacePointOrdering(std::string inputfile="/uboone/data/users/kduffy/CC1pi/book/v06_26_01_13_PionReco_v1/pionrecovalid_merged.root", std::string outputdir="./", int nplots=20){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile *fin = new TFile(inputfile.c_str(),"OPEN");
  TTree *tr = (TTree*)fin->Get("pirecovalid/PiRecoTree");

  PiRecoTree PiRecoTree;
  SetPiRecoTreeVariables(&PiRecoTree,tr);

  TCanvas *c1 = new TCanvas();
  c1->Divide(2);

  int plots_made = 0;
  for (int i_evt=0; i_evt<tr->GetEntries(); i_evt++){

    tr->GetEntry(i_evt);

    bool goodevt = PiRecoTree.Setup();
    if (!goodevt) continue;

    PiRecoTree.CalcLogicalPionVars();
    PiRecoTree.GetPFPHierarchy();

    // Fill spacepoints plots
    double minx = 9999;
    double maxx = -9999;
    double miny = 9999;
    double maxy = -9999;
    double minz = 9999;
    double maxz = -9999;

     for (size_t i_pfph=0; i_pfph<PiRecoTree.bestmatch_PFP_hierarchy_idxs.size(); i_pfph++){
         if (i_pfph>0) continue;

         int pfp_idx=PiRecoTree.bestmatch_PFP_hierarchy_idxs.at(i_pfph);

         std::vector<std::vector<double>> ordered_spacepoints = OrderSpacepoints(PiRecoTree.PFP_spacepoints_XYZ->at(pfp_idx),PiRecoTree.PFP_track_start->at(pfp_idx));

         for (size_t i_sp=0; i_sp<ordered_spacepoints.size(); i_sp++){
             double x = ordered_spacepoints.at(i_sp).at(0);
             double y = ordered_spacepoints.at(i_sp).at(1);
             double z = ordered_spacepoints.at(i_sp).at(2);

              if (x<minx) minx = x;
              if (x>maxx) maxx = x;
              if (y<miny) miny = y;
              if (y>maxy) maxy = y;
              if (z<minz) minz = z;
              if (z>maxz) maxz = z;
          }
    }

    // Make TH2Ds
    TH2D *spacepoints_yz = new TH2D("spacepoints_yz","Y-Z spacepoints (PFP 1);Z;Y",100,minz,maxz,100,miny,maxy);
    TH2D *spacepoints_xz = new TH2D("spacepoints_xz","X-Z spacepoints (PFP 1);Z;X",100,minz,maxz,100,minx,maxx);


    // and do it again...

     for (size_t i_pfph=0; i_pfph<PiRecoTree.bestmatch_PFP_hierarchy_idxs.size(); i_pfph++){
         if (i_pfph>0) continue;

         int pfp_idx=PiRecoTree.bestmatch_PFP_hierarchy_idxs.at(i_pfph);

         std::vector<std::vector<double>> ordered_spacepoints = OrderSpacepoints(PiRecoTree.PFP_spacepoints_XYZ->at(pfp_idx),PiRecoTree.PFP_track_start->at(pfp_idx));

         int i_pt = 0;
         for (size_t i_sp=0; i_sp<ordered_spacepoints.size(); i_sp++){
             double x = ordered_spacepoints.at(i_sp).at(0);
             double y = ordered_spacepoints.at(i_sp).at(1);
             double z = ordered_spacepoints.at(i_sp).at(2);

             int binyz = spacepoints_yz->FindBin(z,y);
             spacepoints_yz->SetBinContent(binyz,i_pt);
             int binxz = spacepoints_xz->FindBin(z,x);
             spacepoints_xz->SetBinContent(binxz,i_pt);

             i_pt++;

          }
    }

    if (spacepoints_yz->GetEntries()==0) continue;

    c1->cd(1);
    spacepoints_yz->Draw("colz");
    c1->cd(2);
    spacepoints_xz->Draw("colz");
    TString savename = TString::Format("%s/pion_%d_checksp.pdf",outputdir.c_str(),plots_made);
    c1->Print(savename.Data());

    plots_made++;
    if (plots_made>=nplots) break;
  } // end loop over events in tree (i_evt)
}
*/
