#include "PiRecoTree.h"
#include "CalcLocalLinearity.h"

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


void MakeLocalLinearityPlots(std::string inputfile="/uboone/data/users/kduffy/CC1pi/book/v06_26_01_13_PionReco_v1/pionrecovalid_merged.root", std::string outputdir="./", int nplots=50, int _slider_window=20){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile *fin = new TFile(inputfile.c_str(),"OPEN");
  TTree *tr = (TTree*)fin->Get("pirecovalid/PiRecoTree");

  PiRecoTree PiRecoTree;
  SetPiRecoTreeVariables(&PiRecoTree,tr);

  // Define colours for plotting PFPs
  // kGreen(416)+3, kBlue(600), kViolet(880), kPink(900), kOrange(800)-3
  int pfp_cols[5] = {419,600,880,900,797};
  // TGraph2D *spacepointsgr[5];
  // for (int i=0; i<5; i++){
  //     spacepointsgr[i] = new TGraph2D();
  //     spacepointsgr[i]->SetMarkerColor(pfp_cols[i]);
  //     spacepointsgr[i]->SetMarkerStyle(7);
  // }
  TGraph *spacepoints_yz[5];
  TGraph *spacepoints_xz[5];
  for (int i=0; i<5; i++){
      spacepoints_yz[i] = new TGraph();
      spacepoints_yz[i]->SetMarkerColor(pfp_cols[i]);
      spacepoints_yz[i]->SetMarkerStyle(7);
      spacepoints_xz[i] = new TGraph();
      spacepoints_xz[i]->SetMarkerColor(pfp_cols[i]);
      spacepoints_xz[i]->SetMarkerStyle(7);
  }

  // TGraph2D *startpoint = new TGraph2D(1);
  // startpoint->SetMarkerColor(kRed);
  // startpoint->SetMarkerStyle(29);
  TGraph *startpoint_yz = new TGraph(1);
  startpoint_yz->SetMarkerColor(kRed);
  startpoint_yz->SetMarkerStyle(29);
  TGraph *startpoint_xz = new TGraph(1);
  startpoint_xz->SetMarkerColor(kRed);
  startpoint_xz->SetMarkerStyle(29);

  TGraph *linearity = new TGraph();

  TCanvas *c1 = new TCanvas("c1","",400,500);
  c1->Divide(1,3);

  int plots_made = 0;
  for (int i_evt=0; i_evt<tr->GetEntries(); i_evt++){

    // Clear TGraph2D
    // for (int i=0; i<5; i++){
    //     spacepointsgr[i]->Clear();
    // }
    // startpoint->Clear();
    for (int i=0; i<5; i++){
        spacepoints_yz[i]->Set(0);
        spacepoints_xz[i]->Set(0);
    }
    startpoint_yz->Set(0);
    startpoint_xz->Set(0);
    linearity->Set(0);
    // Reset canvas
    c1->Clear();
    c1->Divide(1,3);

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

    std::vector<double> LocalLin;
    std::vector<int> PFPboundaries;

     std::cout << "Reconstructed " << PiRecoTree.bestmatch_PFP_hierarchy_idxs.size() << " PFPs in hierarchy" << std::endl;

     for (size_t i_pfph=0; i_pfph<PiRecoTree.bestmatch_PFP_hierarchy_idxs.size(); i_pfph++){
       std::cout << i_pfph << std::endl;
         if (i_pfph>4) {
             std::cout << "ERROR more than 5 PFPs! Can't fill plot. Skipping PFPs with index " << i_pfph << std::endl;
             continue;
         }

         int pfp_idx=PiRecoTree.bestmatch_PFP_hierarchy_idxs.at(i_pfph);
         int i_pt = 0;
         for (size_t i_sp=0; i_sp<PiRecoTree.PFP_spacepoints_XYZ->at(pfp_idx).size(); i_sp++){
             double x = PiRecoTree.PFP_spacepoints_XYZ->at(pfp_idx).at(i_sp).at(0);
             double y = PiRecoTree.PFP_spacepoints_XYZ->at(pfp_idx).at(i_sp).at(1);
             double z = PiRecoTree.PFP_spacepoints_XYZ->at(pfp_idx).at(i_sp).at(2);
             // spacepointsgr[i_pfph]->SetPoint(i_pt,x,y,z);
             spacepoints_yz[i_pfph]->SetPoint(i_pt+1,z,y);
             spacepoints_xz[i_pfph]->SetPoint(i_pt+1,z,x);
             i_pt++;

              if (x<minx) minx = x;
              if (x>maxx) maxx = x;
              if (y<miny) miny = y;
              if (y>maxy) maxy = y;
              if (z<minz) minz = z;
              if (z>maxz) maxz = z;
          }
        if (i_pfph==0){
          // std::cout << "Setting start point" << std::endl;
          // std::cout << "pfp_idx = " << pfp_idx << std::endl;
          // std::cout << "i_evt = " << i_evt << std::endl;
          double x = PiRecoTree.PFP_track_start->at(pfp_idx).at(0);
          double y = PiRecoTree.PFP_track_start->at(pfp_idx).at(1);
          double z = PiRecoTree.PFP_track_start->at(pfp_idx).at(2);
          // startpoint->SetPoint(1,x,y,z);
          startpoint_yz->SetPoint(1,z,y);
          startpoint_xz->SetPoint(1,z,x);

          if (x<minx) minx = x;
          if (x>maxx) maxx = x;
          if (y<miny) miny = y;
          if (y>maxy) maxy = y;
          if (z<minz) minz = z;
          if (z>maxz) maxz = z;

          // std::cout << x << ", " << y << ", " << z << std::endl;

          // Order spacepoints
          std::vector<std::vector<double>> ordered_spacepoints = OrderSpacepoints(PiRecoTree.PFP_spacepoints_XYZ->at(pfp_idx),PiRecoTree.PFP_track_start->at(pfp_idx));

          // Get local linearity
          std::vector<double> tmp = GetLocalLinearityVec(ordered_spacepoints, _slider_window);
          LocalLin.insert(LocalLin.end(), tmp.begin(), tmp.end());
          PFPboundaries.push_back(LocalLin.size()-1);
        }

    }

    // Fill local linearity plot
    for (size_t i_pt=0; i_pt < LocalLin.size(); i_pt++){
      linearity->SetPoint(i_pt+1,i_pt,LocalLin.at(i_pt));
    }

    // Formatting
    // spacepointsgr[0]->GetXaxis()->SetRangeUser(minx-5,maxx+5);
    // spacepointsgr[0]->GetYaxis()->SetRangeUser(miny-5,maxy+5);
    // spacepointsgr[0]->GetZaxis()->SetRangeUser(minz-5,maxz+5);
    // spacepointsgr[0]->GetXaxis()->SetTitle("X");
    // spacepointsgr[0]->GetYaxis()->SetTitle("Y");
    // spacepointsgr[0]->GetZaxis()->SetTitle("Z");
    spacepoints_yz[0]->GetYaxis()->SetRangeUser(miny-5,maxy+5);
    spacepoints_yz[0]->GetXaxis()->SetRangeUser(minz-5,maxz+5);
    spacepoints_yz[0]->GetYaxis()->SetTitle("Y (cm)");
    spacepoints_yz[0]->GetXaxis()->SetTitle("Z (cm)");
    spacepoints_xz[0]->GetYaxis()->SetRangeUser(minx-5,maxx+5);
    spacepoints_xz[0]->GetXaxis()->SetRangeUser(minz-5,maxz+5);
    spacepoints_xz[0]->GetYaxis()->SetTitle("X (cm)");
    spacepoints_xz[0]->GetXaxis()->SetTitle("Z (cm)");



    // Draw plot
    c1->cd(2);
    gPad->Divide(2);
    if (spacepoints_yz[0]->GetN()>0/*spacepointsgr[0]->GetN()>0*/){
      // spacepointsgr[0]->Draw("Ap");
      c1->cd(2); gPad->cd(1);
      spacepoints_yz[0]->Draw("Ap");
      c1->cd(2); gPad->cd(2);
      spacepoints_xz[0]->Draw("Ap");
    }
    else{
      std::cout << "No PFPs to draw for event " << i_evt << std::endl;
      continue;
    }
    for (int i=1;i<5;i++){
        // if (spacepointsgr[i]->GetN()>0) spacepointsgr[i]->Draw("same p");
        c1->cd(2); gPad->cd(1);
        if (spacepoints_yz[i]->GetN()>0) spacepoints_yz[i]->Draw("same p");
        c1->cd(2); gPad->cd(2);
        if (spacepoints_xz[i]->GetN()>0) spacepoints_xz[i]->Draw("same p");
    }
    // startpoint->Draw("same p");
    // c1->cd(2); gPad->cd(1);
    // startpoint_yz->Draw("same p");
    // c1->cd(2); gPad->cd(2);
    // startpoint_xz->Draw("same p");

    // Draw local linearity plot
    c1->cd(3);
    linearity->GetXaxis()->SetTitle("Spacepoint index in PFP (1st PFP only)");
    linearity->GetYaxis()->SetTitle("Local linearity");
    linearity->GetYaxis()->SetRangeUser(0,1.05);
    linearity->Draw("APL");

    // Add lines to separate out PFPs
    // for (size_t i_pfp=1; i_pfp<PiRecoTree.bestmatch_PFP_hierarchy_idxs.size(); i_pfp++){
    //   TLine *l = new TLine(PFPboundaries.at(i_pfp),0,PFPboundaries.at(i_pfp),1.05);
    //   l->SetLineStyle(2);
    //   l->SetLineColor(kGreen+2);
    //   l->SetLineWidth(2);
    //   l->Draw("same");
    // }

    // Make legend
    c1->cd(1);
    TPaveText *pt = new TPaveText(0.05,0.41,0.95,0.95,"NB");
    pt->SetLineColor(kWhite);
    pt->SetFillStyle(0);
    pt->SetTextFont(132);
    pt->AddText(TString::Format("%d #pi^{+} in 'logical pion'",(int)PiRecoTree.PiPlusHierarchy_LogicalPion_MCPids->size()).Data());
    pt->AddText(TString::Format("Daughters: %d protons (%d with p>300 MeV),  %d #mu^{+}",PiRecoTree.n_p,PiRecoTree.n_visp,PiRecoTree.n_mu).Data());
    pt->AddText(TString::Format("Primary PFP matched to PDG code %d",PiRecoTree.primarypfp_pdgcode).Data());
    pt->Draw();

    TLegend *l = new TLegend(0.05,0.05,0.95,0.39);
    l->SetNColumns(3);
    l->SetLineColor(kWhite);
    l->SetFillStyle(0);
    l->SetTextFont(132);
    // l->AddEntry(startpoint,"Reco track start","p");
    // for (int i=0;i<5;i++){
    //     if (spacepointsgr[i]->GetN()>0) l->AddEntry(spacepointsgr[i],TString::Format("PFP %d",i+1).Data(),"p");
    //     else l->AddEntry((TObject*)0,"","");
    // }
    // l->AddEntry(startpoint_yz,"Reco track start","p");
    for (int i=0;i<5;i++){
        if (spacepoints_yz[i]->GetN()>0) l->AddEntry(spacepoints_yz[i],TString::Format("PFP %d",i+1).Data(),"p");
        else l->AddEntry((TObject*)0,"","");
    }
    l->Draw();

    c1->Update();
    c1->Draw();
    TString savename = TString::Format("%s/pion_%d.pdf",outputdir.c_str(),plots_made);
    std::cout << savename << std::endl;
    c1->Print(savename.Data());

    plots_made++;
    if (plots_made>=nplots) break;
  } // end loop over events in tree (i_evt)
}







void MakeLocalLinearitySummaryPlots(std::string inputfile="/uboone/data/users/kduffy/CC1pi/book/v06_26_01_13_PionReco_v1/pionrecovalid_merged.root", std::string outputdir="./", int _slider_window=20){
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
