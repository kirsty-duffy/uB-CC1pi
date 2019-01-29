// #include "KinkFindingTree.h"
#include "KinkFindingTree.h"
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

#include <fstream>

void EvalKinkFinding(std::string inputfile="/uboone/data/users/kduffy/CC1pi/book/v06_26_01_13_PionReco_v1/pionrecovalid_merged.root", std::string outputdir="./"){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile *fin = new TFile(inputfile.c_str(),"OPEN");
  TTree *tr_haskink = (TTree*)fin->Get("kinkfinding/KinkFindingTree");
  // TTree *tr_nokink = (TTree*)fin->Get("nokink_tree");

  TH1F *nkinks_nokink = new TH1F("nkinks_nokink","No kink (truth);Found one or more reco kink(s)?;Area normalised",2,-0.5,1.5);
  TH1F *nkinks_manymcp = new TH1F("nkinks_manymcp","Has kink (truth);Found one or more reco kink(s)?;Area normalised",2,-0.5,1.5);

  TH1F *kinkdistance = new TH1F("kinkdistance",";Min. reco-true kink distance;Area normalised",100,0,500);

  TH1F *nkinks_nokink_byangle = new TH1F("nkinks_nokink_byangle","No kink (truth);Found one or more reco kink(s)?;Area normalised",2,-0.5,1.5);
  TH1F *nkinks_manymcp_byangle = new TH1F("nkinks_manymcp_byangle","Has kink (truth);Found one or more reco kink(s)?;Area normalised",2,-0.5,1.5);

  TH1F *kinkdistance_byangle = new TH1F("kinkdistance_byangle",";Min. reco-true kink distance;Area normalised",100,0,500);

  KinkFindingTree KinkFindingTree_haskink;
  SetKinkFindingTreeVariables(&KinkFindingTree_haskink,tr_haskink);

  KinkFindingTree KinkFindingTree_nokink;
  // SetKinkFindingTreeVariables(&KinkFindingTree_nokink,tr_nokink);


  for (int i_evt=0; i_evt<tr_haskink->GetEntries(); i_evt++){
    if (i_evt%50==0) std::cout << "tree haskink: " << i_evt << "/" << tr_haskink->GetEntries() << std::endl;

    tr_haskink->GetEntry(i_evt);

    bool goodevt = KinkFindingTree_haskink.Setup();
    if (!goodevt) continue;

    std::vector<double> px_0 = KinkFindingTree_haskink.MCP_Px_eachpoint->at(0);
    std::vector<double> py_0 = KinkFindingTree_haskink.MCP_Py_eachpoint->at(0);
    std::vector<double> pz_0 = KinkFindingTree_haskink.MCP_Pz_eachpoint->at(0);
    std::vector<double> px_1 = KinkFindingTree_haskink.MCP_Px_eachpoint->at(1);
    std::vector<double> py_1 = KinkFindingTree_haskink.MCP_Py_eachpoint->at(1);
    std::vector<double> pz_1 = KinkFindingTree_haskink.MCP_Pz_eachpoint->at(1);

    // Find last non-(0,0,0) momentum in MCP 0
    double vec0_x=0, vec0_y=0, vec0_z=0;
    for (int i_p=0; i_p<px_0.size(); i_p++){
      if (!(px_0.at(i_p)==0 && py_0.at(i_p) == 0 && pz_0.at(i_p) ==0)){
        vec0_x = px_0.at(i_p);
        vec0_y = py_0.at(i_p);
        vec0_z = pz_0.at(i_p);
      }
    }

    TVector3 vec0(vec0_x,vec0_y,vec0_z);
    TVector3 vec1(px_1.at(0),py_1.at(0),pz_1.at(0));
    vec0.SetMag(1);
    vec1.SetMag(1);

    std::cout << "Event " << i_evt << ": true cos(th) = " << vec0.Dot(vec1) << std::endl;

    for (size_t i_pfph=0; i_pfph<KinkFindingTree_haskink.n_PFPs; i_pfph++){
      int pfp_idx=i_pfph;
      int pfpid = KinkFindingTree_haskink.PFP_ID->at(pfp_idx);

      std::vector<std::vector<double>> truekinks;
      truekinks.push_back(KinkFindingTree_haskink.PFP_trueMCPend->at(pfp_idx));

      // Now calculate local linearity variables to look for kinks
      std::vector<int> kinkidxs = KinkFindingTree_haskink.FindRecoKinks_orderedsp(pfp_idx,vec0.Dot(vec1),i_evt,truekinks);

      // Now look at kinks by angle betwen groups of spacepoints
      std::vector<int> kinkidxs_byangle = KinkFindingTree_haskink.FindRecoKinks_byangle_orderedsp(pfp_idx,vec0.Dot(vec1),i_evt,0,20,truekinks);


      // Fill plots!
      if (kinkidxs.size()>0) nkinks_manymcp->Fill(1);
      else nkinks_manymcp->Fill(0);

      if (kinkidxs_byangle.size()>0) nkinks_manymcp_byangle->Fill(1);
      else nkinks_manymcp_byangle->Fill(0);


      for (int i_truekink=0; i_truekink<truekinks.size(); i_truekink++){
        std::vector<double> truekink_xyz = truekinks.at(i_truekink);

        // Now loop through all reco kinks and find the closest one
        double dist = 9999999;
        for (int i_r=0; i_r<kinkidxs.size(); i_r++){
          std::vector<double> recokink_xyz = KinkFindingTree_haskink.PFP_ordered_spacepoints->at(pfp_idx).at(kinkidxs.at(i_r));
          double dx = recokink_xyz.at(0)-truekink_xyz.at(0);
          double dy = recokink_xyz.at(1)-truekink_xyz.at(1);
          double dz = recokink_xyz.at(2)-truekink_xyz.at(2);
          double distance = TMath::Sqrt(dx*dx+dy*dy+dz*dz);
          if (distance < dist){
            dist = distance;
          }
        } // end loop over reco indices

        // Fill plot
        if (dist>kinkdistance->GetXaxis()->GetXmax() && dist!=9999999) dist = kinkdistance->GetXaxis()->GetXmax()-0.5;
        kinkdistance->Fill(dist);

        // Repeat for angle local linearity method
        dist = 9999999;
        for (int i_r=0; i_r<kinkidxs_byangle.size(); i_r++){
          std::vector<double> recokink_xyz = KinkFindingTree_haskink.PFP_ordered_spacepoints->at(pfp_idx).at(kinkidxs_byangle.at(i_r));
          double dx = recokink_xyz.at(0)-truekink_xyz.at(0);
          double dy = recokink_xyz.at(1)-truekink_xyz.at(1);
          double dz = recokink_xyz.at(2)-truekink_xyz.at(2);
          double distance = TMath::Sqrt(dx*dx+dy*dy+dz*dz);
          if (distance < dist){
            dist = distance;
          }
        } // end loop over reco indices

        // Fill plot
        if (dist>kinkdistance_byangle->GetXaxis()->GetXmax() && dist!=9999999) dist = kinkdistance_byangle->GetXaxis()->GetXmax()-0.5;
        kinkdistance_byangle->Fill(dist);
      }

    } // end loop i_pfph

  } // end loop i_evt (tree_haskink)

  // Now do the same for the nokink tree
  // for (int i_evt=0; i_evt<tr_nokink->GetEntries(); i_evt++){
  //   if (i_evt%50==0) std::cout << "tree nokink: " << i_evt << "/" << tr_nokink->GetEntries() << std::endl;
  //
  //   tr_nokink->GetEntry(i_evt);
  //
  //   bool goodevt = KinkFindingTree_nokink.Setup();
  //   if (!goodevt) continue;
  //
  //   for (size_t i_pfph=0; i_pfph<KinkFindingTree_nokink.n_PFPs; i_pfph++){
  //     int pfp_idx=i_pfph;
  //
  //     // Now calculate local linearity variables to look for kinks
  //     std::vector<int> kinkidxs = KinkFindingTree_nokink.FindRecoKinks_orderedsp(pfp_idx,-999,i_evt*-1);
  //
  //     // Now look at kinks by angle betwen groups of spacepoints
  //     std::vector<int> kinkidxs_byangle = KinkFindingTree_nokink.FindRecoKinks_byangle_orderedsp(pfp_idx,-999,i_evt*-1,0,20);
  //
  //
  //     // Fill plots!
  //     if (kinkidxs.size()>0) nkinks_nokink->Fill(1);
  //     else nkinks_nokink->Fill(0);
  //
  //     if (kinkidxs_byangle.size()>0) nkinks_nokink_byangle->Fill(1);
  //     else nkinks_nokink_byangle->Fill(0);
  //
  //   } // end loop i_pfph
  //
  // } // end loop i_evt (tree_haskink)

  TFile *fout = new TFile("kinkfindingplots.root","recreate");

  TCanvas *c1 = new TCanvas();

  // nkinks_nokink->Scale(1./nkinks_nokink->Integral());
  // nkinks_nokink->GetYaxis()->SetRangeUser(0,1.1);
  // nkinks_nokink->SetLineWidth(2);
  // nkinks_nokink->SetLineColor(kBlack);
  // nkinks_nokink->Draw();
  nkinks_manymcp->Scale(1./nkinks_manymcp->Integral());
  nkinks_manymcp->SetLineWidth(2);
  nkinks_manymcp->SetLineStyle(2);
  nkinks_manymcp->SetLineColor(kRed+2);
  nkinks_manymcp->Draw("");

  TLegend *l = new TLegend(0.7,0.84,0.93,0.97);
  l->SetFillColor(kWhite);
  l->AddEntry(nkinks_nokink,"No kink","l");
  l->AddEntry(nkinks_manymcp,"Multiple MCPs","l");
  l->Draw();

  c1->Print("nkinks_sorted.pdf");

  nkinks_nokink->Write("nkinks_nokink");
  nkinks_manymcp->Write("nkinks_manymcp");
  c1->Write("nkinks_handscanned_canv");


  // nkinks_nokink_byangle->Scale(1./nkinks_nokink_byangle->Integral());
  // nkinks_nokink_byangle->GetYaxis()->SetRangeUser(0,1.1);
  // nkinks_nokink_byangle->SetLineWidth(2);
  // nkinks_nokink_byangle->SetLineColor(kBlack);
  // nkinks_nokink_byangle->Draw();
  nkinks_manymcp_byangle->Scale(1./nkinks_manymcp_byangle->Integral());
  nkinks_manymcp_byangle->SetLineWidth(2);
  nkinks_manymcp_byangle->SetLineStyle(2);
  nkinks_manymcp_byangle->SetLineColor(kRed+2);
  nkinks_manymcp_byangle->Draw("");
  l->Draw();
  c1->Print("nkinks_handscanned_byangle.pdf");

  nkinks_nokink->Write("nkinks_nokink_byangle");
  nkinks_manymcp->Write("nkinks_manymcp_byangle");
  c1->Write("nkinks_handscanned_byangle_canv");

  kinkdistance->SetLineWidth(2);
  kinkdistance->SetLineStyle(2);
  kinkdistance->SetLineColor(kRed+2);
  kinkdistance->DrawNormalized();
  c1->Print("kinkdist_manymcp_haskink.pdf");
  kinkdistance->Write("kinkdistance");
  c1->Write("kinkdistance_canv");

  kinkdistance_byangle->SetLineWidth(2);
  kinkdistance_byangle->SetLineStyle(2);
  kinkdistance_byangle->SetLineColor(kRed+2);
  kinkdistance_byangle->DrawNormalized();
  c1->Print("kinkdist_byangle_manymcp_haskink.pdf");
  kinkdistance->Write("kinkdistance_byangle");
  c1->Write("kinkdistance_byangle_canv");

  fout->Close();
}
