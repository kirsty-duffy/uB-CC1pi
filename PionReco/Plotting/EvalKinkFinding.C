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
  TTree *tr_haskink = (TTree*)fin->Get("haskink_tree");
  TTree *tr_nokink = (TTree*)fin->Get("nokink_tree");

  TH1F *nkinks_nokink = new TH1F("nkinks_nokink","No kink (truth);Found one or more reco kink(s)?;Area normalised",2,-0.5,1.5);
  TH1F *nkinks_manymcp = new TH1F("nkinks_manymcp","Has kink (truth);Found one or more reco kink(s)?;Area normalised",2,-0.5,1.5);

  TH1F *kinkdistance = new TH1F("kinkdistance",";Min. reco-true kink distance;Area normalised",100,0,500);

  TH1F *nkinks_nokink_byangle = new TH1F("nkinks_nokink_byangle","No kink (truth);Found one or more reco kink(s)?;Area normalised",2,-0.5,1.5);
  TH1F *nkinks_manymcp_byangle = new TH1F("nkinks_manymcp_byangle","Has kink (truth);Found one or more reco kink(s)?;Area normalised",2,-0.5,1.5);

  TH1F *kinkdistance_byangle = new TH1F("kinkdistance_byangle",";Min. reco-true kink distance;Area normalised",100,0,500);

  KinkFindingTree KinkFindingTree_haskink;
  SetKinkFindingTreeVariables(&KinkFindingTree_haskink,tr_haskink);

  KinkFindingTree KinkFindingTree_nokink;
  SetKinkFindingTreeVariables(&KinkFindingTree_nokink,tr_nokink);


  for (int i_evt=0; i_evt<tr_haskink->GetEntries(); i_evt++){
    if (i_evt%50==0) std::cout << "tree haskink: " << i_evt << "/" << tr_haskink->GetEntries() << std::endl;

    tr_haskink->GetEntry(i_evt);

    bool goodevt = KinkFindingTree_haskink.Setup();
    if (!goodevt) continue;

    for (size_t i_pfph=0; i_pfph<KinkFindingTree_haskink.n_PFPs; i_pfph++){
      int pfp_idx=i_pfph;
      int pfpid = KinkFindingTree_haskink.PFP_ID->at(pfp_idx);

      std::vector<int> MCP_vectorposes = KinkFindingTree_haskink.GetGoodRecoMCPsfromPFP(pfpid);


      // For events with multiple MCPs matched to a single PFP, that are in the "has kink" list, and has at least one reco vertex found, how far is the reco vertex from the true one?
      std::vector<std::vector<double>> truekinks;
      for (int i_mcp=0; i_mcp < MCP_vectorposes.size(); i_mcp++){
        int mcp_vectorpos_1 = MCP_vectorposes.at(i_mcp);

        // std::cout << "MCP: " << std::endl;
        // std::cout << " -- start " << "(" << KinkFindingTree_haskink.MCP_StartXYZ->at(mcp_vectorpos_1).at(0) << "," << KinkFindingTree_haskink.MCP_StartXYZ->at(mcp_vectorpos_1).at(1) << "," << KinkFindingTree_haskink.MCP_StartXYZ->at(mcp_vectorpos_1).at(2) << ")" << std::endl;
        // std::cout << " -- end " << "(" << KinkFindingTree_haskink.MCP_EndXYZ->at(mcp_vectorpos_1).at(0) << "," << KinkFindingTree_haskink.MCP_EndXYZ->at(mcp_vectorpos_1).at(1) << "," << KinkFindingTree_haskink.MCP_EndXYZ->at(mcp_vectorpos_1).at(2) << ")" << std::endl;

        // if (KinkFindingTree_haskink.MCP_MotherID->at(mcp_vectorpos_1)==0) {
        //   std::cout << "MCP is neutrino daughter: not using start point" << std::endl;
        //   continue;
        // }

        if (i_mcp == MCP_vectorposes.size()-1) continue;

        // If MCP is not a daughter of the neutrino (i.e. it is the daughter of another MCP), take the MCP start point and assume that is the location of a kink
        TVector3 start_1_xyz(KinkFindingTree_haskink.MCP_StartXYZ->at(mcp_vectorpos_1).at(0),KinkFindingTree_haskink.MCP_StartXYZ->at(mcp_vectorpos_1).at(1),KinkFindingTree_haskink.MCP_StartXYZ->at(mcp_vectorpos_1).at(2));
        TVector3 end_1_xyz(KinkFindingTree_haskink.MCP_EndXYZ->at(mcp_vectorpos_1).at(0),KinkFindingTree_haskink.MCP_EndXYZ->at(mcp_vectorpos_1).at(1),KinkFindingTree_haskink.MCP_EndXYZ->at(mcp_vectorpos_1).at(2));
        TVector3 start_2_xyz(KinkFindingTree_haskink.MCP_StartXYZ->at(mcp_vectorpos_1+1).at(0),KinkFindingTree_haskink.MCP_StartXYZ->at(mcp_vectorpos_1+1).at(1),KinkFindingTree_haskink.MCP_StartXYZ->at(mcp_vectorpos_1+1).at(2));
        TVector3 end_2_xyz(KinkFindingTree_haskink.MCP_EndXYZ->at(mcp_vectorpos_1+1).at(0),KinkFindingTree_haskink.MCP_EndXYZ->at(mcp_vectorpos_1+1).at(1),KinkFindingTree_haskink.MCP_EndXYZ->at(mcp_vectorpos_1+1).at(2));

        if ((start_1_xyz-start_2_xyz).Mag()<1.0 || (start_1_xyz-end_2_xyz).Mag()<1.0){
          std::vector<double> dummy{start_1_xyz.X(),start_1_xyz.Y(),start_1_xyz.Z()};
          truekinks.push_back(dummy);
        }
        if ((end_1_xyz-start_2_xyz).Mag()<1.0 || (end_1_xyz-end_2_xyz).Mag()<1.0){
          std::vector<double> dummy{end_1_xyz.X(),end_1_xyz.Y(),end_1_xyz.Z()};
          truekinks.push_back(dummy);
        }

      }

      // Now calculate local linearity variables to look for kinks
      std::vector<int> kinkidxs = KinkFindingTree_haskink.FindRecoKinks_orderedsp(pfp_idx,-999,i_evt,truekinks);

      // Now look at kinks by angle betwen groups of spacepoints
      std::vector<int> kinkidxs_byangle = KinkFindingTree_haskink.FindRecoKinks_byangle_orderedsp(pfp_idx,-999,i_evt,0,20,truekinks);


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
  for (int i_evt=0; i_evt<tr_nokink->GetEntries(); i_evt++){
    if (i_evt%50==0) std::cout << "tree nokink: " << i_evt << "/" << tr_nokink->GetEntries() << std::endl;

    tr_nokink->GetEntry(i_evt);

    bool goodevt = KinkFindingTree_nokink.Setup();
    if (!goodevt) continue;

    for (size_t i_pfph=0; i_pfph<KinkFindingTree_nokink.n_PFPs; i_pfph++){
      int pfp_idx=i_pfph;

      // Now calculate local linearity variables to look for kinks
      std::vector<int> kinkidxs = KinkFindingTree_nokink.FindRecoKinks_orderedsp(pfp_idx,-999,i_evt*-1);

      // Now look at kinks by angle betwen groups of spacepoints
      std::vector<int> kinkidxs_byangle = KinkFindingTree_nokink.FindRecoKinks_byangle_orderedsp(pfp_idx,-999,i_evt*-1,0,20);


      // Fill plots!
      if (kinkidxs.size()>0) nkinks_nokink->Fill(1);
      else nkinks_nokink->Fill(0);

      if (kinkidxs_byangle.size()>0) nkinks_nokink_byangle->Fill(1);
      else nkinks_nokink_byangle->Fill(0);

    } // end loop i_pfph

  } // end loop i_evt (tree_haskink)

  TFile *fout = new TFile("kinkfindingplots.root","recreate");

  TCanvas *c1 = new TCanvas();

  nkinks_nokink->Scale(1./nkinks_nokink->Integral());
  nkinks_nokink->GetYaxis()->SetRangeUser(0,1.1);
  nkinks_nokink->SetLineWidth(2);
  nkinks_nokink->SetLineColor(kBlack);
  nkinks_nokink->Draw();
  nkinks_manymcp->Scale(1./nkinks_manymcp->Integral());
  nkinks_manymcp->SetLineWidth(2);
  nkinks_manymcp->SetLineStyle(2);
  nkinks_manymcp->SetLineColor(kRed+2);
  nkinks_manymcp->Draw("same");

  TLegend *l = new TLegend(0.7,0.84,0.93,0.97);
  l->SetFillColor(kWhite);
  l->AddEntry(nkinks_nokink,"No kink","l");
  l->AddEntry(nkinks_manymcp,"Multiple MCPs","l");
  l->Draw();

  c1->Print("nkinks_sorted.pdf");

  nkinks_nokink->Write("nkinks_nokink");
  nkinks_manymcp->Write("nkinks_manymcp");
  c1->Write("nkinks_handscanned_canv");


  nkinks_nokink_byangle->Scale(1./nkinks_nokink_byangle->Integral());
  nkinks_nokink_byangle->GetYaxis()->SetRangeUser(0,1.1);
  nkinks_nokink_byangle->SetLineWidth(2);
  nkinks_nokink_byangle->SetLineColor(kBlack);
  nkinks_nokink_byangle->Draw();
  nkinks_manymcp_byangle->Scale(1./nkinks_manymcp_byangle->Integral());
  nkinks_manymcp_byangle->SetLineWidth(2);
  nkinks_manymcp_byangle->SetLineStyle(2);
  nkinks_manymcp_byangle->SetLineColor(kRed+2);
  nkinks_manymcp_byangle->Draw("same");
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

void MakeKinkTrees(std::string inputfile="/uboone/data/users/kduffy/CC1pi/book/v06_26_01_13_PionReco_v1/pionrecovalid_merged.root", std::string outputfile="kinkfindingtrees_out.root"){

  TFile *fin = new TFile(inputfile.c_str(),"OPEN");
  TTree *tr = (TTree*)fin->Get("kinkfinding/KinkFindingTree");

  KinkFindingTree KinkFindingTree;
  //SetupKinkFindingTreeForReading(&KinkFindingTree);
  SetKinkFindingTreeVariables(&KinkFindingTree,tr);

  TFile* outfile = new TFile(outputfile.c_str(),"RECREATE");

  TTree *tree_kink = new TTree();
  struct KinkFindingTree savevars_kink;
  SetKinkFindingTreeBranches(&savevars_kink,tree_kink);
  std::vector<std::vector<std::vector<double>>> *PFP_ordered_spacepoints_kink=nullptr;
  std::vector<std::vector<int>> *PFP_ordered_spacepoints_origidx_kink=nullptr;
  tree_kink->Branch("PFP_ordered_spacepoints",&PFP_ordered_spacepoints_kink);
  tree_kink->Branch("PFP_ordered_spacepoints_origidx",&PFP_ordered_spacepoints_origidx_kink);

  TTree *tree_nokink = new TTree();
  struct KinkFindingTree savevars_nokink;
  SetKinkFindingTreeBranches(&savevars_nokink,tree_nokink);
  std::vector<std::vector<std::vector<double>>> *PFP_ordered_spacepoints_nokink=nullptr;
  std::vector<std::vector<int>> *PFP_ordered_spacepoints_origidx_nokink=nullptr;
  tree_nokink->Branch("PFP_ordered_spacepoints",&PFP_ordered_spacepoints_nokink);
  tree_nokink->Branch("PFP_ordered_spacepoints_origidx",&PFP_ordered_spacepoints_origidx_nokink);

  int n_nokink = 0;
  int n_kink = 0;

  for (int i_evt=300; i_evt</*tr->GetEntries()*/400; i_evt++){
    if (i_evt%20==0) std::cout << i_evt << "/" << tr->GetEntries() << std::endl;

    if (n_nokink>200 && n_kink>200) return;

    tr->GetEntry(i_evt);
    ClearKinkFindingTree(&savevars_nokink);
    ClearKinkFindingTree(&savevars_kink);
    PFP_ordered_spacepoints_kink->clear();
    PFP_ordered_spacepoints_origidx_kink->clear();
    PFP_ordered_spacepoints_nokink->clear();
    PFP_ordered_spacepoints_origidx_nokink->clear();
    CopyMCPs(&KinkFindingTree,&savevars_nokink);
    CopyMCPs(&KinkFindingTree,&savevars_kink);

    bool goodevt = KinkFindingTree.Setup();
    if (!goodevt) continue;

    KinkFindingTree.GetPFPHierarchy();

    for (size_t i_pfph=0; i_pfph<KinkFindingTree.bestmatch_PFP_hierarchy_idxs.size(); i_pfph++){
      int pfp_idx=KinkFindingTree.bestmatch_PFP_hierarchy_idxs.at(i_pfph);
      int pfpid = KinkFindingTree.PFP_ID->at(pfp_idx);

      std::vector<int> MCP_vectorposes = KinkFindingTree.GetGoodRecoMCPsfromPFP(pfpid);

      // If MCP_vectorposes has size 0, move on
      if (MCP_vectorposes.size()==0){
        continue;
      }
      else if (MCP_vectorposes.size()==1 && n_nokink<=200){
        CopyPFP(pfp_idx,&KinkFindingTree,&savevars_nokink);

        std::vector<std::pair<std::vector<double>,int>> ordered_spacepoints_pair_origidx =  OrderSpacepoints(KinkFindingTree.PFP_spacepoints_XYZ->at(pfp_idx));
        std::vector<std::vector<double>> ordered_spacepoints;
        std::vector<int> ordered_spacepoints_origidx;
        for (int i_sp=0; i_sp<ordered_spacepoints_pair_origidx.size(); i_sp++){
          ordered_spacepoints.push_back(ordered_spacepoints_pair_origidx.at(i_sp).first);
          ordered_spacepoints_origidx.push_back(ordered_spacepoints_pair_origidx.at(i_sp).second);
        }
        PFP_ordered_spacepoints_nokink->push_back(ordered_spacepoints);
        PFP_ordered_spacepoints_origidx_nokink->push_back(ordered_spacepoints_origidx);

        n_nokink++;
      }
      else if (MCP_vectorposes.size()>1 && n_kink<=200){
        CopyPFP(pfp_idx,&KinkFindingTree,&savevars_kink);

        std::vector<std::pair<std::vector<double>,int>> ordered_spacepoints_pair_origidx =  OrderSpacepoints(KinkFindingTree.PFP_spacepoints_XYZ->at(pfp_idx));
        std::vector<std::vector<double>> ordered_spacepoints;
        std::vector<int> ordered_spacepoints_origidx;
        for (int i_sp=0; i_sp<ordered_spacepoints_pair_origidx.size(); i_sp++){
          ordered_spacepoints.push_back(ordered_spacepoints_pair_origidx.at(i_sp).first);
          ordered_spacepoints_origidx.push_back(ordered_spacepoints_pair_origidx.at(i_sp).second);
        }
        PFP_ordered_spacepoints_kink->push_back(ordered_spacepoints);
        PFP_ordered_spacepoints_origidx_kink->push_back(ordered_spacepoints_origidx);

        n_kink++;
      }
    } // end loop i_pfph

    // Calculate nPFPs and fill trees
    savevars_nokink.n_PFPs = savevars_nokink.PFP_ID->size();
    savevars_kink.n_PFPs = savevars_kink.PFP_ID->size();

    if (i_evt%100==0) std::cout << n_nokink << " no kinks, " << n_kink << " with kinks" << std::endl;

    // std::cout << savevars_nokink.n_PFPs << ", " << savevars_kink.n_PFPs << std::endl;

    if (savevars_nokink.n_PFPs>0) tree_nokink->Fill();
    if (savevars_kink.n_PFPs>0) tree_kink->Fill();

  } // end loop i_evt

  // std::cout << n_nokink << " no kinks, " << n_kink << " with kinks" << std::endl;

  outfile->cd();
  tree_kink->Write("haskink_tree");
  tree_nokink->Write("nokink_tree");
  outfile->Close();
}

void CheckSpacePointOrdering(std::string inputfile="/uboone/data/users/kduffy/CC1pi/book/v06_26_01_13_PionReco_v1/pionrecovalid_merged.root", std::string outputdir="./", int nplots=20){
  gStyle->SetOptStat(0);
  gStyle->SetOptTitle(0);

  TFile *fin = new TFile(inputfile.c_str(),"OPEN");
  TTree *tr = (TTree*)fin->Get("haskink_tree");

  KinkFindingTree KinkFindingTree;
  SetKinkFindingTreeVariables(&KinkFindingTree,tr);

  TCanvas *c1 = new TCanvas();
  c1->Divide(2);

  int plots_made = 0;
  for (int i_evt=0; i_evt<tr->GetEntries(); i_evt++){

    tr->GetEntry(i_evt);

    // Fill spacepoints plots
    double minx = 9999;
    double maxx = -9999;
    double miny = 9999;
    double maxy = -9999;
    double minz = 9999;
    double maxz = -9999;

     for (size_t i_pfph=0; i_pfph<KinkFindingTree.n_PFPs; i_pfph++){

         std::vector<std::vector<double>> ordered_spacepoints = KinkFindingTree.PFP_ordered_spacepoints->at(i_pfph);

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

     for (size_t i_pfph=0; i_pfph<KinkFindingTree.n_PFPs; i_pfph++){

         std::vector<std::vector<double>> ordered_spacepoints = KinkFindingTree.PFP_ordered_spacepoints->at(i_pfph);

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

      if (spacepoints_yz->GetEntries()==0) continue;

      c1->cd(1);
      spacepoints_yz->Draw("colz");
      c1->cd(2);
      spacepoints_xz->Draw("colz");
      TString savename = TString::Format("%s/evt_%d_pfp_%d_checksp_pandorasort.pdf",outputdir.c_str(),i_evt,(int)i_pfph);
      c1->Print(savename.Data());

      plots_made++;
      if (plots_made>=nplots) break;

    }
  } // end loop over events in tree (i_evt)
}
