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

  std::string line;
  std::vector<int> nokink_evts;
  std::vector<int> nokink_pfpidxs;
  std::vector<int> greyarea_evts;
  std::vector<int> greyarea_pfpidxs;
  std::vector<int> haskink_evts;
  std::vector<int> haskink_pfpidxs;

  ifstream infile_nokink;
	infile_nokink.open ("/uboone/app/users/kduffy/CC1pi/CC1pi_uboonecode/srcs/uboonecode/uboone/CC1pi/ignore/2018-11-20/nokink.txt");
	if (infile_nokink.is_open())
	{
		while (getline(infile_nokink,line)){
      if (line == std::string("Event PFP")) continue;
      std::istringstream iss (line);
      int val;
      iss >> val;
      nokink_evts.push_back(val);
      iss >> val;
      nokink_pfpidxs.push_back(val);
    }
		infile_nokink.close();
	}
	else
	{
		cout << "Error opening file infile_nokink";
	}

  ifstream infile_greyarea;
	infile_greyarea.open ("/uboone/app/users/kduffy/CC1pi/CC1pi_uboonecode/srcs/uboonecode/uboone/CC1pi/ignore/2018-11-20/greyarea.txt");
	if (infile_greyarea.is_open())
	{
		while (getline(infile_greyarea,line)){
      if (line == std::string("Event PFP")) continue;
      std::istringstream iss (line);
      int val;
      iss >> val;
      greyarea_evts.push_back(val);
      iss >> val;
      greyarea_pfpidxs.push_back(val);
    }
		infile_greyarea.close();
	}
	else
	{
		cout << "Error opening file infile_nokink";
	}

  ifstream infile_haskink;
	infile_haskink.open ("/uboone/app/users/kduffy/CC1pi/CC1pi_uboonecode/srcs/uboonecode/uboone/CC1pi/ignore/2018-11-20/haskink.txt");
	if (infile_haskink.is_open())
	{
		while (getline(infile_haskink,line)){
      if (line == std::string("Event PFP")) continue;
      std::istringstream iss (line);
      int val;
      iss >> val;
      haskink_evts.push_back(val);
      iss >> val;
      haskink_pfpidxs.push_back(val);
    }
		infile_haskink.close();
	}
	else
	{
		cout << "Error opening file infile_haskink";
	}

  TFile *fin = new TFile(inputfile.c_str(),"OPEN");
  TTree *tr = (TTree*)fin->Get("kinkfinding/KinkFindingTree");

  TH2F *nkinks_vs_mincosth = new TH2F("nkinks_vs_mincosth",";No. reco kinks;True min(cos(#theta))",5,0,5,100,-1,1);
  TH2F *nkinks_vs_mincosth_manymcp = new TH2F("nkinks_vs_mincosth_manymcp",";No. reco kinks;True min(cos(#theta))",5,0,5,100,-1,1);

  TH1F *nkinks_nokink = new TH1F("nkinks_nokink","No kink (truth);Found one or more reco kink(s)?;Area normalised",2,-0.5,1.5);
  TH1F *nkinks_haskink = new TH1F("nkinks_haskink","Has kink (truth);Found one or more reco kink(s)?;Area normalised",2,-0.5,1.5);
  TH1F *nkinks_greyarea = new TH1F("nkinks_greyarea","Grey area (truth);Found one or more reco kink(s)?;Area normalised",2,-0.5,1.5);
  TH1F *nkinks_manymcp = new TH1F("nkinks_manymcp","Multiple GEANT particles (truth);Found one or more reco kink(s)?;Area normalised",2,-0.5,1.5);

  TH1F *kinkdistance = new TH1F("kinkdistance",";Min. reco-true kink distance;Area normalised",100,0,500);

  TH1F *nkinks_nokink_byangle = new TH1F("nkinks_nokink_byangle","No kink (truth);Found one or more reco kink(s)?;Area normalised",2,-0.5,1.5);
  TH1F *nkinks_haskink_byangle = new TH1F("nkinks_haskink_byangle","Has kink (truth);Found one or more reco kink(s)?;Area normalised",2,-0.5,1.5);
  TH1F *nkinks_greyarea_byangle = new TH1F("nkinks_greyarea_byangle","Grey area (truth);Found one or more reco kink(s)?;Area normalised",2,-0.5,1.5);
  TH1F *nkinks_manymcp_byangle = new TH1F("nkinks_manymcp_byangle","Multiple GEANT particles (truth);Found one or more reco kink(s)?;Area normalised",2,-0.5,1.5);

  TH1F *kinkdistance_byangle = new TH1F("kinkdistance_byangle",";Min. reco-true kink distance;Area normalised",100,0,500);

  KinkFindingTree KinkFindingTree;
  SetKinkFindingTreeVariables(&KinkFindingTree,tr);


  for (int i_evt=0; i_evt<tr->GetEntries(); i_evt++){

    tr->GetEntry(i_evt);

    bool goodevt = KinkFindingTree.Setup();
    if (!goodevt) continue;

    //KinkFindingTree.CalcLogicalPionVars();
    KinkFindingTree.GetPFPHierarchy();

    // std::cout << "Reconstructed " << KinkFindingTree.bestmatch_PFP_hierarchy_idxs.size() << " PFPs in hierarchy" << std::endl;

    for (size_t i_pfph=0; i_pfph<KinkFindingTree.bestmatch_PFP_hierarchy_idxs.size(); i_pfph++){
      int pfp_idx=KinkFindingTree.bestmatch_PFP_hierarchy_idxs.at(i_pfph);
      int pfpid = KinkFindingTree.PFP_ID->at(pfp_idx);

      std::vector<int> MCP_vectorposes = KinkFindingTree.GetGoodRecoMCPsfromPFP(pfpid);

      // std::cout << "- PFP " << i_pfph << ": ";
      // if (MCP_vectorposes.size()==0) std::cout << "Not well reconstructed" << std::endl;
      // else {
      //   std::cout << MCP_vectorposes.size() << " matched MCP(s)" << std::endl;
      //   for (int i_m=0; i_m<MCP_vectorposes.size(); i_m++){
      //     std::cout << KinkFindingTree.MCP_PDGcode->at(MCP_vectorposes.at(i_m)) << "(" <<  KinkFindingTree.MCP_totaldepE->at(MCP_vectorposes.at(i_m)) << ") ";
      //   }
      //   std::cout << std::endl;
      // }

      // If MCP_vectorposes has size 0, move on
      if (MCP_vectorposes.size()==0) continue;

      // Now ask if the track has a kink
      // std::cout << "Track has true kink?" << std::endl;
      double mincosth = KinkFindingTree.IsTrueKink(MCP_vectorposes,0.1);
      // std::cout << istruekink << std::endl;

      // Now calculate local linearity variables to look for kinks
      std::vector<int> kinkidxs = KinkFindingTree.FindRecoKinks(pfp_idx,mincosth,i_evt);

      // Now look at kinks by angle betwen groups of spacepoints
      std::vector<int> kinkidxs_byangle = KinkFindingTree.FindRecoKinks_byangle(pfp_idx,mincosth,i_evt,0,20);


      // Fill plots!

      if (MCP_vectorposes.size()==1) nkinks_vs_mincosth->Fill(kinkidxs.size(),mincosth);
      else nkinks_vs_mincosth_manymcp->Fill(kinkidxs.size(),mincosth+2.0);

      if (MCP_vectorposes.size()==1/*std::find(nokink_evts.begin(), nokink_evts.end(),i_evt)!=nokink_evts.end() && std::find(nokink_pfpidxs.begin(), nokink_pfpidxs.end(),pfp_idx)!=nokink_pfpidxs.end()*/){
        // std::cout << "No kink: " << kinkidxs.size() << std::endl;
        if (kinkidxs.size()>0) nkinks_nokink->Fill(1);
        else nkinks_nokink->Fill(0);

        if (kinkidxs_byangle.size()>0) nkinks_nokink_byangle->Fill(1);
        else nkinks_nokink_byangle->Fill(0);
      }
      else if (std::find(haskink_evts.begin(), haskink_evts.end(),i_evt)!=haskink_evts.end() && std::find(haskink_pfpidxs.begin(), haskink_pfpidxs.end(),pfp_idx)!=haskink_pfpidxs.end()){
        // std::cout << "Has kink: " << kinkidxs.size() << std::endl;
        if (kinkidxs.size()>0) nkinks_haskink->Fill(1);
        else nkinks_haskink->Fill(0);

        if (kinkidxs_byangle.size()>0) nkinks_haskink_byangle->Fill(1);
        else nkinks_haskink_byangle->Fill(0);
      }
      else if (std::find(greyarea_evts.begin(), greyarea_evts.end(),i_evt)!=greyarea_evts.end() && std::find(greyarea_pfpidxs.begin(), greyarea_pfpidxs.end(),pfp_idx)!=greyarea_pfpidxs.end()){
        // std::cout << "Grey area kink: " << kinkidxs.size() << std::endl;
        if (kinkidxs.size()>0) nkinks_greyarea->Fill(1);
        else nkinks_greyarea->Fill(0);

        if (kinkidxs_byangle.size()>0) nkinks_greyarea_byangle->Fill(1);
        else nkinks_greyarea_byangle->Fill(0);
      }

      if (MCP_vectorposes.size()>1 /*&& std::find(haskink_evts.begin(), haskink_evts.end(),i_evt)!=haskink_evts.end() && std::find(haskink_pfpidxs.begin(), haskink_pfpidxs.end(),pfp_idx)!=haskink_pfpidxs.end()*/){
        // std::cout << "Many MCPs and has kink: " << kinkidxs.size() << std::endl;
        if (kinkidxs.size()>0) nkinks_manymcp->Fill(1);
        else nkinks_manymcp->Fill(0);

        if (kinkidxs_byangle.size()>0) nkinks_manymcp_byangle->Fill(1);
        else nkinks_manymcp_byangle->Fill(0);
      }


      // For events with multiple MCPs matched to a single PFP, that are in the "has kink" list, and has at least one reco vertex found, how far is the reco vertex from the true one?
      if (MCP_vectorposes.size()>1 /*&& std::find(haskink_evts.begin(), haskink_evts.end(),i_evt)!=haskink_evts.end() && std::find(haskink_pfpidxs.begin(), haskink_pfpidxs.end(),pfp_idx)!=haskink_pfpidxs.end()*/){
        // Find true vertex (loop through all MCPs and find all potential true vertices)
        std::vector<std::vector<double>> truekinks;
        for (int i_mcp=0; i_mcp < MCP_vectorposes.size(); i_mcp++){
          int mcp_vectorpos_1 = MCP_vectorposes.at(i_mcp);

          // std::cout << "MCP: " << std::endl;
          // std::cout << " -- start " << "(" << KinkFindingTree.MCP_StartXYZ->at(mcp_vectorpos_1).at(0) << "," << KinkFindingTree.MCP_StartXYZ->at(mcp_vectorpos_1).at(1) << "," << KinkFindingTree.MCP_StartXYZ->at(mcp_vectorpos_1).at(2) << ")" << std::endl;
          // std::cout << " -- end " << "(" << KinkFindingTree.MCP_EndXYZ->at(mcp_vectorpos_1).at(0) << "," << KinkFindingTree.MCP_EndXYZ->at(mcp_vectorpos_1).at(1) << "," << KinkFindingTree.MCP_EndXYZ->at(mcp_vectorpos_1).at(2) << ")" << std::endl;

          // if (KinkFindingTree.MCP_MotherID->at(mcp_vectorpos_1)==0) {
          //   std::cout << "MCP is neutrino daughter: not using start point" << std::endl;
          //   continue;
          // }

          if (i_mcp == MCP_vectorposes.size()-1) continue;

          // If MCP is not a daughter of the neutrino (i.e. it is the daughter of another MCP), take the MCP start point and assume that is the location of a kink
          TVector3 start_1_xyz(KinkFindingTree.MCP_StartXYZ->at(mcp_vectorpos_1).at(0),KinkFindingTree.MCP_StartXYZ->at(mcp_vectorpos_1).at(1),KinkFindingTree.MCP_StartXYZ->at(mcp_vectorpos_1).at(2));
          TVector3 end_1_xyz(KinkFindingTree.MCP_EndXYZ->at(mcp_vectorpos_1).at(0),KinkFindingTree.MCP_EndXYZ->at(mcp_vectorpos_1).at(1),KinkFindingTree.MCP_EndXYZ->at(mcp_vectorpos_1).at(2));
          TVector3 start_2_xyz(KinkFindingTree.MCP_StartXYZ->at(mcp_vectorpos_1+1).at(0),KinkFindingTree.MCP_StartXYZ->at(mcp_vectorpos_1+1).at(1),KinkFindingTree.MCP_StartXYZ->at(mcp_vectorpos_1+1).at(2));
          TVector3 end_2_xyz(KinkFindingTree.MCP_EndXYZ->at(mcp_vectorpos_1+1).at(0),KinkFindingTree.MCP_EndXYZ->at(mcp_vectorpos_1+1).at(1),KinkFindingTree.MCP_EndXYZ->at(mcp_vectorpos_1+1).at(2));

          if ((start_1_xyz-start_2_xyz).Mag()<1.0 || (start_1_xyz-end_2_xyz).Mag()<1.0){
            std::vector<double> dummy{start_1_xyz.X(),start_1_xyz.Y(),start_1_xyz.Z()};
            truekinks.push_back(dummy);
          }
          if ((end_1_xyz-start_2_xyz).Mag()<1.0 || (end_1_xyz-end_2_xyz).Mag()<1.0){
            std::vector<double> dummy{end_1_xyz.X(),end_1_xyz.Y(),end_1_xyz.Z()};
            truekinks.push_back(dummy);
          }

        }

        for (int i_truekink=0; i_truekink<truekinks.size(); i_truekink++){
          std::vector<double> truekink_xyz = truekinks.at(i_truekink);

          // Now loop through all reco kinks and find the closest one
          double dist = 9999999;
          for (int i_r=0; i_r<kinkidxs.size(); i_r++){
            std::vector<double> recokink_xyz = KinkFindingTree.PFP_spacepoints_XYZ->at(pfp_idx).at(kinkidxs.at(i_r));
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
            std::vector<double> recokink_xyz = KinkFindingTree.PFP_spacepoints_XYZ->at(pfp_idx).at(kinkidxs_byangle.at(i_r));
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

      }

    } // end loop i_pfph

  } // end loop i_evt

  TCanvas *c1 = new TCanvas();
  nkinks_vs_mincosth->Draw("colz");
  c1->Print("nkinks_vs_mincosth.pdf");

  nkinks_vs_mincosth_manymcp->Draw("colz");
  c1->Print("nkinks_vs_mincosth_manymcp.pdf");

  nkinks_nokink->Scale(1./nkinks_nokink->Integral());
  nkinks_nokink->GetYaxis()->SetRangeUser(0,1.1);
  nkinks_nokink->SetLineWidth(2);
  nkinks_nokink->SetLineColor(kBlack);
  nkinks_nokink->Draw();
  nkinks_greyarea->Scale(1./nkinks_greyarea->Integral());
  nkinks_greyarea->SetLineWidth(2);
  nkinks_greyarea->SetLineColor(kGray+1);
  nkinks_greyarea->SetLineStyle(2);
  nkinks_greyarea->Draw("same");
  nkinks_haskink->Scale(1./nkinks_haskink->Integral());
  nkinks_haskink->SetLineWidth(2);
  nkinks_haskink->SetLineColor(kRed);
  nkinks_haskink->Draw("same");
  nkinks_manymcp->Scale(1./nkinks_manymcp->Integral());
  nkinks_manymcp->SetLineWidth(2);
  nkinks_manymcp->SetLineStyle(2);
  nkinks_manymcp->SetLineColor(kRed+2);
  nkinks_manymcp->Draw("same");

  TLegend *l = new TLegend(0.7,0.84,0.93,0.97);
  l->SetFillColor(kWhite);
  l->AddEntry(nkinks_nokink,"No kink","l");
  l->AddEntry(nkinks_haskink,"Has kink","l");
  l->AddEntry(nkinks_greyarea,"Grey area","l");
  l->AddEntry(nkinks_manymcp,"Multiple MCPs","l");
  l->Draw();

  c1->Print("nkinks_handscanned.pdf");



  nkinks_nokink_byangle->Scale(1./nkinks_nokink_byangle->Integral());
  nkinks_nokink_byangle->GetYaxis()->SetRangeUser(0,1.1);
  nkinks_nokink_byangle->SetLineWidth(2);
  nkinks_nokink_byangle->SetLineColor(kBlack);
  nkinks_nokink_byangle->Draw();
  nkinks_greyarea_byangle->Scale(1./nkinks_greyarea_byangle->Integral());
  nkinks_greyarea_byangle->SetLineWidth(2);
  nkinks_greyarea_byangle->SetLineColor(kGray+1);
  nkinks_greyarea_byangle->SetLineStyle(2);
  nkinks_greyarea_byangle->Draw("same");
  nkinks_haskink_byangle->Scale(1./nkinks_haskink_byangle->Integral());
  nkinks_haskink_byangle->SetLineWidth(2);
  nkinks_haskink_byangle->SetLineColor(kRed);
  nkinks_haskink_byangle->Draw("same");
  nkinks_manymcp_byangle->Scale(1./nkinks_manymcp_byangle->Integral());
  nkinks_manymcp_byangle->SetLineWidth(2);
  nkinks_manymcp_byangle->SetLineStyle(2);
  nkinks_manymcp_byangle->SetLineColor(kRed+2);
  nkinks_manymcp_byangle->Draw("same");
  l->Draw();
  c1->Print("nkinks_handscanned_byangle.pdf");

  kinkdistance->SetLineWidth(2);
  kinkdistance->SetLineStyle(2);
  kinkdistance->SetLineColor(kRed+2);
  kinkdistance->DrawNormalized();
  c1->Print("kinkdist_manymcp_haskink.pdf");

  kinkdistance_byangle->SetLineWidth(2);
  kinkdistance_byangle->SetLineStyle(2);
  kinkdistance_byangle->SetLineColor(kRed+2);
  kinkdistance_byangle->DrawNormalized();
  c1->Print("kinkdist_byangle_manymcp_haskink.pdf");
}
