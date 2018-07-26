#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TLegend.h"
#include "TAxis.h"
#include "TStyle.h"

#include <iostream>

#include "Marco_dQdxCut.h"

void PlotdQdxLength(std::string infile)
{
  // Get tree from input file
  // (Assumes you've already hadded the files -- much faster than TChaining anyway
  TFile *fin = new TFile(infile.c_str(),"read");
  TTree *cc1pitree = (TTree*)fin->Get("cc1piselec/outtree");
  
  std::vector<double> *track_length = nullptr;
  std::vector<double> *dqdx_trunc_uncalib = nullptr;
  std::vector<bool> *MIPConsistency = nullptr;
  bool *isSelected;

  cc1pitree->SetBranchAddress("track_length", &track_length);
  cc1pitree->SetBranchAddress("dqdx_trunc_uncalib", &dqdx_trunc_uncalib);
  cc1pitree->SetBranchAddress("MIPConsistency", &MIPConsistency);
  cc1pitree->SetBranchAddress("isSelected",&isSelected);
  
  // Vectors to store things for plotting (since we don't know the size a priori)
  std::vector<double> trlen_MIP, trlen_notMIP, dqdx_MIP, dqdx_notMIP;
  
  for (int i_event=0; i_event<cc1pitree->GetEntries(); i_event++){
    //std::cout << "Testtree: Entry " << i_event << std::endl;
    // If the event is not selected, nothing else will be filled and we'll get segfaults
    // Skip these events
    if (!isSelected) continue;
    
    // If the event is selected (as tested by testtree, continue on and evaluate things)
    cc1pitree->GetEntry(i_event);
    //std::cout << "CC1pitree: Entry " << i_event << std::endl;
    
    // Check the three vectors are the same size (i.e. each has one entry per track)
    if (!((track_length->size() == dqdx_trunc_uncalib->size()) && (track_length->size() == MIPConsistency->size()))){
      std::cout << "[ERROR] track_length->size()       = " << track_length->size() << std::endl
		<< "        dqdx_trunc_uncalib->size() = " << dqdx_trunc_uncalib->size() << std::endl
		<< "        MIPConsistency->size()     = " << MIPConsistency->size() << std::endl;
    }

    // Loop through tracks (entries in vector) and push back points in vectors
    for (unsigned int i_tr=0; i_tr < track_length->size(); i_tr++){
      double tracklength_itr = track_length->at(i_tr);
      double dqdx_itr = dqdx_trunc_uncalib->at(i_tr);
      bool MIPConsistent_itr = MIPConsistency->at(i_tr);

      if (MIPConsistent_itr){
	trlen_MIP.push_back(tracklength_itr);
	dqdx_MIP.push_back(dqdx_itr);
      }
      else if (!MIPConsistent_itr){
	trlen_notMIP.push_back(tracklength_itr);
	dqdx_notMIP.push_back(dqdx_itr);
      }
      else std::cout << "What is happening? MIPConsistent_itr = " << MIPConsistent_itr << std::endl;
    }// loop over tracks (i_tr)
    
  }// loop over events (entries in cc1pitree)

  // Now make TGraphs and plot!
  TGraph *g_MIP = new TGraph(trlen_MIP.size());
  TGraph *g_notMIP = new TGraph(trlen_notMIP.size());

  for (unsigned int i_MIP=0; i_MIP < trlen_MIP.size(); i_MIP++){
    g_MIP->SetPoint(i_MIP, dqdx_MIP.at(i_MIP), trlen_MIP.at(i_MIP));
  }// i_MIP (fill g_MIP)

  for (unsigned int i_notMIP=0; i_notMIP < trlen_notMIP.size(); i_notMIP++){
    g_notMIP->SetPoint(i_notMIP, dqdx_notMIP.at(i_notMIP), trlen_notMIP.at(i_notMIP));
  }// i_notMIP (fill g_notMIP)

  gStyle->SetOptTitle(0);
  
  g_MIP->SetMarkerStyle(20);
  g_MIP->SetMarkerSize(0.4);
  g_MIP->SetMarkerColor(kRed);
  g_notMIP->SetMarkerStyle(20);
  g_notMIP->SetMarkerSize(0.4);
  g_notMIP->SetMarkerColor(kBlue);

  TLegend *leg = new TLegend(0.6,0.8,0.95,0.95);
  leg->SetFillColor(kWhite);
  leg->SetTextFont(132);
  leg->AddEntry(g_MIP,"Classified as MIP-like","p");
  leg->AddEntry(g_notMIP,"Classified as not MIP-like","p");
  
  g_MIP->GetXaxis()->SetRangeUser(0,75e6);
  g_MIP->GetXaxis()->SetTitle("Track <dQ/dx>(trunc.)");
  g_MIP->GetYaxis()->SetTitle("Track Length (cm)");

  // TGraph to show Marco's cut
  TGraph *g_cut = new TGraph(_dqdx_cutvals.size());
  for (unsigned int i_cutval=0; i_cutval < _dqdx_cutvals.size(); i_cutval++){
    g_cut->SetPoint(i_cutval,_dqdx_cutvals.at(i_cutval),i_cutval);
  }

  g_MIP->Draw("AP");
  g_notMIP->Draw("P same");
  g_cut->Draw("L same");
  leg->Draw();

}
