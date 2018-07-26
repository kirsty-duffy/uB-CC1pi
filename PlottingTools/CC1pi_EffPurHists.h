#ifndef __CC1PI_EFFPUR_H__
#define __CC1PI_EFFPUR_H__

#include "../Algorithms/TopologyEnums.h"
#include "TH1.h"
#include "TColor.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TStyle.h"
// #include <vector>

#include "CC1pi_treevars.h"
#include "CC1pi_cuts.h"


// --------------------------------------------------- //
// This struct contains signal vs background histograms and all functions related to them

struct histCC1piselEffPur{
   TH1D *h_cc1pi_sel;
   TH1D *h_bg_sel;
   TH1D *h_cc1pi_notsel;
   // TH1D *h_bg_notsel;

   TLegend *l;

   // Constructor for this struct of hists
   histCC1piselEffPur(std::string name, std::string title, double nbins, double binlow, double binhigh){
      h_cc1pi_sel = new TH1D(std::string(name+"_cc1pi_sel").c_str(),title.c_str(),nbins,binlow,binhigh);
      h_bg_sel = new TH1D(std::string(name+"_bg_sel").c_str(),title.c_str(),nbins,binlow,binhigh);
      h_cc1pi_notsel = new TH1D(std::string(name+"_cc1pi_notsel").c_str(),title.c_str(),nbins,binlow,binhigh);
      // h_bg_notsel = new TH1D(std::string(name+"_bg_notsel").c_str(),title.c_str(),nbins,binlow,binhigh);


      // h_cc1pi_sel->SetFillColor(TColor::GetColor(8,64,129));
      // h_cc1pi_all->SetFillColor(TColor::GetColor(8,64,129));
      // h_bg_sel->SetFillColor(TColor::GetColor(197,197,197));
      // h_bg_all->SetFillColor(TColor::GetColor(197,197,197));
      //
      // l = new TLegend(0.59,0.64,0.81,0.87);
      // l->AddEntry(h_p,"True cc1#pi^{+}","f");
      // l->AddEntry(h_mu,"True other","f");
   }

   void SetBinLabel(int ibin, std::string label){
      h_cc1pi_sel -> GetXaxis() -> SetBinLabel(ibin, label.c_str());
      h_bg_sel -> GetXaxis() -> SetBinLabel(ibin, label.c_str());
      h_cc1pi_notsel -> GetXaxis() -> SetBinLabel(ibin, label.c_str());
   }
};

void FillCC1piEffPurHist(histCC1piselEffPur *hists, treevars *vars, int i_cut){
   // Loop through all bins in the histograms, and evaluate a cut value at the centre of each bin
   for (int i_bin=1; i_bin < hists->h_cc1pi_sel->GetXaxis()->GetNbins()+1; i_bin++){
      double cutval = hists->h_cc1pi_sel->GetXaxis()->GetBinCenter(i_bin);

      bool isSelected = IsEventSelected_SingleCut(cutval, vars, i_cut);

      // Fill selection info into the relevant histograms
      if (isSelected){
         if (vars->Truth_topology == kCC1piplus0p || vars->Truth_topology == kCC1piplus1p || vars->Truth_topology == kCC1piplusNp){
            hists->h_cc1pi_sel->Fill(cutval);
         }
         else {
            hists->h_bg_sel->Fill(cutval);
         }
      }
      else { // if not selected
         if (vars->Truth_topology == kCC1piplus0p || vars->Truth_topology == kCC1piplus1p || vars->Truth_topology == kCC1piplusNp){
            hists->h_cc1pi_notsel->Fill(cutval);
         }
         // else {
         //   hists->h_bg_notsel->Fill(cutval);
         // }
      }
   } // End loop over bins in the histograms
}

void FillNminus1EffPurHist(histCC1piselEffPur *hists, treevars *vars){
   int ncuts = GetCutVars(vars).size();
   // Loop through all bins in the histogram and evaluate the cuts
   // Each bin represents one cut (note: that means that for a given bin, we apply all cuts EXCEPT that one), with the final bin showing what we get when we apply all cuts.
   // "Applying" cuts means using the cut values given in the vector
   for (int i_bin=1; i_bin < hists->h_cc1pi_sel->GetXaxis()->GetNbins()+1; i_bin++){
      // Now loop over cuts and apply them
      bool isSelected = true;
      if (i_bin<ncuts+2){
         for (int i_cut=0; i_cut < ncuts; i_cut++){
            // Skip cut associated to current bin
            // For bin with i_bin=ncuts+2, it will not be equal to i_cut+1 for any cut, so it should just evaluate all of them
            if (i_cut+1 == i_bin) {
               // std::cout << "skipping cut" << std::endl;
               continue;
            }
            // Evaluate whether event passes this cut. If it doesn't, set isSelected to false for the event
            bool isSelected_i = IsEventSelected_SingleCut(GetCutVars(vars).at(i_cut).CutValue, vars, i_cut);
            if (isSelected_i == false) isSelected = false;

            // if (isSelected) std::cout << "bin: " << i_bin << ", cut: " << i_cut << "isSelected = " << isSelected << std::endl;
         } // end loop over cuts
      }
      // Bin ncuts+2: only evaluate whether the event passed Marco's selection and has 2 tracks
      else if (i_bin==ncuts+2){
         int n_tracks = 0;
         for (size_t i_track=0; i_track<GetCutVars(vars).at(0).Var->size(); i_track++){
            double value = GetCutVars(vars).at(0).Var->at(i_track);

            // Ignore PFPs that failed to reco as tracks. The neutrino PFP will also have a bogus value.
            if (value == -9999 || value == -999) {
               continue;
            }

            n_tracks++;
         } // end loop over tracks in TPCObject
         if (n_tracks>=2 && vars->Marco_selected){
            isSelected = true;
         }
         else{
            isSelected = false;
         }
      }
      // Bin ncuts+3: only evaluate whether the event passed Marco's selection
      if (i_bin==ncuts+3){
         if (vars->Marco_selected){
            isSelected=true;
         }
         else{
            isSelected=false;
         }
      }

      // Now fill selection info into the relevant histograms
      double fillval = hists->h_cc1pi_sel->GetBinCenter(i_bin);
      //if (isSelected) std::cout << "Bin " << i_bin << ", isSelected = " << isSelected << ", fillval = " << fillval << std::endl;
      if (isSelected){
         if (vars->Truth_topology == kCC1piplus0p || vars->Truth_topology == kCC1piplus1p || vars->Truth_topology == kCC1piplusNp){
            hists->h_cc1pi_sel->Fill(fillval);
         }
         else {
            hists->h_bg_sel->Fill(fillval);
         }
      }
      else { // if not selected
         if (vars->Truth_topology == kCC1piplus0p || vars->Truth_topology == kCC1piplus1p || vars->Truth_topology == kCC1piplusNp){
            hists->h_cc1pi_notsel->Fill(fillval);
         }
         // else {
         //   hists->h_bg_notsel->Fill(fillval);
         // }
      }

   } // end loop over bins (i.e. cuts)

}

void DrawCC1piMCEffPur(TCanvas *c, histCC1piselEffPur *hists, std::string drawopt="p",bool isNminus1=false){
   TH1D *heff = (TH1D*)hists->h_cc1pi_sel->Clone("heff");
   TH1D *hpur = (TH1D*)hists->h_cc1pi_sel->Clone("hpur");
   TH1D *heffpur = (TH1D*)hists->h_cc1pi_sel->Clone("heffpur");

   heff->Clear();
   hpur->Clear();
   heffpur->Clear();

   TLegend *l = new TLegend(0.59,0.64,0.81,0.87);
   l->SetTextFont(132);
   l->SetLineColor(kWhite);
   l->SetFillColor(kWhite);

   for (int i_bin=1; i_bin < heff->GetXaxis()->GetNbins()+1; i_bin++){
      double selected_cc1pi = hists->h_cc1pi_sel->GetBinContent(i_bin);
      double total_cc1pi = hists->h_cc1pi_sel->GetBinContent(i_bin)+hists->h_cc1pi_notsel->GetBinContent(i_bin);
      double selected_all = hists->h_cc1pi_sel->GetBinContent(i_bin)+hists->h_bg_sel->GetBinContent(i_bin);

      double eff = selected_cc1pi/total_cc1pi;
      double pur = selected_cc1pi/selected_all;
      if (selected_all==0) pur=0;
      //
      // std::cout << "eff = " << eff << std::endl;
      // std::cout << "pur = " << pur << std::endl;

      heff->SetBinContent(i_bin,eff);
      hpur->SetBinContent(i_bin,pur);
      heffpur->SetBinContent(i_bin,eff*pur);
   }

   heff->SetLineColor(kRed);
   heff->SetMarkerColor(kRed);
   heff->SetMarkerStyle(20);
   heff->SetMarkerSize(.3);

   hpur->SetLineColor(kBlue);
   hpur->SetMarkerColor(kBlue);
   hpur->SetMarkerStyle(20);
   hpur->SetMarkerSize(.3);

   heffpur->SetLineColor(kBlack);
   heffpur->SetMarkerColor(kBlack);
   heffpur->SetMarkerStyle(20);
   heffpur->SetMarkerSize(.3);

   heff->GetYaxis()->SetRangeUser(0,1);

   gStyle->SetOptStat(0); // No stats box

   c->cd();
   heff->Draw(drawopt.c_str());
   hpur->Draw((std::string("same")+drawopt).c_str());
   heffpur->Draw((std::string("same")+drawopt).c_str());

   l->AddEntry(heff,"Efficiency","p");
   l->AddEntry(hpur,"Purity","p");
   l->AddEntry(heffpur,"Efficiency #times Purity","p");
   l->Draw();

   // For N-1 plots only, we want to overlay a histogram that has the "all cuts" values in all bins, for easy comparison. Hardcode that "all cuts" will be in bin 2.
   if (isNminus1){
      TH1D *heff_allcuts = (TH1D*)heff->Clone("heff_allcuts");
      TH1D *hpur_allcuts = (TH1D*)hpur->Clone("hpur_allcuts");
      TH1D *heffpur_allcuts = (TH1D*)heffpur->Clone("heffpur_allcuts");

      heff_allcuts->SetLineStyle(2);
      hpur_allcuts->SetLineStyle(2);
      heffpur_allcuts->SetLineStyle(2);

      int allcuts_bin = heff_allcuts->GetXaxis()->GetNbins()-2;

      for (int i_bin=1; i_bin<heff_allcuts->GetXaxis()->GetNbins()+1; i_bin++){
         heff_allcuts->SetBinContent(i_bin,heff->GetBinContent(allcuts_bin));
         hpur_allcuts->SetBinContent(i_bin,hpur->GetBinContent(allcuts_bin));
         heffpur_allcuts->SetBinContent(i_bin,heffpur->GetBinContent(allcuts_bin));
      }

      std::cout << "Efficiency (all cuts): " << heff->GetBinContent(allcuts_bin) << std::endl;
      std::cout << "Purity (all cuts): " << hpur->GetBinContent(allcuts_bin) << std::endl;
      std::cout << "Efficiency x Purity (all cuts): " << heffpur->GetBinContent(allcuts_bin) << std::endl;

      heff_allcuts->Draw((std::string("same")+drawopt).c_str());
      hpur_allcuts->Draw((std::string("same")+drawopt).c_str());
      heffpur_allcuts->Draw((std::string("same")+drawopt).c_str());

      double xlow = heff_allcuts->GetBinLowEdge(allcuts_bin);
      TLine *line = new TLine(xlow,0,xlow,1); // Hack: this is just to give a vertical line separating the true N-1 bins from the other type of bins we have in the plot
      line->SetLineColor(kBlack);
      line->SetLineWidth(2);
      line->Draw("same");

      c->SetBottomMargin(0.13);
      //c->SetRightMargin(0.15);

      l->SetX1(0.15);
      l->SetX2(0.37);
   }

   c->Draw();
   c->Modified();
   c->Update();

}

#endif
