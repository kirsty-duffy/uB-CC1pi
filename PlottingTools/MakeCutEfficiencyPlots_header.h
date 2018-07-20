#include "../../../../larana/larana/TruncatedMean/Algorithm/TruncMean.cxx"
#include "../Algorithms/TopologyEnums.h"
#include "StackedHistPDGCode.h"
#include "StackedHistTopology.h"
#include "TTree.h"
#include "TH1.h"
#include "TColor.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TStyle.h"
#include "TFile.h"
#include "TVector3.h"
// #include <vector>

#include "treevars_header_structdef.h"
#include "CutValues_header.h"
#include "treevars_header.h"


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

bool IsEventSelected(double cutval, std::vector<double> value_vec, bool KeepBelowCut, bool Marco_selected, std::vector<bool> TPCObj_PFP_isDaughter, bool OnlyDaughters, std::string TracksNeeded, bool isMIPCut){
   // value_vec is a vector corresponding to the TPCObject. It has one entry per track in the TPCObject.

   int n_tracks = 0;
   int n_failed = 0;
   for (size_t i_track=0; i_track<value_vec.size(); i_track++){
      double value = value_vec.at(i_track);

      // Ignore PFPs that failed to reco as tracks. The neutrino PFP will also have a bogus value.
      if (value == -9999 || value == -999) {
         n_failed++;
         continue;
      }

      bool track_passes_cut=false;
      // (Optionally) only conisder direct daughters of the neutrino
      if (OnlyDaughters && !TPCObj_PFP_isDaughter.at(i_track)){
         track_passes_cut=false;
         continue;
      }

      if (KeepBelowCut && value < cutval){
         track_passes_cut=true;
      }
      else if (!KeepBelowCut && value > cutval){
         track_passes_cut=true;
      }

      if (track_passes_cut){
         n_tracks++;
      }
   } // end loop over tracks in TPCObject

   // Throw away events that failed the CC inclusive selection or that have PFPs that failed to reco as tracks.
   // Then check if the correct number of tracks pass the cut, given the chosen options
   bool isSelected = false;
   if(n_failed < 2 && Marco_selected) { // The neutrino PFP will look "failed", so we expect every event to have 1
      if(TracksNeeded == "atleasttwo" && n_tracks >= 2) isSelected = true;
      else if (TracksNeeded == "exactlytwo" && n_tracks == 2) isSelected = true;
      else if (TracksNeeded == "all") {
         if (OnlyDaughters && n_tracks == std::count(TPCObj_PFP_isDaughter.begin(),TPCObj_PFP_isDaughter.end(),true)) isSelected = true;
         else if(!OnlyDaughters && n_tracks == value_vec.size()-1) isSelected = true; // Again, 1 less because of the neutrino
      }
   }
   return isSelected;
}

void FillCC1piEffPurHist(histCC1piselEffPur *hists, std::vector<double> value_vec, NuIntTopology topology, bool KeepBelowCut, bool Marco_selected, std::vector<bool> TPCObj_PFP_isDaughter, bool OnlyDaughters, std::string TracksNeeded, bool isMIPCut){
   // Loop through all bins in the histograms, and evaluate a cut value at the centre of each bin
   for (int i_bin=1; i_bin < hists->h_cc1pi_sel->GetXaxis()->GetNbins()+1; i_bin++){
      double cutval = hists->h_cc1pi_sel->GetXaxis()->GetBinCenter(i_bin);

      bool isSelected = IsEventSelected(cutval, value_vec, KeepBelowCut, Marco_selected, TPCObj_PFP_isDaughter, OnlyDaughters, TracksNeeded, isMIPCut);

      // Fill selection info into the relevant histograms
      if (isSelected){
         if (topology == kCC1piplus0p || topology == kCC1piplus1p || topology == kCC1piplusNp){
            hists->h_cc1pi_sel->Fill(cutval);
         }
         else {
            hists->h_bg_sel->Fill(cutval);
         }
      }
      else { // if not selected
         if (topology == kCC1piplus0p || topology == kCC1piplus1p || topology == kCC1piplusNp){
            hists->h_cc1pi_notsel->Fill(cutval);
         }
         // else {
         //   hists->h_bg_notsel->Fill(cutval);
         // }
      }
   } // End loop over bins in the histograms
}

void FillNminus1EffPurHist(histCC1piselEffPur *hists, std::vector<std::vector<double>> value_vec, NuIntTopology topology, std::vector<bool> KeepBelowCut, bool Marco_selected, std::vector<bool> TPCObj_PFP_isDaughter, std::vector<bool> OnlyDaughters, std::vector<std::string> TracksNeeded, std::vector<double> cutvalues){
   int ncuts = value_vec.size();
   // Loop through all bins in the histogram and evaluate the cuts
   // Each bin represents one cut (note: that means that for a given bin, we apply all cuts EXCEPT that one), with the final bin showing what we get when we apply all cuts.
   // "Applying" cuts means using the cut values given in the vector
   for (int i_bin=1; i_bin < hists->h_cc1pi_sel->GetXaxis()->GetNbins()+1; i_bin++){
      // Now loop over cuts and apply them
      bool isSelected = true;
      if (i_bin<ncuts+2){
         for (int i_cut=0; i_cut < cutvalues.size(); i_cut++){
            // Skip cut associated to current bin
            // For bin with i_bin=ncuts+2, it will not be equal to i_cut+1 for any cut, so it should just evaluate all of them
            if (i_cut+1 == i_bin) {
               // std::cout << "skipping cut" << std::endl;
               continue;
            }
            // Evaluate whether event passes this cut. If it doesn't, set isSelected to false for the event
            bool isSelected_i = IsEventSelected(cutvalues.at(i_cut), value_vec.at(i_cut), KeepBelowCut.at(i_cut), Marco_selected, TPCObj_PFP_isDaughter, OnlyDaughters.at(i_cut), TracksNeeded.at(i_cut), false);
            if (isSelected_i == false) isSelected = false;

            // if (isSelected) std::cout << "bin: " << i_bin << ", cut: " << i_cut << "isSelected = " << isSelected << std::endl;
         } // end loop over cuts
      }
      // Bin ncuts+2: only evaluate whether the event passed Marco's selection and has 2 tracks
      else if (i_bin==ncuts+2){
         int n_tracks = 0;
         for (size_t i_track=0; i_track<value_vec.at(0).size(); i_track++){
            double value = value_vec.at(0).at(i_track);

            // Ignore PFPs that failed to reco as tracks. The neutrino PFP will also have a bogus value.
            if (value == -9999 || value == -999) {
               continue;
            }

            n_tracks++;
         } // end loop over tracks in TPCObject
         if (n_tracks>=2 && Marco_selected){
            isSelected = true;
         }
         else{
            isSelected = false;
         }
      }
      // Bin ncuts+3: only evaluate whether the event passed Marco's selection
      if (i_bin==ncuts+3){
         if (Marco_selected){
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
         if (topology == kCC1piplus0p || topology == kCC1piplus1p || topology == kCC1piplusNp){
            hists->h_cc1pi_sel->Fill(fillval);
         }
         else {
            hists->h_bg_sel->Fill(fillval);
         }
      }
      else { // if not selected
         if (topology == kCC1piplus0p || topology == kCC1piplus1p || topology == kCC1piplusNp){
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

struct histCC1piselEffPur2D{
   TH2D *h_cc1pi_sel;
   TH2D *h_bg_sel;
   TH2D *h_cc1pi_notsel;
   // TH2D *h_bg_notsel;

   TLegend *l;

   // Constructor for this struct of hists
   histCC1piselEffPur2D(std::string name, std::string title, double nbinsx, double binlowx, double binhighx, double nbinsy, double binlowy, double binhighy){
      h_cc1pi_sel = new TH2D(std::string(name+"_cc1pi_sel").c_str(),title.c_str(),nbinsx,binlowx,binhighx,nbinsy,binlowy,binhighy);
      h_bg_sel = new TH2D(std::string(name+"_bg_sel").c_str(),title.c_str(),nbinsx,binlowx,binhighx,nbinsy,binlowy,binhighy);
      h_cc1pi_notsel = new TH2D(std::string(name+"_cc1pi_notsel").c_str(),title.c_str(),nbinsx,binlowx,binhighx,nbinsy,binlowy,binhighy);
      // h_bg_notsel = new TH2D(std::string(name+"_bg_notsel").c_str(),title.c_str(),nbinsx,binlowx,binhighx,nbinsy,binlowy,binhighy);


      // h_cc1pi_sel->SetFillColor(TColor::GetColor(8,64,129));
      // h_cc1pi_all->SetFillColor(TColor::GetColor(8,64,129));
      // h_bg_sel->SetFillColor(TColor::GetColor(197,197,197));
      // h_bg_all->SetFillColor(TColor::GetColor(197,197,197));
      //
      // l = new TLegend(0.59,0.64,0.81,0.87);
      // l->AddEntry(h_p,"True cc1#pi^{+}","f");
      // l->AddEntry(h_mu,"True other","f");
   }

};

void FillCC1piEffPurHist2D(histCC1piselEffPur2D *hists, std::vector<std::vector<double>> value_vec, NuIntTopology topology, std::vector<bool> KeepBelowCut, bool Marco_selected, std::vector<bool> TPCObj_PFP_isDaughter, std::vector<bool> OnlyDaughters, std::vector<std::string> TracksNeeded, std::vector<bool> isMIPCut){

   // Loop through all bins in the histograms, and evaluate a cut value at the centre of each bin
   for (int i_binx=1; i_binx < hists->h_cc1pi_sel->GetXaxis()->GetNbins()+1; i_binx++){
      double cutvalx = hists->h_cc1pi_sel->GetXaxis()->GetBinCenter(i_binx);
      for (int i_biny=1; i_biny < hists->h_cc1pi_sel->GetYaxis()->GetNbins()+1; i_biny++){
         double cutvaly = hists->h_cc1pi_sel->GetYaxis()->GetBinCenter(i_biny);

         std::vector<double> cutvalues = {cutvalx, cutvaly};

         // isSelected is now determined by a loop
         bool isSelected = true;
         for (int i_cut=0; i_cut < cutvalues.size(); i_cut++){
            // Evaluate whether event passes this cut. If it doesn't, set isSelected to false for the event
            bool isSelected_i = IsEventSelected(cutvalues.at(i_cut), value_vec.at(i_cut), KeepBelowCut.at(i_cut), Marco_selected, TPCObj_PFP_isDaughter, OnlyDaughters.at(i_cut), TracksNeeded.at(i_cut), isMIPCut.at(i_cut));
            if (isSelected_i == false) {
               isSelected = false;
               continue;
            }

            // if (isSelected) std::cout << "bin: " << i_bin << ", cut: " << i_cut << "isSelected = " << isSelected << std::endl;
         } // end loop over cuts

         // Fill selection info into the relevant histograms
         if (isSelected){
            if (topology == kCC1piplus0p || topology == kCC1piplus1p || topology == kCC1piplusNp){
               hists->h_cc1pi_sel->Fill(cutvalx,cutvaly);
            }
            else {
               hists->h_bg_sel->Fill(cutvalx,cutvaly);
            }
         }
         else { // if not selected
            if (topology == kCC1piplus0p || topology == kCC1piplus1p || topology == kCC1piplusNp){
               hists->h_cc1pi_notsel->Fill(cutvalx,cutvaly);
            }
            // else {
            //   hists->h_bg_notsel->Fill(cutval);
            // }
         }
      } // End loop over y bins in the histograms
   } // end loop over x bins
}

void DrawCC1piMCEffPur2D(TCanvas *c, histCC1piselEffPur2D *hists){
   TH2D *heff = (TH2D*)hists->h_cc1pi_sel->Clone("heff");
   TH2D *hpur = (TH2D*)hists->h_cc1pi_sel->Clone("hpur");
   TH2D *heffpur = (TH2D*)hists->h_cc1pi_sel->Clone("heffpur");

   heff->Clear();
   hpur->Clear();
   heffpur->Clear();

   for (int i_binx=1; i_binx < heff->GetXaxis()->GetNbins()+1; i_binx++){
      for (int i_biny=1; i_biny < heff->GetXaxis()->GetNbins()+1; i_biny++){
         double selected_cc1pi = hists->h_cc1pi_sel->GetBinContent(i_binx,i_biny);
         double total_cc1pi = hists->h_cc1pi_sel->GetBinContent(i_binx,i_biny)+hists->h_cc1pi_notsel->GetBinContent(i_binx,i_biny);
         double selected_all = hists->h_cc1pi_sel->GetBinContent(i_binx,i_biny)+hists->h_bg_sel->GetBinContent(i_binx,i_biny);

         double eff = selected_cc1pi/total_cc1pi;
         double pur = selected_cc1pi/selected_all;

         heff->SetBinContent(i_binx,i_biny,eff);
         hpur->SetBinContent(i_binx,i_biny,pur);
         heffpur->SetBinContent(i_binx,i_biny,eff*pur);

         heff->SetTitle("Efficiency");
         hpur->SetTitle("Purity");
         heffpur->SetTitle("Efficiency #times Purity");
      }
   }

   gStyle->SetOptStat(0); // No stats box

   c->Clear();
   c->Divide(2,2,0.0005,0.0005);

   std::vector<TH2D*> histstoeval = {
      heff,
      hpur,
      heffpur
   };

   for (int i_h=0; i_h<histstoeval.size(); i_h++){
      c->cd(i_h+1);
//      histstoeval.at(i_h)->GetZaxis()->SetRangeUser(0,1);
      histstoeval.at(i_h)->Draw("colz");
   }

   c->Modified();
   c->Update();

}
