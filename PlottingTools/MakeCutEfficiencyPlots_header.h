#include "../Algorithms/TopologyEnums.h"
#include "TTree.h"
#include "TH1.h"
#include "TColor.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TStyle.h"
#include "TFile.h"
#include "TVector3.h"
#include <vector>

struct treevars{
   // These are the variables that are filled directly from the tree
   std::vector<std::vector<double>> *TPCObj_PFP_LH_fwd_mu = nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_LH_fwd_p = nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_LH_fwd_pi = nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_LH_bwd_mu = nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_LH_bwd_p = nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_LH_bwd_pi = nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_LH_MIP = nullptr;
   NuIntTopology Truth_topology = kUnknown;
   bool Marco_selected = false;
   std::vector<bool> *TPCObj_PFP_isDaughter = nullptr;
   std::vector<double> *TPCObj_PFP_track_theta = nullptr;
   std::vector<double> *TPCObj_PFP_track_phi = nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_track_start = nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_track_end = nullptr;
   std::vector<double> *TPCObj_PFP_track_residual_mean = nullptr;
   std::vector<double> *TPCObj_PFP_track_residual_std = nullptr;
   std::vector<double> *TPCObj_PFP_track_perc_used_hits = nullptr;

   // These are derived quantities - derived from the values above in Calcvars
   std::vector<double> *TPCObj_PFP_LH_p;
   std::vector<double> *TPCObj_PFP_LH_mu;
   std::vector<double> *TPCObj_PFP_LH_pi;
   std::vector<double> *TPCObj_PFP_LH_mip;
   std::vector<double> *TPCObj_PFP_Lmipoverp;
   std::vector<double> *TPCObj_PFP_Lmumipovermumipp;
   std::vector<double> *TPCObj_PFP_BrokenTrackAngle;
};

void settreevars(TTree *intree, treevars *varstoset){
   intree->SetBranchStatus("*",0);
   intree->SetBranchStatus("TPCObj_PFP_LH_fwd_mu",1);
   intree->SetBranchAddress("TPCObj_PFP_LH_fwd_mu", &(varstoset->TPCObj_PFP_LH_fwd_mu));
   intree->SetBranchStatus("TPCObj_PFP_LH_fwd_p",1);
   intree->SetBranchAddress("TPCObj_PFP_LH_fwd_p", &(varstoset->TPCObj_PFP_LH_fwd_p));
   intree->SetBranchStatus("TPCObj_PFP_LH_fwd_pi",1);
   intree->SetBranchAddress("TPCObj_PFP_LH_fwd_pi", &(varstoset->TPCObj_PFP_LH_fwd_pi));
   intree->SetBranchStatus("TPCObj_PFP_LH_bwd_mu",1);
   intree->SetBranchAddress("TPCObj_PFP_LH_bwd_mu", &(varstoset->TPCObj_PFP_LH_bwd_mu));
   intree->SetBranchStatus("TPCObj_PFP_LH_bwd_p",1);
   intree->SetBranchAddress("TPCObj_PFP_LH_bwd_p", &(varstoset->TPCObj_PFP_LH_bwd_p));
   intree->SetBranchStatus("TPCObj_PFP_LH_bwd_pi",1);
   intree->SetBranchAddress("TPCObj_PFP_LH_bwd_pi", &(varstoset->TPCObj_PFP_LH_bwd_pi));
   intree->SetBranchStatus("TPCObj_PFP_LH_MIP",1);
   intree->SetBranchAddress("TPCObj_PFP_LH_MIP", &(varstoset->TPCObj_PFP_LH_MIP));
   intree->SetBranchStatus("Truth_topology",1);
   intree->SetBranchAddress("Truth_topology", &(varstoset->Truth_topology));
   intree->SetBranchStatus("Marco_selected",1);
   intree->SetBranchAddress("Marco_selected", &(varstoset->Marco_selected));
   intree->SetBranchStatus("TPCObj_PFP_isDaughter",1);
   intree->SetBranchAddress("TPCObj_PFP_isDaughter", &(varstoset->TPCObj_PFP_isDaughter));
   intree->SetBranchStatus("TPCObj_PFP_track_theta",1);
   intree->SetBranchAddress("TPCObj_PFP_track_theta", &(varstoset->TPCObj_PFP_track_theta));
   intree->SetBranchStatus("TPCObj_PFP_track_phi",1);
   intree->SetBranchAddress("TPCObj_PFP_track_phi", &(varstoset->TPCObj_PFP_track_phi));
   intree->SetBranchStatus("TPCObj_PFP_track_start",1);
   intree->SetBranchAddress("TPCObj_PFP_track_start", &(varstoset->TPCObj_PFP_track_start));
   intree->SetBranchStatus("TPCObj_PFP_track_end",1);
   intree->SetBranchAddress("TPCObj_PFP_track_end", &(varstoset->TPCObj_PFP_track_end));
   intree->SetBranchStatus("TPCObj_PFP_track_residual_mean",1);
   intree->SetBranchAddress("TPCObj_PFP_track_residual_mean", &(varstoset->TPCObj_PFP_track_residual_mean));
   intree->SetBranchStatus("TPCObj_PFP_track_residual_std",1);
   intree->SetBranchAddress("TPCObj_PFP_track_residual_std", &(varstoset->TPCObj_PFP_track_residual_std));
   intree->SetBranchStatus("TPCObj_PFP_track_perc_used_hits",1);
   intree->SetBranchAddress("TPCObj_PFP_track_perc_used_hits", &(varstoset->TPCObj_PFP_track_perc_used_hits));
}

void Calcvars(treevars *vars){

   // Initialise vectors that we are going to fill with calculated values
   int vecsize = vars->TPCObj_PFP_LH_fwd_p->size();

   vars->TPCObj_PFP_LH_p = new std::vector<double>(vecsize);
   vars->TPCObj_PFP_LH_mu = new std::vector<double>(vecsize);
   vars->TPCObj_PFP_LH_pi = new std::vector<double>(vecsize);
   vars->TPCObj_PFP_LH_mip = new std::vector<double>(vecsize);
   vars->TPCObj_PFP_Lmipoverp = new std::vector<double>(vecsize);
   vars->TPCObj_PFP_Lmumipovermumipp = new std::vector<double>(vecsize);
   vars->TPCObj_PFP_BrokenTrackAngle = new std::vector<double>(vecsize);

   // Just use collection plane for now
   int i_pl = 2;

   // Now calculate the values for all variables
   for (int i_track=0; i_track < vecsize; i_track++){

      vars->TPCObj_PFP_LH_p->at(i_track) = std::max(vars->TPCObj_PFP_LH_fwd_p->at(i_track).at(i_pl), vars->TPCObj_PFP_LH_bwd_p->at(i_track).at(i_pl));
      vars->TPCObj_PFP_LH_mu->at(i_track) = std::max(vars->TPCObj_PFP_LH_fwd_mu->at(i_track).at(i_pl), vars->TPCObj_PFP_LH_bwd_mu->at(i_track).at(i_pl));
      vars->TPCObj_PFP_LH_pi->at(i_track) = std::max(vars->TPCObj_PFP_LH_fwd_pi->at(i_track).at(i_pl), vars->TPCObj_PFP_LH_bwd_pi->at(i_track).at(i_pl));
      vars->TPCObj_PFP_LH_mip->at(i_track) = vars->TPCObj_PFP_LH_MIP->at(i_track).at(i_pl);
      vars->TPCObj_PFP_Lmipoverp->at(i_track) = vars->TPCObj_PFP_LH_mip->at(i_track) / vars->TPCObj_PFP_LH_p->at(i_track);
      vars->TPCObj_PFP_Lmumipovermumipp->at(i_track) = (vars->TPCObj_PFP_LH_mu->at(i_track)+vars->TPCObj_PFP_LH_mip->at(i_track))/(vars->TPCObj_PFP_LH_mu->at(i_track)+vars->TPCObj_PFP_LH_mip->at(i_track)+vars->TPCObj_PFP_LH_p->at(i_track));

      // Angle between tracks (using theta and phi)...
      int track1 = i_track;

      // We want non-tracks to get bogus values (-9999)
      // But we want tracks that just don't have other tracks close by to get 0
      double maxangle = 0;
      if(vars->TPCObj_PFP_track_theta->at(track1) == -9999) maxangle = -9999;

      for (int track2 = track1+1; track2 < vecsize; track2++) {
         if(vars->TPCObj_PFP_track_theta->at(track1) == -9999) continue;
         if(vars->TPCObj_PFP_track_theta->at(track2) == -9999) continue;

         TVector3 v1(1,1,1);
         TVector3 v2(1,1,1);
         v1.SetMag(1);
         v2.SetMag(1);
         v1.SetTheta(vars->TPCObj_PFP_track_theta->at(track1));
         v2.SetTheta(vars->TPCObj_PFP_track_theta->at(track2));
         v1.SetPhi(vars->TPCObj_PFP_track_phi->at(track1));
         v2.SetPhi(vars->TPCObj_PFP_track_phi->at(track2));

         TVector3 start1(vars->TPCObj_PFP_track_start->at(track1).at(0),vars->TPCObj_PFP_track_start->at(track1).at(1),vars->TPCObj_PFP_track_start->at(track1).at(2));
         TVector3 end1(vars->TPCObj_PFP_track_end->at(track1).at(0),vars->TPCObj_PFP_track_end->at(track1).at(1),vars->TPCObj_PFP_track_end->at(track1).at(2));
         TVector3 start2(vars->TPCObj_PFP_track_start->at(track2).at(0),vars->TPCObj_PFP_track_start->at(track2).at(1),vars->TPCObj_PFP_track_start->at(track2).at(2));
         TVector3 end2(vars->TPCObj_PFP_track_end->at(track2).at(0),vars->TPCObj_PFP_track_end->at(track2).at(1),vars->TPCObj_PFP_track_end->at(track2).at(2));
         std::vector<double> distances = { (start1-start2).Mag(),(start1-end2).Mag(), (end1-start2).Mag(), (end1-end2).Mag() };
         double mindist = *std::min_element(distances.begin(),distances.end());

         // Flip tracks if the closest point isn't start-start
         if(mindist==distances.at(1) || mindist==distances.at(3)) {
            v2=-1*v2;
            TVector3 oldstart2=start2;
            start2=end2;
            end2=oldstart2;
         }
         if(mindist==distances.at(2) || mindist==distances.at(3)) {
            v1=-1*v1;
            TVector3 oldstart1=start1;
            start1=end1;
            end1=oldstart1;
         }

         double angle = TMath::ACos(v1.Dot(v2));

         if(mindist < 3 && angle > maxangle) maxangle = angle;
      }

      vars->TPCObj_PFP_BrokenTrackAngle->at(i_track) = maxangle;
   }
}


// --------------------------------------------------- //
// Function to clear memory from calculated variables


void Clearvars(treevars *vars){

   delete vars->TPCObj_PFP_LH_p;
   delete vars->TPCObj_PFP_LH_mu;
   delete vars->TPCObj_PFP_LH_pi;
   delete vars->TPCObj_PFP_LH_mip;
   delete vars->TPCObj_PFP_Lmipoverp;
   delete vars->TPCObj_PFP_Lmumipovermumipp;
   delete vars->TPCObj_PFP_BrokenTrackAngle;
}


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

bool IsEventSelected(double cutval, std::vector<double> value_vec, bool KeepBelowCut, bool Marco_selected, std::vector<bool> TPCObj_PFP_isDaughter, bool OnlyDaughters, std::string TracksNeeded){
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

      // (Optionally) only conisder direct daughters of the neutrino
      if (OnlyDaughters && !TPCObj_PFP_isDaughter.at(i_track)) continue;

      if (KeepBelowCut && value < cutval){
         n_tracks++;
      }
      else if (!KeepBelowCut && value > cutval){
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

void FillCC1piEffPurHist(histCC1piselEffPur *hists, std::vector<double> value_vec, NuIntTopology topology, bool KeepBelowCut, bool Marco_selected, std::vector<bool> TPCObj_PFP_isDaughter, bool OnlyDaughters, std::string TracksNeeded){
   // Loop through all bins in the histograms, and evaluate a cut value at the centre of each bin
   for (int i_bin=1; i_bin < hists->h_cc1pi_sel->GetXaxis()->GetNbins()+1; i_bin++){
      double cutval = hists->h_cc1pi_sel->GetXaxis()->GetBinCenter(i_bin);

      bool isSelected = IsEventSelected(cutval, value_vec, KeepBelowCut, Marco_selected, TPCObj_PFP_isDaughter, OnlyDaughters, TracksNeeded);

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
   // Loop through all bins in the histogram and evaluate the cuts
   // Each bin represents one cut (note: that means that for a given bin, we apply all cuts EXCEPT that one), with the final bin showing what we get when we apply all cuts.
   // "Applying" cuts means using the cut values given in the vector
   for (int i_bin=1; i_bin < hists->h_cc1pi_sel->GetXaxis()->GetNbins()+1; i_bin++){
      // Now loop over cuts and apply them
      bool isSelected = true;
      for (int i_cut=0; i_cut < cutvalues.size(); i_cut++){
         // Skip cut associated to current bin
         // For final bin, it will not be equal to i_cut for any cut, so it should just evaluate all of them
         if (i_cut+1 == i_bin) {
            std::cout << "skipping cut" << std::endl;
            continue;
         }
         // Evaluate whether event passes this cut. If it doesn't, set isSelected to false for the event
         bool isSelected_i = IsEventSelected(cutvalues.at(i_cut), value_vec.at(i_cut), KeepBelowCut.at(i_cut), Marco_selected, TPCObj_PFP_isDaughter, OnlyDaughters.at(i_cut), TracksNeeded.at(i_cut));
         if (isSelected_i == false) isSelected = false;

         std::cout << "bin: " << i_bin << ", cut: " << i_cut << "isSelected = " << isSelected << std::endl;
      } // end loop over cuts

      // Now fill selection info into the relevant histograms
      double fillval = hists->h_cc1pi_sel->GetBinCenter(i_bin);
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

void DrawCC1piMCEffPur(TCanvas *c, histCC1piselEffPur *hists){
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
   heff->Draw("p");
   hpur->Draw("same p");
   heffpur->Draw("same p");

   l->AddEntry(heff,"Efficiency","p");
   l->AddEntry(hpur,"Purity","p");
   l->AddEntry(heffpur,"Efficiency #times Purity","p");
   l->Draw();

   c->Draw();
   c->Modified();
   c->Update();

}
