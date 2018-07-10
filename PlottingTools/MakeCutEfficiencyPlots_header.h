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
   std::vector<double> *TPCObj_reco_vtx = nullptr;
   std::vector<bool> *TPCObj_PFP_track_isContained = nullptr;
   std::vector<int> *TPCObj_PFP_truePDG = nullptr;

   // These are derived quantities - derived from the values above in Calcvars
   std::vector<double> *TPCObj_PFP_LH_p;
   std::vector<double> *TPCObj_PFP_LH_mu;
   std::vector<double> *TPCObj_PFP_LH_pi;
   std::vector<double> *TPCObj_PFP_LH_mip;
   std::vector<double> *TPCObj_PFP_Lmipoverp;
   std::vector<double> *TPCObj_PFP_Lmumipovermumipp;
   std::vector<double> *TPCObj_PFP_BrokenTrackAngle;
   std::vector<double> *TPCObj_PFP_VtxTrackDist;
   std::vector<double> *TPCObj_PFP_isContained_double;
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
   intree->SetBranchStatus("TPCObj_reco_vtx",1);
   intree->SetBranchAddress("TPCObj_reco_vtx", &(varstoset->TPCObj_reco_vtx));
   intree->SetBranchStatus("TPCObj_PFP_track_isContained",1);
   intree->SetBranchAddress("TPCObj_PFP_track_isContained", &(varstoset->TPCObj_PFP_track_isContained));
   intree->SetBranchStatus("TPCObj_PFP_truePDG",1);
   intree->SetBranchAddress("TPCObj_PFP_truePDG", &(varstoset->TPCObj_PFP_truePDG));
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
   vars->TPCObj_PFP_VtxTrackDist = new std::vector<double>(vecsize);
   vars->TPCObj_PFP_isContained_double = new std::vector<double>(vecsize);

   // Just use collection plane for now
   int i_pl = 2;

   // Now calculate the values for all variables
   for (int i_track=0; i_track < vecsize; i_track++){

      // The neutrino and any other PFPs that didn't get reco'd as tracks should get bogus values
      if(vars->TPCObj_PFP_track_theta->at(i_track) == -9999) {

         vars->TPCObj_PFP_LH_p->at(i_track) = -9999;
         vars->TPCObj_PFP_LH_mu->at(i_track) = -9999;
         vars->TPCObj_PFP_LH_pi->at(i_track) = -9999;
         vars->TPCObj_PFP_LH_mip->at(i_track) = -9999;
         vars->TPCObj_PFP_Lmipoverp->at(i_track) = -9999;
         vars->TPCObj_PFP_Lmumipovermumipp->at(i_track) = -9999;
         vars->TPCObj_PFP_VtxTrackDist->at(i_track) = -9999;
         vars->TPCObj_PFP_BrokenTrackAngle->at(i_track) = -9999;
         vars->TPCObj_PFP_isContained_double->at(i_track) = -9999;

         continue;
      }


      vars->TPCObj_PFP_LH_p->at(i_track) = std::max(vars->TPCObj_PFP_LH_fwd_p->at(i_track).at(i_pl), vars->TPCObj_PFP_LH_bwd_p->at(i_track).at(i_pl));
      vars->TPCObj_PFP_LH_mu->at(i_track) = std::max(vars->TPCObj_PFP_LH_fwd_mu->at(i_track).at(i_pl), vars->TPCObj_PFP_LH_bwd_mu->at(i_track).at(i_pl));
      vars->TPCObj_PFP_LH_pi->at(i_track) = std::max(vars->TPCObj_PFP_LH_fwd_pi->at(i_track).at(i_pl), vars->TPCObj_PFP_LH_bwd_pi->at(i_track).at(i_pl));
      vars->TPCObj_PFP_LH_mip->at(i_track) = vars->TPCObj_PFP_LH_MIP->at(i_track).at(i_pl);
      vars->TPCObj_PFP_Lmipoverp->at(i_track) = vars->TPCObj_PFP_LH_mip->at(i_track) / vars->TPCObj_PFP_LH_p->at(i_track);
      vars->TPCObj_PFP_Lmumipovermumipp->at(i_track) = (vars->TPCObj_PFP_LH_mu->at(i_track)+vars->TPCObj_PFP_LH_mip->at(i_track))/(vars->TPCObj_PFP_LH_mu->at(i_track)+vars->TPCObj_PFP_LH_mip->at(i_track)+vars->TPCObj_PFP_LH_p->at(i_track));
      vars->TPCObj_PFP_VtxTrackDist->at(i_track) = std::hypot(std::hypot(vars->TPCObj_reco_vtx->at(0) - vars->TPCObj_PFP_track_start->at(i_track).at(0), vars->TPCObj_reco_vtx->at(1) - vars->TPCObj_PFP_track_start->at(i_track).at(1)), vars->TPCObj_reco_vtx->at(2) - vars->TPCObj_PFP_track_start->at(i_track).at(2));
      vars->TPCObj_PFP_isContained_double->at(i_track) = (double)(vars->TPCObj_PFP_track_isContained->at(i_track));

      // Angle between tracks (using theta and phi)...
      int track1 = i_track;

      // We want non-tracks to get bogus values (-9999)
      // But we want tracks that just don't have other tracks close by to get 0
      double maxangle = -1;

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
         // Not needed: cos is symmetric
         // if(mindist==distances.at(1) || mindist==distances.at(3)) {
         //    v2=-1*v2;
         //    TVector3 oldstart2=start2;
         //    start2=end2;
         //    end2=oldstart2;
         // }
         // if(mindist==distances.at(2) || mindist==distances.at(3)) {
         //    v1=-1*v1;
         //    TVector3 oldstart1=start1;
         //    start1=end1;
         //    end1=oldstart1;
         // }

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
   delete vars->TPCObj_PFP_VtxTrackDist;
   delete vars->TPCObj_PFP_isContained_double;

}

// What variables do we want the plots as a function of?
std::vector<std::vector<double>> GetCutvarstoplot(treevars *vars){
   std::vector<std::vector<double>> varstoplot = {
      *(vars->TPCObj_PFP_Lmipoverp),
//      *(vars->TPCObj_PFP_Lmumipovermumipp),
      *(vars->TPCObj_PFP_BrokenTrackAngle),
//      *(vars->TPCObj_PFP_track_residual_mean),
//      *(vars->TPCObj_PFP_track_residual_mean),
//      *(vars->TPCObj_PFP_track_residual_std),
//      *(vars->TPCObj_PFP_track_residual_std),
      *(vars->TPCObj_PFP_track_perc_used_hits),
      *(vars->TPCObj_PFP_VtxTrackDist)
//      *(vars->TPCObj_PFP_isContained_double)
   };
   return varstoplot;
};

// For efficiency/purity we need to know whether we want to be keep tracks that have values above or below the cut value
std::vector<bool> KeepBelowCut = {
   false, // Lmipoverp
//   false, // Lmumipovermumipp
   true,  // BrokenTrackAngle
//   true,  // residual_mean_up
//   false, // residual_mean_down
//   true,  // residual_std_up
//   false, // residual_std_down
   false, // perc_used_hits
   true  // VtxTrackDist
//   false  // isContained
};

// Do we want to consider just the direct daughters of the neutrino?
std::vector<bool> OnlyDaughters = {
   true,  // Lmipoverp
//   true,  // Lmumipovermumipp
   false, // BrokenTrackAngle
//   true,  // residual_mean_up
//   true,  // residual_mean_down
//   true,  // residual_std_up
//   true,  // residual_std_down
   true,  // perc_used_hits
   false // VtxTrackDist
//   false  // isContained
};

// How many tracks do we want to pass the cut? (Options are atleasttwo, exactlytwo, all)
std::vector<std::string> TracksNeeded = {
   "exactlytwo",  // Lmipoverp
//   "exactlytwo",  // Lmumipovermumipp
   "all",         // BrokenTrackAngle
//   "atleasttwo",  // residual_mean_up
//   "atleasttwo",  // residual_mean_down
//   "atleasttwo",  // residual_std_up
//   "atleasttwo",  // residual_std_down
   "atleasttwo",  // perc_used_hits
   "atleasttwo"  // VtxTrackDist
//   "all"          // isConainted
};

// Cut values for N-1 plot
std::vector<double> CutValues = {
   1.,   // Lmipoverp
//   0.66, // Lmumipovermumipp
   3.05, // BrokenTrackAngle
//   0.7,  // residual_mean_up
//   -0.7, // residual_mean_down
//   2.5,  // residual_std_up
//   0.,   // residual_std_down
   0.7,  // perc_used_hits
   15.  // VtxTrackDist
//   0.5   // isContained
};

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
   // if (isSelected) std::cout << "isSelected = " << isSelected << std::endl;
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
            bool isSelected_i = IsEventSelected(cutvalues.at(i_cut), value_vec.at(i_cut), KeepBelowCut.at(i_cut), Marco_selected, TPCObj_PFP_isDaughter, OnlyDaughters.at(i_cut), TracksNeeded.at(i_cut));
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

void FillCC1piEffPurHist2D(histCC1piselEffPur2D *hists, std::vector<std::vector<double>> value_vec, NuIntTopology topology, std::vector<bool> KeepBelowCut, bool Marco_selected, std::vector<bool> TPCObj_PFP_isDaughter, std::vector<bool> OnlyDaughters, std::vector<std::string> TracksNeeded){

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
            bool isSelected_i = IsEventSelected(cutvalues.at(i_cut), value_vec.at(i_cut), KeepBelowCut.at(i_cut), Marco_selected, TPCObj_PFP_isDaughter, OnlyDaughters.at(i_cut), TracksNeeded.at(i_cut));
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
