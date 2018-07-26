#ifndef __CC1PI_EFFPUR_2D_H__
#define __CC1PI_EFFPUR_2D_H__

#include "../Algorithms/TopologyEnums.h"
#include "TColor.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "THStack.h"
#include "TStyle.h"
// #include <vector>

#include "CC1pi_treevars.h"
#include "CC1pi_cuts.h"


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

void FillCC1piEffPurHist2D(histCC1piselEffPur2D *hists, treevars *vars, int i_cutx, int i_cuty){

   // Loop through all bins in the histograms, and evaluate a cut value at the centre of each bin
   for (int i_binx=1; i_binx < hists->h_cc1pi_sel->GetXaxis()->GetNbins()+1; i_binx++){
      double cutvalx = hists->h_cc1pi_sel->GetXaxis()->GetBinCenter(i_binx);
      for (int i_biny=1; i_biny < hists->h_cc1pi_sel->GetYaxis()->GetNbins()+1; i_biny++){
         double cutvaly = hists->h_cc1pi_sel->GetYaxis()->GetBinCenter(i_biny);

         // isSelected is now determined by evaluating x and y cuts
         bool isSelected = true;
         bool isSelected_x = IsEventSelected_SingleCut(cutvalx, vars, i_cutx);
         if (isSelected_x==false) isSelected=false;
         bool isSelected_y = IsEventSelected_SingleCut(cutvaly, vars, i_cuty);
         if (isSelected_y==false) isSelected=false;


         // Fill selection info into the relevant histograms
         if (isSelected){
            if (vars->Truth_topology == kCC1piplus0p || vars->Truth_topology == kCC1piplus1p || vars->Truth_topology == kCC1piplusNp){
               hists->h_cc1pi_sel->Fill(cutvalx,cutvaly);
            }
            else {
               hists->h_bg_sel->Fill(cutvalx,cutvaly);
            }
         }
         else { // if not selected
            if (vars->Truth_topology == kCC1piplus0p || vars->Truth_topology == kCC1piplus1p || vars->Truth_topology == kCC1piplusNp){
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

#endif
