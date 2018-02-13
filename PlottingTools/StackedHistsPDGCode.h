#ifndef STACKEDHISTS_PDGCODE
#define STACKEDHISTS_PDGCODE

#include "uboone/CC1pi/Algorithms/PDGEnums.h"

#include "TH1F.h"
#include "THStack.h"

#include <vector>
#include <string>

class StackedHistPDGCode{

 public:
  StackedHistPDGCode(std::string histname, std::string title, int nbins, double lowlimit, double highlimit);// Constructor
  void Fill(PDGCode particle_pdg, double value);
  void Fill(PDGCode particle_pdg, double value, double weight);
  void DrawStack(double norm, TCanvas *c1);

 protected:
  int nHists;
  std::vector<PDGCode> hist_order; // For keeping track of which hist goes with which PDG code
  
  THStack *stack;
  TH1F *hists[25];

  void StyleHists();
  unsigned int GetHistN(PDGCode particle_pdg);
  PDGCode GetPDGFromHistN(unsigned int hist_n);
};

// Define functions here instead of a .C file so we don't need too many includes

// -------------------------- Constructor -------------------------- //
// Intended to be very general so you can make histograms of whatever you like
StackedHistPDGCode::StackedHistPDGCode(std::string histname, std::string title, int nbins, double lowlimit, double highlimit)
{
  // PDG categories for the histograms
  hist_order.push_back(kNuMu);
  hist_order.push_back(kNuMuBar);
  hist_order.push_back(kNuE);
  hist_order.push_back(kNuEBar);
  hist_order.push_back(kNuTau);
  hist_order.push_back(kNuTauBar);
  hist_order.push_back(kMuMinus);
  hist_order.push_back(kMuPlus);
  hist_order.push_back(kPiMinus);
  hist_order.push_back(kPiPlus);
  hist_order.push_back(kPiZero);
  hist_order.push_back(kElectron);
  hist_order.push_back(kPositron);
  hist_order.push_back(kTauMinus);
  hist_order.push_back(kTauPlus);
  hist_order.push_back(kPhoton);
  hist_order.push_back(kProton);
  hist_order.push_back(kNeutron);
  hist_order.push_back(kPDGUnknown);
  
  nHists = hist_order.size();
  
  stack = new THStack(histname.c_str(),title.c_str());
  for (int i_hist=0; i_hist < nHists; i_hist++){
    std::string histname_i = std::string(histname)+std::string("_")+std::to_string(i_hist);
    hists[i_hist] = new TH1F(histname_i.c_str(),"",nbins,lowlimit,highlimit);
  }
}

// -------------------------- Function to fill the correct histogram -------------------------- //
void StackedHistPDGCode::Fill(PDGCode particle_pdg, double value)
{
  unsigned int n_hist = StackedHistPDGCode::GetHistN(particle_pdg);
  hists[n_hist]->Fill(value);
}

// -------------------------- Function to fill the correct histogram -------------------------- //
// Overloaded to allow for weights
void StackedHistPDGCode::Fill(PDGCode particle_pdg, double value, double weight)
{
  unsigned int n_hist = StackedHistPDGCode::GetHistN(particle_pdg);
  hists[n_hist]->Fill(value, weight);
}

// -------------------------- Function to draw the histograms -------------------------- //
void StackedHistPDGCode::DrawStack(double norm, TCanvas *c1)
{
  // First: style the histograms
  StackedHistPDGCode::StyleHists();

  // Next: add histogramst to the stack and make TLegend
  // Only do this for histograms that have entries
  TLegend *leg = new TLegend(0.6,0.7,0.95,0.95);
  
  for (int i_hist=0; i_hist < nHists; i_hist++){
    if (hists[i_hist]->GetEntries() == 0) continue;
    
    hists[i_hist]->Scale(norm);
    stack->Add(hists[i_hist]);
    PDGCode pdg_for_legend = StackedHistPDGCode::GetPDGFromHistN((unsigned int)i_hist);
    leg->AddEntry(hists[i_hist],PDGenum2str(pdg_for_legend).c_str(),"f");
    
  }

  c1->cd();
  stack->Draw("hist");
  leg->Draw();
}

// -------------------------- Function to style the histograms -------------------------- //
// Private: only called by DrawStack function in this file
void StackedHistPDGCode::StyleHists()
{
  // Set fill color for all histograms
  hists[0] ->SetFillColor(kOrange);
  hists[1] ->SetFillColor(kOrange-3);
  hists[2] ->SetFillColor(kOrange+2);
  hists[3] ->SetFillColor(kRed);
  hists[4] ->SetFillColor(kRed+2);
  hists[5] ->SetFillColor(kPink-7);
  hists[6] ->SetFillColor(kPink+10);
  hists[7] ->SetFillColor(kMagenta+1);
  hists[8] ->SetFillColor(kViolet+1);
  hists[9] ->SetFillColor(kBlue+2);
  hists[10]->SetFillColor(kBlue);
  hists[11]->SetFillColor(kAzure+1);
  hists[12]->SetFillColor(kCyan+2);
  hists[13]->SetFillColor(kCyan);
  hists[14]->SetFillColor(kGreen+1);
  hists[15]->SetFillColor(kGreen+3);
  hists[16]->SetFillColor(kGray);
  hists[17]->SetFillColor(kGray+2);
  hists[18]->SetFillColor(kBlack);
}


// ---------------------- Function to get histogram number for given PDG code ---------------------- //
// Private: only called by functions in this class
unsigned int StackedHistPDGCode::GetHistN(PDGCode particle_pdg)
{
  unsigned int HistN;
  bool found_hist=false;

  for (unsigned int i=0; i<nHists; i++){
    if (hist_order.at(i) == particle_pdg){
      HistN = i;
      found_hist = true;
      break;
    }
  }

  if (!found_hist){
    std::cout << "[ERROR: StackedHistsPDGCode.h] Could not find histogram for PDG code " << particle_pdg << std::endl;
    HistN = 9999;
  }

  return HistN;
}


// ---------------------- Function to get PDG code for given histogram number ---------------------- //
// Private: only called by functions in this class
PDGCode StackedHistPDGCode::GetPDGFromHistN(unsigned int hist_n)
{
  return hist_order.at(hist_n);
}


#endif
