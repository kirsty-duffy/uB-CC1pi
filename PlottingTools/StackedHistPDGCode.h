#ifndef STACKEDHIST_PDGCODE
#define STACKEDHIST_PDGCODE

//#include "uboone/CC1pi/Algorithms/PDGEnums.h"
#include "../Algorithms/PDGEnums.h"

#include "TH1F.h"
#include "THStack.h"

#include <vector>
#include <string>

class StackedHistPDGCode{

 public:
  StackedHistPDGCode(std::string histname, std::string title, int nbins, double lowlimit, double highlimit);// Constructor
  void Fill(PDGCode particle_pdg, double value);
  void Fill(PDGCode particle_pdg, double value, double weight);
  void DrawStack(double norm, TCanvas *c1, TString option);
  void DrawOverlay(double norm, TCanvas *c1, TString option);
  void DrawOverlayMuPi(double norm, TCanvas *c1, TString option);

 protected:
  int nHists;
  std::vector<PDGCode> hist_order; // For keeping track of which hist goes with which PDG code

  THStack *stack;
  TH1F *hists[25];

  void StyleHistsStack();
  void StyleHistsOverlay();
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
  hist_order.push_back(kPiPlus);
  hist_order.push_back(kPiMinus);
  hist_order.push_back(kPiZero);
  hist_order.push_back(kElectron);
  hist_order.push_back(kPositron);
  hist_order.push_back(kTauMinus);
  hist_order.push_back(kTauPlus);
  hist_order.push_back(kPhoton);
  hist_order.push_back(kProton);
  hist_order.push_back(kNeutron);
  hist_order.push_back(kKaPlus);
  hist_order.push_back(kKaMinus);
  hist_order.push_back(kPDGUnknown);

  nHists = hist_order.size();

  stack = new THStack(histname.c_str(),title.c_str());
  for (int i_hist=0; i_hist < nHists; i_hist++){
    std::string histname_i = std::string(histname)+std::string("_")+std::to_string(i_hist);
    hists[i_hist] = new TH1F(histname_i.c_str(),title.c_str(),nbins,lowlimit,highlimit);
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
void StackedHistPDGCode::DrawStack(double norm, TCanvas *c1, TString option="")
{
  // First: style the histograms
  StackedHistPDGCode::StyleHistsStack();

  // Next: add histogramst to the stack and make TLegend
  // Only do this for histograms that have entries
  TLegend *leg = new TLegend(0.1,0.88,0.9,0.99);
  // leg->SetTextFont(132);
  leg->SetLineColor(kWhite);
  leg->SetTextAlign(12);
  leg->SetNColumns(6);

  for (int i_hist=0; i_hist < nHists; i_hist++){
    if (hists[i_hist]->GetEntries() == 0) continue;

    hists[i_hist]->Scale(norm);
    stack->Add(hists[i_hist]);
    PDGCode pdg_for_legend = StackedHistPDGCode::GetPDGFromHistN((unsigned int)i_hist);
    leg->AddEntry(hists[i_hist],PDGenum2str(pdg_for_legend).c_str(),"f");

  }

  c1->cd();
  c1->SetTopMargin(0.13);
  stack->Draw("hist"+option);
  leg->Draw();
}

// -------------------------- Function to draw the histograms -------------------------- //
void StackedHistPDGCode::DrawOverlay(double norm, TCanvas *c1, TString option="")
{
  // First: style the histograms
  StackedHistPDGCode::StyleHistsOverlay();

  // Next: add histogramst to the stack and make TLegend
  // Only do this for histograms that have entries
  TLegend *leg = new TLegend(0.1,0.88,0.9,0.99);
  // leg->SetTextFont(132);
  leg->SetLineColor(kWhite);
  leg->SetTextAlign(12);
  leg->SetNColumns(6);

  leg -> SetNColumns(2);

  bool drawn_first = false;

  for (int i_hist=0; i_hist < nHists; i_hist++){
    if (hists[i_hist]->GetEntries() == 0) continue;

    if (norm == 0.) hists[i_hist]->Scale(1.0/hists[i_hist]->Integral());
    else hists[i_hist]->Scale(norm);

    PDGCode pdg_for_legend = StackedHistPDGCode::GetPDGFromHistN((unsigned int)i_hist);
    leg->AddEntry(hists[i_hist],PDGenum2str(pdg_for_legend).c_str(),"f");

    hists[i_hist]->SetFillStyle(3004);
    hists[i_hist]->SetLineWidth(2);

    c1->cd();
    c1->SetTopMargin(0.13);
    if (!drawn_first){
      if (norm == 0.) hists[i_hist]->GetYaxis()->SetRangeUser(0,1.);
      hists[i_hist]->Draw("hist"+option);
      drawn_first = true;
    }
    else{
      hists[i_hist]->Draw("hist same"+option);
    }
  }

  c1->cd();
  leg->Draw();
}

// -------------------------- Function to draw the histograms -------------------------- //
void StackedHistPDGCode::DrawOverlayMuPi(double norm, TCanvas *c1, TString option="")
{
  // First: style the histograms
  StackedHistPDGCode::StyleHistsOverlay();

  // Next: make TLegend
  // Only do this for histograms that have entries
  TLegend *leg = new TLegend(0.1,0.88,0.9,0.99);
  // leg->SetTextFont(132);
  leg->SetLineColor(kWhite);
  leg->SetTextAlign(12);
  leg->SetNColumns(2);

  bool drawn_first = false;

  TH1D *clone_hist_muplus = nullptr;
  TH1D *clone_hist_piplus = nullptr;

  for (int i_hist=0; i_hist < nHists; i_hist++){
    if (hists[i_hist]->GetEntries() == 0) continue;

    TH1D *hist = (TH1D*)hists[i_hist]->Clone("hist");

    std::string legend_entry = "";
    PDGCode pdg_for_legend = StackedHistPDGCode::GetPDGFromHistN((unsigned int)i_hist);

    if (pdg_for_legend == kMuPlus) {
      clone_hist_muplus = (TH1D*)hist->Clone("clone_hist_muplus");
      continue;
    }
    else if (pdg_for_legend == kPiPlus){
      clone_hist_piplus = (TH1D*)hist->Clone("clone_hist_piplus");
      continue;
    }
    else if (pdg_for_legend == kMuMinus){
      legend_entry = "#mu^{+/-}";
      if (clone_hist_muplus){
        hist->Add(clone_hist_muplus);
      }
    }
    else if (pdg_for_legend == kPiMinus){
      legend_entry = "#pi^{+/-}";
      if (clone_hist_piplus){
        hist->Add(clone_hist_piplus);
      }
    }
    else continue; // Only plot muons and pions here

    leg->AddEntry(hist,legend_entry.c_str(),"f");

    if (norm == 0.) hist->Scale(1.0/hist->Integral());
    else hist->Scale(norm);

    hist->SetFillStyle(3004);
    hist->SetLineWidth(2);
    hist->Rebin(2);

    c1->cd();
    c1->SetTopMargin(0.13);
    if (!drawn_first){
      if (norm == 0.) hist->GetYaxis()->SetRangeUser(0,1.);
      hist->Draw("hist"+option);
      drawn_first = true;
    }
    else{
      hist->Draw("hist same"+option);
    }
  }

  c1->cd();
  leg->Draw();
}

// -------------------------- Function to style the histograms -------------------------- //
// Private: only called by DrawStack function in this file
void StackedHistPDGCode::StyleHistsStack()
{
  // Set fill color for all histograms
  hists[0] ->SetFillColor(kOrange); // kNuMu
  hists[1] ->SetFillColor(kOrange-3); // kNuMuBar
  hists[2] ->SetFillColor(kOrange+2); // kNuE
  hists[3] ->SetFillColor(kRed); // kNuEBar
  hists[4] ->SetFillColor(kRed+2); // kNuTau
  hists[5] ->SetFillColor(kPink-7); // kNuTauBar
  hists[6] ->SetFillColor(kPink+10); // kMuMinus
  hists[8] ->SetFillColor(kViolet+1); // kPiMinus
  hists[7] ->SetFillColor(kMagenta+1); // kMuPlus
  hists[9] ->SetFillColor(kBlue+2); // kPiPlus
  hists[10]->SetFillColor(kBlue); // kPiZero
  hists[11]->SetFillColor(kAzure+1); // kElectron
  hists[12]->SetFillColor(kCyan+2); // kPositron
  hists[13]->SetFillColor(kCyan); // kTauMinus
  hists[14]->SetFillColor(kGreen+1); // kTauPlus
  hists[15]->SetFillColor(kOrange-2); // kPhoton
  hists[16]->SetFillColor(kGray); // kProton
  hists[17]->SetFillColor(kGray+2); // kNeutron
  hists[18]->SetFillColor(kGreen+3); // kKaonPlus
  hists[19]->SetFillColor(kGreen+2); // kKaonMinus
  hists[20]->SetFillColor(kBlack); // kPDGUnknown
}

// -------------------------- Function to style the histograms -------------------------- //
// Private: only called by DrawOverlay function in this file
void StackedHistPDGCode::StyleHistsOverlay()
{
  // Set line color for all histograms
  hists[0] ->SetLineColor(kOrange); // kNuMu
  hists[1] ->SetLineColor(kOrange-3); // kNuMuBar
  hists[2] ->SetLineColor(kOrange+2); // kNuE
  hists[3] ->SetLineColor(kRed); // kNuEBar
  hists[4] ->SetLineColor(kRed+2); // kNuTau
  hists[5] ->SetLineColor(kPink-7); // kNuTauBar
  hists[6] ->SetLineColor(kPink+10); // kMuMinus
  hists[8] ->SetLineColor(kViolet+1); // kPiMinus
  hists[7] ->SetLineColor(kMagenta+1); // kMuPlus
  hists[9] ->SetLineColor(kBlue+2); // kPiPlus
  hists[10]->SetLineColor(kBlue); // kPiZero
  hists[11]->SetLineColor(kAzure+1); // kElectron
  hists[12]->SetLineColor(kCyan+2); // kPositron
  hists[13]->SetLineColor(kCyan); // kTauMinus
  hists[14]->SetLineColor(kGreen+1); // kTauPlus
  hists[15]->SetLineColor(kOrange-2); // kPhoton
  hists[16]->SetLineColor(kGray); // kProton
  hists[17]->SetLineColor(kGray+2); // kNeutron
  hists[18]->SetLineColor(kGreen+3); // kKaonPlus
  hists[19]->SetLineColor(kGreen+2); // kKaonMinus
  hists[20]->SetLineColor(kBlack); // kPDGUnknown
}


// ---------------------- Function to get histogram number for given PDG code ---------------------- //
// Private: only called by functions in this class
unsigned int StackedHistPDGCode::GetHistN(PDGCode particle_pdg)
{
  unsigned int HistN;
  bool found_hist=false;

  for (int i=0; i<nHists; i++){
    if (hist_order.at(i) == particle_pdg){
      HistN = i;
      found_hist = true;
      break;
    }
  }

  if (!found_hist){
    //std::cout << "[ERROR: StackedHistPDGCode.h] Could not find histogram for PDG code " << particle_pdg << std::endl;
    HistN = nHists-1;
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
