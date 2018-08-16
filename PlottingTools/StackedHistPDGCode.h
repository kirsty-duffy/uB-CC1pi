#ifndef STACKEDHIST_PDGCODE
#define STACKEDHIST_PDGCODE

//#include "uboone/CC1pi/Algorithms/PDGEnums.h"
#include "../Algorithms/PDGEnums.h"

#include "TH1F.h"
#include "TH2F.h"
#include "THStack.h"

#include <vector>
#include <string>

class StackedHistPDGCode{

 public:
  StackedHistPDGCode(std::string histname, std::string title, int nbins, double lowlimit, double highlimit);// Constructor
  StackedHistPDGCode(std::string histname, std::string title, int nbinsx, double lowx, double highx, int nbinsy, double lowy, double highy); // Overloaded constructor for 2D histograms
  void Fill(PDGCode particle_pdg, double value);
  void Fill(PDGCode particle_pdg, double value, double weight);
  void Fill2D(PDGCode particle_pdg, double value_x, double value_y);
  void DrawStack(double norm, TCanvas *c1, TString option);
  void DrawOverlay(double norm, TCanvas *c1, TString option);
  void DrawOverlayMuPi(double norm, TCanvas *c1, TString option);
  void Draw2D(TCanvas *c1, TString option);

 protected:
  int nHists;
  std::vector<PDGCode> hist_order; // For keeping track of which hist goes with which PDG code

  THStack *stack;
  TH1F *hists[25];
  TH2F *hists2D[25];

  bool is2Dhists;

  double invalid_total_x;
  double invalid_total_y;

  void InstantiateHistOrder();
  void StyleHistsStack();
  void StyleHistsOverlay();
  void StyleHists2D();
  unsigned int GetHistN(PDGCode particle_pdg);
  PDGCode GetPDGFromHistN(unsigned int hist_n);
};

// Define functions here instead of a .C file so we don't need too many includes

// -------------------------- Constructor -------------------------- //
// Intended to be very general so you can make histograms of whatever you like
StackedHistPDGCode::StackedHistPDGCode(std::string histname, std::string title, int nbins, double lowlimit, double highlimit)
{
  InstantiateHistOrder();
  stack = new THStack(histname.c_str(),title.c_str());
  for (int i_hist=0; i_hist < nHists; i_hist++){
    std::string histname_i = std::string(histname)+std::string("_")+std::to_string(i_hist);
    hists[i_hist] = new TH1F(histname_i.c_str(),title.c_str(),nbins,lowlimit,highlimit);
  }
  is2Dhists = false;
  invalid_total_x = 0.;
}

// -------------------------- Constructor for 2D hists -------------------------- //
// Intended to be very general so you can make histograms of whatever you like
StackedHistPDGCode::StackedHistPDGCode(std::string histname, std::string title, int nbinsx, double lowx, double highx, int nbinsy, double lowy, double highy)
{
  InstantiateHistOrder();
  for (int i_hist=0; i_hist < nHists; i_hist++){
    std::string histname_i = std::string(histname)+std::string("_")+std::to_string(i_hist);
    hists2D[i_hist] = new TH2F(histname_i.c_str(),title.c_str(),nbinsx,lowx,highx,nbinsy,lowy,highy);
  }
  is2Dhists = true;
  invalid_total_x = 0.;
  invalid_total_y = 0.;
}

// Set histogram order (moved to be in a separate function from the constructor to allow for overload of constructor for 1D or 2D histograms)
void StackedHistPDGCode::InstantiateHistOrder(){
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
}

// -------------------------- Function to fill the correct histogram -------------------------- //
void StackedHistPDGCode::Fill(PDGCode particle_pdg, double value)
{
  if (is2Dhists){
    std::cout << "[StackedHistPDGCode] ERROR: cannot call Fill for 2D hists. Call Fill2D instead. Exiting..." << std::endl;
    throw;
  }

  unsigned int n_hist = StackedHistPDGCode::GetHistN(particle_pdg);
  hists[n_hist]->Fill(value);

  if (value==-999 || value==-9999) invalid_total_x++;
}

// -------------------------- Function to fill the correct histogram -------------------------- //
// Overloaded to allow for weights
void StackedHistPDGCode::Fill2D(PDGCode particle_pdg, double value_x, double value_y)
{
  if (!is2Dhists){
    std::cout << "[StackedHistPDGCode] ERROR: cannot call Fill2D for 21D hists. Call Fill instead. Exiting..." << std::endl;
    throw;
  }

  unsigned int n_hist = StackedHistPDGCode::GetHistN(particle_pdg);
  hists2D[n_hist]->Fill(value_x, value_y);

  if (value_x == -999 || value_x == -9999) invalid_total_x++;
  if (value_y == -999 || value_y == -9999) invalid_total_y++;
}

// -------------------------- Function to draw the histograms -------------------------- //
void StackedHistPDGCode::DrawStack(double norm, TCanvas *c1, TString option="")
{
  if (is2Dhists){
    std::cout << "[StackedHistPDGCode] ERROR: cannot call DrawStack for 2D hists. Call Draw2D instead. Exiting..." << std::endl;
    throw;
  }

  // First: style the histograms
  StackedHistPDGCode::StyleHistsStack();

  // Next: add histogramst to the stack and make TLegend
  // Only do this for histograms that have entries
  TLegend *leg = new TLegend(0.1,0.88,0.65,0.99);
  // leg->SetTextFont(132);
  leg->SetLineColor(kWhite);
  leg->SetTextAlign(12);
  leg->SetNColumns(6);

  TPaveText *pt = new TPaveText(0.65,0.88,0.9,0.99,"NDC NB");
  pt->SetLineColor(kWhite);
  pt->SetFillColor(kWhite);
  pt->SetTextAlign(12);

  double underflow_total = 0.;
  double overflow_total = 0.;

  for (int i_hist=0; i_hist < nHists; i_hist++){
    if (hists[i_hist]->GetEntries() == 0) continue;

    hists[i_hist]->Scale(norm);
    stack->Add(hists[i_hist]);
    PDGCode pdg_for_legend = StackedHistPDGCode::GetPDGFromHistN((unsigned int)i_hist);
    leg->AddEntry(hists[i_hist],PDGenum2str(pdg_for_legend).c_str(),"f");

    underflow_total += hists[i_hist]->GetBinContent(0);
    overflow_total += hists[i_hist]->GetBinContent(hists[i_hist]->GetXaxis()->GetNbins()+1);

  }

  pt->AddText(TString::Format("Underflow (Invalid): %.2f (%.2f)",underflow_total,invalid_total_x).Data());
  pt->AddText(TString::Format("Overflow: %.2f",overflow_total).Data());

  c1->cd();
  c1->SetTopMargin(0.13);
  stack->Draw("hist"+option);
  leg->Draw();
  pt->Draw();
}

// -------------------------- Function to draw the histograms -------------------------- //
void StackedHistPDGCode::DrawOverlay(double norm, TCanvas *c1, TString option="")
{
  if (is2Dhists){
    std::cout << "[StackedHistPDGCode] ERROR: cannot call DrawOverlay for 2D hists. Call Draw2D instead. Exiting..." << std::endl;
    throw;
  }

  // First: style the histograms
  StackedHistPDGCode::StyleHistsOverlay();

  // Next: add histogramst to the stack and make TLegend
  // Only do this for histograms that have entries
  TLegend *leg = new TLegend(0.1,0.88,0.65,0.99);
  // leg->SetTextFont(132);
  leg->SetLineColor(kWhite);
  leg->SetTextAlign(12);
  leg->SetNColumns(6);

  TPaveText *pt = new TPaveText(0.65,0.88,0.9,0.99,"NDC NB");
  pt->SetLineColor(kWhite);
  pt->SetFillColor(kWhite);
  pt->SetTextAlign(12);

  double underflow_total = 0.;
  double overflow_total = 0.;

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

    underflow_total += hists[i_hist]->GetBinContent(0);
    overflow_total += hists[i_hist]->GetBinContent(hists[i_hist]->GetXaxis()->GetNbins()+1);

  }

  pt->AddText(TString::Format("Underflow (Invalid): %.2f (%.2f)",underflow_total,invalid_total_x).Data());
  pt->AddText(TString::Format("Overflow: %.2f",overflow_total).Data());

  c1->cd();
  leg->Draw();
  pt->Draw();
}

// -------------------------- Function to draw the histograms -------------------------- //
void StackedHistPDGCode::DrawOverlayMuPi(double norm, TCanvas *c1, TString option="")
{
  if (is2Dhists){
    std::cout << "[StackedHistPDGCode] ERROR: cannot call DrawOverlayMuPi for 2D hists. Call Draw2D instead. Exiting..." << std::endl;
    throw;
  }

  // First: style the histograms
  StackedHistPDGCode::StyleHistsOverlay();

  // Next: make TLegend
  // Only do this for histograms that have entries
  TLegend *leg = new TLegend(0.1,0.88,0.65,0.99);
  // leg->SetTextFont(132);
  leg->SetLineColor(kWhite);
  leg->SetTextAlign(12);
  leg->SetNColumns(2);

  TPaveText *pt = new TPaveText(0.65,0.88,0.9,0.99,"NDC NB");
  pt->SetLineColor(kWhite);
  pt->SetFillColor(kWhite);
  pt->SetTextAlign(12);

  double underflow_total = 0.;
  double overflow_total = 0.;

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

    underflow_total += hists[i_hist]->GetBinContent(0);
    overflow_total += hists[i_hist]->GetBinContent(hists[i_hist]->GetXaxis()->GetNbins()+1);

  }

  pt->AddText(TString::Format("Underflow (Invalid): %.2f (%.2f)",underflow_total,invalid_total_x).Data());
  pt->AddText(TString::Format("Overflow: %.2f",overflow_total).Data());

  c1->cd();
  leg->Draw();
  pt->Draw();

}

// -------------------------- Function to draw the histograms -------------------------- //
void StackedHistPDGCode::Draw2D(TCanvas *c1, TString option="")
{
  if (!is2Dhists){
    std::cout << "[StackedHistPDGCode] ERROR: cannot call Draw2D for 1D hists. Call DrawStack, DrawOverlay, or DrawOverlayMuPi instead. Exiting..." << std::endl;
    throw;
  }

  // First: style the histograms
  StackedHistPDGCode::StyleHists2D();

  // Next: make TLegend
  // Only do this for histograms that have entries
  TLegend *leg = new TLegend(0.1,0.88,0.65,0.99);
  // leg->SetTextFont(132);
  leg->SetLineColor(kWhite);
  leg->SetTextAlign(12);
  leg->SetNColumns(6);

  TPaveText *pt = new TPaveText(0.65,0.88,0.9,0.99,"NDC NB");
  pt->SetLineColor(kWhite);
  pt->SetFillColor(kWhite);
  pt->SetTextAlign(12);

  double underflow_total_x = 0.;
  double overflow_total_x = 0.;
  double underflow_total_y = 0.;
  double overflow_total_y = 0.;

  bool drawn_first = false;

  for (int i_hist=0; i_hist < nHists; i_hist++){
    if (hists2D[i_hist]->GetEntries() == 0) continue;

    PDGCode pdg_for_legend = StackedHistPDGCode::GetPDGFromHistN((unsigned int)i_hist);
    leg->AddEntry(hists2D[i_hist],PDGenum2str(pdg_for_legend).c_str(),"f");

    c1->cd();
    c1->SetTopMargin(0.13);
    if (!drawn_first){
      hists2D[i_hist]->Draw("scat"+option);
      drawn_first = true;
    }
    else{
      hists2D[i_hist]->Draw("scat same"+option);
    }

    underflow_total_x += hists2D[i_hist]->Integral(0,0,0,hists2D[i_hist]->GetYaxis()->GetNbins()+1);
    overflow_total_x += hists2D[i_hist]->Integral(hists2D[i_hist]->GetXaxis()->GetNbins()+1,hists2D[i_hist]->GetXaxis()->GetNbins()+1,0,hists2D[i_hist]->GetYaxis()->GetNbins()+1);
    underflow_total_y += hists2D[i_hist]->Integral(0,hists2D[i_hist]->GetXaxis()->GetNbins()+1,0,0);
    overflow_total_y += hists2D[i_hist]->Integral(0,hists2D[i_hist]->GetXaxis()->GetNbins()+1,hists2D[i_hist]->GetYaxis()->GetNbins()+1,hists2D[i_hist]->GetYaxis()->GetNbins()+1);

  }

  pt->AddText(TString::Format("U/f x (Inv): %.2f (%.2f)",underflow_total_x,invalid_total_x).Data());
  pt->AddText(TString::Format("O/f x: %.2f",overflow_total_x).Data());
  pt->AddText(TString::Format("U/f y (Inv): %.2f (%.2f)",underflow_total_y,invalid_total_y).Data());
  pt->AddText(TString::Format("O/f y: %.2f",overflow_total_y).Data());

  c1->cd();
  leg->Draw();
  pt->Draw();
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

// -------------------------- Function to style the histograms -------------------------- //
// Private: only called by Draw2D function in this file
void StackedHistPDGCode::StyleHists2D()
{
  // Set marker color for all histograms
  hists2D[0] ->SetMarkerColor(kOrange); // kNuMu
  hists2D[1] ->SetMarkerColor(kOrange-3); // kNuMuBar
  hists2D[2] ->SetMarkerColor(kOrange+2); // kNuE
  hists2D[3] ->SetMarkerColor(kRed); // kNuEBar
  hists2D[4] ->SetMarkerColor(kRed+2); // kNuTau
  hists2D[5] ->SetMarkerColor(kPink-7); // kNuTauBar
  hists2D[6] ->SetMarkerColor(kPink+10); // kMuMinus
  hists2D[8] ->SetMarkerColor(kViolet+1); // kPiMinus
  hists2D[7] ->SetMarkerColor(kMagenta+1); // kMuPlus
  hists2D[9] ->SetMarkerColor(kBlue+2); // kPiPlus
  hists2D[10]->SetMarkerColor(kBlue); // kPiZero
  hists2D[11]->SetMarkerColor(kAzure+1); // kElectron
  hists2D[12]->SetMarkerColor(kCyan+2); // kPositron
  hists2D[13]->SetMarkerColor(kCyan); // kTauMinus
  hists2D[14]->SetMarkerColor(kGreen+1); // kTauPlus
  hists2D[15]->SetMarkerColor(kOrange-2); // kPhoton
  hists2D[16]->SetMarkerColor(kGray); // kProton
  hists2D[17]->SetMarkerColor(kGray+2); // kNeutron
  hists2D[18]->SetMarkerColor(kGreen+3); // kKaonPlus
  hists2D[19]->SetMarkerColor(kGreen+2); // kKaonMinus
  hists2D[20]->SetMarkerColor(kBlack); // kPDGUnknown
  // Set marker size for all histograms
  for (size_t i=0; i<nHists; i++){
    hists2D[i]->SetMarkerStyle(20);
    hists2D[i]->SetMarkerSize(.5);
  }
  // Also set fill color for legend
  hists2D[0] ->SetFillColor(kOrange); // kNuMu
  hists2D[1] ->SetFillColor(kOrange-3); // kNuMuBar
  hists2D[2] ->SetFillColor(kOrange+2); // kNuE
  hists2D[3] ->SetFillColor(kRed); // kNuEBar
  hists2D[4] ->SetFillColor(kRed+2); // kNuTau
  hists2D[5] ->SetFillColor(kPink-7); // kNuTauBar
  hists2D[6] ->SetFillColor(kPink+10); // kMuMinus
  hists2D[8] ->SetFillColor(kViolet+1); // kPiMinus
  hists2D[7] ->SetFillColor(kMagenta+1); // kMuPlus
  hists2D[9] ->SetFillColor(kBlue+2); // kPiPlus
  hists2D[10]->SetFillColor(kBlue); // kPiZero
  hists2D[11]->SetFillColor(kAzure+1); // kElectron
  hists2D[12]->SetFillColor(kCyan+2); // kPositron
  hists2D[13]->SetFillColor(kCyan); // kTauMinus
  hists2D[14]->SetFillColor(kGreen+1); // kTauPlus
  hists2D[15]->SetFillColor(kOrange-2); // kPhoton
  hists2D[16]->SetFillColor(kGray); // kProton
  hists2D[17]->SetFillColor(kGray+2); // kNeutron
  hists2D[18]->SetFillColor(kGreen+3); // kKaonPlus
  hists2D[19]->SetFillColor(kGreen+2); // kKaonMinus
  hists2D[20]->SetFillColor(kBlack); // kPDGUnknown
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
