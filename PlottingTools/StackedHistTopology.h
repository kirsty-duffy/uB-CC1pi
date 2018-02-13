#ifndef STACKEDHIST_TOPOLOGY
#define STACKEDHIST_TOPOLOGY

#include "uboone/CC1pi/Algorithms/TopologyEnums.h"

#include "TH1F.h"
#include "THStack.h"

#include <vector>
#include <string>

class StackedHistTopology{

 public:
  StackedHistTopology(std::string histname, std::string title, int nbins, double lowlimit, double highlimit);// Constructor
  void Fill(NuIntTopology topology, double value);
  void Fill(NuIntTopology topology, double value, double weight);
  void DrawStack(double norm, TCanvas *c1);

 protected:
  int nHists;
  std::vector<NuIntTopology> hist_order; // For keeping track of which hist goes with which topology
  
  THStack *stack;
  TH1F *hists[25];

  void StyleHists();
  unsigned int GetHistN(NuIntTopology topology);
  NuIntTopology GetTopologyFromHistN(unsigned int hist_n);
};

// Define functions here instead of a .C file so we don't need too many includes

// -------------------------- Constructor -------------------------- //
// Intended to be very general so you can make histograms of whatever you like
StackedHistTopology::StackedHistTopology(std::string histname, std::string title, int nbins, double lowlimit, double highlimit)
{
  // Topology categories for the histograms
  hist_order.push_back(kCC0pi0p);
  hist_order.push_back(kCC0pi1p);
  hist_order.push_back(kCC0piNp);
  hist_order.push_back(kCC1piplus0p);
  hist_order.push_back(kCC1piplus1p);
  hist_order.push_back(kCC1piplusNp);
  hist_order.push_back(kCC1piminus0p);
  hist_order.push_back(kCC1piminus1p);
  hist_order.push_back(kCC1piminusNp);
  hist_order.push_back(kCC1pizero0p);
  hist_order.push_back(kCC1pizero1p);
  hist_order.push_back(kCC1pizeroNp);
  hist_order.push_back(kCCmultipi0p);
  hist_order.push_back(kCCmultipi1p);
  hist_order.push_back(kCCmultipiNp);
  hist_order.push_back(kCCother);
  hist_order.push_back(kNC);
  hist_order.push_back(kUnknown);

  nHists = hist_order.size();
  
  stack = new THStack(histname.c_str(),title.c_str());
  for (int i_hist=0; i_hist < nHists; i_hist++){
    std::string histname_i = std::string(histname)+std::string("_")+std::to_string(i_hist);
    hists[i_hist] = new TH1F(histname_i.c_str(),"",nbins,lowlimit,highlimit);
  }
}

// -------------------------- Function to fill the correct histogram -------------------------- //
void StackedHistTopology::Fill(NuIntTopology topology, double value)
{
  unsigned int n_hist = StackedHistTopology::GetHistN(topology);
  hists[n_hist]->Fill(value);
}

// -------------------------- Function to fill the correct histogram -------------------------- //
// Overloaded to allow for weights
void StackedHistTopology::Fill(NuIntTopology topology, double value, double weight)
{
  unsigned int n_hist = StackedHistTopology::GetHistN(topology);
  hists[n_hist]->Fill(value, weight);
}

// -------------------------- Function to draw the histograms -------------------------- //
void StackedHistTopology::DrawStack(double norm, TCanvas *c1)
{
  // First: style the histograms
  StackedHistTopology::StyleHists();

  // Next: add histogramst to the stack and make TLegend
  // Only do this for histograms that have entries
  TLegend *leg = new TLegend(0.6,0.7,0.95,0.95);
  
  for (int i_hist=0; i_hist < nHists; i_hist++){
    if (hists[i_hist]->GetEntries() == 0) continue;
    
    hists[i_hist]->Scale(norm);
    stack->Add(hists[i_hist]);
    NuIntTopology topology_for_legend = StackedHistTopology::GetTopologyFromHistN((unsigned int)i_hist);
    leg->AddEntry(hists[i_hist],topologyenum2str(topology_for_legend).c_str(),"f");
    
  }

  c1->cd();
  stack->Draw("hist");
  leg->Draw();
}

// -------------------------- Function to style the histograms -------------------------- //
// Private: only called by DrawStack function in this file
void StackedHistTopology::StyleHists()
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


// ---------------------- Function to get histogram number for given topology ---------------------- //
// Private: only called by functions in this class
unsigned int StackedHistTopology::GetHistN(NuIntTopology topology)
{
  unsigned int HistN;
  bool found_hist=false;

  for (unsigned int i=0; i<nHists; i++){
    if (hist_order.at(i) == topology){
      HistN = i;
      found_hist = true;
      break;
    }
  }

  if (!found_hist){
    std::cout << "[ERROR: StackedHistTopology.h] Could not find histogram for topology " << topology << std::endl;
    HistN = 9999;
  }

  return HistN;
}


// ---------------------- Function to get topology code for given histogram number ---------------------- //
// Private: only called by functions in this class
NuIntTopology StackedHistTopology::GetTopologyFromHistN(unsigned int hist_n)
{
  return hist_order.at(hist_n);
}


#endif
