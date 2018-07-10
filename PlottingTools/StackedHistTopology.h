#ifndef STACKEDHIST_TOPOLOGY
#define STACKEDHIST_TOPOLOGY

//#include "uboone/CC1pi/Algorithms/TopologyEnums.h"
#include "../Algorithms/TopologyEnums.h"

#include "TH1F.h"
#include "THStack.h"

#include <vector>
#include <string>

class StackedHistTopology{

 public:
  StackedHistTopology(std::string histname, std::string title, int nbins, double lowlimit, double highlimit);// Constructor
  void Fill(NuIntTopology topology, double value);
  void Fill(NuIntTopology topology, double value, double weight);
  void DrawStack(double norm, TCanvas *c1, bool coarse=false);

 protected:
  int nHists;
  int nHists_coarse;
  std::vector<NuIntTopology> hist_order; // For keeping track of which hist goes with which topology
  std::vector<NuIntTopology> hist_order_coarse; // For keeping track of which hist goes with which topology

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
  hist_order.push_back(kCCNue);
  hist_order.push_back(kNC);
  hist_order.push_back(kOutFV);
  hist_order.push_back(kCosmic);
  hist_order.push_back(kMixed);
  hist_order.push_back(kUnknown);

  nHists = hist_order.size();

  // Topology categories for the histograms
  hist_order_coarse.push_back(kCC0pi0p);
  hist_order_coarse.push_back(kCC0pi1p);
  hist_order_coarse.push_back(kCC0piNp);
  hist_order_coarse.push_back(kCC1piplus0p);
  hist_order_coarse.push_back(kCC1piplus1p);
  hist_order_coarse.push_back(kCC1piplusNp);
  hist_order_coarse.push_back(kCC1piminus0p);
  hist_order_coarse.push_back(kCC1piminus1p);
  hist_order_coarse.push_back(kCC1piminusNp);
  hist_order_coarse.push_back(kCC1pizero0p);
  hist_order_coarse.push_back(kCC1pizero1p);
  hist_order_coarse.push_back(kCC1pizeroNp);
  hist_order_coarse.push_back(kCCmultipi0p);
  hist_order_coarse.push_back(kCCmultipi1p);
  hist_order_coarse.push_back(kCCmultipiNp);
  hist_order_coarse.push_back(kCCother);
  hist_order_coarse.push_back(kCCNue);
  hist_order_coarse.push_back(kNC);
  hist_order_coarse.push_back(kOutFV);
  hist_order_coarse.push_back(kCosmic);
  hist_order_coarse.push_back(kMixed);
  hist_order_coarse.push_back(kUnknown);

  nHists = hist_order_coarse.size();

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
void StackedHistTopology::DrawStack(double norm, TCanvas *c1, bool coarse)
{
  // First: style the histograms
  StackedHistTopology::StyleHists();

  // Next: add histogramst to the stack and make TLegend
  // Only do this for histograms that have entries
  TLegend *leg = new TLegend(0.55,0.7,0.95,0.95);

  leg -> SetNColumns(2);

  if (coarse) {
  std::map<std::string, TH1F *> coarse_histos;

  // Collapse fine topo enums onto coarse topo names by using
  // topologyenum2str_coarse
  for (int i_hist = nHists-1; i_hist >=0; i_hist--) {
    hists[i_hist]->Scale(norm);
    NuIntTopology topology_for_legend =
        StackedHistTopology::GetTopologyFromHistN((unsigned int)i_hist);
    std::string coarse_topo_legend_title =
        topologyenum2str_coarse(topology_for_legend);

    if (coarse_histos.find(coarse_topo_legend_title) ==
        coarse_histos.end()) { // Haven't already got one in this coarse topo
      coarse_histos[coarse_topo_legend_title] =
          static_cast<TH1F *>(hists[i_hist]->Clone());
    } else { // Already have one of this topo, so just add it
      coarse_histos[coarse_topo_legend_title]->Add(hists[i_hist]);
    }
  }

  for (std::pair<std::string, TH1F *> ch : coarse_histos) {
    stack->Add(ch.second);
    leg->AddEntry(ch.second, ch.first.c_str(), "f");
  }

}
else { // fine
  for (int i_hist = 0; i_hist < nHists; i_hist++) {
    if (hists[i_hist]->GetEntries() == 0)
      continue;

    hists[i_hist]->Scale(norm);
    stack->Add(hists[i_hist]);
    NuIntTopology topology_for_legend =
        StackedHistTopology::GetTopologyFromHistN((unsigned int)i_hist);
    leg->AddEntry(hists[i_hist],
                  topologyenum2str(topology_for_legend).c_str(), "f");
  }
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
  hists[0] ->SetFillColor(kOrange); // CC0pi0p
  hists[1] ->SetFillColor(kOrange-3); // CC0pi1p
  hists[2] ->SetFillColor(kOrange+2); // CC0piNp
  hists[3] ->SetFillColor(kRed); // CC1piplus0p
  hists[4] ->SetFillColor(kRed+2); // CC1piplus1p
  hists[5] ->SetFillColor(kPink-7); // CC1piplusNp
  hists[6] ->SetFillColor(kPink+10); // CC1piminus0p
  hists[7] ->SetFillColor(kMagenta+1); // CC1piminus1p
  hists[8] ->SetFillColor(kViolet+1); // CC1piminusNp
  hists[9] ->SetFillColor(kBlue+2); // CC1pizero0p
  hists[10]->SetFillColor(kBlue); // CC1pizero1p
  hists[11]->SetFillColor(kAzure+1); // CC1pizeroNp
  hists[12]->SetFillColor(kCyan+2); // CCmultipi0p
  hists[13]->SetFillColor(kCyan); // CCmultipi1p
  hists[14]->SetFillColor(kGreen+1); // CCmultipiNp
  hists[15]->SetFillColor(kGreen+3); // CCother
  hists[16]->SetFillColor(kYellow+1); // CCnue
  hists[17]->SetFillColor(kBlue-10); // NC
  hists[18]->SetFillColor(kBlue-5); // outFV
  hists[19]->SetFillColor(kGray); // Cosmic
  hists[20]->SetFillColor(kGray+2); // Mixed
  hists[21]->SetFillColor(kBlack); // Unknown
}


// ---------------------- Function to get histogram number for given topology ---------------------- //
// Private: only called by functions in this class
unsigned int StackedHistTopology::GetHistN(NuIntTopology topology)
{
  unsigned int HistN;
  bool found_hist=false;

  for (int i=0; i<nHists; i++){
    if (hist_order.at(i) == topology){
      HistN = i;
      found_hist = true;
      break;
    }
  }

  if (!found_hist){
    std::cout << "[ERROR: StackedHistTopology.h] Could not find histogram for topology " << topology << std::endl;
    HistN = nHists-1;
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
