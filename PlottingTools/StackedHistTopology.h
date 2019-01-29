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
  StackedHistTopology(std::string histname, std::string title, int nbinsx, double lowx, double highx, int nbinsy, double lowy, double highy); // Overloaded constructor for 2D histograms
  void Fill(NuIntTopology topology, double value);
  void Fill(NuIntTopology topology, double value, double weight);
  void Fill2D(NuIntTopology topology, double value_x, double value_y);
  void DrawStack(double mc_scaling, TCanvas *c1, bool coarse=true, TH1F *onbeam_h=nullptr, TH1F *offbeam_h=nullptr, double offbeam_scaling=1.0, bool onminusoffbeam=false);
  void PrintHistIntegrals(bool coarse=true);
  double GetCC1piIntegral();
  double GetTotalIntegral();
  void Draw2D(TCanvas *c1, bool coarse=true, TString option="");

 protected:
  int nHists;
  int nHists_coarse;
  std::vector<NuIntTopology> hist_order; // For keeping track of which hist goes with which topology

  std::map<std::string, TH1F *> coarse_histos; // For keeping track of "coarse" topology histograms (e.g. CC1pi instead of CC1pi1p)

  std::map<std::string, TH2F *> coarse_histos_2D; // For keeping track of "coarse" topology histograms (e.g. CC1pi instead of CC1pi1p)

  THStack *stack;
  TH1F *hists[25];
  TH2F *hists2D[25];

  bool is2Dhists;

  double invalid_total_x;
  double invalid_total_y;

  void InstantiateHistOrder();
  void StyleHistsStack();
  void StyleHists2D();
  unsigned int GetHistN(NuIntTopology topology);
  NuIntTopology GetTopologyFromHistN(unsigned int hist_n);
  void GenerateCoarseHistos();
};

// Define functions here instead of a .C file so we don't need too many includes

// -------------------------- Constructor -------------------------- //
// Intended to be very general so you can make histograms of whatever you like
StackedHistTopology::StackedHistTopology(std::string histname, std::string title, int nbins, double lowlimit, double highlimit)
{
  InstantiateHistOrder();
  stack = new THStack(histname.c_str(),title.c_str());
  for (int i_hist=0; i_hist < nHists; i_hist++){
    std::string histname_i = std::string(histname)+std::string("_")+std::to_string(i_hist);
    hists[i_hist] = new TH1F(histname_i.c_str(),"",nbins,lowlimit,highlimit);
  }

  // Style the histograms
  StyleHistsStack();

  is2Dhists = false;
  invalid_total_x = 0.;
}

// -------------------------- Constructor for 2D hists -------------------------- //
// Intended to be very general so you can make histograms of whatever you like
StackedHistTopology::StackedHistTopology(std::string histname, std::string title, int nbinsx, double lowlimitx, double highlimitx, int nbinsy, double lowlimity, double highlimity)
{
  InstantiateHistOrder();
  for (int i_hist=0; i_hist < nHists; i_hist++){
    std::string histname_i = std::string(histname)+std::string("_")+std::to_string(i_hist);
    hists2D[i_hist] = new TH2F(histname_i.c_str(),title.c_str(),nbinsx,lowlimitx,highlimitx,nbinsy,lowlimity,highlimity);
  }

  // Style the histograms
  StyleHists2D();

  is2Dhists = true;
  invalid_total_x = 0.;
  invalid_total_y = 0.;
}

// Set histogram order (moved to be in a separate function from the constructor to allow for overload of constructor for 1D or 2D histograms)
void StackedHistTopology::InstantiateHistOrder(){
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
}

// -------------------------- Function to fill the correct histogram -------------------------- //
void StackedHistTopology::Fill(NuIntTopology topology, double value)
{
  if (is2Dhists){
    std::cout << "[StackedHistTopology] ERROR: cannot call Fill for 2D hists. Call Fill2D instead. Exiting..." << std::endl;
    throw;
  }

  unsigned int n_hist = StackedHistTopology::GetHistN(topology);
  hists[n_hist]->Fill(value);

  if (value==-999 || value==-9999) invalid_total_x++;
}

// -------------------------- Function to fill the correct histogram -------------------------- //
// Overloaded to allow for weights
void StackedHistTopology::Fill(NuIntTopology topology, double value, double weight)
{
  unsigned int n_hist = StackedHistTopology::GetHistN(topology);
  hists[n_hist]->Fill(value, weight);

  if (value==-999 || value==-9999) invalid_total_x++;
}

// -------------------------- Function to fill the correct histogram -------------------------- //
// For 2D histograms
void StackedHistTopology::Fill2D(NuIntTopology topology, double value_x, double value_y)
{
  if (!is2Dhists){
    std::cout << "[StackedHistTopology] ERROR: cannot call Fill2D for 21D hists. Call Fill instead. Exiting..." << std::endl;
    throw;
  }

  unsigned int n_hist = StackedHistTopology::GetHistN(topology);
  hists2D[n_hist]->Fill(value_x, value_y);

  if (value_x == -999 || value_x == -9999) invalid_total_x++;
  if (value_y == -999 || value_y == -9999) invalid_total_y++;
}

// -------------------------- Function to collapse fine topo enums onto coarse topo names -------------------------- //
void StackedHistTopology::GenerateCoarseHistos()
{
  if (!is2Dhists){
    // Clear coarse histos and remake
    coarse_histos.clear();


    for (int i_hist = nHists-1; i_hist >=0; i_hist--) {
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
  }
  else{
    // Clear coarse histos and remake
    coarse_histos_2D.clear();

    for (int i_hist = nHists-1; i_hist >=0; i_hist--) {
      NuIntTopology topology_for_legend =
          StackedHistTopology::GetTopologyFromHistN((unsigned int)i_hist);
      std::string coarse_topo_legend_title =
          topologyenum2str_coarse(topology_for_legend);

      if (coarse_histos_2D.find(coarse_topo_legend_title) ==
          coarse_histos_2D.end()) { // Haven't already got one in this coarse topo
        coarse_histos_2D[coarse_topo_legend_title] =
            static_cast<TH2F *>(hists2D[i_hist]->Clone());
      } else { // Already have one of this topo, so just add it
        coarse_histos_2D[coarse_topo_legend_title]->Add(hists2D[i_hist]);
      }
    }
  }
};

// -------------------------- Function to draw the histograms -------------------------- //
void StackedHistTopology::DrawStack(double mc_scaling, TCanvas *c1, bool coarse=true, TH1F *onbeam_h=nullptr, TH1F *offbeam_h=nullptr, double offbeam_scaling=1.0, bool onminusoffbeam=false)
{

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

  if (mc_scaling==0) mc_scaling = 1.0;

  for (int i_hist = 0; i_hist < nHists; i_hist++){
    hists[i_hist]->Scale(mc_scaling);
  }

  if (coarse) {

  StackedHistTopology::GenerateCoarseHistos();

  for (std::pair<std::string, TH1F *> ch : coarse_histos) {
    stack->Add(ch.second);
    leg->AddEntry(ch.second, ch.first.c_str(), "f");
    //std::cout << "Integral for topology " << ch.first.c_str() << ": " << ch.second->Integral() << std::endl;

    underflow_total += ch.second->GetBinContent(0);
    overflow_total += ch.second->GetBinContent(ch.second->GetXaxis()->GetNbins()+1);
  }
}
else { // fine
  for (int i_hist = 0; i_hist < nHists; i_hist++) {
    if (hists[i_hist]->GetEntries() == 0)
      continue;

    stack->Add(hists[i_hist]);
    NuIntTopology topology_for_legend =
        StackedHistTopology::GetTopologyFromHistN((unsigned int)i_hist);
    leg->AddEntry(hists[i_hist],
                  topologyenum2str(topology_for_legend).c_str(), "f");
    //std::cout << "Integral for topology " << topologyenum2str(topology_for_legend).c_str() << ": " << hists[i_hist]->Integral() << std::endl;

    underflow_total += hists[i_hist]->GetBinContent(0);
    overflow_total += hists[i_hist]->GetBinContent(hists[i_hist]->GetXaxis()->GetNbins()+1);
  }
} // end if coarse/fine

  pt->AddText(TString::Format("Underflow (Invalid): %.2f (%.2f)",underflow_total,invalid_total_x).Data());
  pt->AddText(TString::Format("Overflow: %.2f",overflow_total).Data());

  // If data histograms are given, deal with them here
  TH1F *datatodraw;
  if (onbeam_h && offbeam_h){
    offbeam_h->Sumw2();
    offbeam_h->Scale(offbeam_scaling);
    datatodraw = (TH1F*)onbeam_h->Clone();
    datatodraw->Sumw2();
    datatodraw->SetMarkerStyle(20);
    datatodraw->SetMarkerSize(0.6);
    datatodraw->SetFillStyle(0);
    if (onminusoffbeam){
      datatodraw->Add(offbeam_h,-1);
      leg->AddEntry(datatodraw,"Data (on-off beam)","lp");
    }
    else{
      offbeam_h->SetFillStyle(3345);
      offbeam_h->SetFillColor(kBlack);
      offbeam_h->SetLineColor(kBlack);
      stack->Add(offbeam_h);
      leg->AddEntry(offbeam_h,"Beam-off data","f");
      leg->AddEntry(datatodraw,"Beam-on data","lp");
    }
    pt->AddText(TString::Format("Beam-off Data Underflow/Overflow: %.2f/%.2f",offbeam_h->GetBinContent(0),offbeam_h->GetBinContent(offbeam_h->GetXaxis()->GetNbins()+1)).Data());
    pt->AddText(TString::Format("Beam-on Data Underflow/Overflow: %.2f/%.2f",onbeam_h->GetBinContent(0),onbeam_h->GetBinContent(onbeam_h->GetXaxis()->GetNbins()+1)).Data());


    c1->cd();
    TPad *topPad = new TPad("topPad", "", 0.005, 0.3, 0.995, 0.995);
    TPad *bottomPad = new TPad("bottomPad", "", 0.005, 0.005, 0.995, 0.3);
    topPad->SetTopMargin(0.13);
    topPad->SetBottomMargin(0.02);
    bottomPad->SetTopMargin(0.0);
    bottomPad->SetBottomMargin(0.20);
    bottomPad->SetGridy();
    topPad->Draw();
    bottomPad->Draw();

    topPad->cd();
    stack->Draw("hist");
    stack->GetXaxis()->SetTitleSize(0);
    stack->GetXaxis()->SetLabelSize(0);
    datatodraw->Draw("same p E1");
    leg->Draw();
    pt->Draw();

    bottomPad->cd();
    TH1* MCtotal = (TH1*)stack->GetStack()->Last()->Clone("MCtotal");
    TH1* dataratio = (TH1*)datatodraw->Clone("dataratio");
    dataratio->Divide(MCtotal);
    dataratio->GetXaxis()->SetTitle(stack->GetXaxis()->GetTitle());
    dataratio->GetYaxis()->SetTitle("Data/MC ratio");
    dataratio->GetYaxis()->SetTitleOffset(0.57);
    dataratio->GetYaxis()->SetTitleSize(0.07);
    dataratio->GetYaxis()->SetLabelSize(0.075);
    dataratio->GetXaxis()->SetTitleSize(0.09);
    dataratio->GetXaxis()->SetLabelSize(0.075);
    TLine *l1 = new TLine(dataratio->GetXaxis()->GetBinLowEdge(1),1.0,dataratio->GetXaxis()->GetBinUpEdge(dataratio->GetXaxis()->GetNbins()),1.0);
    l1->SetLineStyle(2);
    l1->SetLineColor(kBlack);
    dataratio->Draw("p E1");
    l1->Draw("same");

  } // end if (data histograms)
else{
    c1->cd();
    c1->SetTopMargin(0.13);
    stack->Draw("hist");
    leg->Draw();
    pt->Draw();
  }

  // TList * histKeys = stack->GetHists();
  // TIter next(histKeys);
  // TObject* object = 0;
  // double total_integral = 0.;
  //
  // while ((object = next()))
  // {
  //   total_integral += ((TH1*)object)->Integral();
  // }
  // std::cout << "------------------------------------------" << std::endl;
  // std::cout << "Total Integral over all topologies: " << total_integral << std::endl;
  // std::cout << "------------------------------------------" << std::endl;
}



// -------------------------- Function to draw the histograms -------------------------- //
void StackedHistTopology::Draw2D(TCanvas *c1, bool coarse=true, TString option="")
{
  if (!is2Dhists){
    std::cout << "[StackedHistTopology] ERROR: cannot call Draw2D for 1D hists. Call DrawStack instead. Exiting..." << std::endl;
    throw;
  }

  // First: style the histograms
  StyleHists2D();

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

  if (coarse) {

  StackedHistTopology::GenerateCoarseHistos();

  for (std::pair<std::string, TH2F *> ch : coarse_histos_2D) {
      if (ch.second->GetEntries() == 0) continue;

      c1->cd();
      c1->SetTopMargin(0.13);
      if (!drawn_first){
        ch.second->Draw("scat"+option);
        drawn_first = true;
      }
      else{
        ch.second->Draw("scat same"+option);
      }

    leg->AddEntry(ch.second, ch.first.c_str(), "f");
    //std::cout << "Integral for topology " << ch.first.c_str() << ": " << ch.second->Integral() << std::endl;

    underflow_total_x += ch.second->Integral(0,0,0,ch.second->GetYaxis()->GetNbins()+1);
    overflow_total_x += ch.second->Integral(ch.second->GetXaxis()->GetNbins()+1,ch.second->GetXaxis()->GetNbins()+1,0,ch.second->GetYaxis()->GetNbins()+1);
    underflow_total_y += ch.second->Integral(0,ch.second->GetXaxis()->GetNbins()+1,0,0);
    overflow_total_y += ch.second->Integral(0,ch.second->GetXaxis()->GetNbins()+1,ch.second->GetYaxis()->GetNbins()+1,ch.second->GetYaxis()->GetNbins()+1);
    }
  }
  else { // fine
    for (int i_hist = 0; i_hist < nHists; i_hist++) {
      if (hists2D[i_hist]->GetEntries() == 0)
        continue;

        c1->cd();
        c1->SetTopMargin(0.13);
        if (!drawn_first){
          hists2D[i_hist]->Draw("scat"+option);
          drawn_first = true;
        }
        else{
          hists2D[i_hist]->Draw("scat same"+option);
        }

      NuIntTopology topology_for_legend =
          StackedHistTopology::GetTopologyFromHistN((unsigned int)i_hist);
      leg->AddEntry(hists2D[i_hist],
                    topologyenum2str(topology_for_legend).c_str(), "f");
      //std::cout << "Integral for topology " << topologyenum2str(topology_for_legend).c_str() << ": " << hists[i_hist]->Integral() << std::endl;

      underflow_total_x += hists2D[i_hist]->Integral(0,0,0,hists2D[i_hist]->GetYaxis()->GetNbins()+1);
      overflow_total_x += hists2D[i_hist]->Integral(hists2D[i_hist]->GetXaxis()->GetNbins()+1,hists2D[i_hist]->GetXaxis()->GetNbins()+1,0,hists2D[i_hist]->GetYaxis()->GetNbins()+1);
      underflow_total_y += hists2D[i_hist]->Integral(0,hists2D[i_hist]->GetXaxis()->GetNbins()+1,0,0);
      overflow_total_y += hists2D[i_hist]->Integral(0,hists2D[i_hist]->GetXaxis()->GetNbins()+1,hists2D[i_hist]->GetYaxis()->GetNbins()+1,hists2D[i_hist]->GetYaxis()->GetNbins()+1);
    }
  } // end if coarse/fine

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
void StackedHistTopology::StyleHistsStack()
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

// -------------------------- Function to style the histograms -------------------------- //
// Private: only called by DrawStack function in this file
void StackedHistTopology::StyleHists2D()
{
  // Set fill color for all histograms
  hists2D[0] ->SetMarkerColor(kOrange); // CC0pi0p
  hists2D[1] ->SetMarkerColor(kOrange-3); // CC0pi1p
  hists2D[2] ->SetMarkerColor(kOrange+2); // CC0piNp
  hists2D[3] ->SetMarkerColor(kRed); // CC1piplus0p
  hists2D[4] ->SetMarkerColor(kRed+2); // CC1piplus1p
  hists2D[5] ->SetMarkerColor(kPink-7); // CC1piplusNp
  hists2D[6] ->SetMarkerColor(kPink+10); // CC1piminus0p
  hists2D[7] ->SetMarkerColor(kMagenta+1); // CC1piminus1p
  hists2D[8] ->SetMarkerColor(kViolet+1); // CC1piminusNp
  hists2D[9] ->SetMarkerColor(kBlue+2); // CC1pizero0p
  hists2D[10]->SetMarkerColor(kBlue); // CC1pizero1p
  hists2D[11]->SetMarkerColor(kAzure+1); // CC1pizeroNp
  hists2D[12]->SetMarkerColor(kCyan+2); // CCmultipi0p
  hists2D[13]->SetMarkerColor(kCyan); // CCmultipi1p
  hists2D[14]->SetMarkerColor(kGreen+1); // CCmultipiNp
  hists2D[15]->SetMarkerColor(kGreen+3); // CCother
  hists2D[16]->SetMarkerColor(kYellow+1); // CCnue
  hists2D[17]->SetMarkerColor(kBlue-10); // NC
  hists2D[18]->SetMarkerColor(kBlue-5); // outFV
  hists2D[19]->SetMarkerColor(kGray); // Cosmic
  hists2D[20]->SetMarkerColor(kGray+2); // Mixed
  hists2D[21]->SetMarkerColor(kBlack); // Unknown
  // Set fill color for all histograms
  hists2D[0] ->SetFillColor(kOrange); // CC0pi0p
  hists2D[1] ->SetFillColor(kOrange-3); // CC0pi1p
  hists2D[2] ->SetFillColor(kOrange+2); // CC0piNp
  hists2D[3] ->SetFillColor(kRed); // CC1piplus0p
  hists2D[4] ->SetFillColor(kRed+2); // CC1piplus1p
  hists2D[5] ->SetFillColor(kPink-7); // CC1piplusNp
  hists2D[6] ->SetFillColor(kPink+10); // CC1piminus0p
  hists2D[7] ->SetFillColor(kMagenta+1); // CC1piminus1p
  hists2D[8] ->SetFillColor(kViolet+1); // CC1piminusNp
  hists2D[9] ->SetFillColor(kBlue+2); // CC1pizero0p
  hists2D[10]->SetFillColor(kBlue); // CC1pizero1p
  hists2D[11]->SetFillColor(kAzure+1); // CC1pizeroNp
  hists2D[12]->SetFillColor(kCyan+2); // CCmultipi0p
  hists2D[13]->SetFillColor(kCyan); // CCmultipi1p
  hists2D[14]->SetFillColor(kGreen+1); // CCmultipiNp
  hists2D[15]->SetFillColor(kGreen+3); // CCother
  hists2D[16]->SetFillColor(kYellow+1); // CCnue
  hists2D[17]->SetFillColor(kBlue-10); // NC
  hists2D[18]->SetFillColor(kBlue-5); // outFV
  hists2D[19]->SetFillColor(kGray); // Cosmic
  hists2D[20]->SetFillColor(kGray+2); // Mixed
  hists2D[21]->SetFillColor(kBlack); // Unknown
  // Set marker size for all histograms
  for (size_t i=0; i<nHists; i++){
    hists2D[i]->SetMarkerStyle(20);
    hists2D[i]->SetMarkerSize(.5);
  }
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

// ---------------------- Function to get histogram integrals ---------------------- //
void StackedHistTopology::PrintHistIntegrals(bool coarse)
{
  // Calculate and print out relative integrals (percentage of events that are each topology)
   if (coarse) {
     StackedHistTopology::GenerateCoarseHistos();

    double total_integral = 0.;
   for (std::pair<std::string, TH1F *> ch : coarse_histos) {
     total_integral += ch.second->Integral();
   }

   for (std::pair<std::string, TH1F *> ch : coarse_histos) {
     std::cout << "Integral for topology " << ch.first.c_str() << ": " << ch.second->Integral()/total_integral << std::endl;
   }
 }
 else{
   // Compute total integral
   double total_integral = 0;
   for (int i_hist=0; i_hist < nHists; i_hist++){
      if (hists[i_hist]->GetEntries() == 0) continue;

      double integral = hists[i_hist]->Integral();
      total_integral += integral;
   }
   for (int i_hist=0; i_hist < nHists; i_hist++){
      if (hists[i_hist]->GetEntries() == 0) continue;

      NuIntTopology topology = StackedHistTopology::GetTopologyFromHistN((unsigned int)i_hist);

      double integral = hists[i_hist]->Integral();

      std::cout << "Integral for toplogy " << topologyenum2str(topology) << ": " << integral/total_integral << std::endl;

   }
 }
}

// ---------------------- Function to get histogram integrals ---------------------- //
double StackedHistTopology::GetTotalIntegral()
{
   // Compute total integral
   double total_integral = 0;
   for (int i_hist=0; i_hist < nHists; i_hist++){
      if (hists[i_hist]->GetEntries() == 0) continue;

      double integral = hists[i_hist]->Integral();
      total_integral += integral;
   }
   return total_integral;
}

// ---------------------- Function to get total number of CC1pi events ---------------------- //
double StackedHistTopology::GetCC1piIntegral()
{
  // Do "coarse" integral because we want all CC1pi events, and don't care about subcategories
  StackedHistTopology::GenerateCoarseHistos();

   for (std::pair<std::string, TH1F *> ch : coarse_histos) {
     if (ch.first == "#nu_{#mu} CC 1#pi^{+}") return ch.second->Integral();
   }

   return 0.;
}

#endif
