#include "MakePIDPlots_header.h"

// What variables do we want these plots as a function of?
std::vector<std::vector<std::vector<double>>> GetPIDvarstoplot(treevars *vars){
  std::vector<std::vector<std::vector<double>>> varstoplot = {
    *(vars->TPCObj_PFP_max_Lp)
    ,*(vars->TPCObj_PFP_max_Lmu)
    ,*(vars->TPCObj_PFP_max_Lpi)
    ,*(vars->TPCObj_PFP_max_Lmip)
    ,*(vars->TPCObj_PFP_max_Lmaxmumip)
    ,*(vars->TPCObj_PFP_track_Chi2Muon)
    ,*(vars->TPCObj_PFP_track_Chi2Proton)
    ,*(vars->TPCObj_PFP_track_Chi2Pion)
    ,*(vars->TPCObj_PFP_PIDA)
    ,*(vars->TPCObj_PFP_Lmuoverp)
    ,*(vars->TPCObj_PFP_Lmipoverp)
    ,*(vars->TPCObj_PFP_Lmaxmumipoverp)
    ,*(vars->TPCObj_PFP_chi2_muminusp)
    ,*(vars->TPCObj_PFP_Lmu_0to1)
    ,*(vars->TPCObj_PFP_Lmip_0to1)
    ,*(vars->TPCObj_PFP_Lpi_0to1)
    ,*(vars->TPCObj_PFP_Lp_0to1)
    ,*(vars->TPCObj_PFP_Lmumip_0to1)
    ,*(vars->TPCObj_PFP_Lmumippi_0to1)
    ,*(vars->TPCObj_PFP_depE_minus_rangeE_mu)
    ,*(vars->TPCObj_PFP_depE_minus_rangeE_p)
  };
  return varstoplot;
};

// How many planes do we want to consider?
// 4 planes: plane0 is plane0, plane1 is plane1, plane2 is plane2, plane3 is (plane0+plane1+plane2)/3
const size_t nplanes = 4;


// Binning (nbins, binlow, binhigh) in the same order as the vector above
std::vector<std::vector<double>> bins = {
                    {20,0,0.6}, // track_max_neglogl_p
                    {20,0,0.6}, // track_max_neglogl_mu
                    {20,0,0.6}, // track_max_neglogl_pi
                    {20,0,0.6}, // track_max_neglogl_mip
                    {20,0,0.6}, // track_max_neglogl_minmumip
                    {25,0,125}, // track_chi2mu
                    {15,0,300}, // track_chi2p
                    {25,0,125}, // track_chi2pi
                    {20,0,30}, // track_PIDA
                    {30,0,60}, // track_neglogl_muoverp
                    {30,0,60}, // track_neglogl_mipoverp
                    {30,0,60}, // track_neglogl_minmumipoverp
                    {25,-400,100}, // track_chi2_muminusp
                    {25,0,1}, // track_Lmu_0to1
                    {25,0,1}, // track_Lmip_0to1
                    {25,0,1}, // track_Lpi_0to1
                    {25,0,1}, // track_Lp_0to1
                    {25,0,1}, // track_Lmumip_0to1
                    {25,0,1}, // track_Lmumippi_0to1
                    {25,-100,100}, // track_depE_minus_rangeE_mu
                    {30,-200,100} // track_depE_minus_rangeE_p
                    };

// Histogram titles in the same order as the vector above
std::vector<std::string> histtitles = {
                    ";L_{p};",
                    ";L_{#mu};",
                    ";L_{#pi};",
                    ";L_{MIP};",
                    ";L_{#mu/MIP};",
                    ";#chi^{2}_{#mu};",
                    ";#chi^{2}_{p};",
                    ";#chi^{2}_{#pi};",
                    ";PIDa (by median);",
                    ";(L_{#mu})/(L_{p});",
                    ";(L_{MIP})/(L_{p});",
                    ";(L_{#mu/MIP})/(L_{p});",
                    ";#chi^{2}_{#mu}-#chi^{2}_{p};",
                    ";L_{#mu}/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";L_{MIP}/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";L_{#pi}/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";L_{p}/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";(L_{#mu}+L_{MIP})/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";(L_{#mu}+L_{MIP}+L_{#pi})/(L_{#mu}+L_{MIP}+L_{#pi}+L_{p}+L_{K});",
                    ";Dep. E - E. by range (muon assumption) [MeV];",
                    ";Dep. E - E. by range (proton assumption) [MeV];"
                  };

// What to call saved plots in the same order as the vector above
std::vector<std::string> histnames = {
                  "effpur_Lp",
                  "effpur_Lmu",
                  "effpur_Lpi",
                  "effpur_Lmip",
                  "effpur_Lmumip",
                  "effpur_chi2mu",
                  "effpur_chi2p",
                  "effpur_chi2pi",
                  "effpur_pida_median",
                  "effpur_Lmuoverp",
                  "effpur_Lmipoverp",
                  "effpur_Lmumipoverp",
                  "effpur_chi2muminusp",
                  "effpur_Lmu0to1",
                  "effpur_Lmip0to1",
                  "effpur_Lpi0to1",
                  "effpur_Lp0to1",
                  "effpur_Lmumip0to1",
                  "effpur_Lmumippi0to1",
                  "effpur_depErangeEmu",
                  "effpur_depErangeEp"
                };

// For efficiency/purity we need to know whether MIPs are supposed to be low or high. In the same order as the vector above
std::vector<bool> MIPlow = {
                    true, // L_p
                    false, // L_mu
                    false, // L_pi
                    false, // L_mip
                    false, // L_minmumip
                    true, // chi2mu
                    false, // chi2p
                    true, // chi2pi
                    true, // PIDA_median
                    false, // L_muoverp
                    false, // L_mipoverp
                    false, // L_mumipoverp
                    true, // chi2_muminusp
                    false, // Lmu_0to1
                    false, // Lmip_0to1
                    false, // Lpi_0to1
                    true, // Lp_0to1
                    false, // Lmumip_0to1
                    false, // Lmumippi_0to1
                    true, // depE_minus_rangeE_mu
                    true // depE_minus_rangeE_p
                    };

// ---------------------------------------------------- //
//  Now the function starts
// ---------------------------------------------------- //

void MakePIDEfficiencyPlots(std::string mcfile){

  gStyle->SetTitleX(0.5);
  gStyle->SetTitleAlign(23);
  gStyle->SetOptTitle(1);
  gStyle->SetTitleBorderSize(0.);

  TFile *f_bnbcos = new TFile(mcfile.c_str(), "read");
  TTree *t_bnbcos = (TTree*)f_bnbcos->Get("cc1piselec/outtree");
  treevars mc_vars;
  settreevars(t_bnbcos,&mc_vars);

  // Sanity check: the plot vectors should be the same size
  t_bnbcos->GetEntry(0);
  CalcPIDvars(&mc_vars);
  std::vector<std::vector<std::vector<double>>> PIDvarstoplot_dummy = GetPIDvarstoplot(&mc_vars);
  // if (PIDvarstoplot_dummy.size() != bins.size()) std::cout << "WARNING PIDvarstoplot_dummy.size() = " << PIDvarstoplot_dummy.size() << "and bins.size() = " << bins.size() << ". This is going to cause you problems!" << std::endl;
  std::cout << "PIDvarstoplot.size() = " << PIDvarstoplot_dummy.size() << std::endl;
  std::cout << "bins.size() = " << bins.size() << std::endl;
  std::cout << "histtitles.size() = " << histtitles.size() << std::endl;
  std::cout << "histnames.size() = " << histnames.size() << std::endl;
  std::cout << "MIPlow.size() = " << MIPlow.size() << std::endl;

  // ----------------- MC

  // Make histograms to fill
  const size_t nplots = PIDvarstoplot_dummy.size();
  hist1D *mc_hists[nplanes][nplots];
  histCC1piselEffPur *mc_hists_cc1pieffpur[nplanes][nplots];
  for (size_t i_pl=0; i_pl<nplanes; i_pl++){
    for (size_t i_h=0; i_h<nplots; i_h++){
      mc_hists[i_pl][i_h] = new hist1D(std::string("h_")+histnames.at(i_h)+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_pl)+histtitles.at(i_h),bins.at(i_h).at(0),bins.at(i_h).at(1),bins.at(i_h).at(2));

      mc_hists_cc1pieffpur[i_pl][i_h] = new histCC1piselEffPur(std::string("hCC1pieffpur_")+histnames.at(i_h)+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_pl)+histtitles.at(i_h),bins.at(i_h).at(0),bins.at(i_h).at(1),bins.at(i_h).at(2));
    }
  }

  // Loop through MC tree and fill plots
  for (int i = 0; i < t_bnbcos->GetEntries(); i++){
    t_bnbcos->GetEntry(i);
    CalcPIDvars(&mc_vars);
    std::vector<std::vector<std::vector<double>>> PIDvarstoplot = GetPIDvarstoplot(&mc_vars);

    //for (size_t i_pl=0; i_pl < nplanes; i_pl++){
    for (size_t i_pl=2; i_pl < 3; i_pl++){
      for (size_t i_h = 0; i_h < nplots; i_h++){
        FillHist(mc_hists[i_pl][i_h],PIDvarstoplot.at(i_h),i_pl,mc_vars.TPCObj_PFP_PandoraClassedAsTrack,mc_vars.TPCObj_PFP_truePDG);

        FillCC1piEffPurHist(mc_hists_cc1pieffpur[i_pl][i_h],PIDvarstoplot.at(i_h),i_pl,mc_vars.TPCObj_PFP_PandoraClassedAsTrack,mc_vars.Truth_topology,MIPlow.at(i_h));
      }
    }


  } // end loop over entries in tree


  // -------------------- Now make all the plots

  // for (size_t i_pl=0; i_pl < nplanes; i_pl++){
  for (size_t i_pl=2; i_pl < 3; i_pl++){
    for (size_t i_h=0; i_h < nplots; i_h++){
      TCanvas *c1 = new TCanvas("c1","c1");
      DrawMCEffPur(c1, mc_hists[i_pl][i_h],MIPlow.at(i_h));
      std::string printname = std::string(histnames[i_h]+std::string("_plane")+std::to_string(i_pl)+".png");
      c1->Print(printname.c_str());
      c1->Clear();

      DrawCC1piMCEffPur(c1, mc_hists_cc1pieffpur[i_pl][i_h]);
      c1->Print(std::string(std::string("CC1pi_")+printname).c_str());
      delete c1;
    }
  }

  // for (size_t i_pl=0; i_pl < nplanes; i_pl++){
  for (size_t i_pl=2; i_pl < 3; i_pl++){
    for (size_t i_h=0; i_h < nplots; i_h++){

    }
  }
}
