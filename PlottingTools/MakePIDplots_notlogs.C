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
                  "Lp",
                  "Lmu",
                  "Lpi",
                  "Lmip",
                  "Lmumip",
                  "chi2mu",
                  "chi2p",
                  "chi2pi",
                  "pida_median",
                  "Lmuoverp",
                  "Lmipoverp",
                  "Lmumipoverp",
                  "chi2muminusp",
                  "Lmu0to1",
                  "Lmip0to1",
                  "Lpi0to1",
                  "Lp0to1",
                  "Lmumip0to1",
                  "Lmumippi0to1",
                  "depErangeEmu",
                  "depErangeEp"
                };

// ---------------------------------------------------- //
//  Now the function starts
// ---------------------------------------------------- //

void MakePIDplots(std::string mcfile, double POTscaling=0., std::string onbeamdatafile="", std::string offbeamdatafile="", double offbeamscaling=0., bool onminusoffbeam=true){

  gStyle->SetTitleX(0.1f);
  gStyle->SetTitleW(0.8f);
  gStyle->SetTitleBorderSize(0.);

  TFile *f_bnbcos = new TFile(mcfile.c_str(), "read");
  TTree *t_bnbcos = (TTree*)f_bnbcos->Get("cc1piselec/outtree");
  treevars mc_vars;
  settreevars(t_bnbcos,&mc_vars);

  TFile *f_onbeam=nullptr;
  TTree *t_onbeam=nullptr;
  treevars onbeam_vars;
  if (onbeamdatafile!=""){
    std::cout << "Making data-MC comparisons" << std::endl;
    f_onbeam = new TFile(onbeamdatafile.c_str(), "read");
    t_onbeam = (TTree*)f_onbeam->Get("cc1piselec/outtree");
    settreevars(t_onbeam,&onbeam_vars);
  }

  TFile *f_offbeam=nullptr;
  TTree *t_offbeam=nullptr;
  treevars offbeam_vars;
  if (offbeamdatafile!=""){
    f_offbeam = new TFile(offbeamdatafile.c_str(), "read");
    t_offbeam = (TTree*)f_offbeam->Get("cc1piselec/outtree");
    settreevars(t_offbeam,&offbeam_vars);
  }


  // Sanity check: the plot vectors should be the same size
  t_bnbcos->GetEntry(0);
  CalcPIDvars(&mc_vars);
  std::vector<std::vector<std::vector<double>>> PIDvarstoplot_dummy = GetPIDvarstoplot(&mc_vars);

  std::cout << "PIDvarstoplot.size() = " << PIDvarstoplot_dummy.size() << std::endl;
  std::cout << "bins.size() = " << bins.size() << std::endl;
  std::cout << "histtitles.size() = " << histtitles.size() << std::endl;
  std::cout << "histnames.size() = " << histnames.size() << std::endl;

  // ----------------- MC

  // Make histograms to fill
  const size_t nplots = PIDvarstoplot_dummy.size();
  hist1D *mc_hists[nplanes][nplots];
  for (int i_pl=0; i_pl<nplanes; i_pl++){
    for (int i_h=0; i_h<nplots; i_h++){
      mc_hists[i_pl][i_h] = new hist1D(std::string("h_")+histnames.at(i_h)+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_pl)+histtitles.at(i_h),bins.at(i_h).at(0),bins.at(i_h).at(1),bins.at(i_h).at(2));
    }
  }


  // Loop through MC tree and fill plots
  for (int i = 0; i < t_bnbcos->GetEntries(); i++){
    t_bnbcos->GetEntry(i);
    CalcPIDvars(&mc_vars);
    std::vector<std::vector<std::vector<double>>> PIDvarstoplot = GetPIDvarstoplot(&mc_vars);

    for (size_t i_pl=0; i_pl < nplanes; i_pl++){
      for (size_t i_h = 0; i_h < nplots; i_h++){
        FillHist(mc_hists[i_pl][i_h],PIDvarstoplot.at(i_h),i_pl,mc_vars.TPCObj_PFP_PandoraClassedAsTrack,mc_vars.TPCObj_PFP_truePDG);
      }
    }


  } // end loop over entries in tree

  // ----------------- On-beam data
  hist1D *onb_hists[nplanes][nplots];
  if (t_onbeam){
    // Make histograms to fill
    for (size_t i_pl=0; i_pl < nplanes; i_pl++){
      for (int i_h=0; i_h<nplots; i_h++){
        onb_hists[i_pl][i_h] = new hist1D(std::string("h_ondat_")+histnames.at(i_h)+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_pl)+histtitles.at(i_h),bins.at(i_h).at(0),bins.at(i_h).at(1),bins.at(i_h).at(2));
      }
    }


    // Loop through on-beam data tree and fill plots
    for (int i = 0; i < t_onbeam->GetEntries(); i++){
      t_onbeam->GetEntry(i);
      CalcPIDvars(&onbeam_vars);
      std::vector<std::vector<std::vector<double>>> PIDvarstoplot = GetPIDvarstoplot(&onbeam_vars);


      for (size_t i_pl=0; i_pl < nplanes; i_pl++){
        for (size_t i_h = 0; i_h < nplots; i_h++){
          FillHist(onb_hists[i_pl][i_h],PIDvarstoplot.at(i_h),i_pl,onbeam_vars.TPCObj_PFP_PandoraClassedAsTrack); // 0 because there is no "true PDG" for data
        }
      }

    }
  }

    // ----------------- Off-beam data
  hist1D *offb_hists[nplanes][nplots];
  if (t_offbeam){
    // Make histograms to fill
    for (size_t i_pl=0; i_pl < nplanes; i_pl++){
      for (int i_h=0; i_h<nplots; i_h++){
        offb_hists[i_pl][i_h] = new hist1D(std::string("h_offdat_")+histnames.at(i_h)+std::string("_plane")+std::to_string(i_pl),std::string("Plane ")+std::to_string(i_pl)+histtitles.at(i_h),bins.at(i_h).at(0),bins.at(i_h).at(1),bins.at(i_h).at(2));
      }
    }


    // Loop through tree and fill plots
    for (int i = 0; i < t_offbeam->GetEntries(); i++){
      t_offbeam->GetEntry(i);
      CalcPIDvars(&offbeam_vars);
      std::vector<std::vector<std::vector<double>>> PIDvarstoplot = GetPIDvarstoplot(&offbeam_vars);

      for (size_t i_pl=0; i_pl < nplanes; i_pl++){
        for (size_t i_h = 0; i_h < nplots; i_h++){
          FillHist(offb_hists[i_pl][i_h],PIDvarstoplot.at(i_h),i_pl,offbeam_vars.TPCObj_PFP_PandoraClassedAsTrack); // 0 because there is no "true PDG" for data
        }
      }
    } // end loop over entries in tree
  }


  // -------------------- Now make all the plots

  for (size_t i_pl=0; i_pl < nplanes; i_pl++){
    for (size_t i_h=0; i_h < nplots; i_h++){
      TCanvas *c1 = new TCanvas();
      if (onminusoffbeam){
        DrawMC(mc_hists[i_pl][i_h],POTscaling);
        if (f_onbeam && f_offbeam){
          OverlayOnMinusOffData(c1,onb_hists[i_pl][i_h],offb_hists[i_pl][i_h],offbeamscaling,POTscaling);
        }
      }
      else{
        if (f_onbeam && f_offbeam){
          DrawMCPlusOffbeam(mc_hists[i_pl][i_h], offb_hists[i_pl][i_h], POTscaling, offbeamscaling);
          OverlayOnBeamData(c1, onb_hists[i_pl][i_h]);
        }
        else{
          DrawMC(mc_hists[i_pl][i_h],POTscaling);
        }
      }
      c1->Print(std::string(histnames[i_h]+std::string("_plane")+std::to_string(i_pl)+".png").c_str());
    }
  }
}
