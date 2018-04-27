void MakePIDPlots(std::string filename)
{
  // Get input file and TTree
  TFile *infile = new TFile(filename.c_str(),"open");
  TTree *tree   = (TTree*)infile->Get("cc1piselec/outtree");

  // Make output file
  TFile *outfile = new TFile("CC1piselec_PIDplots_out.root","recreate");

  // Make histograms
  TH1F *hdEdx_all = new TH1F("hdEdx_all",";dE/dx;",100,0,10);
  TH1F *hdEdx_Muon = new TH1F("hdEdx_Muon",";dE/dx;",100,0,10);
  TH1F *hdEdx_Proton = new TH1F("hdEdx_Proton",";dE/dx;",100,0,10);

  TH2F *hdEdx_rr_all = new TH2F("hdEdx_rr_all",";Residual Range (cm);dE/dx",1000,0,200,1000,0,50);
  TH2F *hdEdx_rr_Muon = new TH2F("hdEdx_rr_Muon",";Residual Range (cm);dE/dx",1000,0,200,1000,0,50);
  TH2F *hdEdx_rr_Proton = new TH2F("hdEdx_rr_Proton",";Residual Range (cm);dE/dx",1000,0,200,1000,0,50);

  TH1F *hnegLLmu_all = new TH1F("neg2llhmu_all",";-2ln(L_{#mu});",60,0,30);
  TH1F *hnegLLmu_Muon = new TH1F("neg2llhmu_Muon","",60,0,30);
  TH1F *hnegLLmu_Proton = new TH1F("neg2llhmu_Proton","",60,0,30);
  TH1F *hnegLLmu_Pion = new TH1F("neg2llhmu_Pion","",60,0,30);

  TH1F *hnegLLMIP_all = new TH1F("neg2llhMIP_all",";-2ln(L_{MIP});",60,0,30);
  TH1F *hnegLLMIP_Muon = new TH1F("neg2llhMIP_Muon","",60,0,30);
  TH1F *hnegLLMIP_Proton = new TH1F("neg2llhMIP_Proton","",60,0,30);
  TH1F *hnegLLMIP_Pion = new TH1F("neg2llhMIP_Pion","",60,0,30);

  TH1F *hnegLLmuMIP_all = new TH1F("neg2llhmuMIP_all",";-2ln(L_{#mu/MIP});",60,0,30);
  TH1F *hnegLLmuMIP_Muon = new TH1F("neg2llhmuMIP_Muon","",60,0,30);
  TH1F *hnegLLmuMIP_Proton = new TH1F("neg2llhmuMIP_Proton","",60,0,30);
  TH1F *hnegLLmuMIP_Pion = new TH1F("neg2llhmuMIP_Pion","",60,0,30);

  TH1F *hnegLLp_all = new TH1F("neg2llhp_all",";-2ln(L_{p});",60,0,30);
  TH1F *hnegLLp_Muon = new TH1F("neg2llhp_Muon","",60,0,30);
  TH1F *hnegLLp_Proton = new TH1F("neg2llhp_Proton","",60,0,30);
  TH1F *hnegLLp_Pion = new TH1F("neg2llhp_Pion","",60,0,30);

  TH1F *hnegLLmuminusp_all = new TH1F("neg2llhmuminusp_all",";(-2ln(L_{#mu}))-(-2ln(L_{p}));",40,-20,7);
  TH1F *hnegLLmuminusp_Muon = new TH1F("neg2llhmuminusp_Muon","",40,-20,7);
  TH1F *hnegLLmuminusp_Proton = new TH1F("neg2llhmuminusp_Proton","",40,-20,7);
  TH1F *hnegLLmuminusp_Pion = new TH1F("neg2llhmuminusp_Pion","",40,-20,7);

  TH1F *hnegLLmuMIPminusp_all = new TH1F("neg2llhmuMIPminusp_all",";(-2ln(L_{#mu/MIP}))-(-2ln(L_{p}));",50,-20,7);
  TH1F *hnegLLmuMIPminusp_Muon = new TH1F("neg2llhmuMIPminusp_Muon","",50,-20,7);
  TH1F *hnegLLmuMIPminusp_Proton = new TH1F("neg2llhmuMIPminusp_Proton","",50,-20,7);
  TH1F *hnegLLmuMIPminusp_Pion = new TH1F("neg2llhmuMIPminusp_Pion","",50,-20,7);

  TH2F *hdepErangeEmu_all = new TH2F("hdepErangeEmu_all",";Deposited Energy;Energy by range (muon assumption)",100,0,800,100,0,800);
  TH2F *hdepErangeEmu_Muon = new TH2F("hdepErangeEmu_Muon",";Deposited Energy;Energy by range (muon assumption)",100,0,800,100,0,800);
  TH2F *hdepErangeEmu_Proton = new TH2F("hdepErangeEmu_Proton",";Deposited Energy;Energy by range (muon assumption)",100,0,800,100,0,800);
  TH2F *hdepErangeEmu_Pion = new TH2F("hdepErangeEmu_Pion",";Deposited Energy;Energy by range (muon assumption)",100,0,800,100,0,800);

  TH2F *hdepErangeEp_all = new TH2F("hdepErangeEp_all",";Deposited Energy;Energy by range (proton assumption)",100,0,800,100,0,800);
  TH2F *hdepErangeEp_Muon = new TH2F("hdepErangeEp_Muon",";Deposited Energy;Energy by range (proton assumption)",100,0,800,100,0,800);
  TH2F *hdepErangeEp_Proton = new TH2F("hdepErangeEp_Proton",";Deposited Energy;Energy by range (proton assumption)",100,0,800,100,0,800);
  TH2F *hdepErangeEp_Pion = new TH2F("hdepErangeEp_Pion",";Deposited Energy;Energy by range (proton assumption)",100,0,800,100,0,800);

  TH1F *hdepEminusrangeEmu_all = new TH1F("hdepEminusrangeEmu_all",";Deposited energy - energy by range (muon assumption)",100,-800,800);
  TH1F *hdepEminusrangeEmu_Muon = new TH1F("hdepEminusrangeEmu_Muon",";Deposited energy - energy by range (muon assumption)",100,-800,800);
  TH1F *hdepEminusrangeEmu_Proton = new TH1F("hdepEminusrangeEmu_Proton",";Deposited energy - energy by range (muon assumption)",100,-800,800);
  TH1F *hdepEminusrangeEmu_Pion = new TH1F("hdepEminusrangeEmu_Pion",";Deposited energy - energy by range (muon assumption)",100,-800,800);

  TH1F *hdepEminusrangeEp_all = new TH1F("hdepEminusrangeEp_all",";Deposited energy - energy by range (proton assumption)",100,-800,800);
  TH1F *hdepEminusrangeEp_Muon = new TH1F("hdepEminusrangeEp_Muon",";Deposited energy - energy by range (proton assumption)",100,-800,800);
  TH1F *hdepEminusrangeEp_Proton = new TH1F("hdepEminusrangeEp_Proton",";Deposited energy - energy by range (proton assumption)",100,-800,800);
  TH1F *hdepEminusrangeEp_Pion = new TH1F("hdepEminusrangeEp_Pion",";Deposited energy - energy by range (proton assumption)",100,-800,800);

  TH2F *hdepEmrangeEmuvsp_all = new TH2F("hdepEmrangeEmuvsp_all",";Deposited energy - energy by range (muon assumption);Deposited energy - energy by range (proton assumption)",100,-800,800,100,-800,800);
  TH2F *hdepEmrangeEmuvsp_Muon = new TH2F("hdepEmrangeEmuvsp_Muon",";Deposited energy - energy by range (muon assumption);Deposited energy - energy by range (proton assumption)",100,-800,800,100,-800,800);
  TH2F *hdepEmrangeEmuvsp_Proton = new TH2F("hdepEmrangeEmuvsp_Proton",";Deposited energy - energy by range (muon assumption);Deposited energy - energy by range (proton assumption)",100,-800,800,100,-800,800);
  TH2F *hdepEmrangeEmuvsp_Pion = new TH2F("hdepEmrangeEmuvsp_Pion",";Deposited energy - energy by range (muon assumption);Deposited energy - energy by range (proton assumption)",100,-800,800,100,-800,800);

  TH1F *hdepEminusrangeEp_cut_all = new TH1F("hdepEminusrangeEp_cut_all","Events with neg2LLmuMIPminusp < -2;Deposited energy - energy by range (proton assumption)",100,-800,800);
  TH1F *hdepEminusrangeEp_cut_Muon = new TH1F("hdepEminusrangeEp_cut_Muon","Events with neg2LLmuMIPminusp < -2;Deposited energy - energy by range (proton assumption)",100,-800,800);
  TH1F *hdepEminusrangeEp_cut_Proton = new TH1F("hdepEminusrangeEp_cut_Proton","Events with neg2LLmuMIPminusp < -2;Deposited energy - energy by range (proton assumption)",100,-800,800);
  TH1F *hdepEminusrangeEp_cut_Pion = new TH1F("hdepEminusrangeEp_cut_Pion","Events with neg2LLmuMIPminusp < -2;Deposited energy - energy by range (proton assumption)",100,-800,800);

  // Set histogram styles
  hnegLLmu_all->SetLineColor(kBlack);
  hnegLLmu_all->SetFillColor(TColor::GetColor(197,197,197));
  hnegLLMIP_all->SetLineColor(kBlack);
  hnegLLMIP_all->SetFillColor(TColor::GetColor(197,197,197));
  hnegLLmuMIP_all->SetLineColor(kBlack);
  hnegLLmuMIP_all->SetFillColor(TColor::GetColor(197,197,197));
  hnegLLp_all->SetLineColor(kBlack);
  hnegLLp_all->SetFillColor(TColor::GetColor(197,197,197));
  hnegLLmuminusp_all->SetLineColor(kBlack);
  hnegLLmuminusp_all->SetFillColor(TColor::GetColor(197,197,197));
  hnegLLmuMIPminusp_all->SetLineColor(kBlack);
  hnegLLmuMIPminusp_all->SetFillColor(TColor::GetColor(197,197,197));
  hdepErangeEmu_all->SetMarkerColor(TColor::GetColor(197,197,197));
  hdepErangeEp_all->SetMarkerColor(TColor::GetColor(197,197,197));
  hdepEmrangeEmuvsp_all->SetMarkerColor(TColor::GetColor(197,197,197));
  hdepEminusrangeEmu_all->SetLineColor(TColor::GetColor(197,197,197));
  hdepEminusrangeEmu_all->SetFillColor(TColor::GetColor(197,197,197));
  hdepEminusrangeEp_all->SetLineColor(TColor::GetColor(197,197,197));
  hdepEminusrangeEp_all->SetFillColor(TColor::GetColor(197,197,197));
  hdepEminusrangeEp_cut_all->SetLineColor(TColor::GetColor(197,197,197));
  hdepEminusrangeEp_cut_all->SetFillColor(TColor::GetColor(197,197,197));

  hnegLLmu_Muon->SetLineColor(TColor::GetColor(8,64,129));
  hnegLLmu_Muon->SetFillColor(TColor::GetColor(8,64,129));
  hnegLLMIP_Muon->SetLineColor(TColor::GetColor(8,64,129));
  hnegLLMIP_Muon->SetFillColor(TColor::GetColor(8,64,129));
  hnegLLmuMIP_Muon->SetLineColor(TColor::GetColor(8,64,129));
  hnegLLmuMIP_Muon->SetFillColor(TColor::GetColor(8,64,129));
  hnegLLp_Muon->SetLineColor(TColor::GetColor(8,64,129));
  hnegLLp_Muon->SetFillColor(TColor::GetColor(8,64,129));
  hnegLLmuminusp_Muon->SetLineColor(TColor::GetColor(8,64,129));
  hnegLLmuminusp_Muon->SetFillColor(TColor::GetColor(8,64,129));
  hnegLLmuMIPminusp_Muon->SetLineColor(TColor::GetColor(8,64,129));
  hnegLLmuMIPminusp_Muon->SetFillColor(TColor::GetColor(8,64,129));
  hdepErangeEmu_Muon->SetMarkerColor(TColor::GetColor(8,64,129));
  hdepErangeEp_Muon->SetMarkerColor(TColor::GetColor(8,64,129));
  hdepEmrangeEmuvsp_Muon->SetMarkerColor(TColor::GetColor(8,64,129));
  hdepEminusrangeEmu_Muon->SetLineColor(TColor::GetColor(8,64,129));
  hdepEminusrangeEmu_Muon->SetFillColor(TColor::GetColor(8,64,129));
  hdepEminusrangeEp_Muon->SetLineColor(TColor::GetColor(8,64,129));
  hdepEminusrangeEp_Muon->SetFillColor(TColor::GetColor(8,64,129));
  hdepEminusrangeEp_cut_Muon->SetLineColor(TColor::GetColor(8,64,129));
  hdepEminusrangeEp_cut_Muon->SetFillColor(TColor::GetColor(8,64,129));

  hnegLLmu_Proton->SetLineColor(TColor::GetColor(215, 48, 39));
  hnegLLmu_Proton->SetFillColor(TColor::GetColor(215, 48, 39));
  hnegLLMIP_Proton->SetLineColor(TColor::GetColor(215, 48, 39));
  hnegLLMIP_Proton->SetFillColor(TColor::GetColor(215, 48, 39));
  hnegLLmuMIP_Proton->SetLineColor(TColor::GetColor(215, 48, 39));
  hnegLLmuMIP_Proton->SetFillColor(TColor::GetColor(215, 48, 39));
  hnegLLp_Proton->SetLineColor(TColor::GetColor(215, 48, 39));
  hnegLLp_Proton->SetFillColor(TColor::GetColor(215, 48, 39));
  hnegLLmuminusp_Proton->SetLineColor(TColor::GetColor(215, 48, 39));
  hnegLLmuminusp_Proton->SetFillColor(TColor::GetColor(215, 48, 39));
  hnegLLmuMIPminusp_Proton->SetLineColor(TColor::GetColor(215, 48, 39));
  hnegLLmuMIPminusp_Proton->SetFillColor(TColor::GetColor(215, 48, 39));
  hdepErangeEmu_Proton->SetMarkerColor(TColor::GetColor(215, 48, 39));
  hdepErangeEp_Proton->SetMarkerColor(TColor::GetColor(215, 48, 39));
  hdepEmrangeEmuvsp_Proton->SetMarkerColor(TColor::GetColor(215, 48, 39));
  hdepEminusrangeEmu_Proton->SetLineColor(TColor::GetColor(215, 48, 39));
  hdepEminusrangeEmu_Proton->SetFillColor(TColor::GetColor(215, 48, 39));
  hdepEminusrangeEp_Proton->SetLineColor(TColor::GetColor(215, 48, 39));
  hdepEminusrangeEp_Proton->SetFillColor(TColor::GetColor(215, 48, 39));
  hdepEminusrangeEp_cut_Proton->SetLineColor(TColor::GetColor(215, 48, 39));
  hdepEminusrangeEp_cut_Proton->SetFillColor(TColor::GetColor(215, 48, 39));

  hnegLLmu_Pion->SetLineColor(TColor::GetColor(166,217,106));
  hnegLLmu_Pion->SetFillColor(TColor::GetColor(166,217,106));
  hnegLLMIP_Pion->SetLineColor(TColor::GetColor(166,217,106));
  hnegLLMIP_Pion->SetFillColor(TColor::GetColor(166,217,106));
  hnegLLmuMIP_Pion->SetLineColor(TColor::GetColor(166,217,106));
  hnegLLmuMIP_Pion->SetFillColor(TColor::GetColor(166,217,106));
  hnegLLp_Pion->SetLineColor(TColor::GetColor(166,217,106));
  hnegLLp_Pion->SetFillColor(TColor::GetColor(166,217,106));
  hnegLLmuminusp_Pion->SetLineColor(TColor::GetColor(166,217,106));
  hnegLLmuminusp_Pion->SetFillColor(TColor::GetColor(166,217,106));
  hnegLLmuMIPminusp_Pion->SetLineColor(TColor::GetColor(166,217,106));
  hnegLLmuMIPminusp_Pion->SetFillColor(TColor::GetColor(166,217,106));
  hdepErangeEmu_Pion->SetMarkerColor(TColor::GetColor(166,217,106));
  hdepErangeEp_Pion->SetMarkerColor(TColor::GetColor(166,217,106));
  hdepEmrangeEmuvsp_Pion->SetMarkerColor(TColor::GetColor(166,217,106));
  hdepEminusrangeEmu_Pion->SetLineColor(TColor::GetColor(166,217,106));
  hdepEminusrangeEmu_Pion->SetFillColor(TColor::GetColor(166,217,106));
  hdepEminusrangeEp_Pion->SetLineColor(TColor::GetColor(166,217,106));
  hdepEminusrangeEp_Pion->SetFillColor(TColor::GetColor(166,217,106));
  hdepEminusrangeEp_cut_Pion->SetLineColor(TColor::GetColor(166,217,106));
  hdepEminusrangeEp_cut_Pion->SetFillColor(TColor::GetColor(166,217,106));

  // Variables to read from tree
  std::vector<std::vector<double>> *dEdxvec = nullptr;
  std::vector<std::vector<double>> *resrangevec = nullptr;
  std::vector<double> *n2LLH_fwd_mu = nullptr;
  std::vector<double> *n2LLH_bwd_mu = nullptr;
  std::vector<double> *n2LLH_fwd_p = nullptr;
  std::vector<double> *n2LLH_bwd_p = nullptr;
  std::vector<double> *n2LLH_MIP = nullptr;
  std::vector<double> *depE = nullptr;
  std::vector<double> *rangeE_mu = nullptr;
  std::vector<double> *rangeE_p = nullptr;
  std::vector<double> *truePDG = nullptr;
  std::vector<bool> *isTrack = nullptr;

  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("TPCObj_PFP_track_dedx_perhit",1);
  tree->SetBranchAddress("TPCObj_PFP_track_dedx_perhit",&dEdxvec);
  tree->SetBranchStatus("TPCObj_PFP_track_resrange_perhit",1);
  tree->SetBranchAddress("TPCObj_PFP_track_resrange_perhit",&resrangevec);
  tree->SetBranchStatus("TPCObj_PFP_n2LLH_fwd_mu",1);
  tree->SetBranchAddress("TPCObj_PFP_n2LLH_fwd_mu",&n2LLH_fwd_mu);
  tree->SetBranchStatus("TPCObj_PFP_n2LLH_bwd_mu",1);
  tree->SetBranchAddress("TPCObj_PFP_n2LLH_bwd_mu",&n2LLH_bwd_mu);
  tree->SetBranchStatus("TPCObj_PFP_n2LLH_fwd_p",1);
  tree->SetBranchAddress("TPCObj_PFP_n2LLH_fwd_p",&n2LLH_fwd_p);
  tree->SetBranchStatus("TPCObj_PFP_n2LLH_bwd_p",1);
  tree->SetBranchAddress("TPCObj_PFP_n2LLH_bwd_p",&n2LLH_bwd_p);
  tree->SetBranchStatus("TPCObj_PFP_n2LLH_MIP",1);
  tree->SetBranchAddress("TPCObj_PFP_n2LLH_MIP",&n2LLH_MIP);
  tree->SetBranchStatus("TPCObj_PFP_track_depE",1);
  tree->SetBranchAddress("TPCObj_PFP_track_depE",&depE);
  tree->SetBranchStatus("TPCObj_PFP_track_rangeE_mu",1);
  tree->SetBranchAddress("TPCObj_PFP_track_rangeE_mu",&rangeE_mu);
  tree->SetBranchStatus("TPCObj_PFP_track_rangeE_p",1);
  tree->SetBranchAddress("TPCObj_PFP_track_rangeE_p",&rangeE_p);
  tree->SetBranchStatus("TPCObj_PFP_truePDG",1);
  tree->SetBranchAddress("TPCObj_PFP_truePDG",&truePDG);
  tree->SetBranchStatus("TPCObj_PFP_isTrack",1);
  tree->SetBranchAddress("TPCObj_PFP_isTrack",&isTrack);


  // Loop through tree and plot things
  for (size_t i_tree=0; i_tree<tree->GetEntries(); i_tree++){
    tree->GetEntry(i_tree);

    for (size_t i_pfp=0; i_pfp<truePDG->size(); i_pfp++){

      // PID variables are only filled for tracks. If one isn't filled, probably none of them are, so just skip the entire track.
      if (n2LLH_fwd_mu->at(i_pfp) == -999) continue;
      //if (!rangeE_mu) continue;

      // Skip PFPs that Pandora doesn't reconstruct as tracks
      // This is probably not what we want to do for the CC1pi analysis, but is useful for making plots for general consumption/PID presentations
      if (!isTrack->at(i_pfp)) continue;

      // Make dE/dx plots
      std::vector<double> dEdx = dEdxvec->at(i_pfp);
      std::vector<double> resrange = resrangevec->at(i_pfp);

      for (size_t i_dEdx=0; i_dEdx<dEdx.size(); i_dEdx++){
        if (resrange.at(i_dEdx) > 100 && resrange.at(i_dEdx) < 150){
          hdEdx_all->Fill(dEdx.at(i_dEdx));
          hdEdx_rr_all->Fill(resrange.at(i_dEdx),dEdx.at(i_dEdx));
        }
        if (resrange.at(i_dEdx) > 100 && resrange.at(i_dEdx) < 150 && TMath::Abs(truePDG->at(i_pfp))==13){
          hdEdx_Muon->Fill(dEdx.at(i_dEdx));
          hdEdx_rr_Muon->Fill(resrange.at(i_dEdx),dEdx.at(i_dEdx));
        }
        if (resrange.at(i_dEdx) > 100 && resrange.at(i_dEdx) < 150 && TMath::Abs(truePDG->at(i_pfp))==2212){
          hdEdx_Proton->Fill(dEdx.at(i_dEdx));
          hdEdx_rr_Proton->Fill(resrange.at(i_dEdx),dEdx.at(i_dEdx));
        }
      } // loop over hits (dEdx and resrange)


      // Fill n2llh histograms
      double tr_n2llh_mu = std::min(n2LLH_fwd_mu->at(i_pfp),n2LLH_bwd_mu->at(i_pfp));
      double tr_n2llh_MIP = n2LLH_MIP->at(i_pfp);
      double tr_n2llh_muMIP = std::min({n2LLH_fwd_mu->at(i_pfp),n2LLH_bwd_mu->at(i_pfp),n2LLH_MIP->at(i_pfp)});
      double tr_n2llh_p = std::min(n2LLH_fwd_p->at(i_pfp),n2LLH_bwd_p->at(i_pfp));
      double tr_n2llh_muminusp = tr_n2llh_mu - tr_n2llh_p;
      double tr_n2llh_muMIPminusp = tr_n2llh_muMIP - tr_n2llh_p;

      hnegLLmu_all->Fill(tr_n2llh_mu);
      hnegLLMIP_all->Fill(tr_n2llh_MIP);
      hnegLLmuMIP_all->Fill(tr_n2llh_muMIP);
      hnegLLp_all->Fill(tr_n2llh_p);
      hnegLLmuminusp_all->Fill(tr_n2llh_muminusp);
      hnegLLmuMIPminusp_all->Fill(tr_n2llh_muMIPminusp);
      hdepErangeEmu_all->Fill(depE->at(i_pfp),rangeE_mu->at(i_pfp));
      hdepErangeEp_all->Fill(depE->at(i_pfp),rangeE_p->at(i_pfp));
      hdepEminusrangeEmu_all->Fill(depE->at(i_pfp)-rangeE_mu->at(i_pfp));
      hdepEminusrangeEp_all->Fill(depE->at(i_pfp)-rangeE_p->at(i_pfp));
      hdepEmrangeEmuvsp_all->Fill(depE->at(i_pfp)-rangeE_mu->at(i_pfp),depE->at(i_pfp)-rangeE_p->at(i_pfp));

      if (tr_n2llh_muMIPminusp < -2){
        hdepEminusrangeEp_cut_all->Fill(depE->at(i_pfp)-rangeE_p->at(i_pfp));
      }

      if (TMath::Abs(truePDG->at(i_pfp))==13){
        hnegLLmu_Muon->Fill(tr_n2llh_mu);
        hnegLLMIP_Muon->Fill(tr_n2llh_MIP);
        hnegLLmuMIP_Muon->Fill(tr_n2llh_muMIP);
        hnegLLp_Muon->Fill(tr_n2llh_p);
        hnegLLmuminusp_Muon->Fill(tr_n2llh_muminusp);
        hnegLLmuMIPminusp_Muon->Fill(tr_n2llh_muMIPminusp);
        hdepErangeEmu_Muon->Fill(depE->at(i_pfp),rangeE_mu->at(i_pfp));
        hdepErangeEp_Muon->Fill(depE->at(i_pfp),rangeE_p->at(i_pfp));
        hdepEminusrangeEmu_Muon->Fill(depE->at(i_pfp)-rangeE_mu->at(i_pfp));
        hdepEminusrangeEp_Muon->Fill(depE->at(i_pfp)-rangeE_p->at(i_pfp));
        hdepEmrangeEmuvsp_Muon->Fill(depE->at(i_pfp)-rangeE_mu->at(i_pfp),depE->at(i_pfp)-rangeE_p->at(i_pfp));

        if (tr_n2llh_muMIPminusp < -2){
          hdepEminusrangeEp_cut_Muon->Fill(depE->at(i_pfp)-rangeE_p->at(i_pfp));
        }
      }
      else if (TMath::Abs(truePDG->at(i_pfp))==2212){
        hnegLLmu_Proton->Fill(tr_n2llh_mu);
        hnegLLMIP_Proton->Fill(tr_n2llh_MIP);
        hnegLLmuMIP_Proton->Fill(tr_n2llh_muMIP);
        hnegLLp_Proton->Fill(tr_n2llh_p);
        hnegLLmuminusp_Proton->Fill(tr_n2llh_muminusp);
        hnegLLmuMIPminusp_Proton->Fill(tr_n2llh_muMIPminusp);
        hdepErangeEmu_Proton->Fill(depE->at(i_pfp),rangeE_mu->at(i_pfp));
        hdepErangeEp_Proton->Fill(depE->at(i_pfp),rangeE_p->at(i_pfp));
        hdepEminusrangeEmu_Proton->Fill(depE->at(i_pfp)-rangeE_mu->at(i_pfp));
        hdepEminusrangeEp_Proton->Fill(depE->at(i_pfp)-rangeE_p->at(i_pfp));
        hdepEmrangeEmuvsp_Proton->Fill(depE->at(i_pfp)-rangeE_mu->at(i_pfp),depE->at(i_pfp)-rangeE_p->at(i_pfp));

        if (tr_n2llh_muMIPminusp < -2){
          hdepEminusrangeEp_cut_Proton->Fill(depE->at(i_pfp)-rangeE_p->at(i_pfp));
        }
      }
      else if (TMath::Abs(truePDG->at(i_pfp))==211){
        hnegLLmu_Pion->Fill(tr_n2llh_mu);
        hnegLLMIP_Pion->Fill(tr_n2llh_MIP);
        hnegLLmuMIP_Pion->Fill(tr_n2llh_muMIP);
        hnegLLp_Pion->Fill(tr_n2llh_p);
        hnegLLmuminusp_Pion->Fill(tr_n2llh_muminusp);
        hnegLLmuMIPminusp_Pion->Fill(tr_n2llh_muMIPminusp);
        hdepErangeEmu_Pion->Fill(depE->at(i_pfp),rangeE_mu->at(i_pfp));
        hdepErangeEp_Pion->Fill(depE->at(i_pfp),rangeE_p->at(i_pfp));
        hdepEminusrangeEmu_Pion->Fill(depE->at(i_pfp)-rangeE_mu->at(i_pfp));
        hdepEminusrangeEp_Pion->Fill(depE->at(i_pfp)-rangeE_p->at(i_pfp));
        hdepEmrangeEmuvsp_Pion->Fill(depE->at(i_pfp)-rangeE_mu->at(i_pfp),depE->at(i_pfp)-rangeE_p->at(i_pfp));

        if (tr_n2llh_muMIPminusp < -2){
          hdepEminusrangeEp_cut_Pion->Fill(depE->at(i_pfp)-rangeE_p->at(i_pfp));
        }
      }

    } // loop over PFPs

  } // loop over entries in the tree

  TLegend *l = new TLegend(0.65,0.7,0.8,0.85);
  l->SetFillColor(kWhite);
  l->SetLineColor(kWhite);
  l->SetTextFont(132);
  l->AddEntry(hnegLLmu_Muon,"True muons","f");
  l->AddEntry(hnegLLmu_Proton,"True protons","f");
  l->AddEntry(hnegLLmu_Pion,"True pions","f");
  l->AddEntry(hnegLLmu_all,"MC other","f");

  // Plot things!
  outfile->cd();
  hdEdx_all->Write("hdEdx_all");
  if (hdEdx_Muon->Integral() > 0) hdEdx_Muon->Write("hdEdx_Muon");
  if (hdEdx_Proton->Integral() > 0) hdEdx_Proton->Write("hdEdx_Proton");
  hdEdx_rr_all->Write("hdEdx_rr_all");
  if (hdEdx_rr_Muon->Integral() > 0) hdEdx_rr_Muon->Write("hdEdx_rr_Muon");
  if (hdEdx_rr_Proton->Integral() > 0) hdEdx_rr_Proton->Write("hdEdx_rr_Proton");

  TCanvas *c1 = new TCanvas();

  hnegLLmu_all->Write();
  hnegLLmu_Muon->Write();
  hnegLLmu_Proton->Write();
  hnegLLmu_Pion->Write();

  hnegLLmu_all->Draw();
  THStack *st_negLLmu = new THStack("st_negLLmu","");
  st_negLLmu->Add(hnegLLmu_Proton);
  st_negLLmu->Add(hnegLLmu_Muon);
  st_negLLmu->Add(hnegLLmu_Pion);
  st_negLLmu->Draw("same");
  l->Draw();
  c1->Print("neg2LLmu.pdf");
  c1->Write("neg2LLmucanv");

  hnegLLMIP_all->Write();
  hnegLLMIP_Muon->Write();
  hnegLLMIP_Proton->Write();
  hnegLLMIP_Pion->Write();

  hnegLLMIP_all->Draw();
  THStack *st_negLLMIP = new THStack("st_negLLMIP","");
  st_negLLMIP->Add(hnegLLMIP_Proton);
  st_negLLMIP->Add(hnegLLMIP_Muon);
  st_negLLMIP->Add(hnegLLMIP_Pion);
  st_negLLMIP->Draw("same");
  l->Draw();
  c1->Print("neg2LLMIP.pdf");
  c1->Write("neg2LLMIPcanv");

  hnegLLmuMIP_all->Write();
  hnegLLmuMIP_Muon->Write();
  hnegLLmuMIP_Proton->Write();
  hnegLLmuMIP_Pion->Write();


  hnegLLmuMIP_all->Draw();
  THStack *st_negLLmuMIP = new THStack("st_negLLmuMIP","");
  st_negLLmuMIP->Add(hnegLLmuMIP_Proton);
  st_negLLmuMIP->Add(hnegLLmuMIP_Muon);
  st_negLLmuMIP->Add(hnegLLmuMIP_Pion);
  st_negLLmuMIP->Draw("same");
  l->Draw();
  c1->Print("neg2LLmuMIP.pdf");
  c1->Write("neg2LLmuMIPcanv");

  hnegLLp_all->Write();
  hnegLLp_Muon->Write();
  hnegLLp_Proton->Write();
  hnegLLp_Pion->Write();

  hnegLLp_all->Draw();
  THStack *st_negLLp = new THStack("st_negLLp","");
  st_negLLp->Add(hnegLLp_Proton);
  st_negLLp->Add(hnegLLp_Muon);
  st_negLLp->Add(hnegLLp_Pion);
  st_negLLp->Draw("same");
  l->Draw();
  c1->Print("neg2LLp.pdf");
  c1->Write("neg2LLpcanv");

  hnegLLmuminusp_all->Write();
  hnegLLmuminusp_Muon->Write();
  hnegLLmuminusp_Proton->Write();
  hnegLLmuminusp_Pion->Write();

  hnegLLmuminusp_all->Draw();
  THStack *st_negLLmuminusp = new THStack("st_negLLmuminusp","");
  st_negLLmuminusp->Add(hnegLLmuminusp_Proton);
  st_negLLmuminusp->Add(hnegLLmuminusp_Muon);
  st_negLLmuminusp->Add(hnegLLmuminusp_Pion);
  st_negLLmuminusp->Draw("same");
  l->Draw();
  c1->Print("neg2LLmuminusp.pdf");
  c1->Write("neg2LLmuminuspcanv");

  hnegLLmuMIPminusp_all->Write();
  hnegLLmuMIPminusp_Muon->Write();
  hnegLLmuMIPminusp_Proton->Write();
  hnegLLmuMIPminusp_Pion->Write();

  hnegLLmuMIPminusp_all->Draw();
  THStack *st_negLLmuMIPminusp = new THStack("st_negLLmuMIPminusp","");
  st_negLLmuMIPminusp->Add(hnegLLmuMIPminusp_Proton);
  st_negLLmuMIPminusp->Add(hnegLLmuMIPminusp_Muon);
  st_negLLmuMIPminusp->Add(hnegLLmuMIPminusp_Pion);
  st_negLLmuMIPminusp->Draw("same");
  l->Draw();
  c1->Print("neg2LLmuMIPminusp.pdf");
  c1->Write("neg2LLmuMIPminuspcanv");

  hdepErangeEmu_all->Write();
  hdepErangeEmu_Muon->Write();
  hdepErangeEmu_Proton->Write();
  hdepErangeEmu_Pion->Write();

  hdepErangeEmu_all->Draw();
  hdepErangeEmu_Muon->Draw("same");
  hdepErangeEmu_Pion->Draw("same");
  hdepErangeEmu_Proton->Draw("same");
  l->Draw();
  c1->Print("depErangeEmu.pdf");
  c1->Write("depErangeEmucanv");

  hdepErangeEp_all->Write();
  hdepErangeEp_Muon->Write();
  hdepErangeEp_Proton->Write();
  hdepErangeEp_Pion->Write();

  hdepErangeEp_all->Draw();
  hdepErangeEp_Muon->Draw("same");
  hdepErangeEp_Pion->Draw("same");
  hdepErangeEp_Proton->Draw("same");
  l->Draw();
  c1->Print("depErangeEp.pdf");
  c1->Write("depErangeEpcanv");

  hdepEmrangeEmuvsp_all->Draw();
  hdepEmrangeEmuvsp_Muon->Draw("same");
  hdepEmrangeEmuvsp_Pion->Draw("same");
  hdepEmrangeEmuvsp_Proton->Draw("same");
  l->Draw();
  c1->Print("depEmrangeEmuvsp.pdf");
  c1->Write("depEmrangeEmuvspcanv");

  hdepEminusrangeEmu_all->Write();
  hdepEminusrangeEmu_Muon->Write();
  hdepEminusrangeEmu_Proton->Write();
  hdepEminusrangeEmu_Pion->Write();

  hdepEminusrangeEmu_all->Draw();
  THStack *st_depEminusrangeEmu = new THStack("st_depEminusrangeEmu","");
  st_depEminusrangeEmu->Add(hdepEminusrangeEmu_Muon);
  st_depEminusrangeEmu->Add(hdepEminusrangeEmu_Pion);
  st_depEminusrangeEmu->Add(hdepEminusrangeEmu_Proton);
  st_depEminusrangeEmu->Draw("same");
  l->Draw();
  c1->Print("depEminusrangeEmu.pdf");
  c1->Write("depEminusrangeEmucanv");

  hdepEminusrangeEp_all->Write();
  hdepEminusrangeEp_Muon->Write();
  hdepEminusrangeEp_Proton->Write();
  hdepEminusrangeEp_Pion->Write();

  hdepEminusrangeEp_all->Draw();
  THStack *st_depEminusrangeEp = new THStack("st_depEminusrangeEp","");
  st_depEminusrangeEp->Add(hdepEminusrangeEp_Proton);
  st_depEminusrangeEp->Add(hdepEminusrangeEp_Muon);
  st_depEminusrangeEp->Add(hdepEminusrangeEp_Pion);
  st_depEminusrangeEp->Draw("same");
  l->Draw();
  c1->Print("depEminusrangeEp.pdf");
  c1->Write("depEminusrangeEpcanv");

  hdepEminusrangeEp_cut_all->Write();
  hdepEminusrangeEp_cut_Muon->Write();
  hdepEminusrangeEp_cut_Proton->Write();
  hdepEminusrangeEp_cut_Pion->Write();

  hdepEminusrangeEp_cut_all->Draw();
  THStack *st_depEminusrangeEp_cut = new THStack("st_depEminusrangeEp_cut","");
  st_depEminusrangeEp_cut->Add(hdepEminusrangeEp_cut_Proton);
  st_depEminusrangeEp_cut->Add(hdepEminusrangeEp_cut_Muon);
  st_depEminusrangeEp_cut->Add(hdepEminusrangeEp_cut_Pion);
  st_depEminusrangeEp_cut->Draw("same");
  l->Draw();
  c1->Print("depEminusrangeEp_cut.pdf");
  c1->Write("depEminusrangeEp_cut_canv");

}
