struct treevars{
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
};

// ----------------------------------------------------------- //

void settreevars(TTree *intree, treevars *varstoset){
  intree->SetBranchStatus("*",0);
  intree->SetBranchStatus("TPCObj_PFP_track_dedx_perhit",1);
  intree->SetBranchAddress("TPCObj_PFP_track_dedx_perhit",&(varstoset->dEdxvec));
  intree->SetBranchStatus("TPCObj_PFP_track_resrange_perhit",1);
  intree->SetBranchAddress("TPCObj_PFP_track_resrange_perhit",&(varstoset->resrangevec));
  intree->SetBranchStatus("TPCObj_PFP_n2LLH_fwd_mu",1);
  intree->SetBranchAddress("TPCObj_PFP_n2LLH_fwd_mu",&(varstoset->n2LLH_fwd_mu));
  intree->SetBranchStatus("TPCObj_PFP_n2LLH_bwd_mu",1);
  intree->SetBranchAddress("TPCObj_PFP_n2LLH_bwd_mu",&(varstoset->n2LLH_bwd_mu));
  intree->SetBranchStatus("TPCObj_PFP_n2LLH_fwd_p",1);
  intree->SetBranchAddress("TPCObj_PFP_n2LLH_fwd_p",&(varstoset->n2LLH_fwd_p));
  intree->SetBranchStatus("TPCObj_PFP_n2LLH_bwd_p",1);
  intree->SetBranchAddress("TPCObj_PFP_n2LLH_bwd_p",&(varstoset->n2LLH_bwd_p));
  intree->SetBranchStatus("TPCObj_PFP_n2LLH_MIP",1);
  intree->SetBranchAddress("TPCObj_PFP_n2LLH_MIP",&(varstoset->n2LLH_MIP));
  intree->SetBranchStatus("TPCObj_PFP_track_depE",1);
  intree->SetBranchAddress("TPCObj_PFP_track_depE",&(varstoset->depE));
  intree->SetBranchStatus("TPCObj_PFP_track_rangeE_mu",1);
  intree->SetBranchAddress("TPCObj_PFP_track_rangeE_mu",&(varstoset->rangeE_mu));
  intree->SetBranchStatus("TPCObj_PFP_track_rangeE_p",1);
  intree->SetBranchAddress("TPCObj_PFP_track_rangeE_p",&(varstoset->rangeE_p));
  intree->SetBranchStatus("TPCObj_PFP_truePDG",1);
  intree->SetBranchAddress("TPCObj_PFP_truePDG",&(varstoset->truePDG));
  intree->SetBranchStatus("TPCObj_PFP_isTrack",1);
  intree->SetBranchAddress("TPCObj_PFP_isTrack",&(varstoset->isTrack));
}

// ----------------------------------------------------------- //

struct histograms{
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
};

// ----------------------------------------------------------- //

void hists_setstyleMC(histograms *hists){
  // Set histogram styles
  hists->hnegLLmu_all->SetLineColor(kBlack);
  hists->hnegLLmu_all->SetFillColor(TColor::GetColor(197,197,197));
  hists->hnegLLMIP_all->SetLineColor(kBlack);
  hists->hnegLLMIP_all->SetFillColor(TColor::GetColor(197,197,197));
  hists->hnegLLmuMIP_all->SetLineColor(kBlack);
  hists->hnegLLmuMIP_all->SetFillColor(TColor::GetColor(197,197,197));
  hists->hnegLLp_all->SetLineColor(kBlack);
  hists->hnegLLp_all->SetFillColor(TColor::GetColor(197,197,197));
  hists->hnegLLmuminusp_all->SetLineColor(kBlack);
  hists->hnegLLmuminusp_all->SetFillColor(TColor::GetColor(197,197,197));
  hists->hnegLLmuMIPminusp_all->SetLineColor(kBlack);
  hists->hnegLLmuMIPminusp_all->SetFillColor(TColor::GetColor(197,197,197));
  hists->hdepErangeEmu_all->SetMarkerColor(TColor::GetColor(197,197,197));
  hists->hdepErangeEp_all->SetMarkerColor(TColor::GetColor(197,197,197));
  hists->hdepEmrangeEmuvsp_all->SetMarkerColor(TColor::GetColor(197,197,197));
  hists->hdepEminusrangeEmu_all->SetLineColor(TColor::GetColor(197,197,197));
  hists->hdepEminusrangeEmu_all->SetFillColor(TColor::GetColor(197,197,197));
  hists->hdepEminusrangeEp_all->SetLineColor(TColor::GetColor(197,197,197));
  hists->hdepEminusrangeEp_all->SetFillColor(TColor::GetColor(197,197,197));
  hists->hdepEminusrangeEp_cut_all->SetLineColor(TColor::GetColor(197,197,197));
  hists->hdepEminusrangeEp_cut_all->SetFillColor(TColor::GetColor(197,197,197));

  hists->hnegLLmu_Muon->SetLineColor(TColor::GetColor(8,64,129));
  hists->hnegLLmu_Muon->SetFillColor(TColor::GetColor(8,64,129));
  hists->hnegLLMIP_Muon->SetLineColor(TColor::GetColor(8,64,129));
  hists->hnegLLMIP_Muon->SetFillColor(TColor::GetColor(8,64,129));
  hists->hnegLLmuMIP_Muon->SetLineColor(TColor::GetColor(8,64,129));
  hists->hnegLLmuMIP_Muon->SetFillColor(TColor::GetColor(8,64,129));
  hists->hnegLLp_Muon->SetLineColor(TColor::GetColor(8,64,129));
  hists->hnegLLp_Muon->SetFillColor(TColor::GetColor(8,64,129));
  hists->hnegLLmuminusp_Muon->SetLineColor(TColor::GetColor(8,64,129));
  hists->hnegLLmuminusp_Muon->SetFillColor(TColor::GetColor(8,64,129));
  hists->hnegLLmuMIPminusp_Muon->SetLineColor(TColor::GetColor(8,64,129));
  hists->hnegLLmuMIPminusp_Muon->SetFillColor(TColor::GetColor(8,64,129));
  hists->hdepErangeEmu_Muon->SetMarkerColor(TColor::GetColor(8,64,129));
  hists->hdepErangeEp_Muon->SetMarkerColor(TColor::GetColor(8,64,129));
  hists->hdepEmrangeEmuvsp_Muon->SetMarkerColor(TColor::GetColor(8,64,129));
  hists->hdepEminusrangeEmu_Muon->SetLineColor(TColor::GetColor(8,64,129));
  hists->hdepEminusrangeEmu_Muon->SetFillColor(TColor::GetColor(8,64,129));
  hists->hdepEminusrangeEp_Muon->SetLineColor(TColor::GetColor(8,64,129));
  hists->hdepEminusrangeEp_Muon->SetFillColor(TColor::GetColor(8,64,129));
  hists->hdepEminusrangeEp_cut_Muon->SetLineColor(TColor::GetColor(8,64,129));
  hists->hdepEminusrangeEp_cut_Muon->SetFillColor(TColor::GetColor(8,64,129));

  hists->hnegLLmu_Proton->SetLineColor(TColor::GetColor(215, 48, 39));
  hists->hnegLLmu_Proton->SetFillColor(TColor::GetColor(215, 48, 39));
  hists->hnegLLMIP_Proton->SetLineColor(TColor::GetColor(215, 48, 39));
  hists->hnegLLMIP_Proton->SetFillColor(TColor::GetColor(215, 48, 39));
  hists->hnegLLmuMIP_Proton->SetLineColor(TColor::GetColor(215, 48, 39));
  hists->hnegLLmuMIP_Proton->SetFillColor(TColor::GetColor(215, 48, 39));
  hists->hnegLLp_Proton->SetLineColor(TColor::GetColor(215, 48, 39));
  hists->hnegLLp_Proton->SetFillColor(TColor::GetColor(215, 48, 39));
  hists->hnegLLmuminusp_Proton->SetLineColor(TColor::GetColor(215, 48, 39));
  hists->hnegLLmuminusp_Proton->SetFillColor(TColor::GetColor(215, 48, 39));
  hists->hnegLLmuMIPminusp_Proton->SetLineColor(TColor::GetColor(215, 48, 39));
  hists->hnegLLmuMIPminusp_Proton->SetFillColor(TColor::GetColor(215, 48, 39));
  hists->hdepErangeEmu_Proton->SetMarkerColor(TColor::GetColor(215, 48, 39));
  hists->hdepErangeEp_Proton->SetMarkerColor(TColor::GetColor(215, 48, 39));
  hists->hdepEmrangeEmuvsp_Proton->SetMarkerColor(TColor::GetColor(215, 48, 39));
  hists->hdepEminusrangeEmu_Proton->SetLineColor(TColor::GetColor(215, 48, 39));
  hists->hdepEminusrangeEmu_Proton->SetFillColor(TColor::GetColor(215, 48, 39));
  hists->hdepEminusrangeEp_Proton->SetLineColor(TColor::GetColor(215, 48, 39));
  hists->hdepEminusrangeEp_Proton->SetFillColor(TColor::GetColor(215, 48, 39));
  hists->hdepEminusrangeEp_cut_Proton->SetLineColor(TColor::GetColor(215, 48, 39));
  hists->hdepEminusrangeEp_cut_Proton->SetFillColor(TColor::GetColor(215, 48, 39));

  hists->hnegLLmu_Pion->SetLineColor(TColor::GetColor(166,217,106));
  hists->hnegLLmu_Pion->SetFillColor(TColor::GetColor(166,217,106));
  hists->hnegLLMIP_Pion->SetLineColor(TColor::GetColor(166,217,106));
  hists->hnegLLMIP_Pion->SetFillColor(TColor::GetColor(166,217,106));
  hists->hnegLLmuMIP_Pion->SetLineColor(TColor::GetColor(166,217,106));
  hists->hnegLLmuMIP_Pion->SetFillColor(TColor::GetColor(166,217,106));
  hists->hnegLLp_Pion->SetLineColor(TColor::GetColor(166,217,106));
  hists->hnegLLp_Pion->SetFillColor(TColor::GetColor(166,217,106));
  hists->hnegLLmuminusp_Pion->SetLineColor(TColor::GetColor(166,217,106));
  hists->hnegLLmuminusp_Pion->SetFillColor(TColor::GetColor(166,217,106));
  hists->hnegLLmuMIPminusp_Pion->SetLineColor(TColor::GetColor(166,217,106));
  hists->hnegLLmuMIPminusp_Pion->SetFillColor(TColor::GetColor(166,217,106));
  hists->hdepErangeEmu_Pion->SetMarkerColor(TColor::GetColor(166,217,106));
  hists->hdepErangeEp_Pion->SetMarkerColor(TColor::GetColor(166,217,106));
  hists->hdepEmrangeEmuvsp_Pion->SetMarkerColor(TColor::GetColor(166,217,106));
  hists->hdepEminusrangeEmu_Pion->SetLineColor(TColor::GetColor(166,217,106));
  hists->hdepEminusrangeEmu_Pion->SetFillColor(TColor::GetColor(166,217,106));
  hists->hdepEminusrangeEp_Pion->SetLineColor(TColor::GetColor(166,217,106));
  hists->hdepEminusrangeEp_Pion->SetFillColor(TColor::GetColor(166,217,106));
  hists->hdepEminusrangeEp_cut_Pion->SetLineColor(TColor::GetColor(166,217,106));
  hists->hdepEminusrangeEp_cut_Pion->SetFillColor(TColor::GetColor(166,217,106));
}

// ----------------------------------------------------------- //

void hists_areanormMC(histograms *hists){
  hists->hnegLLmu_Muon->Scale(1.0/hists->hnegLLmu_all->Integral());
  hists->hnegLLMIP_Muon->Scale(1.0/hists->hnegLLMIP_all->Integral());
  hists->hnegLLmuMIP_Muon->Scale(1.0/hists->hnegLLmuMIP_all->Integral());
  hists->hnegLLp_Muon->Scale(1.0/hists->hnegLLp_all->Integral());
  hists->hnegLLmuminusp_Muon->Scale(1.0/hists->hnegLLmuminusp_all->Integral());
  hists->hnegLLmuMIPminusp_Muon->Scale(1.0/hists->hnegLLmuMIPminusp_all->Integral());
  hists->hdepErangeEmu_Muon->Scale(1.0/hists->hdepErangeEmu_all->Integral());
  hists->hdepErangeEp_Muon->Scale(1.0/hists->hdepErangeEp_all->Integral());
  hists->hdepEmrangeEmuvsp_Muon->Scale(1.0/hists->hdepEmrangeEmuvsp_all->Integral());
  hists->hdepEminusrangeEmu_Muon->Scale(1.0/hists->hdepEminusrangeEmu_all->Integral());
  hists->hdepEminusrangeEp_Muon->Scale(1.0/hists->hdepEminusrangeEp_all->Integral());
  hists->hdepEminusrangeEp_cut_Muon->Scale(1.0/hists->hdepEminusrangeEp_cut_all->Integral());

  hists->hnegLLmu_Proton->Scale(1.0/hists->hnegLLmu_all->Integral());
  hists->hnegLLMIP_Proton->Scale(1.0/hists->hnegLLMIP_all->Integral());
  hists->hnegLLmuMIP_Proton->Scale(1.0/hists->hnegLLmuMIP_all->Integral());
  hists->hnegLLp_Proton->Scale(1.0/hists->hnegLLp_all->Integral());
  hists->hnegLLmuminusp_Proton->Scale(1.0/hists->hnegLLmuminusp_all->Integral());
  hists->hnegLLmuMIPminusp_Proton->Scale(1.0/hists->hnegLLmuMIPminusp_all->Integral());
  hists->hdepErangeEmu_Proton->Scale(1.0/hists->hdepErangeEmu_all->Integral());
  hists->hdepErangeEp_Proton->Scale(1.0/hists->hdepErangeEp_all->Integral());
  hists->hdepEmrangeEmuvsp_Proton->Scale(1.0/hists->hdepEmrangeEmuvsp_all->Integral());
  hists->hdepEminusrangeEmu_Proton->Scale(1.0/hists->hdepEminusrangeEmu_all->Integral());
  hists->hdepEminusrangeEp_Proton->Scale(1.0/hists->hdepEminusrangeEp_all->Integral());
  hists->hdepEminusrangeEp_cut_Proton->Scale(1.0/hists->hdepEminusrangeEp_cut_all->Integral());

  hists->hnegLLmu_Pion->Scale(1.0/hists->hnegLLmu_all->Integral());
  hists->hnegLLMIP_Pion->Scale(1.0/hists->hnegLLMIP_all->Integral());
  hists->hnegLLmuMIP_Pion->Scale(1.0/hists->hnegLLmuMIP_all->Integral());
  hists->hnegLLp_Pion->Scale(1.0/hists->hnegLLp_all->Integral());
  hists->hnegLLmuminusp_Pion->Scale(1.0/hists->hnegLLmuminusp_all->Integral());
  hists->hnegLLmuMIPminusp_Pion->Scale(1.0/hists->hnegLLmuMIPminusp_all->Integral());
  hists->hdepErangeEmu_Pion->Scale(1.0/hists->hdepErangeEmu_all->Integral());
  hists->hdepErangeEp_Pion->Scale(1.0/hists->hdepErangeEp_all->Integral());
  hists->hdepEmrangeEmuvsp_Pion->Scale(1.0/hists->hdepEmrangeEmuvsp_all->Integral());
  hists->hdepEminusrangeEmu_Pion->Scale(1.0/hists->hdepEminusrangeEmu_all->Integral());
  hists->hdepEminusrangeEp_Pion->Scale(1.0/hists->hdepEminusrangeEp_all->Integral());
  hists->hdepEminusrangeEp_cut_Pion->Scale(1.0/hists->hdepEminusrangeEp_cut_all->Integral());

  hists->hnegLLmu_all->Scale(1.0/hists->hnegLLmu_all->Integral());
  hists->hnegLLMIP_all->Scale(1.0/hists->hnegLLMIP_all->Integral());
  hists->hnegLLmuMIP_all->Scale(1.0/hists->hnegLLmuMIP_all->Integral());
  hists->hnegLLp_all->Scale(1.0/hists->hnegLLp_all->Integral());
  hists->hnegLLmuminusp_all->Scale(1.0/hists->hnegLLmuminusp_all->Integral());
  hists->hnegLLmuMIPminusp_all->Scale(1.0/hists->hnegLLmuMIPminusp_all->Integral());
  hists->hdepErangeEmu_all->Scale(1.0/hists->hdepErangeEmu_all->Integral());
  hists->hdepErangeEp_all->Scale(1.0/hists->hdepErangeEp_all->Integral());
  hists->hdepEmrangeEmuvsp_all->Scale(1.0/hists->hdepEmrangeEmuvsp_all->Integral());
  hists->hdepEminusrangeEmu_all->Scale(1.0/hists->hdepEminusrangeEmu_all->Integral());
  hists->hdepEminusrangeEp_all->Scale(1.0/hists->hdepEminusrangeEp_all->Integral());
  hists->hdepEminusrangeEp_cut_all->Scale(1.0/hists->hdepEminusrangeEp_cut_all->Integral());
}

// ----------------------------------------------------------- //

void hists_setstyleData(histograms *hists){
  // Set histogram styles
  // Only bother to style the _all histograms because this is all we plot for data
  hists->hnegLLmu_all->SetLineColor(kBlack);
  hists->hnegLLmu_all->SetMarkerColor(kBlack);
  hists->hnegLLMIP_all->SetLineColor(kBlack);
  hists->hnegLLMIP_all->SetMarkerColor(kBlack);
  hists->hnegLLmuMIP_all->SetLineColor(kBlack);
  hists->hnegLLmuMIP_all->SetMarkerColor(kBlack);
  hists->hnegLLp_all->SetLineColor(kBlack);
  hists->hnegLLp_all->SetMarkerColor(kBlack);
  hists->hnegLLmuminusp_all->SetLineColor(kBlack);
  hists->hnegLLmuminusp_all->SetMarkerColor(kBlack);
  hists->hnegLLmuMIPminusp_all->SetLineColor(kBlack);
  hists->hnegLLmuMIPminusp_all->SetMarkerColor(kBlack);
  hists->hdepErangeEmu_all->SetMarkerColor(kBlack);
  hists->hdepErangeEp_all->SetMarkerColor(kBlack);
  hists->hdepEmrangeEmuvsp_all->SetMarkerColor(kBlack);
  hists->hdepEminusrangeEmu_all->SetLineColor(kBlack);
  hists->hdepEminusrangeEmu_all->SetMarkerColor(kBlack);
  hists->hdepEminusrangeEp_all->SetLineColor(kBlack);
  hists->hdepEminusrangeEp_all->SetMarkerColor(kBlack);
  hists->hdepEminusrangeEp_cut_all->SetLineColor(kBlack);
  hists->hdepEminusrangeEp_cut_all->SetMarkerColor(kBlack);
}

// ----------------------------------------------------------- //

void FillHists(treevars *vars, histograms *hists){
  for (size_t i_pfp=0; i_pfp<vars->truePDG->size(); i_pfp++){

    // PID variables are only filled for tracks. If one isn't filled, probably none of them are, so just skip the entire track.
    if (vars->n2LLH_fwd_mu->at(i_pfp) == -999) continue;
    //if (!rangeE_mu) continue;

    // Skip PFPs that Pandora doesn't reconstruct as tracks
    // This is probably not what we want to do for the CC1pi analysis, but is useful for making plots for general consumption/PID presentations
    if (!vars->isTrack->at(i_pfp)) continue;

    // Make dE/dx plots
    std::vector<double> dEdx = vars->dEdxvec->at(i_pfp);
    std::vector<double> resrange = vars->resrangevec->at(i_pfp);

    for (size_t i_dEdx=0; i_dEdx<dEdx.size(); i_dEdx++){
      if (resrange.at(i_dEdx) > 100 && resrange.at(i_dEdx) < 150){
        hists->hdEdx_all->Fill(dEdx.at(i_dEdx));
        hists->hdEdx_rr_all->Fill(resrange.at(i_dEdx),dEdx.at(i_dEdx));
      }
      if (resrange.at(i_dEdx) > 100 && resrange.at(i_dEdx) < 150 && TMath::Abs(vars->truePDG->at(i_pfp))==13){
        hists->hdEdx_Muon->Fill(dEdx.at(i_dEdx));
        hists->hdEdx_rr_Muon->Fill(resrange.at(i_dEdx),dEdx.at(i_dEdx));
      }
      if (resrange.at(i_dEdx) > 100 && resrange.at(i_dEdx) < 150 && TMath::Abs(vars->truePDG->at(i_pfp))==2212){
        hists->hdEdx_Proton->Fill(dEdx.at(i_dEdx));
        hists->hdEdx_rr_Proton->Fill(resrange.at(i_dEdx),dEdx.at(i_dEdx));
      }
    } // loop over hits (dEdx and resrange)


    // Fill n2llh histograms
    double tr_n2llh_mu = std::min(vars->n2LLH_fwd_mu->at(i_pfp),vars->n2LLH_bwd_mu->at(i_pfp));
    double tr_n2llh_MIP = vars->n2LLH_MIP->at(i_pfp);
    double tr_n2llh_muMIP = std::min({vars->n2LLH_fwd_mu->at(i_pfp),vars->n2LLH_bwd_mu->at(i_pfp),vars->n2LLH_MIP->at(i_pfp)});
    double tr_n2llh_p = std::min(vars->n2LLH_fwd_p->at(i_pfp),vars->n2LLH_bwd_p->at(i_pfp));
    double tr_n2llh_muminusp = tr_n2llh_mu - tr_n2llh_p;
    double tr_n2llh_muMIPminusp = tr_n2llh_muMIP - tr_n2llh_p;

    hists->hnegLLmu_all->Fill(tr_n2llh_mu);
    hists->hnegLLMIP_all->Fill(tr_n2llh_MIP);
    hists->hnegLLmuMIP_all->Fill(tr_n2llh_muMIP);
    hists->hnegLLp_all->Fill(tr_n2llh_p);
    hists->hnegLLmuminusp_all->Fill(tr_n2llh_muminusp);
    hists->hnegLLmuMIPminusp_all->Fill(tr_n2llh_muMIPminusp);
    hists->hdepErangeEmu_all->Fill(vars->depE->at(i_pfp),vars->rangeE_mu->at(i_pfp));
    hists->hdepErangeEp_all->Fill(vars->depE->at(i_pfp),vars->rangeE_p->at(i_pfp));
    hists->hdepEminusrangeEmu_all->Fill(vars->depE->at(i_pfp)-vars->rangeE_mu->at(i_pfp));
    hists->hdepEminusrangeEp_all->Fill(vars->depE->at(i_pfp)-vars->rangeE_p->at(i_pfp));
    hists->hdepEmrangeEmuvsp_all->Fill(vars->depE->at(i_pfp)-vars->rangeE_mu->at(i_pfp),vars->depE->at(i_pfp)-vars->rangeE_p->at(i_pfp));

    if (tr_n2llh_muMIPminusp < -2){
      hists->hdepEminusrangeEp_cut_all->Fill(vars->depE->at(i_pfp)-vars->rangeE_p->at(i_pfp));
    }

    if (TMath::Abs(vars->truePDG->at(i_pfp))==13){
      hists->hnegLLmu_Muon->Fill(tr_n2llh_mu);
      hists->hnegLLMIP_Muon->Fill(tr_n2llh_MIP);
      hists->hnegLLmuMIP_Muon->Fill(tr_n2llh_muMIP);
      hists->hnegLLp_Muon->Fill(tr_n2llh_p);
      hists->hnegLLmuminusp_Muon->Fill(tr_n2llh_muminusp);
      hists->hnegLLmuMIPminusp_Muon->Fill(tr_n2llh_muMIPminusp);
      hists->hdepErangeEmu_Muon->Fill(vars->depE->at(i_pfp),vars->rangeE_mu->at(i_pfp));
      hists->hdepErangeEp_Muon->Fill(vars->depE->at(i_pfp),vars->rangeE_p->at(i_pfp));
      hists->hdepEminusrangeEmu_Muon->Fill(vars->depE->at(i_pfp)-vars->rangeE_mu->at(i_pfp));
      hists->hdepEminusrangeEp_Muon->Fill(vars->depE->at(i_pfp)-vars->rangeE_p->at(i_pfp));
      hists->hdepEmrangeEmuvsp_Muon->Fill(vars->depE->at(i_pfp)-vars->rangeE_mu->at(i_pfp),vars->depE->at(i_pfp)-vars->rangeE_p->at(i_pfp));

      if (tr_n2llh_muMIPminusp < -2){
        hists->hdepEminusrangeEp_cut_Muon->Fill(vars->depE->at(i_pfp)-vars->rangeE_p->at(i_pfp));
      }
    }
    else if (TMath::Abs(vars->truePDG->at(i_pfp))==2212){
      hists->hnegLLmu_Proton->Fill(tr_n2llh_mu);
      hists->hnegLLMIP_Proton->Fill(tr_n2llh_MIP);
      hists->hnegLLmuMIP_Proton->Fill(tr_n2llh_muMIP);
      hists->hnegLLp_Proton->Fill(tr_n2llh_p);
      hists->hnegLLmuminusp_Proton->Fill(tr_n2llh_muminusp);
      hists->hnegLLmuMIPminusp_Proton->Fill(tr_n2llh_muMIPminusp);
      hists->hdepErangeEmu_Proton->Fill(vars->depE->at(i_pfp),vars->rangeE_mu->at(i_pfp));
      hists->hdepErangeEp_Proton->Fill(vars->depE->at(i_pfp),vars->rangeE_p->at(i_pfp));
      hists->hdepEminusrangeEmu_Proton->Fill(vars->depE->at(i_pfp)-vars->rangeE_mu->at(i_pfp));
      hists->hdepEminusrangeEp_Proton->Fill(vars->depE->at(i_pfp)-vars->rangeE_p->at(i_pfp));
      hists->hdepEmrangeEmuvsp_Proton->Fill(vars->depE->at(i_pfp)-vars->rangeE_mu->at(i_pfp),vars->depE->at(i_pfp)-vars->rangeE_p->at(i_pfp));

      if (tr_n2llh_muMIPminusp < -2){
        hists->hdepEminusrangeEp_cut_Proton->Fill(vars->depE->at(i_pfp)-vars->rangeE_p->at(i_pfp));
      }
    }
    else if (TMath::Abs(vars->truePDG->at(i_pfp))==211){
      hists->hnegLLmu_Pion->Fill(tr_n2llh_mu);
      hists->hnegLLMIP_Pion->Fill(tr_n2llh_MIP);
      hists->hnegLLmuMIP_Pion->Fill(tr_n2llh_muMIP);
      hists->hnegLLp_Pion->Fill(tr_n2llh_p);
      hists->hnegLLmuminusp_Pion->Fill(tr_n2llh_muminusp);
      hists->hnegLLmuMIPminusp_Pion->Fill(tr_n2llh_muMIPminusp);
      hists->hdepErangeEmu_Pion->Fill(vars->depE->at(i_pfp),vars->rangeE_mu->at(i_pfp));
      hists->hdepErangeEp_Pion->Fill(vars->depE->at(i_pfp),vars->rangeE_p->at(i_pfp));
      hists->hdepEminusrangeEmu_Pion->Fill(vars->depE->at(i_pfp)-vars->rangeE_mu->at(i_pfp));
      hists->hdepEminusrangeEp_Pion->Fill(vars->depE->at(i_pfp)-vars->rangeE_p->at(i_pfp));
      hists->hdepEmrangeEmuvsp_Pion->Fill(vars->depE->at(i_pfp)-vars->rangeE_mu->at(i_pfp),vars->depE->at(i_pfp)-vars->rangeE_p->at(i_pfp));

      if (tr_n2llh_muMIPminusp < -2){
        hists->hdepEminusrangeEp_cut_Pion->Fill(vars->depE->at(i_pfp)-vars->rangeE_p->at(i_pfp));
      }
    }

  } // loop over PFPs
}

// ----------------------------------------------------------- //

void MakePlots(TFile *outfile, histograms *MChists, histograms *EXThists=nullptr, histograms *BNBhists=nullptr, double EXTtoBNBscaling=0){
  TLegend *l = new TLegend(0.65,0.7,0.8,0.85);
  l->SetFillColor(kWhite);
  l->SetLineColor(kWhite);
  l->SetTextFont(132);
  l->AddEntry(MChists->hnegLLmu_Muon,"True muons","f");
  l->AddEntry(MChists->hnegLLmu_Proton,"True protons","f");
  l->AddEntry(MChists->hnegLLmu_Pion,"True pions","f");
  l->AddEntry(MChists->hnegLLmu_all,"MC other","f");
  if (EXThists && BNBhists) l->AddEntry(EXThists->hnegLLmu_all,"BNB-EXT data","p");

  // Plot things!
  outfile->cd();
  MChists->hdEdx_all->Write("hdEdx_all");
  if (MChists->hdEdx_Muon->Integral() > 0) MChists->hdEdx_Muon->Write("hdEdx_Muon");
  if (MChists->hdEdx_Proton->Integral() > 0) MChists->hdEdx_Proton->Write("hdEdx_Proton");
  MChists->hdEdx_rr_all->Write("hdEdx_rr_all");
  if (MChists->hdEdx_rr_Muon->Integral() > 0) MChists->hdEdx_rr_Muon->Write("hdEdx_rr_Muon");
  if (MChists->hdEdx_rr_Proton->Integral() > 0) MChists->hdEdx_rr_Proton->Write("hdEdx_rr_Proton");

  if (BNBhists) BNBhists->hdEdx_all->Write("hdEdx_all_BNB");
  if (EXThists) EXThists->hdEdx_all->Write("hdEdx_all_EXT");

  TCanvas *c1 = new TCanvas();

  MChists->hnegLLmu_all->Write();
  MChists->hnegLLmu_Muon->Write();
  MChists->hnegLLmu_Proton->Write();
  MChists->hnegLLmu_Pion->Write();

  if (BNBhists) BNBhists->hnegLLmu_all->Write("hnegLLmu_all_BNB");
  if (EXThists) EXThists->hnegLLmu_all->Write("hnegLLmu_all_EXT");

  MChists->hnegLLmu_all->Draw("hist");
  THStack *st_negLLmu = new THStack("st_negLLmu","");
  st_negLLmu->Add(MChists->hnegLLmu_Proton);
  st_negLLmu->Add(MChists->hnegLLmu_Muon);
  st_negLLmu->Add(MChists->hnegLLmu_Pion);
  st_negLLmu->Draw("same hist");
  if (BNBhists && EXThists){
    BNBhists->hnegLLmu_all->Add(EXThists->hnegLLmu_all,-1.0*EXTtoBNBscaling);
    BNBhists->hnegLLmu_all->Sumw2();
    BNBhists->hnegLLmu_all->Scale(1.0/BNBhists->hnegLLmu_all->Integral());
    BNBhists->hnegLLmu_all->Draw("same p e0");
  }
  l->Draw();
  c1->Print("neg2LLmu.pdf");
  c1->Write("neg2LLmucanv");

  if (BNBhists) BNBhists->hnegLLMIP_all->Write("hnegLLMIP_all_BNB");
  if (EXThists) EXThists->hnegLLMIP_all->Write("hnegLLMIP_all_EXT");

  MChists->hnegLLMIP_all->Write();
  MChists->hnegLLMIP_Muon->Write();
  MChists->hnegLLMIP_Proton->Write();
  MChists->hnegLLMIP_Pion->Write();

  MChists->hnegLLMIP_all->Draw("hist");
  THStack *st_negLLMIP = new THStack("st_negLLMIP","");
  st_negLLMIP->Add(MChists->hnegLLMIP_Proton);
  st_negLLMIP->Add(MChists->hnegLLMIP_Muon);
  st_negLLMIP->Add(MChists->hnegLLMIP_Pion);
  st_negLLMIP->Draw("same hist");
  if (BNBhists && EXThists){
    BNBhists->hnegLLMIP_all->Add(EXThists->hnegLLMIP_all,-1.0*EXTtoBNBscaling);
    BNBhists->hnegLLMIP_all->Sumw2();
    BNBhists->hnegLLMIP_all->Scale(1.0/BNBhists->hnegLLMIP_all->Integral());
    BNBhists->hnegLLMIP_all->Draw("same p e0");
  }
  l->Draw();
  c1->Print("neg2LLMIP.pdf");
  c1->Write("neg2LLMIPcanv");

  if (BNBhists) BNBhists->hnegLLmuMIP_all->Write("hnegLLmuMIP_all_BNB");
  if (EXThists) EXThists->hnegLLmuMIP_all->Write("hnegLLmuMIP_all_EXT");

  MChists->hnegLLmuMIP_all->Write();
  MChists->hnegLLmuMIP_Muon->Write();
  MChists->hnegLLmuMIP_Proton->Write();
  MChists->hnegLLmuMIP_Pion->Write();


  MChists->hnegLLmuMIP_all->Draw("hist");
  THStack *st_negLLmuMIP = new THStack("st_negLLmuMIP","");
  st_negLLmuMIP->Add(MChists->hnegLLmuMIP_Proton);
  st_negLLmuMIP->Add(MChists->hnegLLmuMIP_Muon);
  st_negLLmuMIP->Add(MChists->hnegLLmuMIP_Pion);
  st_negLLmuMIP->Draw("same hist");
  if (BNBhists && EXThists){
    BNBhists->hnegLLmuMIP_all->Add(EXThists->hnegLLmuMIP_all,-1.0*EXTtoBNBscaling);
    BNBhists->hnegLLmuMIP_all->Sumw2();
    BNBhists->hnegLLmuMIP_all->Scale(1.0/BNBhists->hnegLLmuMIP_all->Integral());
    BNBhists->hnegLLmuMIP_all->Draw("same p e0");
  }
  l->Draw();
  c1->Print("neg2LLmuMIP.pdf");
  c1->Write("neg2LLmuMIPcanv");

  if (BNBhists) BNBhists->hnegLLp_all->Write("hnegLLp_all_BNB");
  if (EXThists) EXThists->hnegLLp_all->Write("hnegLLp_all_EXT");

  MChists->hnegLLp_all->Write();
  MChists->hnegLLp_Muon->Write();
  MChists->hnegLLp_Proton->Write();
  MChists->hnegLLp_Pion->Write();

  MChists->hnegLLp_all->Draw("hist");
  THStack *st_negLLp = new THStack("st_negLLp","");
  st_negLLp->Add(MChists->hnegLLp_Proton);
  st_negLLp->Add(MChists->hnegLLp_Muon);
  st_negLLp->Add(MChists->hnegLLp_Pion);
  st_negLLp->Draw("same hist");
  if (BNBhists && EXThists){
    BNBhists->hnegLLp_all->Add(EXThists->hnegLLp_all,-1.0*EXTtoBNBscaling);
    BNBhists->hnegLLp_all->Sumw2();
    BNBhists->hnegLLp_all->Scale(1.0/BNBhists->hnegLLp_all->Integral());
    BNBhists->hnegLLp_all->Draw("same p e0");
  }
  l->Draw();
  c1->Print("neg2LLp.pdf");
  c1->Write("neg2LLpcanv");

  if (BNBhists) BNBhists->hnegLLmuminusp_all->Write("hnegLLmuminusp_all_BNB");
  if (EXThists) EXThists->hnegLLmuminusp_all->Write("hnegLLmuminusp_all_EXT");

  MChists->hnegLLmuminusp_all->Write();
  MChists->hnegLLmuminusp_Muon->Write();
  MChists->hnegLLmuminusp_Proton->Write();
  MChists->hnegLLmuminusp_Pion->Write();

  MChists->hnegLLmuminusp_all->Draw("hist");
  THStack *st_negLLmuminusp = new THStack("st_negLLmuminusp","");
  st_negLLmuminusp->Add(MChists->hnegLLmuminusp_Proton);
  st_negLLmuminusp->Add(MChists->hnegLLmuminusp_Muon);
  st_negLLmuminusp->Add(MChists->hnegLLmuminusp_Pion);
  st_negLLmuminusp->Draw("same hist");
  if (BNBhists && EXThists){
    BNBhists->hnegLLmuminusp_all->Add(EXThists->hnegLLmuminusp_all,-1.0*EXTtoBNBscaling);
    BNBhists->hnegLLmuminusp_all->Sumw2();
    BNBhists->hnegLLmuminusp_all->Scale(1.0/BNBhists->hnegLLmuminusp_all->Integral());
    BNBhists->hnegLLmuminusp_all->Draw("same p e0");
  }
  l->Draw();
  c1->Print("neg2LLmuminusp.pdf");
  c1->Write("neg2LLmuminuspcanv");

  if (BNBhists) BNBhists->hnegLLmuMIPminusp_all->Write("hnegLLmuMIPminusp_all_BNB");
  if (EXThists) EXThists->hnegLLmuMIPminusp_all->Write("hnegLLmuMIPminusp_all_EXT");

  MChists->hnegLLmuMIPminusp_all->Write();
  MChists->hnegLLmuMIPminusp_Muon->Write();
  MChists->hnegLLmuMIPminusp_Proton->Write();
  MChists->hnegLLmuMIPminusp_Pion->Write();

  MChists->hnegLLmuMIPminusp_all->Draw("hist");
  THStack *st_negLLmuMIPminusp = new THStack("st_negLLmuMIPminusp","");
  st_negLLmuMIPminusp->Add(MChists->hnegLLmuMIPminusp_Proton);
  st_negLLmuMIPminusp->Add(MChists->hnegLLmuMIPminusp_Muon);
  st_negLLmuMIPminusp->Add(MChists->hnegLLmuMIPminusp_Pion);
  st_negLLmuMIPminusp->Draw("same hist");
  if (BNBhists && EXThists){
    BNBhists->hnegLLmuMIPminusp_all->Add(EXThists->hnegLLmuMIPminusp_all,-1.0*EXTtoBNBscaling);
    BNBhists->hnegLLmuMIPminusp_all->Sumw2();
    BNBhists->hnegLLmuMIPminusp_all->Scale(1.0/BNBhists->hnegLLmuMIPminusp_all->Integral());
    BNBhists->hnegLLmuMIPminusp_all->Draw("same p e0");
  }
  l->Draw();
  c1->Print("neg2LLmuMIPminusp.pdf");
  c1->Write("neg2LLmuMIPminuspcanv");

  if (BNBhists) BNBhists->hdepErangeEmu_all->Write("hdepErangeEmu_all_BNB");
  if (EXThists) EXThists->hdepErangeEmu_all->Write("hdepErangeEmu_all_EXT");

  MChists->hdepErangeEmu_all->Write();
  MChists->hdepErangeEmu_Muon->Write();
  MChists->hdepErangeEmu_Proton->Write();
  MChists->hdepErangeEmu_Pion->Write();

  MChists->hdepErangeEmu_all->Draw();
  MChists->hdepErangeEmu_Muon->Draw("same");
  MChists->hdepErangeEmu_Pion->Draw("same");
  MChists->hdepErangeEmu_Proton->Draw("same");
  l->Draw();
  c1->Print("depErangeEmu.pdf");
  c1->Write("depErangeEmucanv");

  if (BNBhists) BNBhists->hdepErangeEp_all->Write("hdepErangeEp_all_BNB");
  if (EXThists) EXThists->hdepErangeEp_all->Write("hdepErangeEp_all_EXT");

  MChists->hdepErangeEp_all->Write();
  MChists->hdepErangeEp_Muon->Write();
  MChists->hdepErangeEp_Proton->Write();
  MChists->hdepErangeEp_Pion->Write();

  MChists->hdepErangeEp_all->Draw();
  MChists->hdepErangeEp_Muon->Draw("same");
  MChists->hdepErangeEp_Pion->Draw("same");
  MChists->hdepErangeEp_Proton->Draw("same");
  l->Draw();
  c1->Print("depErangeEp.pdf");
  c1->Write("depErangeEpcanv");

  MChists->hdepEmrangeEmuvsp_all->Draw();
  MChists->hdepEmrangeEmuvsp_Muon->Draw("same");
  MChists->hdepEmrangeEmuvsp_Pion->Draw("same");
  MChists->hdepEmrangeEmuvsp_Proton->Draw("same");
  l->Draw();
  c1->Print("depEmrangeEmuvsp.pdf");
  c1->Write("depEmrangeEmuvspcanv");

  if (BNBhists) BNBhists->hdepEminusrangeEmu_all->Write("hdepEminusrangeEmu_all_BNB");
  if (EXThists) EXThists->hdepEminusrangeEmu_all->Write("hdepEminusrangeEmu_all_EXT");

  MChists->hdepEminusrangeEmu_all->Write();
  MChists->hdepEminusrangeEmu_Muon->Write();
  MChists->hdepEminusrangeEmu_Proton->Write();
  MChists->hdepEminusrangeEmu_Pion->Write();

  MChists->hdepEminusrangeEmu_all->Draw("hist");
  THStack *st_depEminusrangeEmu = new THStack("st_depEminusrangeEmu","");
  st_depEminusrangeEmu->Add(MChists->hdepEminusrangeEmu_Muon);
  st_depEminusrangeEmu->Add(MChists->hdepEminusrangeEmu_Pion);
  st_depEminusrangeEmu->Add(MChists->hdepEminusrangeEmu_Proton);
  st_depEminusrangeEmu->Draw("same hist");
  if (BNBhists && EXThists){
    BNBhists->hdepEminusrangeEmu_all->Add(EXThists->hdepEminusrangeEmu_all,-1.0*EXTtoBNBscaling);
    BNBhists->hdepEminusrangeEmu_all->Sumw2();
    BNBhists->hdepEminusrangeEmu_all->Scale(1.0/BNBhists->hdepEminusrangeEmu_all->Integral());
    BNBhists->hdepEminusrangeEmu_all->Draw("same p e0");
  }
  l->Draw();
  c1->Print("depEminusrangeEmu.pdf");
  c1->Write("depEminusrangeEmucanv");

  if (BNBhists) BNBhists->hdepEminusrangeEp_all->Write("hdepEminusrangeEp_all_BNB");
  if (EXThists) EXThists->hdepEminusrangeEp_all->Write("hdepEminusrangeEp_all_EXT");

  MChists->hdepEminusrangeEp_all->Write();
  MChists->hdepEminusrangeEp_Muon->Write();
  MChists->hdepEminusrangeEp_Proton->Write();
  MChists->hdepEminusrangeEp_Pion->Write();

  MChists->hdepEminusrangeEp_all->Draw("hist");
  THStack *st_depEminusrangeEp = new THStack("st_depEminusrangeEp","");
  st_depEminusrangeEp->Add(MChists->hdepEminusrangeEp_Proton);
  st_depEminusrangeEp->Add(MChists->hdepEminusrangeEp_Muon);
  st_depEminusrangeEp->Add(MChists->hdepEminusrangeEp_Pion);
  st_depEminusrangeEp->Draw("same hist");
  if (BNBhists && EXThists){
    BNBhists->hdepEminusrangeEp_all->Add(EXThists->hdepEminusrangeEp_all,-1.0*EXTtoBNBscaling);
    BNBhists->hdepEminusrangeEp_all->Sumw2();
    BNBhists->hdepEminusrangeEp_all->Scale(1.0/BNBhists->hdepEminusrangeEp_all->Integral());
    BNBhists->hdepEminusrangeEp_all->Draw("same p e0");
  }
  l->Draw();
  c1->Print("depEminusrangeEp.pdf");
  c1->Write("depEminusrangeEpcanv");

  if (BNBhists) BNBhists->hdepEminusrangeEp_cut_all->Write("hdepEminusrangeEp_cut_all_BNB");
  if (EXThists) EXThists->hdepEminusrangeEp_cut_all->Write("hdepEminusrangeEp_cut_all_EXT");

  MChists->hdepEminusrangeEp_cut_all->Write();
  MChists->hdepEminusrangeEp_cut_Muon->Write();
  MChists->hdepEminusrangeEp_cut_Proton->Write();
  MChists->hdepEminusrangeEp_cut_Pion->Write();

  MChists->hdepEminusrangeEp_cut_all->Draw("hist");
  THStack *st_depEminusrangeEp_cut = new THStack("st_depEminusrangeEp_cut","");
  st_depEminusrangeEp_cut->Add(MChists->hdepEminusrangeEp_cut_Proton);
  st_depEminusrangeEp_cut->Add(MChists->hdepEminusrangeEp_cut_Muon);
  st_depEminusrangeEp_cut->Add(MChists->hdepEminusrangeEp_cut_Pion);
  st_depEminusrangeEp_cut->Draw("same hist");
  if (BNBhists && EXThists){
    BNBhists->hdepEminusrangeEp_cut_all->Add(EXThists->hdepEminusrangeEp_cut_all,-1.0*EXTtoBNBscaling);
    BNBhists->hdepEminusrangeEp_cut_all->Sumw2();
    BNBhists->hdepEminusrangeEp_cut_all->Scale(1.0/BNBhists->hdepEminusrangeEp_cut_all->Integral());
    BNBhists->hdepEminusrangeEp_cut_all->Draw("same p e0");
  }
  l->Draw();
  c1->Print("depEminusrangeEp_cut.pdf");
  c1->Write("depEminusrangeEp_cut_canv");
}

// ----------------------------------------------------------- //
// ----------------------------------------------------------- //

void MakePIDPlots(std::string filename, std::string filename_EXT="", std::string filename_BNB="", double EXTtoBNBscaling=0.)
{
  // Make output file
  TFile *outfile = new TFile("CC1piselec_PIDplots_out.root","recreate");

  // --- MC --- //
  // Get input file and TTree
  TFile *infile_MC = new TFile(filename.c_str(),"open");
  TTree *tree_MC   = (TTree*)infile_MC->Get("cc1piselec/outtree");

  // Make histograms
  histograms MChists;
  hists_setstyleMC(&MChists);

  // Variables to read from tree
  treevars MCvars;
  settreevars(tree_MC,&MCvars);

  // Loop through MC tree and fill variables/histograms
  for (size_t i_tree_MC=0; i_tree_MC<tree_MC->GetEntries(); i_tree_MC++){
    tree_MC->GetEntry(i_tree_MC);
    FillHists(&MCvars, &MChists);
  } // loop over entries in MC tree

  // Area-normalise MC histogrms
  hists_areanormMC(&MChists);

  // --- Data --- //
  TFile *infile_EXT;
  TTree *tree_EXT;
  TFile *infile_BNB;
  TTree *tree_BNB;
  if (filename_EXT != "" && filename_BNB != ""){
    std::cout << "Making Data/MC comparisons" << std::endl;
    // --- EXT --- //
    infile_EXT = new TFile(filename_EXT.c_str(),"open");
    tree_EXT   = (TTree*)infile_EXT->Get("cc1piselec/outtree");

    // Make histograms
    histograms EXThists;
    hists_setstyleData(&EXThists);

    // Variables to read from tree
    treevars EXTvars;
    settreevars(tree_EXT,&EXTvars);

    // Loop through EXT tree and fill variables/histograms
    if (filename_EXT != ""){
      for (size_t i_tree_EXT=0; i_tree_EXT<tree_EXT->GetEntries(); i_tree_EXT++){
        tree_EXT->GetEntry(i_tree_EXT);
        FillHists(&EXTvars, &EXThists);
      } // loop over entries in EXT tree
    }

    // --- BNB --- //
    infile_BNB = new TFile(filename_BNB.c_str(),"open");
    tree_BNB   = (TTree*)infile_BNB->Get("cc1piselec/outtree");

    // Make histograms
    histograms BNBhists;
    hists_setstyleData(&BNBhists);

    // Variables to read from tree
    treevars BNBvars;
    settreevars(tree_BNB,&BNBvars);

    // Loop through BNB tree and fill variables/histograms
    if (filename_BNB != ""){
      for (size_t i_tree_BNB=0; i_tree_BNB<tree_BNB->GetEntries(); i_tree_BNB++){
        tree_BNB->GetEntry(i_tree_BNB);
        FillHists(&BNBvars, &BNBhists);
      } // loop over entries in BNB tree
    }

    // Make Plots
     MakePlots(outfile, &MChists, &EXThists, &BNBhists, EXTtoBNBscaling);
  } // if EXT and BNB
  else{
    std::cout << "Making plots for MC only" << std::endl;
     MakePlots(outfile, &MChists);
  }


}
