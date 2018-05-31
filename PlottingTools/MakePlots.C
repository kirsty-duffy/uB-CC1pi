#include "MakePlots.h"

double GetPOT(TString FileName) {

   TChain *t_pot = new TChain("cc1piselec/pottree");
   t_pot -> Add(FileName);

   double pot;
   t_pot -> SetBranchStatus("*",0);
   t_pot -> SetBranchStatus("pot",1);
   t_pot -> SetBranchAddress("pot", &pot);

   double totalPOT = 0;

   const int potentries = t_pot -> GetEntries();
   for (int i = 0; i < potentries; i++) {
      t_pot -> GetEntry(i);
      totalPOT += pot;
   }

   return totalPOT;

}

//TPC boundary + FV cut
const double FVxmin = 0 + 12;
const double FVxmax = 256.35 - 12;
const double FVymin = -115.53 + 35;
const double FVymax = 117.47 - 35;
const double FVzmin = 0.1 + 25;
const double FVzmax = 1036.9 - 85;

//Cut out additional section of mostly dead wires in z
const double deadzmin = 675.1;
const double deadzmax = 775.1;

bool inFV(double x, double y, double z){
   if (x > FVxmin && x < FVxmax && y > FVymin && y < FVymax && z > FVzmin && z < FVzmax && (z < deadzmin || z > deadzmax)) return true;
   else return false;
}


void MakePlots(std::map<std::string,bool> SelectionCutflow, std::string SaveString, TString FileName) {

   //Setup tree...
   TChain *t = new TChain("cc1piselec/outtree");
   t -> Add(FileName);

   // Event variables
   bool isData;
   unsigned int run_num;
   unsigned int subrun_num;
   unsigned int event_num;

   NuIntTopology Truth_topology = kUnknown;

   // Selection result variables
   std::map<std::string,bool> *Marco_cutflow = nullptr;
   bool Marco_selected = 0;
   std::map<std::string,bool> *CC1picutflow = nullptr;
   
   // Track variables
   std::vector<double> *TPCObj_PFP_track_length = nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_track_start = nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_track_end = nullptr;
   std::vector<double> *TPCObj_PFP_track_theta = nullptr;
   std::vector<double> *TPCObj_PFP_track_phi = nullptr;
   std::vector<double> *TPCObj_PFP_track_mom = nullptr;
   std::vector<bool> *TPCObj_PFP_isMIP = nullptr;
   std::vector<bool> *TPCObj_PFP_track_isContained = nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_track_AngleBetweenTracks = nullptr;

   // Shower variables
   std::vector<double> *TPCObj_PFP_shower_length = nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_shower_start = nullptr;

   // PFP reco variables
   int TPCObj_NPFPs = 0;
   int TPCObj_NTracks = 0;
   int TPCObj_NShowers = 0;
//   std::vector<bool> *TPCObj_PFP_PandoraClassedAsTrack = nullptr;
//   std::vector<bool> *TPCObj_PFP_PandoraClassedAsShower = nullptr;
   std::vector<bool> *TPCObj_PFP_isDaughter = nullptr;
   std::vector<std::vector<int>> *TPCObj_PFP_daughterids = nullptr;
   std::vector<int> *TPCObj_PFP_id = nullptr;
   std::vector<double> *TPCObj_reco_vtx = nullptr;

   // PFP true variables
   std::vector<int> *TPCObj_PFP_MCPid = nullptr;
   std::vector<int> *TPCObj_PFP_truePDG = nullptr;
   std::vector<double> *TPCObj_PFP_trueE = nullptr;
   std::vector<double> *TPCObj_PFP_trueKE = nullptr;
   int TPCObj_origin = 0;
   int TPCObj_origin_extra = 0;

   // MCParticle variables 
   std::vector<int> *MCP_PDG = nullptr;
   std::vector<double> *MCP_length = nullptr;
   std::vector<std::string> *MCP_process = nullptr;
   std::vector<std::string> *MCP_endprocess = nullptr;
   std::vector<int> *MCP_numdaughters = nullptr;
   std::vector<double> *MCP_P = nullptr;
   std::vector<double> *MCP_Px = nullptr;
   std::vector<double> *MCP_Py = nullptr;
   std::vector<double> *MCP_Pz = nullptr;
   std::vector<double> *MCP_E = nullptr;
   std::vector<double> *MCP_KE = nullptr;
   std::vector<bool> *MCP_isContained = nullptr;

   // True neutrino variables
   std::vector<double> *nu_vtx = nullptr;
   std::vector<double> *nu_vtx_spacecharge = nullptr;
   bool nu_isCC = 0;
   int nu_PDG = 0;
   double nu_E = 0;

   // PID variables
   std::vector<std::vector<double>> *n2LLH_fwd_mu = nullptr;
   std::vector<std::vector<double>> *n2LLH_bwd_mu = nullptr;
   std::vector<std::vector<double>> *n2LLH_fwd_p = nullptr;
   std::vector<std::vector<double>> *n2LLH_bwd_p = nullptr;
   std::vector<std::vector<double>> *n2LLH_MIP = nullptr;

   std::vector<std::vector<double>> *TPCObj_PFP_track_depE = nullptr;
   std::vector<double> *TPCObj_PFP_track_rangeE_mu = nullptr;
   std::vector<double> *TPCObj_PFP_track_rangeE_p = nullptr;


   t -> SetBranchStatus("*",0);

   t -> SetBranchStatus("isData",1);
   t -> SetBranchAddress("isData", &isData);
   t -> SetBranchStatus("run_num",1);
   t -> SetBranchAddress("run_num", &run_num);
   t -> SetBranchStatus("subrun_num",1);
   t -> SetBranchAddress("subrun_num", &subrun_num);
   t -> SetBranchStatus("event_num",1);
   t -> SetBranchAddress("event_num", &event_num);

   t -> SetBranchStatus("Truth_topology",1);
   t -> SetBranchAddress("Truth_topology", &Truth_topology);

   t -> SetBranchStatus("Marco_cutflow",1);
   t -> SetBranchAddress("Marco_cutflow", &Marco_cutflow);
   t -> SetBranchStatus("Marco_selected",1);
   t -> SetBranchAddress("Marco_selected", &Marco_selected);
   t -> SetBranchStatus("CC1picutflow",1);
   t -> SetBranchAddress("CC1picutflow", &CC1picutflow);

   t -> SetBranchStatus("TPCObj_PFP_track_length",1);
   t -> SetBranchAddress("TPCObj_PFP_track_length", &TPCObj_PFP_track_length);
   t -> SetBranchStatus("TPCObj_PFP_track_start",1);
   t -> SetBranchAddress("TPCObj_PFP_track_start", &TPCObj_PFP_track_start);
   t -> SetBranchStatus("TPCObj_PFP_track_end",1);
   t -> SetBranchAddress("TPCObj_PFP_track_end", &TPCObj_PFP_track_end);
   t -> SetBranchStatus("TPCObj_PFP_track_theta",1);
   t -> SetBranchAddress("TPCObj_PFP_track_theta", &TPCObj_PFP_track_theta);
   t -> SetBranchStatus("TPCObj_PFP_track_phi",1);
   t -> SetBranchAddress("TPCObj_PFP_track_phi", &TPCObj_PFP_track_phi);
   t -> SetBranchStatus("TPCObj_PFP_track_mom",1);
   t -> SetBranchAddress("TPCObj_PFP_track_mom", &TPCObj_PFP_track_mom);
   t -> SetBranchStatus("TPCObj_PFP_isMIP",1);
   t -> SetBranchAddress("TPCObj_PFP_isMIP", &TPCObj_PFP_isMIP);
   t -> SetBranchStatus("TPCObj_PFP_track_isContained",1);
   t -> SetBranchAddress("TPCObj_PFP_track_isContained", &TPCObj_PFP_track_isContained);
   t -> SetBranchStatus("TPCObj_PFP_track_AngleBetweenTracks",1);
   t -> SetBranchAddress("TPCObj_PFP_track_AngleBetweenTracks", &TPCObj_PFP_track_AngleBetweenTracks);

   t -> SetBranchStatus("TPCObj_PFP_shower_length",1);
   t -> SetBranchAddress("TPCObj_PFP_shower_length", &TPCObj_PFP_shower_length);
   t -> SetBranchStatus("TPCObj_PFP_shower_start",1);
   t -> SetBranchAddress("TPCObj_PFP_shower_start", &TPCObj_PFP_shower_start);

   t -> SetBranchStatus("TPCObj_NPFPs",1);
   t -> SetBranchAddress("TPCObj_NPFPs", &TPCObj_NPFPs);
   t -> SetBranchStatus("TPCObj_NTracks",1);
   t -> SetBranchAddress("TPCObj_NTracks", &TPCObj_NTracks);
   t -> SetBranchStatus("TPCObj_NShowers",1);
   t -> SetBranchAddress("TPCObj_NShowers", &TPCObj_NShowers);
//   t -> SetBranchStatus("TPCObj_PFP_PandoraClassedAsTrack",1);
//   t -> SetBranchAddress("TPCObj_PFP_PandoraClassedAsTrack", &TPCObj_PFP_PandoraClassedAsTrack);
//   t -> SetBranchStatus("TPCObj_PFP_PandoraClassedAsShower",1);
//   t -> SetBranchAddress("TPCObj_PFP_PandoraClassedAsShower", &TPCObj_PFP_PandoraClassedAsShower);
   t -> SetBranchStatus("TPCObj_PFP_isDaughter",1);
   t -> SetBranchAddress("TPCObj_PFP_isDaughter", &TPCObj_PFP_isDaughter);
   t -> SetBranchStatus("TPCObj_PFP_daughterids",1);
   t -> SetBranchAddress("TPCObj_PFP_daughterids", &TPCObj_PFP_daughterids);
   t -> SetBranchStatus("TPCObj_PFP_id",1);
   t -> SetBranchAddress("TPCObj_PFP_id", &TPCObj_PFP_id);
   t -> SetBranchStatus("TPCObj_reco_vtx",1);
   t -> SetBranchAddress("TPCObj_reco_vtx", &TPCObj_reco_vtx);

   t -> SetBranchStatus("TPCObj_PFP_MCPid",1);
   t -> SetBranchAddress("TPCObj_PFP_MCPid", &TPCObj_PFP_MCPid);
   t -> SetBranchStatus("TPCObj_PFP_truePDG",1);
   t -> SetBranchAddress("TPCObj_PFP_truePDG", &TPCObj_PFP_truePDG);
   t -> SetBranchStatus("TPCObj_PFP_trueE",1);
   t -> SetBranchAddress("TPCObj_PFP_trueE", &TPCObj_PFP_trueE);
   t -> SetBranchStatus("TPCObj_PFP_trueKE",1);
   t -> SetBranchAddress("TPCObj_PFP_trueKE", &TPCObj_PFP_trueKE);
   t -> SetBranchStatus("TPCObj_origin",1);
   t -> SetBranchAddress("TPCObj_origin", &TPCObj_origin);
   t -> SetBranchStatus("TPCObj_origin_extra",1);
   t -> SetBranchAddress("TPCObj_origin_extra", &TPCObj_origin_extra);

   t -> SetBranchStatus("MCP_PDG",1);
   t -> SetBranchAddress("MCP_PDG", &MCP_PDG);
   t -> SetBranchStatus("MCP_length",1);
   t -> SetBranchAddress("MCP_length", &MCP_length);
   t -> SetBranchStatus("MCP_process",1);
   t -> SetBranchAddress("MCP_process", &MCP_process);
   t -> SetBranchStatus("MCP_endprocess",1);
   t -> SetBranchAddress("MCP_endprocess", &MCP_endprocess);
   t -> SetBranchStatus("MCP_numdaughters",1);
   t -> SetBranchAddress("MCP_numdaughters", &MCP_numdaughters);
   t -> SetBranchStatus("MCP_P",1);
   t -> SetBranchAddress("MCP_P", &MCP_P);
   t -> SetBranchStatus("MCP_Px",1);
   t -> SetBranchAddress("MCP_Px", &MCP_Px);
   t -> SetBranchStatus("MCP_Py",1);
   t -> SetBranchAddress("MCP_Py", &MCP_Py);
   t -> SetBranchStatus("MCP_Pz",1);
   t -> SetBranchAddress("MCP_Pz", &MCP_Pz);
   t -> SetBranchStatus("MCP_E",1);
   t -> SetBranchAddress("MCP_E", &MCP_E);
   t -> SetBranchStatus("MCP_KE",1);
   t -> SetBranchAddress("MCP_KE", &MCP_KE);
   t -> SetBranchStatus("MCP_isContained",1);
   t -> SetBranchAddress("MCP_isContained", &MCP_isContained);

   t -> SetBranchStatus("nu_vtx",1);
   t -> SetBranchAddress("nu_vtx", &nu_vtx);
   t -> SetBranchStatus("nu_vtx_spacecharge",1);
   t -> SetBranchAddress("nu_vtx_spacecharge", &nu_vtx_spacecharge);
   t -> SetBranchStatus("nu_isCC",1);
   t -> SetBranchAddress("nu_isCC", &nu_isCC);
   t -> SetBranchStatus("nu_PDG",1);
   t -> SetBranchAddress("nu_PDG", &nu_PDG);
   t -> SetBranchStatus("nu_E",1);
   t -> SetBranchAddress("nu_E", &nu_E);

   t -> SetBranchStatus("TPCObj_PFP_n2LLH_fwd_mu",1);
   t -> SetBranchAddress("TPCObj_PFP_n2LLH_fwd_mu",&n2LLH_fwd_mu);
   t -> SetBranchStatus("TPCObj_PFP_n2LLH_bwd_mu",1);
   t -> SetBranchAddress("TPCObj_PFP_n2LLH_bwd_mu",&n2LLH_bwd_mu);
   t -> SetBranchStatus("TPCObj_PFP_n2LLH_fwd_p",1);
   t -> SetBranchAddress("TPCObj_PFP_n2LLH_fwd_p",&n2LLH_fwd_p);
   t -> SetBranchStatus("TPCObj_PFP_n2LLH_bwd_p",1);
   t -> SetBranchAddress("TPCObj_PFP_n2LLH_bwd_p",&n2LLH_bwd_p);
   t -> SetBranchStatus("TPCObj_PFP_n2LLH_MIP",1);
   t -> SetBranchAddress("TPCObj_PFP_n2LLH_MIP",&n2LLH_MIP);

   t -> SetBranchStatus("TPCObj_PFP_track_depE",1);
   t -> SetBranchAddress("TPCObj_PFP_track_depE",&TPCObj_PFP_track_depE);
   t -> SetBranchStatus("TPCObj_PFP_track_rangeE_mu",1);
   t -> SetBranchAddress("TPCObj_PFP_track_rangeE_mu",&TPCObj_PFP_track_rangeE_mu);
   t -> SetBranchStatus("TPCObj_PFP_track_rangeE_p",1);
   t -> SetBranchAddress("TPCObj_PFP_track_rangeE_p",&TPCObj_PFP_track_rangeE_p);


   TEfficiency* eff_nuE = new TEfficiency("eff_nuE",";True Neutrino Energy [GeV];",15,0,3);
   TEfficiency* eff_muP = new TEfficiency("eff_muP",";True Muon Momentum [GeV];",15,0,2);
   TEfficiency* eff_muCosT = new TEfficiency("eff_muCosT",";True Muon cos(#theta);",10,-1,1);
   TEfficiency* eff_muPhi = new TEfficiency("eff_muPhi",";True Muon #phi angle;",10,-3.2,3.2);

   TEfficiency* pur_nuE = new TEfficiency("pur_nuE",";True Neutrino Energy [GeV];",15,0,3);
   TEfficiency* pur_muP = new TEfficiency("pur_muP",";True Muon Momentum [GeV];",15,0,2);
   TEfficiency* pur_muCosT = new TEfficiency("pur_muCosT",";True Muon cos(#theta);",10,-1,1);
   TEfficiency* pur_muPhi = new TEfficiency("pur_muPhi",";True Muon #phi angle;",10,-3.2,3.2);

   StackedHistTopology* len_stack = new StackedHistTopology("len_stack", ";Longest Track Length [cm];Selected Eevents(normalised to 3.782 #times 10^{19} POT)", 30, 0, 700);
   StackedHistTopology* mips_stack = new StackedHistTopology("mips_stack", ";Number of MIP-like #nu daughter tracks in selected TPCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)", 10, 0, 10);
   StackedHistTopology* vtxtrack_stack = new StackedHistTopology("vtxtrack_stack", ";Distance between vertex and start of #nu daughter tracks in selected TPCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)", 20, 5, 25);
   StackedHistTopology* pfps_stack = new StackedHistTopology("pfps_stack", ";Number of #nu daughter PFPs in selected TPCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)", 10, 0, 10);

   StackedHistPDGCode* isMIP_stack_byPDG = new StackedHistPDGCode("isMIP_stack_byPDG", ";#nu daughter track is MIP-like?;Selected Events (normalised to 3.782 #times 10^{19} POT)", 2, 0, 2);
   StackedHistPDGCode* mips_stack_byPDG = new StackedHistPDGCode("mips_stack_byPDG", ";Number of MIP-like #nu daughter tracks in selected TPCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)", 10, 0, 10);
   StackedHistPDGCode* pfps_stack_byPDG = new StackedHistPDGCode("pfps_stack_byPDG", ";Number of #nu daughter PFPs in selected TPCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)", 10, 0, 10);
   StackedHistPDGCode* n2llh_muMIPminusp_stack_byPDG = new StackedHistPDGCode("n2llh_muMIPminusp_stack_byPDG", ";n2llh_muMIPminusp for #nu daughter PFPs in selected TPCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)", 54, -20, 7);
   StackedHistPDGCode* n2llh_p_stack_byPDG = new StackedHistPDGCode("n2llh_p_stack_byPDG", ";n2llh_p for #nu daughter PFPs in selected TPCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)", 54, 0, 27);
   StackedHistPDGCode* depE_rangeE_mu_stack_byPDG = new StackedHistPDGCode("depE_rangeE_mu_stack_byPDG", ";Deposited energy - energy by range (muon assumption) [MeV];Selected Events (normalised to 3.782 #times 10^{19} POT)", 50, -100, 100);
   StackedHistPDGCode* depE_rangeE_p_stack_byPDG = new StackedHistPDGCode("depE_rangeE_p_stack_byPDG", ";Deposited energy - energy by range (proton assumption) [MeV];Selected Events (normalised to 3.782 #times 10^{19} POT)", 50, -100, 100);

   // Stacked plots for longest/second-longest/third-longest etc. track
   // For now, only consider up to 10 longest tracks but that can be easily changed
   StackedHistPDGCode* len_stack_byPDG[10];
   for (int j=0; j<10; j++){
      std::string histname = std::string("len_stack_byPDG_track")+std::to_string(j);
      if (j==0){
         len_stack_byPDG[j] = new StackedHistPDGCode(histname.c_str(), ";Longest Track Length [cm];Selected Events (normalised to 3.782 #times 10^{19} POT)", 30, 0, 700);
      }
      else{
         len_stack_byPDG[j] = new StackedHistPDGCode(histname.c_str(), ";Longest Track Length [cm];Selected Events (normalised to 3.782 #times 10^{19} POT)", 30, 0, 300);
      }
   }

   StackedHistPDGCode* len_stack_MIPs_byPDG[5];
   for (int j=0; j<5; j++){
      std::string histname = std::string("len_stack_MIPs_byPDG_track")+std::to_string(j);
      if (j==0){
         len_stack_MIPs_byPDG[j] = new StackedHistPDGCode(histname.c_str(), ";Longest Track Length [cm];Selected Events (normalised to 3.782 #times 10^{19} POT)", 30, 0, 700);
      }
      else{
         len_stack_MIPs_byPDG[j] = new StackedHistPDGCode(histname.c_str(), ";Longest Track Length [cm];Selected Events (normalised to 3.782 #times 10^{19} POT)", 30, 0, 300);
      }
   }

   StackedHistPDGCode* mom_stack_byPDG[10];
   for (int j=0; j<10; j++){
      std::string histname = std::string("mom_stack_byPDG_track")+std::to_string(j);
      if (j==0){
         mom_stack_byPDG[j] = new StackedHistPDGCode(histname.c_str(), ";Longest Track Length [cm];Selected Events (normalised to 3.782 #times 10^{19} POT)", 30, 0, 700);
      }
      else{
         mom_stack_byPDG[j] = new StackedHistPDGCode(histname.c_str(), ";Longest Track Length [cm];Selected Events (normalised to 3.782 #times 10^{19} POT)", 30, 0, 300);
      }
   }

   StackedHistPDGCode* mom_stack_MIPs_byPDG[5];
   for (int j=0; j<5; j++){
      std::string histname = std::string("mom_stack_MIPs_byPDG_track")+std::to_string(j);
      if (j==0){
         mom_stack_MIPs_byPDG[j] = new StackedHistPDGCode(histname.c_str(), ";Longest Track Length [cm];Selected Events (normalised to 3.782 #times 10^{19} POT)", 30, 0, 700);
      }
      else{
         mom_stack_MIPs_byPDG[j] = new StackedHistPDGCode(histname.c_str(), ";Longest Track Length [cm];Selected Events (normalised to 3.782 #times 10^{19} POT)", 30, 0, 300);
      }
   }

   TH2D* longest2tracks_byPDG_signal = new TH2D("longest2tracks_byPDG_signal",";True PDG of longest track;True PDG of second longest track",4,0,4,4,0,4);
   TH2D* longest2tracks_byPDG_bg = new TH2D("longest2tracks_byPDG_bg",";True PDG of longest track;True PDG of second longest track",4,0,4,4,0,4);

//   TH1D* vtx_truereco_dist = new TH1D("vtx_truereco_dist",";Distance between true and reco vertices;Selected Events (normalised to 3.782 #times 10^{19} POT)",50,0,50);
   TH1D* vtx_truereco_spacecharge_dist = new TH1D("vtx_truereco_spacecharge_dist",";Distance between true and reco vertices (including spacecharge);Selected Events (normalised to 3.782 #times 10^{19} POT)",50,0,300);

   double MIPcuts_muMIPminusp[] = {-10, -9, -8, -7, -6, -5, -4, -3, -2.5, -2, -1.5, -1, 0, 1, 2, 3, 4, 5, 6, 7};
   int MIPcuts_muMIPminusp_size = sizeof(MIPcuts_muMIPminusp)/sizeof(MIPcuts_muMIPminusp[0]) - 1;
   TEfficiency* eff_n2llh_muMIPminusp_atleast2MIP = new TEfficiency("eff_n2llh_muMIPminusp_atleast2MIP","Require at least 2 MIPs;Value of n2llh_muMIPminusp used for MIP classification;",MIPcuts_muMIPminusp_size,MIPcuts_muMIPminusp);
   TEfficiency* eff_n2llh_muMIPminusp_exactly2MIP = new TEfficiency("eff_n2llh_muMIPminusp_exactly2MIP","Require exactly 2 MIPs;Value of n2llh_muMIPminusp used for MIP classification;",MIPcuts_muMIPminusp_size,MIPcuts_muMIPminusp);
   TEfficiency* pur_n2llh_muMIPminusp_atleast2MIP = new TEfficiency("pur_n2llh_muMIPminusp_atleast2MIP","Require at least 2 MIPs;Value of n2llh_muMIPminusp used for MIP classification;",MIPcuts_muMIPminusp_size,MIPcuts_muMIPminusp);
   TEfficiency* pur_n2llh_muMIPminusp_exactly2MIP = new TEfficiency("pur_n2llh_muMIPminusp_exactly2MIP","Require exactly 2 MIPs;Value of n2llh_muMIPminusp used for MIP classification;",MIPcuts_muMIPminusp_size,MIPcuts_muMIPminusp);
   TH1D* effpur_n2llh_muMIPminusp_atleast2MIP = new TH1D("effpur_n2llh_muMIPminusp_atleast2MIP",";;",MIPcuts_muMIPminusp_size,MIPcuts_muMIPminusp);
   TH1D* effpur_n2llh_muMIPminusp_exactly2MIP = new TH1D("effpur_n2llh_muMIPminusp_exactly2MIP",";;",MIPcuts_muMIPminusp_size,MIPcuts_muMIPminusp);

   double MIPcuts_p[] = {0, 1, 2, 3, 4, 5, 5.5, 6, 6.5, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17};
   int MIPcuts_p_size = sizeof(MIPcuts_p)/sizeof(MIPcuts_p[0]) - 1;
   TEfficiency* eff_n2llh_p_atleast2MIP = new TEfficiency("eff_n2llh_p_atleast2MIP","Require at least 2 MIPs;Value of n2llh_p used for MIP classification;",MIPcuts_p_size,MIPcuts_p);
   TEfficiency* eff_n2llh_p_exactly2MIP = new TEfficiency("eff_n2llh_p_exactly2MIP","Require exactly 2 MIPs;Value of n2llh_p used for MIP classification;",MIPcuts_p_size,MIPcuts_p);
   TEfficiency* pur_n2llh_p_atleast2MIP = new TEfficiency("pur_n2llh_p_atleast2MIP","Require at least 2 MIPs;Value of n2llh_p used for MIP classification;",MIPcuts_p_size,MIPcuts_p);
   TEfficiency* pur_n2llh_p_exactly2MIP = new TEfficiency("pur_n2llh_p_exactly2MIP","Require exactly 2 MIPs;Value of n2llh_p used for MIP classification;",MIPcuts_p_size,MIPcuts_p);
   TH1D* effpur_n2llh_p_atleast2MIP = new TH1D("effpur_n2llh_p_atleast2MIP",";;",MIPcuts_p_size,MIPcuts_p);
   TH1D* effpur_n2llh_p_exactly2MIP = new TH1D("effpur_n2llh_p_exactly2MIP",";;",MIPcuts_p_size,MIPcuts_p);

   double MIPcuts_depE_rangeE_mu[] = {-20, -10, 0, 10, 20, 30, 40, 50, 60, 70, 80};
   int MIPcuts_depE_rangeE_mu_size = sizeof(MIPcuts_depE_rangeE_mu)/sizeof(MIPcuts_depE_rangeE_mu[0]) - 1;
   TEfficiency* eff_depE_rangeE_mu_atleast2MIP = new TEfficiency("eff_depE_rangeE_mu_atleast2MIP","Require at least 2 MIPs;Value of depE_rangeE_mu used for MIP classification;",MIPcuts_depE_rangeE_mu_size,MIPcuts_depE_rangeE_mu);
   TEfficiency* eff_depE_rangeE_mu_exactly2MIP = new TEfficiency("eff_depE_rangeE_mu_exactly2MIP","Require exactly 2 MIPs;Value of depE_rangeE_mu used for MIP classification;",MIPcuts_depE_rangeE_mu_size,MIPcuts_depE_rangeE_mu);
   TEfficiency* pur_depE_rangeE_mu_atleast2MIP = new TEfficiency("pur_depE_rangeE_mu_atleast2MIP","Require at least 2 MIPs;Value of depE_rangeE_mu used for MIP classification;",MIPcuts_depE_rangeE_mu_size,MIPcuts_depE_rangeE_mu);
   TEfficiency* pur_depE_rangeE_mu_exactly2MIP = new TEfficiency("pur_depE_rangeE_mu_exactly2MIP","Require exactly 2 MIPs;Value of depE_rangeE_mu used for MIP classification;",MIPcuts_depE_rangeE_mu_size,MIPcuts_depE_rangeE_mu);
   TH1D* effpur_depE_rangeE_mu_atleast2MIP = new TH1D("effpur_depE_rangeE_mu_atleast2MIP",";;",MIPcuts_depE_rangeE_mu_size,MIPcuts_depE_rangeE_mu);
   TH1D* effpur_depE_rangeE_mu_exactly2MIP = new TH1D("effpur_depE_rangeE_mu_exactly2MIP",";;",MIPcuts_depE_rangeE_mu_size,MIPcuts_depE_rangeE_mu);

   double MIPcuts_depE_rangeE_p[] = {-100, -90, -80, -70, -65, -60, -55, -50, -40, -30, -20, -10, 0};
   int MIPcuts_depE_rangeE_p_size = sizeof(MIPcuts_depE_rangeE_p)/sizeof(MIPcuts_depE_rangeE_p[0]) - 1;
   TEfficiency* eff_depE_rangeE_p_atleast2MIP = new TEfficiency("eff_depE_rangeE_p_atleast2MIP","Require at least 2 MIPs;Value of depE_rangeE_p used for MIP classification;",MIPcuts_depE_rangeE_p_size,MIPcuts_depE_rangeE_p);
   TEfficiency* eff_depE_rangeE_p_exactly2MIP = new TEfficiency("eff_depE_rangeE_p_exactly2MIP","Require exactly 2 MIPs;Value of depE_rangeE_p used for MIP classification;",MIPcuts_depE_rangeE_p_size,MIPcuts_depE_rangeE_p);
   TEfficiency* pur_depE_rangeE_p_atleast2MIP = new TEfficiency("pur_depE_rangeE_p_atleast2MIP","Require at least 2 MIPs;Value of depE_rangeE_p used for MIP classification;",MIPcuts_depE_rangeE_p_size,MIPcuts_depE_rangeE_p);
   TEfficiency* pur_depE_rangeE_p_exactly2MIP = new TEfficiency("pur_depE_rangeE_p_exactly2MIP","Require exactly 2 MIPs;Value of depE_rangeE_p used for MIP classification;",MIPcuts_depE_rangeE_p_size,MIPcuts_depE_rangeE_p);
   TH1D* effpur_depE_rangeE_p_atleast2MIP = new TH1D("effpur_depE_rangeE_p_atleast2MIP",";;",MIPcuts_depE_rangeE_p_size,MIPcuts_depE_rangeE_p);
   TH1D* effpur_depE_rangeE_p_exactly2MIP = new TH1D("effpur_depE_rangeE_p_exactly2MIP",";;",MIPcuts_depE_rangeE_p_size,MIPcuts_depE_rangeE_p);

   TEfficiency* eff_n2llh_muMIPminusp_depE_rangeE_p_atleast2MIP = new TEfficiency("eff_n2llh_muMIPminusp_depE_rangeE_p_atleast2MIP","Require at least 2 MIPs;Value of n2llh_muMIPminusp used for MIP classification;Value of depE_rangeE_p used for MIP classification",MIPcuts_muMIPminusp_size,MIPcuts_muMIPminusp,MIPcuts_depE_rangeE_p_size,MIPcuts_depE_rangeE_p);
   TEfficiency* eff_n2llh_muMIPminusp_depE_rangeE_p_exactly2MIP = new TEfficiency("eff_n2llh_muMIPminusp_depE_rangeE_p_exactly2MIP","Require exactly 2 MIPs;Value of n2llh_muMIPminusp used for MIP classification;Value of depE_rangeE_p used for MIP classification",MIPcuts_muMIPminusp_size,MIPcuts_muMIPminusp,MIPcuts_depE_rangeE_p_size,MIPcuts_depE_rangeE_p);
   TEfficiency* pur_n2llh_muMIPminusp_depE_rangeE_p_atleast2MIP = new TEfficiency("pur_n2llh_muMIPminusp_depE_rangeE_p_atleast2MIP","Require at least 2 MIPs;Value of n2llh_muMIPminusp used for MIP classification;Value of depE_rangeE_p used for MIP classification",MIPcuts_muMIPminusp_size,MIPcuts_muMIPminusp,MIPcuts_depE_rangeE_p_size,MIPcuts_depE_rangeE_p);
   TEfficiency* pur_n2llh_muMIPminusp_depE_rangeE_p_exactly2MIP = new TEfficiency("pur_n2llh_muMIPminusp_depE_rangeE_p_exactly2MIP","Require exactly 2 MIPs;Value of n2llh_muMIPminusp used for MIP classification;Value of depE_rangeE_p used for MIP classification",MIPcuts_muMIPminusp_size,MIPcuts_muMIPminusp,MIPcuts_depE_rangeE_p_size,MIPcuts_depE_rangeE_p);
   TH2D* effpur_n2llh_muMIPminusp_depE_rangeE_p_atleast2MIP = new TH2D("effpur_n2llh_muMIPminusp_depE_rangeE_p_atleast2MIP","Require at least 2 MIPs;Value of n2llh_muMIPminusp used for MIP classification;Value of depE_rangeE_p used for MIP classification",MIPcuts_muMIPminusp_size,MIPcuts_muMIPminusp,MIPcuts_depE_rangeE_p_size,MIPcuts_depE_rangeE_p);
   TH2D* effpur_n2llh_muMIPminusp_depE_rangeE_p_exactly2MIP = new TH2D("effpur_n2llh_muMIPminusp_depE_rangeE_p_exactly2MIP","Require exactly 2 MIPs;Value of n2llh_muMIPminusp used for MIP classification;Value of depE_rangeE_p used for MIP classification",MIPcuts_muMIPminusp_size,MIPcuts_muMIPminusp,MIPcuts_depE_rangeE_p_size,MIPcuts_depE_rangeE_p);

   double distcuts[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15};
   int distcuts_size = sizeof(distcuts)/sizeof(distcuts[0]) - 1;
   TEfficiency* eff_distcut = new TEfficiency("eff_distcut",";Value used for distance cut;",distcuts_size,distcuts);
   TEfficiency* pur_distcut = new TEfficiency("pur_distcut",";Value used for distance cut;",distcuts_size,distcuts);
   TH1D* effpur_distcut = new TH1D("effpur_distcut",";;",distcuts_size,distcuts);

   TH1D* AngleBetweenTracks = new TH1D("AngleBetweenTracks",";Angle between selected tracks;Entries",30,0,3.14);
   StackedHistTopology* AngleBetweenTracks_stack = new StackedHistTopology("AngleBetweenTracks_stack",";Angle between selected tracks;Entries",30,0,3.14);
   TH1D* AngleBetweenTracks_Close = new TH1D("AngleBetweenTracks_Close",";Angle between selected tracks that are within 3 cm;Entries",30,0,3.14);
   StackedHistTopology* AngleBetweenTracks_Close_stack = new StackedHistTopology("AngleBetweenTracks_Close_stack",";Angle between selected tracks that are within 3 cm;Entries",30,0,3.14);
   StackedHistTopology* AngleBetweenTracks_Close_sameID_stack = new StackedHistTopology("AngleBetweenTracks_Close_sameID_stack",";Angle between selected tracks that are within 3 cm and match the same MCP;Entries",30,0,3.14);
   StackedHistTopology* AngleBetweenTracks_Close_diffID_stack = new StackedHistTopology("AngleBetweenTracks_Close_diffID_stack",";Angle between selected tracks that are within 3 cm and match different MCPs;Entries",30,0,3.14);

   //ofstream evdinfo;
   //evdinfo.open("evdinfo.txt");
   //evdinfo << "run_num subrun_num event_num track_daughters Truth_topology" << std::endl;

   const int nentries = t -> GetEntries();

   for (int i = 0; i < nentries; i++) {
      t -> GetEntry(i);

      // Check whether the event passes the requested cutflow map
      bool SelectedEvent = true;
      for(std::map<std::string,bool>::const_iterator iter = SelectionCutflow.begin(); iter != SelectionCutflow.end(); ++iter) {
         if((CC1picutflow -> find(iter -> first)) -> second != iter -> second) {
            SelectedEvent = false;
            break;
         }
      }


      bool isSignal = false;

      // if signal...
      if (Truth_topology == kCC1piplus0p || Truth_topology == kCC1piplus1p || Truth_topology == kCC1piplusNp) {

         isSignal = true;

         // Check neutrino flavour (just in case)
         if (nu_PDG != 14) std::cout << "WARNING: Signal event has neutrino with PDG code " << nu_PDG << " (should be 14)" << std::endl;

      }

      // Test vertex to track start distance cut...
      double maxdist = 0;
      for(int j = 0; j < TPCObj_NPFPs; j++) {
         if(!(TPCObj_PFP_isDaughter -> at(j))) continue;

         double vtxtrackdist = std::hypot(std::hypot(TPCObj_reco_vtx->at(0) - TPCObj_PFP_track_start->at(j).at(0), TPCObj_reco_vtx->at(1) - TPCObj_PFP_track_start->at(j).at(1)), TPCObj_reco_vtx->at(2) - TPCObj_PFP_track_start->at(j).at(2));
         if(vtxtrackdist > maxdist) maxdist = vtxtrackdist;
/*
         for(int k = j+1; k < TPCObj_NPFPs; k++) {
            if(!(TPCObj_PFP_isDaughter -> at(k))) continue;

            double tracktrackdist = std::hypot(std::hypot(TPCObj_PFP_track_start->at(k).at(0) - TPCObj_PFP_track_start->at(j).at(0), TPCObj_PFP_track_start->at(k).at(1) - TPCObj_PFP_track_start->at(j).at(1)), TPCObj_PFP_track_start->at(k).at(2) - TPCObj_PFP_track_start->at(j).at(2));
            if(tracktrackdist > maxdist) maxdist = tracktrackdist;
         }
*/
      }

      // distance cut eff/pur
      for (const double& distcut : distcuts) {
         if(isSignal) eff_distcut->Fill(SelectedEvent && maxdist < distcut, distcut);
         if(SelectedEvent && maxdist < distcut) pur_distcut->Fill(isSignal, distcut);
      }

      for(int j = 1; j <= distcuts_size; j++) {
         effpur_distcut->SetBinContent(j,10*eff_distcut->GetEfficiency(j)*pur_distcut->GetEfficiency(j));
      }


      // Calculate derived MIP vars
      std::vector<double> n2llh_muMIPminusp_vect;
      std::vector<double> n2llh_p_vect;
      std::vector<double> depE_rangeE_mu_vect;
      std::vector<double> depE_rangeE_p_vect;
      for(int j = 0; j < TPCObj_NPFPs; j++) {
         if(TPCObj_PFP_isDaughter -> at(j) && n2LLH_fwd_mu -> at(j).at(2) != -999 && TPCObj_PFP_track_length -> at(j) >= 0) {
            double tr_n2llh_muMIP = std::min({n2LLH_fwd_mu->at(j).at(2),n2LLH_bwd_mu->at(j).at(2),n2LLH_MIP->at(j).at(2)});
            double tr_n2llh_p = std::min(n2LLH_fwd_p->at(j).at(2),n2LLH_bwd_p->at(j).at(2));
            double tr_n2llh_muMIPminusp = tr_n2llh_muMIP - tr_n2llh_p;

            n2llh_muMIPminusp_vect.emplace_back(tr_n2llh_muMIPminusp);
            n2llh_p_vect.emplace_back(tr_n2llh_p);

            depE_rangeE_mu_vect.emplace_back(TPCObj_PFP_track_depE->at(j).at(2) - TPCObj_PFP_track_rangeE_mu->at(j));
            depE_rangeE_p_vect.emplace_back(TPCObj_PFP_track_depE->at(j).at(2) - TPCObj_PFP_track_rangeE_p->at(j));
         }
         else {
            n2llh_muMIPminusp_vect.emplace_back(-999);
            n2llh_p_vect.emplace_back(-999);
            depE_rangeE_mu_vect.emplace_back(-999);
            depE_rangeE_p_vect.emplace_back(-999);
        }
      }

      // n2llh_muMIPminusp eff/pur
      for (const double& MIPcut : MIPcuts_muMIPminusp) {
         std::vector<bool> TPCObj_PFP_isMIP_LLH;
         for(const double& MIPvalue : n2llh_muMIPminusp_vect) {
            if(MIPvalue < MIPcut && MIPvalue != -999) TPCObj_PFP_isMIP_LLH.emplace_back(true);
            else TPCObj_PFP_isMIP_LLH.emplace_back(false);
         }
         int nMIPs = std::count(TPCObj_PFP_isMIP_LLH.begin(),TPCObj_PFP_isMIP_LLH.end(),true);
         if(isSignal) eff_n2llh_muMIPminusp_atleast2MIP->Fill(SelectedEvent && nMIPs>=2, MIPcut);
         if(isSignal) eff_n2llh_muMIPminusp_exactly2MIP->Fill(SelectedEvent && nMIPs==2, MIPcut);
         if(SelectedEvent && nMIPs>=2) pur_n2llh_muMIPminusp_atleast2MIP->Fill(isSignal, MIPcut);
         if(SelectedEvent && nMIPs==2) pur_n2llh_muMIPminusp_exactly2MIP->Fill(isSignal, MIPcut);
      }

      for(int j = 1; j <= MIPcuts_muMIPminusp_size; j++) {
         effpur_n2llh_muMIPminusp_atleast2MIP->SetBinContent(j,10*eff_n2llh_muMIPminusp_atleast2MIP->GetEfficiency(j)*pur_n2llh_muMIPminusp_atleast2MIP->GetEfficiency(j));
         effpur_n2llh_muMIPminusp_exactly2MIP->SetBinContent(j,10*eff_n2llh_muMIPminusp_exactly2MIP->GetEfficiency(j)*pur_n2llh_muMIPminusp_exactly2MIP->GetEfficiency(j));
      }

      // n2llh_p eff/pur
      for (const double& MIPcut : MIPcuts_p) {
         std::vector<bool> TPCObj_PFP_isMIP_LLH;
         for(const double& MIPvalue : n2llh_p_vect) {
            if(MIPvalue > MIPcut && MIPvalue != -999) TPCObj_PFP_isMIP_LLH.emplace_back(true);
            else TPCObj_PFP_isMIP_LLH.emplace_back(false);
         }
         int nMIPs = std::count(TPCObj_PFP_isMIP_LLH.begin(),TPCObj_PFP_isMIP_LLH.end(),true);
         if(isSignal) eff_n2llh_p_atleast2MIP->Fill(SelectedEvent && nMIPs>=2, MIPcut);
         if(isSignal) eff_n2llh_p_exactly2MIP->Fill(SelectedEvent && nMIPs==2, MIPcut);
         if(SelectedEvent && nMIPs>=2) pur_n2llh_p_atleast2MIP->Fill(isSignal, MIPcut);
         if(SelectedEvent && nMIPs==2) pur_n2llh_p_exactly2MIP->Fill(isSignal, MIPcut);
      }

      for(int j = 1; j <= MIPcuts_p_size; j++) {
         effpur_n2llh_p_atleast2MIP->SetBinContent(j,10*eff_n2llh_p_atleast2MIP->GetEfficiency(j)*pur_n2llh_p_atleast2MIP->GetEfficiency(j));
         effpur_n2llh_p_exactly2MIP->SetBinContent(j,10*eff_n2llh_p_exactly2MIP->GetEfficiency(j)*pur_n2llh_p_exactly2MIP->GetEfficiency(j));
      }

      // depE_rangeE_mu eff/pur
      for (const double& MIPcut : MIPcuts_depE_rangeE_mu) {
         std::vector<bool> TPCObj_PFP_isMIP_E;
         for(const double& MIPvalue : depE_rangeE_mu_vect) {
            if(MIPvalue < MIPcut && MIPvalue != -999) TPCObj_PFP_isMIP_E.emplace_back(true);
            else TPCObj_PFP_isMIP_E.emplace_back(false);
         }
         int nMIPs = std::count(TPCObj_PFP_isMIP_E.begin(),TPCObj_PFP_isMIP_E.end(),true);
         if(isSignal) eff_depE_rangeE_mu_atleast2MIP->Fill(SelectedEvent && nMIPs>=2, MIPcut);
         if(isSignal) eff_depE_rangeE_mu_exactly2MIP->Fill(SelectedEvent && nMIPs==2, MIPcut);
         if(SelectedEvent && nMIPs>=2) pur_depE_rangeE_mu_atleast2MIP->Fill(isSignal, MIPcut);
         if(SelectedEvent && nMIPs==2) pur_depE_rangeE_mu_exactly2MIP->Fill(isSignal, MIPcut);
      }

      for(int j = 1; j <= MIPcuts_depE_rangeE_mu_size; j++) {
         effpur_depE_rangeE_mu_atleast2MIP->SetBinContent(j,10*eff_depE_rangeE_mu_atleast2MIP->GetEfficiency(j)*pur_depE_rangeE_mu_atleast2MIP->GetEfficiency(j));
         effpur_depE_rangeE_mu_exactly2MIP->SetBinContent(j,10*eff_depE_rangeE_mu_exactly2MIP->GetEfficiency(j)*pur_depE_rangeE_mu_exactly2MIP->GetEfficiency(j));
      }

      // depE_rangeE_p eff/pur
      for (const double& MIPcut : MIPcuts_depE_rangeE_p) {
         std::vector<bool> TPCObj_PFP_isMIP_E;
         for(const double& MIPvalue : depE_rangeE_p_vect) {
            if(MIPvalue < MIPcut && MIPvalue != -999) TPCObj_PFP_isMIP_E.emplace_back(true);
            else TPCObj_PFP_isMIP_E.emplace_back(false);
         }
         int nMIPs = std::count(TPCObj_PFP_isMIP_E.begin(),TPCObj_PFP_isMIP_E.end(),true);
         if(isSignal) eff_depE_rangeE_p_atleast2MIP->Fill(SelectedEvent && nMIPs>=2, MIPcut);
         if(isSignal) eff_depE_rangeE_p_exactly2MIP->Fill(SelectedEvent && nMIPs==2, MIPcut);
         if(SelectedEvent && nMIPs>=2) pur_depE_rangeE_p_atleast2MIP->Fill(isSignal, MIPcut);
         if(SelectedEvent && nMIPs==2) pur_depE_rangeE_p_exactly2MIP->Fill(isSignal, MIPcut);
      }

      for(int j = 1; j <= MIPcuts_depE_rangeE_p_size; j++) {
         effpur_depE_rangeE_p_atleast2MIP->SetBinContent(j,10*eff_depE_rangeE_p_atleast2MIP->GetEfficiency(j)*pur_depE_rangeE_p_atleast2MIP->GetEfficiency(j));
         effpur_depE_rangeE_p_exactly2MIP->SetBinContent(j,10*eff_depE_rangeE_p_exactly2MIP->GetEfficiency(j)*pur_depE_rangeE_p_exactly2MIP->GetEfficiency(j));
      }

      // 2D eff/pur
      for (const double& MIPcut_muMIPminusp : MIPcuts_muMIPminusp) {
         for (const double& MIPcut_depE_rangeE_p : MIPcuts_depE_rangeE_p) {
            std::vector<bool> TPCObj_PFP_isMIP;
            for(int j = 0; j < n2llh_muMIPminusp_vect.size(); j++) {
               if(n2llh_muMIPminusp_vect.at(j) < MIPcut_muMIPminusp && depE_rangeE_p_vect.at(j) < MIPcut_depE_rangeE_p && n2llh_muMIPminusp_vect.at(j) != -999) TPCObj_PFP_isMIP.emplace_back(true);
               else TPCObj_PFP_isMIP.emplace_back(false);
            }
            int nMIPs = std::count(TPCObj_PFP_isMIP.begin(),TPCObj_PFP_isMIP.end(),true);
            if(isSignal) eff_n2llh_muMIPminusp_depE_rangeE_p_atleast2MIP->Fill(SelectedEvent && nMIPs>=2, MIPcut_muMIPminusp, MIPcut_depE_rangeE_p);
            if(isSignal) eff_n2llh_muMIPminusp_depE_rangeE_p_exactly2MIP->Fill(SelectedEvent && nMIPs==2, MIPcut_muMIPminusp, MIPcut_depE_rangeE_p);
            if(SelectedEvent && nMIPs>=2) pur_n2llh_muMIPminusp_depE_rangeE_p_atleast2MIP->Fill(isSignal, MIPcut_muMIPminusp, MIPcut_depE_rangeE_p);
            if(SelectedEvent && nMIPs==2) pur_n2llh_muMIPminusp_depE_rangeE_p_exactly2MIP->Fill(isSignal, MIPcut_muMIPminusp, MIPcut_depE_rangeE_p);
         }
      }

      for(int j = 1; j <= MIPcuts_muMIPminusp_size; j++) {
         for(int k = 1; k <= MIPcuts_depE_rangeE_p_size; k++) {
            int bin = eff_n2llh_muMIPminusp_depE_rangeE_p_atleast2MIP->GetGlobalBin(j,k);
            effpur_n2llh_muMIPminusp_depE_rangeE_p_atleast2MIP->SetBinContent(bin,10*eff_n2llh_muMIPminusp_depE_rangeE_p_atleast2MIP->GetEfficiency(bin)*pur_n2llh_muMIPminusp_depE_rangeE_p_atleast2MIP->GetEfficiency(bin));
            effpur_n2llh_muMIPminusp_depE_rangeE_p_exactly2MIP->SetBinContent(bin,10*eff_n2llh_muMIPminusp_depE_rangeE_p_exactly2MIP->GetEfficiency(bin)*pur_n2llh_muMIPminusp_depE_rangeE_p_exactly2MIP->GetEfficiency(bin));
         }
      }

      // Implement cut on n2llh_muMIPminusp
      double MIPcut_muMIPminusp = -2.5;
      std::vector<bool> TPCObj_PFP_isMIP_LLH;
      for(const double& MIPvalue : n2llh_muMIPminusp_vect) {
         if(MIPvalue < MIPcut_muMIPminusp && MIPvalue != -999) TPCObj_PFP_isMIP_LLH.emplace_back(true);
         else TPCObj_PFP_isMIP_LLH.emplace_back(false);
      }
//      if(std::count(TPCObj_PFP_isMIP_LLH.begin(),TPCObj_PFP_isMIP_LLH.end(),true)!=2) SelectedEvent = false;

      // Implement cut on depE_rangeE_p
      double MIPcut_depE_rangeE_p = -60;
      std::vector<bool> TPCObj_PFP_isMIP_E;
      for(const double& MIPvalue : depE_rangeE_p_vect) {
         if(MIPvalue < MIPcut_depE_rangeE_p && MIPvalue != -999) TPCObj_PFP_isMIP_E.emplace_back(true);
         else TPCObj_PFP_isMIP_E.emplace_back(false);
      }
//      if(std::count(TPCObj_PFP_isMIP_E.begin(),TPCObj_PFP_isMIP_E.end(),true)<2) SelectedEvent = false;


      // Efficiency...
      if(isSignal) {
         eff_nuE -> Fill(SelectedEvent, nu_E);

         for (size_t j = 0; j < MCP_PDG -> size(); j++) {
            if (MCP_PDG -> at(j) == 13) {
               eff_muP -> Fill(SelectedEvent, MCP_P -> at(j));
               TVector3 muP(MCP_Px -> at(j), MCP_Py -> at(j), MCP_Pz -> at(j));
               eff_muCosT -> Fill(SelectedEvent, muP.CosTheta());
               eff_muPhi -> Fill(SelectedEvent, muP.Phi());
               break;
            }
         }
      } //end if signal



      if(SelectedEvent) {

         //if (nu_PDG != 14) std::cout << "WARNING: Selected event has neutrino with PDG code " << nu_PDG << " (should be 14)" << std::endl;

         // Purity...
         pur_nuE -> Fill(isSignal, nu_E);

         for (size_t j = 0; j < MCP_PDG -> size(); j++) {
            if (MCP_PDG -> at(j) == 13) {
               pur_muP -> Fill(isSignal, MCP_P -> at(j));
               TVector3 muP(MCP_Px -> at(j), MCP_Py -> at(j), MCP_Pz -> at(j));
               pur_muCosT -> Fill(isSignal, muP.CosTheta());
               pur_muPhi -> Fill(isSignal, muP.Phi());
               break;
            }
         }

         // THStacks...

         //Find order of tracks by length/momentum etc.
         std::vector<std::pair<double, int>> pair_trklen_index;
         std::vector<std::pair<double, int>> pair_MIP_trklen_index;
         std::vector<std::pair<double, int>> pair_trkmom_index;
         std::vector<std::pair<double, int>> pair_MIP_trkmom_index;
         for(int j = 0; j < TPCObj_NPFPs; j++) {
            pair_trklen_index.emplace_back(std::make_pair(TPCObj_PFP_track_length->at(j), j));
            pair_trklen_index.emplace_back(std::make_pair(TPCObj_PFP_track_mom->at(j), j));
            if(TPCObj_PFP_isMIP_LLH.at(j)){
               pair_MIP_trklen_index.emplace_back(std::make_pair(TPCObj_PFP_track_length->at(j), j));
               pair_MIP_trklen_index.emplace_back(std::make_pair(TPCObj_PFP_track_mom->at(j), j));
            }
            //std::pair<double,int> tmppairmom(TPCObj_PFP_)
         }

         // Now actually order vector of pairs by track length/momentum etc.
         std::sort(pair_trklen_index.rbegin(), pair_trklen_index.rend());
         std::sort(pair_MIP_trklen_index.rbegin(), pair_MIP_trklen_index.rend());
         std::sort(pair_trkmom_index.rbegin(), pair_trkmom_index.rend());
         std::sort(pair_MIP_trkmom_index.rbegin(), pair_MIP_trkmom_index.rend());

         double maxlen = 0;
         PDGCode maxPDG = kPDGUnknown;
         double maxlen2 = 0;
         PDGCode maxPDG2 = kPDGUnknown;
         if (TPCObj_NPFPs > 0){
            double maxlen = TPCObj_PFP_track_length -> at(pair_trklen_index.at(0).second);
            PDGCode maxPDG = (PDGCode)TPCObj_PFP_truePDG -> at(pair_trklen_index.at(0).second);
         }
         if (TPCObj_NPFPs > 1){
            maxlen2 = TPCObj_PFP_track_length -> at(pair_trklen_index.at(1).second);
            PDGCode maxPDG2 = (PDGCode)TPCObj_PFP_truePDG -> at(pair_trklen_index.at(1).second);
         }

         int track_daughters = 0;
         int mip_daughters = 0;
         std::vector<double> vtxtrack_daughters;

         for(int j = 0; j < TPCObj_NPFPs; j++) {
            if(TPCObj_PFP_isDaughter -> at(j)){
               if(TPCObj_PFP_track_length -> at(j) < 0) {
                  std::cout << "Daughter PFP has no track." << std::endl;
                  continue;
               }
               track_daughters++;
               vtxtrack_daughters.emplace_back(std::hypot(std::hypot(TPCObj_reco_vtx->at(0) - TPCObj_PFP_track_start->at(j).at(0), TPCObj_reco_vtx->at(1) - TPCObj_PFP_track_start->at(j).at(1)), TPCObj_reco_vtx->at(2) - TPCObj_PFP_track_start->at(j).at(2)));
               if(TPCObj_PFP_isMIP_LLH.at(j)) mip_daughters++;
            }
         }

         len_stack -> Fill(Truth_topology, maxlen);
         mips_stack -> Fill(Truth_topology, mip_daughters);
         pfps_stack -> Fill(Truth_topology, track_daughters);

         // Fill plots for longest/second-longest/third-longest etc. track
         // For now, only consider up to 10 longest tracks but that can be easily changed
         for (int j = 0; j < pair_trklen_index.size(); j++){
            if (j>=10) continue;
            int len_index = pair_trklen_index.at(j).second;
            len_stack_byPDG[j]->Fill((PDGCode)TPCObj_PFP_truePDG -> at(len_index), TPCObj_PFP_track_length -> at(len_index));
         }
         for (int j = 0; j < pair_MIP_trklen_index.size(); j++){
            if (j>=5) continue;
            int len_index_MIP = pair_MIP_trklen_index.at(j).second;
            len_stack_MIPs_byPDG[j]->Fill((PDGCode)TPCObj_PFP_truePDG -> at(len_index_MIP), TPCObj_PFP_track_length -> at(len_index_MIP));
         }
         for (int j = 0; j < pair_trkmom_index.size(); j++){
            if (j>=10) continue;
            int mom_index = pair_trkmom_index.at(j).second;
            mom_stack_byPDG[j]->Fill((PDGCode)TPCObj_PFP_truePDG -> at(mom_index), TPCObj_PFP_track_mom -> at(mom_index));
         }
         for (int j = 0; j < pair_MIP_trkmom_index.size(); j++){
            if (j>=5) continue;
            int mom_index_MIP = pair_MIP_trkmom_index.at(j).second;
            mom_stack_MIPs_byPDG[j]->Fill((PDGCode)TPCObj_PFP_truePDG -> at(mom_index_MIP), TPCObj_PFP_track_mom -> at(mom_index_MIP));
         }

         for (int j = 0; j < track_daughters; j++) {
            vtxtrack_stack -> Fill(Truth_topology, vtxtrack_daughters.at(j));
         }
         for(int j = 0; j < TPCObj_NPFPs; j++) {
            if(TPCObj_PFP_isDaughter -> at(j)) {
               if(TPCObj_PFP_track_length -> at(j) < 0) continue;
               isMIP_stack_byPDG -> Fill((PDGCode)TPCObj_PFP_truePDG->at(j), TPCObj_PFP_isMIP_LLH.at(j));
               pfps_stack_byPDG -> Fill((PDGCode)TPCObj_PFP_truePDG->at(j), track_daughters, 1./track_daughters);
               if(TPCObj_PFP_isMIP_LLH.at(j)) mips_stack_byPDG -> Fill((PDGCode)TPCObj_PFP_truePDG->at(j), mip_daughters, 1./mip_daughters);

               double tr_n2llh_muMIP = std::min({n2LLH_fwd_mu->at(j).at(2),n2LLH_bwd_mu->at(j).at(2),n2LLH_MIP->at(j).at(2)});
               double tr_n2llh_p = std::min(n2LLH_fwd_p->at(j).at(2),n2LLH_bwd_p->at(j).at(2));
               double tr_n2llh_muMIPminusp = tr_n2llh_muMIP - tr_n2llh_p;
               if (tr_n2llh_muMIP!=-999) {
                  n2llh_muMIPminusp_stack_byPDG -> Fill((PDGCode)TPCObj_PFP_truePDG->at(j), tr_n2llh_muMIPminusp);
                  n2llh_p_stack_byPDG -> Fill((PDGCode)TPCObj_PFP_truePDG->at(j), tr_n2llh_p);
               }

               depE_rangeE_mu_stack_byPDG -> Fill((PDGCode)TPCObj_PFP_truePDG->at(j), TPCObj_PFP_track_depE->at(j).at(2) - TPCObj_PFP_track_rangeE_mu->at(j));
               depE_rangeE_p_stack_byPDG -> Fill((PDGCode)TPCObj_PFP_truePDG->at(j), TPCObj_PFP_track_depE->at(j).at(2) - TPCObj_PFP_track_rangeE_p->at(j));

            }
         }

         //Longest 2 daughter tracks by PDG...
         if(track_daughters >= 2) {
            double maxPDG_fill, maxPDG2_fill;
            if(maxPDG == kMuMinus) maxPDG_fill = 0.;
            else if(maxPDG == kPiPlus) maxPDG_fill = 1;
            else if(maxPDG == kProton) maxPDG_fill = 2;
            else maxPDG_fill = 3;
            if(maxPDG2 == kMuMinus) maxPDG2_fill = 0.;
            else if(maxPDG2 == kPiPlus) maxPDG2_fill = 1;
            else if(maxPDG2 == kProton) maxPDG2_fill = 2;
            else maxPDG2_fill = 3;

            if(isSignal) longest2tracks_byPDG_signal -> Fill(maxPDG_fill, maxPDG2_fill);
            else longest2tracks_byPDG_bg -> Fill(maxPDG_fill, maxPDG2_fill);
         } //end longest 2 daughter tracks by PDG

         if (Truth_topology == kOutFV) {
//            vtx_truereco_dist -> Fill(std::hypot(std::hypot(TPCObj_reco_vtx->at(0) - nu_vtx->at(0), TPCObj_reco_vtx->at(1) - nu_vtx->at(1)), TPCObj_reco_vtx->at(2) - nu_vtx->at(2)));
            vtx_truereco_spacecharge_dist -> Fill(std::hypot(std::hypot(TPCObj_reco_vtx->at(0) - nu_vtx_spacecharge->at(0), TPCObj_reco_vtx->at(1) - nu_vtx_spacecharge->at(1)), TPCObj_reco_vtx->at(2) - nu_vtx_spacecharge->at(2)));
         }

/*         
         // Angle between tracks...
         for(size_t track1 = 0; track1 < TPCObj_PFP_track_AngleBetweenTracks->size(); track1++) {
            for(size_t track2 = track1+1; track2 < TPCObj_PFP_track_AngleBetweenTracks->size(); track2++) {
               AngleBetweenTracks -> Fill(TPCObj_PFP_track_AngleBetweenTracks->at(track1).at(track2));
               AngleBetweenTracks_stack -> Fill(Truth_topology, TPCObj_PFP_track_AngleBetweenTracks->at(track1).at(track2));
            }
         }
*/
         // Angle between track (using theta and phi)...
         for (int track1 = 0; track1 < TPCObj_NPFPs; track1++) {
            if(TPCObj_PFP_track_length -> at(track1) < 0) continue;

            for (int track2 = track1+1; track2 < TPCObj_NPFPs; track2++) {
               if(TPCObj_PFP_track_length -> at(track2) < 0) continue;

               TVector3 v0(1,1,1);
               TVector3 v1(1,1,1);
               v0.SetMag(1);
               v1.SetMag(1);
               v0.SetTheta(TPCObj_PFP_track_theta->at(track1));
               v1.SetTheta(TPCObj_PFP_track_theta->at(track2));
               v0.SetPhi(TPCObj_PFP_track_phi->at(track1));
               v1.SetPhi(TPCObj_PFP_track_phi->at(track2));
               double angle = TMath::ACos(v0.Dot(v1));
               AngleBetweenTracks -> Fill(angle);
               AngleBetweenTracks_stack -> Fill(Truth_topology, angle);

               TVector3 start1(TPCObj_PFP_track_start->at(track1).at(0),TPCObj_PFP_track_start->at(track1).at(1),TPCObj_PFP_track_start->at(track1).at(2));
               TVector3 end1(TPCObj_PFP_track_end->at(track1).at(0),TPCObj_PFP_track_end->at(track1).at(1),TPCObj_PFP_track_end->at(track1).at(2));
               TVector3 start2(TPCObj_PFP_track_start->at(track2).at(0),TPCObj_PFP_track_start->at(track2).at(1),TPCObj_PFP_track_start->at(track2).at(2));
               TVector3 end2(TPCObj_PFP_track_end->at(track2).at(0),TPCObj_PFP_track_end->at(track2).at(1),TPCObj_PFP_track_end->at(track2).at(2));
               std::vector<double> distances = { (start1-start2).Mag(),(start1-end2).Mag(), (end1-start2).Mag(), (end1-end2).Mag() };
               double mindist = *std::min_element(distances.begin(),distances.end());
               if(mindist < 3) {
                  AngleBetweenTracks_Close -> Fill(angle);
                  AngleBetweenTracks_Close_stack -> Fill(Truth_topology, angle);
                  if(TPCObj_PFP_MCPid->at(track1)==TPCObj_PFP_MCPid->at(track2)) AngleBetweenTracks_Close_sameID_stack -> Fill(Truth_topology, angle);
                  else AngleBetweenTracks_Close_diffID_stack -> Fill(Truth_topology, angle);
               }
            }
         }

         //evdinfo << run_num << " " << subrun_num << " " << event_num << " " << track_daughters << " " << topologyenum2str(Truth_topology) << std::endl;

      } // end if(SelectedEvent)
   } //end loop over entries

   //evdinfo.close();

   gStyle->SetOptStat(0); //No stats box
//   gStyle -> SetOptStat(111111); //Full stats box, including under/overflow

   TCanvas *c1 = new TCanvas("c1", "c1");

   eff_nuE -> Draw("AP");
   pur_nuE -> Draw("P SAME");
   gPad->Update();
   eff_nuE->SetMarkerColor(kRed);
   eff_nuE->SetMarkerStyle(21);
   eff_nuE->SetMarkerSize(1);
   auto graph = eff_nuE->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   pur_nuE->SetMarkerColor(kBlue);
   pur_nuE->SetMarkerStyle(21);
   pur_nuE->SetMarkerSize(1);
   graph = pur_nuE->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   gPad->Update();
   c1 -> SaveAs(TString::Format("%s_effpur_nuE.eps",SaveString.c_str()));

   eff_muP -> Draw("AP");
   pur_muP -> Draw("P SAME");
   gPad->Update();
   eff_muP->SetMarkerColor(kRed);
   eff_muP->SetMarkerStyle(21);
   eff_muP->SetMarkerSize(1);
   graph = eff_muP->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   pur_muP->SetMarkerColor(kBlue);
   pur_muP->SetMarkerStyle(21);
   pur_muP->SetMarkerSize(1);
   graph = pur_muP->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   gPad->Update();
   c1 -> SaveAs(TString::Format("%s_effpur_muP.eps",SaveString.c_str()));

   eff_muCosT -> Draw("AP");
   pur_muCosT -> Draw("P SAME");
   gPad->Update();
   eff_muCosT->SetMarkerColor(kRed);
   eff_muCosT->SetMarkerStyle(21);
   eff_muCosT->SetMarkerSize(1);
   graph = eff_muCosT->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   pur_muCosT->SetMarkerColor(kBlue);
   pur_muCosT->SetMarkerStyle(21);
   pur_muCosT->SetMarkerSize(1);
   graph = pur_muCosT->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   gPad->Update();
   c1 -> SaveAs(TString::Format("%s_effpur_muCosT.eps",SaveString.c_str()));

   eff_muPhi -> Draw("AP");
   pur_muPhi -> Draw("P SAME");
   gPad->Update();
   eff_muPhi->SetMarkerColor(kRed);
   eff_muPhi->SetMarkerStyle(21);
   eff_muPhi->SetMarkerSize(1);
   graph = eff_muPhi->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   pur_muPhi->SetMarkerColor(kBlue);
   pur_muPhi->SetMarkerStyle(21);
   pur_muPhi->SetMarkerSize(1);
   graph = pur_muPhi->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   gPad->Update();
   c1 -> SaveAs(TString::Format("%s_effpur_muPhi.eps",SaveString.c_str()));

   double pot_marco = 3.782e+19;
   double pot_total = GetPOT(FileName);
   double norm = pot_marco/pot_total;

   len_stack -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THStack_len.eps",SaveString.c_str()));

   mips_stack -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THStack_mips.eps",SaveString.c_str()));

   pfps_stack -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THStack_pfps.eps",SaveString.c_str()));

   vtxtrack_stack -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THStack_vtxtrack.eps",SaveString.c_str()));

   for (int j=0; j<10; j++){
      len_stack_byPDG[j]->DrawStack(norm, c1);
      c1->SaveAs(TString::Format("%s_THStack_trklen_track%d_byPDG.eps",SaveString.c_str(),j));
   }

   for (int j=0; j<5; j++){
      len_stack_MIPs_byPDG[j]->DrawStack(norm, c1);
      c1->SaveAs(TString::Format("%s_THStack_MIPs_trklen_track%d_byPDG.eps",SaveString.c_str(),j));
   }

   for (int j=0; j<10; j++){
      mom_stack_byPDG[j]->DrawStack(norm, c1);
      c1->SaveAs(TString::Format("%s_THStack_trkmom_track%d_byPDG.eps",SaveString.c_str(),j));
   }

   for (int j=0; j<5; j++){
      mom_stack_MIPs_byPDG[j]->DrawStack(norm, c1);
      c1->SaveAs(TString::Format("%s_THStack_MIPs_trkmom_track%d_byPDG.eps",SaveString.c_str(),j));
   }

   isMIP_stack_byPDG -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THStack_isMIP_byPDG.eps",SaveString.c_str()));

   mips_stack_byPDG -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THStack_mips_byPDG.eps",SaveString.c_str()));

   pfps_stack_byPDG -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THStack_pfps_byPDG.eps",SaveString.c_str()));

   n2llh_muMIPminusp_stack_byPDG -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THStack_n2llh_muMIPminusp_byPDG.eps",SaveString.c_str()));

   n2llh_muMIPminusp_stack_byPDG -> DrawStack(norm, c1, "pads");
   c1 -> SaveAs(TString::Format("%s_nostack_n2llh_muMIPminusp_byPDG.eps",SaveString.c_str()));

   n2llh_p_stack_byPDG -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THStack_n2llh_p_byPDG.eps",SaveString.c_str()));

   n2llh_p_stack_byPDG -> DrawStack(norm, c1, "pads");
   c1 -> SaveAs(TString::Format("%s_nostack_n2llh_p_byPDG.eps",SaveString.c_str()));

   depE_rangeE_mu_stack_byPDG -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THStack_depE_rangeE_mu_byPDG.eps",SaveString.c_str()));

   depE_rangeE_p_stack_byPDG -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THStack_depE_rangeE_p_byPDG.eps",SaveString.c_str()));


   longest2tracks_byPDG_signal->GetXaxis()->SetBinLabel(1,PDGenum2str(kMuMinus).c_str());
   longest2tracks_byPDG_signal->GetXaxis()->SetBinLabel(2,PDGenum2str(kPiPlus).c_str());
   longest2tracks_byPDG_signal->GetXaxis()->SetBinLabel(3,PDGenum2str(kProton).c_str());
   longest2tracks_byPDG_signal->GetXaxis()->SetBinLabel(4,"Other");
   longest2tracks_byPDG_signal->GetYaxis()->SetBinLabel(1,PDGenum2str(kMuMinus).c_str());
   longest2tracks_byPDG_signal->GetYaxis()->SetBinLabel(2,PDGenum2str(kPiPlus).c_str());
   longest2tracks_byPDG_signal->GetYaxis()->SetBinLabel(3,PDGenum2str(kProton).c_str());
   longest2tracks_byPDG_signal->GetYaxis()->SetBinLabel(4,"Other");
   longest2tracks_byPDG_signal->Draw("COLZ TEXT");
   c1 -> SaveAs(TString::Format("%s_longest2tracks_byPDG_signal.eps",SaveString.c_str()));

   longest2tracks_byPDG_bg->GetXaxis()->SetBinLabel(1,PDGenum2str(kMuMinus).c_str());
   longest2tracks_byPDG_bg->GetXaxis()->SetBinLabel(2,PDGenum2str(kPiPlus).c_str());
   longest2tracks_byPDG_bg->GetXaxis()->SetBinLabel(3,PDGenum2str(kProton).c_str());
   longest2tracks_byPDG_bg->GetXaxis()->SetBinLabel(4,"Other");
   longest2tracks_byPDG_bg->GetYaxis()->SetBinLabel(1,PDGenum2str(kMuMinus).c_str());
   longest2tracks_byPDG_bg->GetYaxis()->SetBinLabel(2,PDGenum2str(kPiPlus).c_str());
   longest2tracks_byPDG_bg->GetYaxis()->SetBinLabel(3,PDGenum2str(kProton).c_str());
   longest2tracks_byPDG_bg->GetYaxis()->SetBinLabel(4,"Other");
   longest2tracks_byPDG_bg->Draw("COLZ TEXT");
   c1 -> SaveAs(TString::Format("%s_longest2tracks_byPDG_bg.eps",SaveString.c_str()));

//   vtx_truereco_dist -> Scale(norm);
//   vtx_truereco_dist -> Draw();
//   c1 -> SaveAs(TString::Format("%s_vtx_truereco_dist.eps",SaveString.c_str()));
   vtx_truereco_spacecharge_dist -> Scale(norm);
   vtx_truereco_spacecharge_dist -> Draw();
   c1 -> SaveAs(TString::Format("%s_vtx_truereco_spacecharge_dist.eps",SaveString.c_str()));

   eff_n2llh_muMIPminusp_atleast2MIP -> Draw("AP");
   pur_n2llh_muMIPminusp_atleast2MIP -> Draw("P SAME");
   effpur_n2llh_muMIPminusp_atleast2MIP -> Draw("P SAME");
   gPad->Update();
   eff_n2llh_muMIPminusp_atleast2MIP->SetMarkerColor(kRed);
   eff_n2llh_muMIPminusp_atleast2MIP->SetMarkerStyle(21);
   eff_n2llh_muMIPminusp_atleast2MIP->SetMarkerSize(1);
   graph = eff_n2llh_muMIPminusp_atleast2MIP->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   pur_n2llh_muMIPminusp_atleast2MIP->SetMarkerColor(kBlue);
   pur_n2llh_muMIPminusp_atleast2MIP->SetMarkerStyle(21);
   pur_n2llh_muMIPminusp_atleast2MIP->SetMarkerSize(1);
   graph = pur_n2llh_muMIPminusp_atleast2MIP->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   effpur_n2llh_muMIPminusp_atleast2MIP->SetMarkerColor(kBlack);
   effpur_n2llh_muMIPminusp_atleast2MIP->SetMarkerStyle(21);
   effpur_n2llh_muMIPminusp_atleast2MIP->SetMarkerSize(1);
   effpur_n2llh_muMIPminusp_atleast2MIP->SetMinimum(0);
   effpur_n2llh_muMIPminusp_atleast2MIP->SetMaximum(1);
   gPad->Update();
   c1 -> SaveAs(TString::Format("%s_effpur_n2llh_muMIPminusp_atleast2MIP.eps",SaveString.c_str()));

   eff_n2llh_muMIPminusp_exactly2MIP -> Draw("AP");
   pur_n2llh_muMIPminusp_exactly2MIP -> Draw("P SAME");
   effpur_n2llh_muMIPminusp_exactly2MIP -> Draw("P SAME");
   gPad->Update();
   eff_n2llh_muMIPminusp_exactly2MIP->SetMarkerColor(kRed);
   eff_n2llh_muMIPminusp_exactly2MIP->SetMarkerStyle(21);
   eff_n2llh_muMIPminusp_exactly2MIP->SetMarkerSize(1);
   graph = eff_n2llh_muMIPminusp_exactly2MIP->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   pur_n2llh_muMIPminusp_exactly2MIP->SetMarkerColor(kBlue);
   pur_n2llh_muMIPminusp_exactly2MIP->SetMarkerStyle(21);
   pur_n2llh_muMIPminusp_exactly2MIP->SetMarkerSize(1);
   graph = pur_n2llh_muMIPminusp_exactly2MIP->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   effpur_n2llh_muMIPminusp_exactly2MIP->SetMarkerColor(kBlack);
   effpur_n2llh_muMIPminusp_exactly2MIP->SetMarkerStyle(21);
   effpur_n2llh_muMIPminusp_exactly2MIP->SetMarkerSize(1);
   effpur_n2llh_muMIPminusp_exactly2MIP->SetMinimum(0);
   effpur_n2llh_muMIPminusp_exactly2MIP->SetMaximum(1);
   gPad->Update();
   c1 -> SaveAs(TString::Format("%s_effpur_n2llh_muMIPminusp_exactly2MIP.eps",SaveString.c_str()));

   eff_n2llh_p_atleast2MIP -> Draw("AP");
   pur_n2llh_p_atleast2MIP -> Draw("P SAME");
   effpur_n2llh_p_atleast2MIP -> Draw("P SAME");
   gPad->Update();
   eff_n2llh_p_atleast2MIP->SetMarkerColor(kRed);
   eff_n2llh_p_atleast2MIP->SetMarkerStyle(21);
   eff_n2llh_p_atleast2MIP->SetMarkerSize(1);
   graph = eff_n2llh_p_atleast2MIP->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   pur_n2llh_p_atleast2MIP->SetMarkerColor(kBlue);
   pur_n2llh_p_atleast2MIP->SetMarkerStyle(21);
   pur_n2llh_p_atleast2MIP->SetMarkerSize(1);
   graph = pur_n2llh_p_atleast2MIP->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   effpur_n2llh_p_atleast2MIP->SetMarkerColor(kBlack);
   effpur_n2llh_p_atleast2MIP->SetMarkerStyle(21);
   effpur_n2llh_p_atleast2MIP->SetMarkerSize(1);
   effpur_n2llh_p_atleast2MIP->SetMinimum(0);
   effpur_n2llh_p_atleast2MIP->SetMaximum(1);
   gPad->Update();
   c1 -> SaveAs(TString::Format("%s_effpur_n2llh_p_atleast2MIP.eps",SaveString.c_str()));

   eff_n2llh_p_exactly2MIP -> Draw("AP");
   pur_n2llh_p_exactly2MIP -> Draw("P SAME");
   effpur_n2llh_p_exactly2MIP -> Draw("P SAME");
   gPad->Update();
   eff_n2llh_p_exactly2MIP->SetMarkerColor(kRed);
   eff_n2llh_p_exactly2MIP->SetMarkerStyle(21);
   eff_n2llh_p_exactly2MIP->SetMarkerSize(1);
   graph = eff_n2llh_p_exactly2MIP->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   pur_n2llh_p_exactly2MIP->SetMarkerColor(kBlue);
   pur_n2llh_p_exactly2MIP->SetMarkerStyle(21);
   pur_n2llh_p_exactly2MIP->SetMarkerSize(1);
   graph = pur_n2llh_p_exactly2MIP->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   effpur_n2llh_p_exactly2MIP->SetMarkerColor(kBlack);
   effpur_n2llh_p_exactly2MIP->SetMarkerStyle(21);
   effpur_n2llh_p_exactly2MIP->SetMarkerSize(1);
   effpur_n2llh_p_exactly2MIP->SetMinimum(0);
   effpur_n2llh_p_exactly2MIP->SetMaximum(1);
   gPad->Update();
   c1 -> SaveAs(TString::Format("%s_effpur_n2llh_p_exactly2MIP.eps",SaveString.c_str()));

   eff_depE_rangeE_mu_atleast2MIP -> Draw("AP");
   pur_depE_rangeE_mu_atleast2MIP -> Draw("P SAME");
   effpur_depE_rangeE_mu_atleast2MIP -> Draw("P SAME");
   gPad->Update();
   eff_depE_rangeE_mu_atleast2MIP->SetMarkerColor(kRed);
   eff_depE_rangeE_mu_atleast2MIP->SetMarkerStyle(21);
   eff_depE_rangeE_mu_atleast2MIP->SetMarkerSize(1);
   graph = eff_depE_rangeE_mu_atleast2MIP->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   pur_depE_rangeE_mu_atleast2MIP->SetMarkerColor(kBlue);
   pur_depE_rangeE_mu_atleast2MIP->SetMarkerStyle(21);
   pur_depE_rangeE_mu_atleast2MIP->SetMarkerSize(1);
   graph = pur_depE_rangeE_mu_atleast2MIP->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   effpur_depE_rangeE_mu_atleast2MIP->SetMarkerColor(kBlack);
   effpur_depE_rangeE_mu_atleast2MIP->SetMarkerStyle(21);
   effpur_depE_rangeE_mu_atleast2MIP->SetMarkerSize(1);
   effpur_depE_rangeE_mu_atleast2MIP->SetMinimum(0);
   effpur_depE_rangeE_mu_atleast2MIP->SetMaximum(1);
   gPad->Update();
   c1 -> SaveAs(TString::Format("%s_effpur_depE_rangeE_mu_atleast2MIP.eps",SaveString.c_str()));

   eff_depE_rangeE_mu_exactly2MIP -> Draw("AP");
   pur_depE_rangeE_mu_exactly2MIP -> Draw("P SAME");
   effpur_depE_rangeE_mu_exactly2MIP -> Draw("P SAME");
   gPad->Update();
   eff_depE_rangeE_mu_exactly2MIP->SetMarkerColor(kRed);
   eff_depE_rangeE_mu_exactly2MIP->SetMarkerStyle(21);
   eff_depE_rangeE_mu_exactly2MIP->SetMarkerSize(1);
   graph = eff_depE_rangeE_mu_exactly2MIP->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   pur_depE_rangeE_mu_exactly2MIP->SetMarkerColor(kBlue);
   pur_depE_rangeE_mu_exactly2MIP->SetMarkerStyle(21);
   pur_depE_rangeE_mu_exactly2MIP->SetMarkerSize(1);
   graph = pur_depE_rangeE_mu_exactly2MIP->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   effpur_depE_rangeE_mu_exactly2MIP->SetMarkerColor(kBlack);
   effpur_depE_rangeE_mu_exactly2MIP->SetMarkerStyle(21);
   effpur_depE_rangeE_mu_exactly2MIP->SetMarkerSize(1);
   effpur_depE_rangeE_mu_exactly2MIP->SetMinimum(0);
   effpur_depE_rangeE_mu_exactly2MIP->SetMaximum(1);
   gPad->Update();
   c1 -> SaveAs(TString::Format("%s_effpur_depE_rangeE_mu_exactly2MIP.eps",SaveString.c_str()));

   eff_depE_rangeE_p_atleast2MIP -> Draw("AP");
   pur_depE_rangeE_p_atleast2MIP -> Draw("P SAME");
   effpur_depE_rangeE_p_atleast2MIP -> Draw("P SAME");
   gPad->Update();
   eff_depE_rangeE_p_atleast2MIP->SetMarkerColor(kRed);
   eff_depE_rangeE_p_atleast2MIP->SetMarkerStyle(21);
   eff_depE_rangeE_p_atleast2MIP->SetMarkerSize(1);
   graph = eff_depE_rangeE_p_atleast2MIP->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   pur_depE_rangeE_p_atleast2MIP->SetMarkerColor(kBlue);
   pur_depE_rangeE_p_atleast2MIP->SetMarkerStyle(21);
   pur_depE_rangeE_p_atleast2MIP->SetMarkerSize(1);
   graph = pur_depE_rangeE_p_atleast2MIP->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   effpur_depE_rangeE_p_atleast2MIP->SetMarkerColor(kBlack);
   effpur_depE_rangeE_p_atleast2MIP->SetMarkerStyle(21);
   effpur_depE_rangeE_p_atleast2MIP->SetMarkerSize(1);
   effpur_depE_rangeE_p_atleast2MIP->SetMinimum(0);
   effpur_depE_rangeE_p_atleast2MIP->SetMaximum(1);
   gPad->Update();
   c1 -> SaveAs(TString::Format("%s_effpur_depE_rangeE_p_atleast2MIP.eps",SaveString.c_str()));

   eff_depE_rangeE_p_exactly2MIP -> Draw("AP");
   pur_depE_rangeE_p_exactly2MIP -> Draw("P SAME");
   effpur_depE_rangeE_p_exactly2MIP -> Draw("P SAME");
   gPad->Update();
   eff_depE_rangeE_p_exactly2MIP->SetMarkerColor(kRed);
   eff_depE_rangeE_p_exactly2MIP->SetMarkerStyle(21);
   eff_depE_rangeE_p_exactly2MIP->SetMarkerSize(1);
   graph = eff_depE_rangeE_p_exactly2MIP->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   pur_depE_rangeE_p_exactly2MIP->SetMarkerColor(kBlue);
   pur_depE_rangeE_p_exactly2MIP->SetMarkerStyle(21);
   pur_depE_rangeE_p_exactly2MIP->SetMarkerSize(1);
   graph = pur_depE_rangeE_p_exactly2MIP->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   effpur_depE_rangeE_p_exactly2MIP->SetMarkerColor(kBlack);
   effpur_depE_rangeE_p_exactly2MIP->SetMarkerStyle(21);
   effpur_depE_rangeE_p_exactly2MIP->SetMarkerSize(1);
   effpur_depE_rangeE_p_exactly2MIP->SetMinimum(0);
   effpur_depE_rangeE_p_exactly2MIP->SetMaximum(1);
   gPad->Update();
   c1 -> SaveAs(TString::Format("%s_effpur_depE_rangeE_p_exactly2MIP.eps",SaveString.c_str()));

   gStyle->SetPaintTextFormat("4.2f");
   eff_n2llh_muMIPminusp_depE_rangeE_p_atleast2MIP -> Draw("COLZ TEXT");
   c1 -> SaveAs(TString::Format("%s_eff_n2llh_muMIPminusp_depE_rangeE_p_atleast2MIP.eps",SaveString.c_str()));
   pur_n2llh_muMIPminusp_depE_rangeE_p_atleast2MIP -> Draw("COLZ TEXT");
   c1 -> SaveAs(TString::Format("%s_pur_n2llh_muMIPminusp_depE_rangeE_p_atleast2MIP.eps",SaveString.c_str()));
   effpur_n2llh_muMIPminusp_depE_rangeE_p_atleast2MIP -> Draw("COLZ TEXT");
   c1 -> SaveAs(TString::Format("%s_effpur_n2llh_muMIPminusp_depE_rangeE_p_atleast2MIP.eps",SaveString.c_str()));
   eff_n2llh_muMIPminusp_depE_rangeE_p_exactly2MIP -> Draw("COLZ TEXT");
   c1 -> SaveAs(TString::Format("%s_eff_n2llh_muMIPminusp_depE_rangeE_p_exactly2MIP.eps",SaveString.c_str()));
   pur_n2llh_muMIPminusp_depE_rangeE_p_exactly2MIP -> Draw("COLZ TEXT");
   c1 -> SaveAs(TString::Format("%s_pur_n2llh_muMIPminusp_depE_rangeE_p_exactly2MIP.eps",SaveString.c_str()));
   effpur_n2llh_muMIPminusp_depE_rangeE_p_exactly2MIP -> Draw("COLZ TEXT");
   c1 -> SaveAs(TString::Format("%s_effpur_n2llh_muMIPminusp_depE_rangeE_p_exactly2MIP.eps",SaveString.c_str()));

   eff_distcut -> Draw("AP");
   pur_distcut -> Draw("P SAME");
   effpur_distcut -> Draw("P SAME");
   gPad->Update();
   eff_distcut->SetMarkerColor(kRed);
   eff_distcut->SetMarkerStyle(21);
   eff_distcut->SetMarkerSize(1);
   graph = eff_distcut->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   pur_distcut->SetMarkerColor(kBlue);
   pur_distcut->SetMarkerStyle(21);
   pur_distcut->SetMarkerSize(1);
   graph = pur_distcut->GetPaintedGraph();
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   effpur_distcut->SetMarkerColor(kBlack);
   effpur_distcut->SetMarkerStyle(21);
   effpur_distcut->SetMarkerSize(1);
   effpur_distcut->SetMinimum(0);
   effpur_distcut->SetMaximum(1);
   gPad->Update();
   c1 -> SaveAs(TString::Format("%s_effpur_distcut.eps",SaveString.c_str()));

   AngleBetweenTracks -> Draw();
   c1 -> SaveAs(TString::Format("%s_AngleBetweenTracks.eps",SaveString.c_str()));
   AngleBetweenTracks_stack -> DrawStack(1,c1);
   c1 -> SaveAs(TString::Format("%s_AngleBetweenTracks_stack.eps",SaveString.c_str()));
   AngleBetweenTracks_Close -> Draw();
   c1 -> SaveAs(TString::Format("%s_AngleBetweenTracks_Close.eps",SaveString.c_str()));
   AngleBetweenTracks_Close_stack -> DrawStack(1,c1);
   c1 -> SaveAs(TString::Format("%s_AngleBetweenTracks_Close_stack.eps",SaveString.c_str()));
   AngleBetweenTracks_Close_sameID_stack -> DrawStack(1,c1);
   c1 -> SaveAs(TString::Format("%s_AngleBetweenTracks_Close_sameID_stack.eps",SaveString.c_str()));
   AngleBetweenTracks_Close_diffID_stack -> DrawStack(1,c1);
   c1 -> SaveAs(TString::Format("%s_AngleBetweenTracks_Close_diffID_stack.eps",SaveString.c_str()));

   std::cout << "Total Efficiency: " << (1. * (eff_nuE -> GetPassedHistogram() -> GetEntries()))/(eff_nuE -> GetTotalHistogram() -> GetEntries()) << std::endl;
   std::cout << "Total Purity: " << (1. * (pur_nuE -> GetPassedHistogram() -> GetEntries()))/(pur_nuE -> GetTotalHistogram() -> GetEntries()) << std::endl;

   gDirectory->GetList()->Delete();

   return;
}
