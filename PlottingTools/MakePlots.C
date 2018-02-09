#include "MakePlots.h"

double GetPOT(TString FileName) {

   TChain *t_pot = new TChain("cc1piselec/pottree");
   t_pot -> Add(FileName);

   double pot;
   t_pot -> SetBranchAddress("pot", &pot);

   double totalPOT = 0;

   const int potentries = t_pot -> GetEntries();
   for (int i = 0; i < potentries; i++) {
      t_pot -> GetEntry(i);
      totalPOT += pot;
   }

   return totalPOT;

}

// Can swap to using Marco's FV calculator to make this cleaner?
// Can also just put this in the tree once we settle on a particular FV def
const double xmin = 0;
const double xmax = 256.35;
const double ymin = -115.53;
const double ymax = 117.47;
const double zmin = 0.1;
const double zmax = 1036.9;
const double FVxmin = xmin + 12;
const double FVxmax = xmax - 12;
const double FVymin = ymin + 35;
const double FVymax = ymax - 35;
const double FVzmin = zmin + 25;
const double FVzmax = zmax - 85;

bool inFV(double x, double y, double z){

//   if (x > FVxmin && x < FVxmax && y > FVymin && y < FVymax && z > FVzmin && z < FVzmax && (z < 675 || z > 775)) return true;
   if (x > FVxmin && x < FVxmax && y > FVymin && y < FVymax && z > FVzmin && z < FVzmax) return true;
   else return false;

}


void MakePlots(std::map<std::string,bool> SelectionCutflow, std::string SaveString, TString FileName) {

   //Setup tree...
   TChain *t = new TChain("cc1piselec/outtree");
   t -> Add(FileName);

   bool isData;
   unsigned int run_num;
   unsigned int subrun_num;
   unsigned int event_num;

   std::map<std::string,bool> *cutflow = NULL;
   bool isSelected;
   std::vector<double> *track_length = NULL;
   std::vector<double> *shower_length = NULL;
   int NPFPs;
   int NTracks;
   int NShowers;
   std::vector<bool> *Sel_PFP_isTrack = NULL;
   std::vector<bool> *Sel_PFP_isShower = NULL;
   std::vector<bool> *Sel_PFP_isDaughter = NULL;
//   std::vector<bool> *Sel_PFP_isMIP = NULL;
   std::vector<int> *Sel_PFP_ID = NULL;
   std::vector<int> *Sel_MCP_ID = NULL;
   std::vector<int> *Sel_MCP_PDG = NULL;
   std::vector<double> *Sel_MCP_E = NULL;
   int tpcobj_origin;
   int tpcobj_origin_extra;

   std::vector<int> *MCP_PDG = NULL;
   std::vector<double> *MCP_length = NULL;
   std::vector<std::string> *MCP_process = NULL;
   std::vector<std::string> *MCP_endprocess = NULL;
   std::vector<int> *MCP_numdaughters = NULL;
   std::vector<double> *MCP_P = NULL;
   std::vector<double> *MCP_Px = NULL;
   std::vector<double> *MCP_Py = NULL;
   std::vector<double> *MCP_Pz = NULL;
   std::vector<double> *MCP_E = NULL;
//   std::vector<bool> *MCP_isContained = NULL;
   std::vector<bool> *MIPConsistency = NULL;

   std::vector<double> *nu_vtxx = NULL;
   std::vector<double> *nu_vtxy = NULL;
   std::vector<double> *nu_vtxz = NULL;
   std::vector<bool> *nu_isCC = NULL;
   std::vector<int> *nu_PDG = NULL;
   std::vector<double> *nu_E = NULL;

   std::vector<double> *dqdx_trunc_uncalib = NULL;

   std::map<std::string,bool> *CC1picutflow = NULL;
   bool PassesCC1piSelec;
   std::string *CC1piSelecFailureReason = NULL;

   t -> SetBranchAddress("isData", &isData);
   t -> SetBranchAddress("run_num", &run_num);
   t -> SetBranchAddress("subrun_num", &subrun_num);
   t -> SetBranchAddress("event_num", &event_num);

   t -> SetBranchAddress("cutflow", &cutflow);
   t -> SetBranchAddress("isSelected", &isSelected);
   t -> SetBranchAddress("track_length", &track_length);
   t -> SetBranchAddress("shower_length", &shower_length);
   t -> SetBranchAddress("NPFPs", &NPFPs);
   t -> SetBranchAddress("NTracks", &NTracks);
   t -> SetBranchAddress("NShowers", &NShowers);
   t -> SetBranchAddress("Sel_PFP_isTrack", &Sel_PFP_isTrack);
   t -> SetBranchAddress("Sel_PFP_isShower", &Sel_PFP_isShower);
   t -> SetBranchAddress("Sel_PFP_isDaughter", &Sel_PFP_isDaughter);
//   t -> SetBranchAddress("Sel_PFP_isMIP", &Sel_PFP_isMIP);
   t -> SetBranchAddress("Sel_PFP_ID", &Sel_PFP_ID);
   t -> SetBranchAddress("Sel_MCP_ID", &Sel_MCP_ID);
   t -> SetBranchAddress("Sel_MCP_PDG", &Sel_MCP_PDG);
   t -> SetBranchAddress("Sel_MCP_E", &Sel_MCP_E);
   t -> SetBranchAddress("tpcobj_origin", &tpcobj_origin);
   t -> SetBranchAddress("tpcobj_origin_extra", &tpcobj_origin_extra);

   t -> SetBranchAddress("MCP_PDG", &MCP_PDG);
   t -> SetBranchAddress("MCP_length", &MCP_length);
   t -> SetBranchAddress("MCP_process", &MCP_process);
   t -> SetBranchAddress("MCP_endprocess", &MCP_endprocess);
   t -> SetBranchAddress("MCP_numdaughters", &MCP_numdaughters);
   t -> SetBranchAddress("MCP_P", &MCP_P);
   t -> SetBranchAddress("MCP_Px", &MCP_Px);
   t -> SetBranchAddress("MCP_Py", &MCP_Py);
   t -> SetBranchAddress("MCP_Pz", &MCP_Pz);
   t -> SetBranchAddress("MCP_E", &MCP_E);
//   t -> SetBranchAddress("MCP_isContained", &MCP_isContained);
   t -> SetBranchAddress("MIPConsistency", &MIPConsistency);

   t -> SetBranchAddress("nu_vtxx", &nu_vtxx);
   t -> SetBranchAddress("nu_vtxy", &nu_vtxy);
   t -> SetBranchAddress("nu_vtxz", &nu_vtxz);
   t -> SetBranchAddress("nu_isCC", &nu_isCC);
   t -> SetBranchAddress("nu_PDG", &nu_PDG);
   t -> SetBranchAddress("nu_E", &nu_E);

   t -> SetBranchAddress("dqdx_trunc_uncalib", &dqdx_trunc_uncalib);

   t -> SetBranchAddress("CC1picutflow", &CC1picutflow);
   t -> SetBranchAddress("PassesCC1piSelec", &PassesCC1piSelec);
   t -> SetBranchAddress("CC1piSelecFailureReason", &CC1piSelecFailureReason);

   TEfficiency* eff_nuE = new TEfficiency("eff_nuE",";True Neutrino Energy [GeV];",15,0,3);
   TEfficiency* eff_muP = new TEfficiency("eff_muP",";True Muon Momentum [GeV];",15,0,2);
   TEfficiency* eff_muCosT = new TEfficiency("eff_muCosT",";True Muon cos(#theta);",10,-1,1);
   TEfficiency* eff_muPhi = new TEfficiency("eff_muPhi",";True Muon #phi angle;",10,-3.2,3.2);

   TEfficiency* pur_nuE = new TEfficiency("pur_nuE",";True Neutrino Energy [GeV];",15,0,3);
   TEfficiency* pur_muP = new TEfficiency("pur_muP",";True Muon Momentum [GeV];",15,0,2);
   TEfficiency* pur_muCosT = new TEfficiency("pur_muCosT",";True Muon cos(#theta);",10,-1,1);
   TEfficiency* pur_muPhi = new TEfficiency("pur_muPhi",";True Muon #phi angle;",10,-3.2,3.2);

   //TO DO: Make a seperate THStack making function
   THStack* len_stack = new THStack("len_stack", ";Longest Track Length [cm];Selected Events (normalised to 3.782 #times 10^{19} POT)");
   TH1F* len_muCC1pip = new TH1F("len_muCC1pip", ";;", 30, 0, 700);
   TH1F* len_otherCC1pip = new TH1F("len_otherCC1pip", ";;", 30, 0, 700);
   TH1F* len_CCother = new TH1F("len_CCother", ";;", 30, 0, 700);
   TH1F* len_NC = new TH1F("len_NC", ";;", 30, 0, 700);
   TH1F* len_cosmic = new TH1F("len_cosmic", ";;", 30, 0, 700);
   TH1F* len_mixed = new TH1F("len_mixed", ";;", 30, 0, 700);
   TH1F* len_unknown = new TH1F("len_unknown", ";;", 30, 0, 700);
   TH1F* len_outFV = new TH1F("len_outFV", ";;", 30, 0, 700);

   THStack* tracks_stack = new THStack("tracks_stack", ";Number of #nu daughter tracks in selected TPCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)");
   TH1F* tracks_muCC1pip = new TH1F("tracks_muCC1pip", ";;", 10, 0, 10);
   TH1F* tracks_otherCC1pip = new TH1F("tracks_otherCC1pip", ";;", 10, 0, 10);
   TH1F* tracks_CCother = new TH1F("tracks_CCother", ";;", 10, 0, 10);
   TH1F* tracks_NC = new TH1F("tracks_NC", ";;", 10, 0, 10);
   TH1F* tracks_cosmic = new TH1F("tracks_cosmic", ";;", 10, 0, 10);
   TH1F* tracks_mixed = new TH1F("tracks_mixed", ";;", 10, 0, 10);
   TH1F* tracks_unknown = new TH1F("tracks_unknown", ";;", 10, 0, 10);
   TH1F* tracks_outFV = new TH1F("tracks_outFV", ";;", 10, 0, 10);

   THStack* showers_stack = new THStack("showers_stack", ";Number of #nu daughter showers in selected TCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)");
   TH1F* showers_muCC1pip = new TH1F("showers_muCC1pip", ";;", 10, 0, 10);
   TH1F* showers_otherCC1pip = new TH1F("showers_otherCC1pip", ";;", 10, 0, 10);
   TH1F* showers_CCother = new TH1F("showers_CCother", ";;", 10, 0, 10);
   TH1F* showers_NC = new TH1F("showers_NC", ";;", 10, 0, 10);
   TH1F* showers_cosmic = new TH1F("showers_cosmic", ";;", 10, 0, 10);
   TH1F* showers_mixed = new TH1F("showers_mixed", ";;", 10, 0, 10);
   TH1F* showers_unknown = new TH1F("showers_unknown", ";;", 10, 0, 10);
   TH1F* showers_outFV = new TH1F("showers_outFV", ";;", 10, 0, 10);

   THStack* mips_stack = new THStack("mips_stack", ";Number of MIP-like #nu daughter tracks in selected TPCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)");
   TH1F* mips_muCC1pip = new TH1F("mips_muCC1pip", ";;", 10, 0, 10);
   TH1F* mips_otherCC1pip = new TH1F("mips_otherCC1pip", ";;", 10, 0, 10);
   TH1F* mips_CCother = new TH1F("mips_CCother", ";;", 10, 0, 10);
   TH1F* mips_NC = new TH1F("mips_NC", ";;", 10, 0, 10);
   TH1F* mips_cosmic = new TH1F("mips_cosmic", ";;", 10, 0, 10);
   TH1F* mips_mixed = new TH1F("mips_mixed", ";;", 10, 0, 10);
   TH1F* mips_unknown = new TH1F("mips_unknown", ";;", 10, 0, 10);
   TH1F* mips_outFV = new TH1F("mips_outFV", ";;", 10, 0, 10);

   THStack* vtxtrack_stack = new THStack("vtxtrack_stack", ";Distance between vertex and start of #nu daughter tracks in selected TPCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)");
   TH1F* vtxtrack_muCC1pip = new TH1F("vtxtrack_muCC1pip", ";;", 50, 0, 50);
   TH1F* vtxtrack_otherCC1pip = new TH1F("vtxtrack_otherCC1pip", ";;", 50, 0, 50);
   TH1F* vtxtrack_CCother = new TH1F("vtxtrack_CCother", ";;", 50, 0, 50);
   TH1F* vtxtrack_NC = new TH1F("vtxtrack_NC", ";;", 50, 0, 50);
   TH1F* vtxtrack_cosmic = new TH1F("vtxtrack_cosmic", ";;", 50, 0, 50);
   TH1F* vtxtrack_mixed = new TH1F("vtxtrack_mixed", ";;", 50, 0, 50);
   TH1F* vtxtrack_unknown = new TH1F("vtxtrack_unknown", ";;", 50, 0, 50);
   TH1F* vtxtrack_outFV = new TH1F("vtxtrack_outFV", ";;", 50, 0, 50);

   THStack* vtxshwr_stack = new THStack("vtxshwr_stack", ";Distance between vertex and start of #nu daughter showers in selected TPCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)");
   TH1F* vtxshwr_muCC1pip = new TH1F("vtxshwr_muCC1pip", ";;", 50, 0, 50);
   TH1F* vtxshwr_otherCC1pip = new TH1F("vtxshwr_otherCC1pip", ";;", 50, 0, 50);
   TH1F* vtxshwr_CCother = new TH1F("vtxshwr_CCother", ";;", 50, 0, 50);
   TH1F* vtxshwr_NC = new TH1F("vtxshwr_NC", ";;", 50, 0, 50);
   TH1F* vtxshwr_cosmic = new TH1F("vtxshwr_cosmic", ";;", 50, 0, 50);
   TH1F* vtxshwr_mixed = new TH1F("vtxshwr_mixed", ";;", 50, 0, 50);
   TH1F* vtxshwr_unknown = new TH1F("vtxshwr_unknown", ";;", 50, 0, 50);
   TH1F* vtxshwr_outFV = new TH1F("vtxshwr_outFV", ";;", 50, 0, 50);

   THStack* shwrtrue_stack = new THStack("shwrtrue_stack", ";True CC1pi: Number of shower #nu daughter tracks in selected TPCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)");
   TH1F *shwrtrue_muplus = new TH1F("shwrtrue_muplus","",10,0,10);
   TH1F *shwrtrue_muminus = new TH1F("shwrtrue_muminus","",10,0,10);
   TH1F *shwrtrue_piplus = new TH1F("shwrtrue_piplus","",10,0,10);
   TH1F *shwrtrue_piminus = new TH1F("shwrtrue_piminus","",10,0,10);
   TH1F *shwrtrue_pizero = new TH1F("shwrtrue_pizero","",10,0,10);
   TH1F *shwrtrue_proton = new TH1F("shwrtrue_proton","",10,0,10);
   TH1F *shwrtrue_electron = new TH1F("shwrtrue_electron","",10,0,10);
   TH1F *shwrtrue_positron = new TH1F("shwrtrue_positron","",10,0,10);
   TH1F *shwrtrue_neutron = new TH1F("shwrtrue_neutron","",10,0,10);
   TH1F *shwrtrue_photon = new TH1F("shwrtrue_photon","",10,0,10);
   TH1F *shwrtrue_other = new TH1F("shwrtrue_other","",10,0,10);

   THStack* shwrtrue_bg_stack = new THStack("shwrtrue_bg_stack", ";Background: Number of shower #nu daughter tracks in selected TPCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)");
   TH1F *shwrtrue_bg_muplus = new TH1F("shwrtrue_bg_muplus","",10,0,10);
   TH1F *shwrtrue_bg_muminus = new TH1F("shwrtrue_bg_muminus","",10,0,10);
   TH1F *shwrtrue_bg_piplus = new TH1F("shwrtrue_bg_piplus","",10,0,10);
   TH1F *shwrtrue_bg_piminus = new TH1F("shwrtrue_bg_piminus","",10,0,10);
   TH1F *shwrtrue_bg_pizero = new TH1F("shwrtrue_bg_pizero","",10,0,10);
   TH1F *shwrtrue_bg_proton = new TH1F("shwrtrue_bg_proton","",10,0,10);
   TH1F *shwrtrue_bg_electron = new TH1F("shwrtrue_bg_electron","",10,0,10);
   TH1F *shwrtrue_bg_positron = new TH1F("shwrtrue_bg_positron","",10,0,10);
   TH1F *shwrtrue_bg_neutron = new TH1F("shwrtrue_bg_neutron","",10,0,10);
   TH1F *shwrtrue_bg_photon = new TH1F("shwrtrue_bg_photon","",10,0,10);
   TH1F *shwrtrue_bg_other = new TH1F("shwrtrue_bg_other","",10,0,10);
   

   len_muCC1pip -> SetFillColor(kRed);
   len_otherCC1pip -> SetFillColor(kOrange);
   len_CCother -> SetFillColor(kMagenta);
   len_NC -> SetFillColor(kGray);
   len_cosmic -> SetFillColor(kBlue);
   len_mixed -> SetFillColor(kCyan);
   len_unknown -> SetFillColor(kBlack);
   len_outFV -> SetFillColor(kGreen);

   tracks_muCC1pip -> SetFillColor(kRed);
   tracks_otherCC1pip -> SetFillColor(kOrange);
   tracks_CCother -> SetFillColor(kMagenta);
   tracks_NC -> SetFillColor(kGray);
   tracks_cosmic -> SetFillColor(kBlue);
   tracks_mixed -> SetFillColor(kCyan);
   tracks_unknown -> SetFillColor(kBlack);
   tracks_outFV -> SetFillColor(kGreen);

   showers_muCC1pip -> SetFillColor(kRed);
   showers_otherCC1pip -> SetFillColor(kOrange);
   showers_CCother -> SetFillColor(kMagenta);
   showers_NC -> SetFillColor(kGray);
   showers_cosmic -> SetFillColor(kBlue);
   showers_mixed -> SetFillColor(kCyan);
   showers_unknown -> SetFillColor(kBlack);
   showers_outFV -> SetFillColor(kGreen);

   mips_muCC1pip -> SetFillColor(kRed);
   mips_otherCC1pip -> SetFillColor(kOrange);
   mips_CCother -> SetFillColor(kMagenta);
   mips_NC -> SetFillColor(kGray);
   mips_cosmic -> SetFillColor(kBlue);
   mips_mixed -> SetFillColor(kCyan);
   mips_unknown -> SetFillColor(kBlack);
   mips_outFV -> SetFillColor(kGreen);

   vtxtrack_muCC1pip -> SetFillColor(kRed);
   vtxtrack_otherCC1pip -> SetFillColor(kOrange);
   vtxtrack_CCother -> SetFillColor(kMagenta);
   vtxtrack_NC -> SetFillColor(kGray);
   vtxtrack_cosmic -> SetFillColor(kBlue);
   vtxtrack_mixed -> SetFillColor(kCyan);
   vtxtrack_unknown -> SetFillColor(kBlack);
   vtxtrack_outFV -> SetFillColor(kGreen);

   vtxshwr_muCC1pip -> SetFillColor(kRed);
   vtxshwr_otherCC1pip -> SetFillColor(kOrange);
   vtxshwr_CCother -> SetFillColor(kMagenta);
   vtxshwr_NC -> SetFillColor(kGray);
   vtxshwr_cosmic -> SetFillColor(kBlue);
   vtxshwr_mixed -> SetFillColor(kCyan);
   vtxshwr_unknown -> SetFillColor(kBlack);
   vtxshwr_outFV -> SetFillColor(kGreen);

   shwrtrue_muminus -> SetFillColor(kRed);
   shwrtrue_muplus -> SetFillColor(kOrange);
   shwrtrue_piplus -> SetFillColor(kMagenta);
   shwrtrue_piminus -> SetFillColor(kGray);
   shwrtrue_pizero -> SetFillColor(kBlue);
   shwrtrue_proton -> SetFillColor(kCyan);
   shwrtrue_electron -> SetFillColor(kBlack);
   shwrtrue_positron -> SetFillColor(kOrange+1);
   shwrtrue_neutron -> SetFillColor(kViolet);
   shwrtrue_photon -> SetFillColor(kGreen+3);
   shwrtrue_other -> SetFillColor(kGreen);

   shwrtrue_bg_muminus -> SetFillColor(kRed);
   shwrtrue_bg_muplus -> SetFillColor(kOrange);
   shwrtrue_bg_piplus -> SetFillColor(kMagenta);
   shwrtrue_bg_piminus -> SetFillColor(kGray);
   shwrtrue_bg_pizero -> SetFillColor(kBlue);
   shwrtrue_bg_proton -> SetFillColor(kCyan);
   shwrtrue_bg_electron -> SetFillColor(kBlack);
   shwrtrue_bg_positron -> SetFillColor(kOrange+1);
   shwrtrue_bg_neutron -> SetFillColor(kViolet);
   shwrtrue_bg_photon -> SetFillColor(kGreen+3);
   shwrtrue_bg_other -> SetFillColor(kGreen);

   TEfficiency* eff_track_muon = new TEfficiency("eff_track_muon",";True Kinetic Energy [GeV];",10,0,2);
   TEfficiency* eff_track_pion = new TEfficiency("eff_track_pion",";True Kinetic Energy [GeV];",10,0,2);
   TEfficiency* eff_track_proton = new TEfficiency("eff_track_proton",";True Kinetic Energy [GeV];",10,0,2);
/*
   TH2F* mu_pi_length = new TH2F("mu_pi_length",";Muon length [cm];Pion length [cm]", 15, 0, 350, 15, 0, 100);
   TH2F* mu_mom_length = new TH2F("mu_mom_length", "Momentum [GeV/c];Length [cm]", 20, 0, 1, 20, 0, 350);
   TH2F* pi_mom_length = new TH2F("pi_mom_length", "Momentum [GeV/c];Length [cm]", 20, 0, 1, 20, 0, 100);
*/

   const int nentries = t -> GetEntries();

   int nsignal = 0;
   int nsignal_contained = 0;
   int longermuon = 0;

   for (int i = 0; i < nentries; i++) {

      t -> GetEntry(i);
/*
      if(nu_isCC -> at(0) && nu_PDG -> at(0) == 14 && inFV(nu_vtxx -> at(0), nu_vtxy -> at(0), nu_vtxz -> at(0)) && std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 13) == 1 && std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 211) == 1 && MCP_PDG -> size() == 2 + std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 2112) + std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 2212) + std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 2000000101)) {
         nsignal++;
         double muon_mom = 0;
         double muon_len = 0;
         double pion_mom = 0;
         double pion_len = 0;
         for(int j = 0; j < MCP_PDG -> size(); j++) {
            if (MCP_PDG -> at(j) == 13) {
               if(MCP_isContained -> at(j)) {
                  muom_mom = MCP_P -> at(j);
                  muon_len = MCP_length -> at(j);
               }
            }
            else if(MCP_PDG -> at(j) == 211) {
               if(MCP_isContained -> at(j)) {
                  pion_mom = MCP_P -> at(j);
                  pion_len = MCP_length -> at(j);
               }
            }
         }
         if(muon_mom != 0) {
            mu_mom_length -> Fill(muon_mom, muon_len);
         }
         if(pion_mom != 0) {
            pi_mom_length -> Fill(pion_mom, pion_len);
         }
         if(muon_mom != 0 && pion_mom != 0) {
            mu_pi_length -> Fill(muon_len, pion_len);
            nsignal_contained++;
            if(muon_len > pion_len) longermuon++;
         }
      }
*/
      bool SelectedEvent = true;
      for(std::map<std::string,bool>::const_iterator iter = SelectionCutflow.begin(); iter != SelectionCutflow.end(); ++iter) {
         if((CC1picutflow -> find(iter -> first)) -> second != iter -> second) {
            SelectedEvent = false;
            break;
         }
      }

      bool isSignal = false;

      // CC; muon neutrino; vtx in FV; with exactly 1 muon, exactly 1 pi+, any nucleons/Ar-40, and nothing else
      if(nu_isCC -> at(0) && nu_PDG -> at(0) == 14 && inFV(nu_vtxx -> at(0), nu_vtxy -> at(0), nu_vtxz -> at(0)) && std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 13) == 1 && std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 211) == 1 && MCP_PDG -> size() == 2 + std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 2112) + std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 2212) + std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 2000000101)) {

         isSignal = true;

         //Efficiency...
         eff_nuE -> Fill(SelectedEvent, nu_E -> at(0));

         for (int j = 0; j < MCP_PDG -> size(); j++) {
            if (MCP_PDG -> at(j) == 13) {
               eff_muP -> Fill(SelectedEvent, MCP_P -> at(j));
               TVector3 muP(MCP_Px -> at(j), MCP_Py -> at(j), MCP_Pz -> at(j));
               eff_muCosT -> Fill(SelectedEvent, muP.CosTheta());
               eff_muPhi -> Fill(SelectedEvent, muP.Phi());
               break;
            }
         }
      }

      if(SelectedEvent) {
         //Purity...
         pur_nuE -> Fill(isSignal, nu_E -> at(0));

         for (int j = 0; j < MCP_PDG -> size(); j++) {
            if (MCP_PDG -> at(j) == 13) {
               pur_muP -> Fill(isSignal, MCP_P -> at(j));
               TVector3 muP(MCP_Px -> at(j), MCP_Py -> at(j), MCP_Pz -> at(j));
               pur_muCosT -> Fill(isSignal, muP.CosTheta());
               pur_muPhi -> Fill(isSignal, muP.Phi());
               break;
            }
         }

         //THStacks...
         double maxlen = *std::max_element(track_length -> begin(), track_length -> end());

         int shower_daughters = 0;
         int track_daughters = 0;
         int mip_daughters = 0;
         int track_counter = 0;
	 std::vector<double> vtxtrack_daughters = {0.};
	 std::vector<double> vtxshwr_daughters = {0.}; 

         for(int j = 0; j < NPFPs; j++) {
	   if(Sel_PFP_isShower -> at(j) && Sel_PFP_isDaughter -> at(j)){
	     shower_daughters++;
	     vtxshwr_daughters.push_back(0); // Add things in the tree!
	   }
            if(Sel_PFP_isTrack -> at(j)) {
               if (Sel_PFP_isDaughter -> at(j)) {
                  track_daughters++;
                  if(MIPConsistency -> at(track_counter)) mip_daughters++;
		  vtxtrack_daughters.push_back(0); // Add things in the tree!
               }
               track_counter++;
            }
         }

         if(tpcobj_origin==1) {
            len_cosmic -> Fill(maxlen);
            tracks_cosmic -> Fill(track_daughters);
            showers_cosmic -> Fill(shower_daughters);
            mips_cosmic -> Fill(mip_daughters);
	    for (int i_tr=0; i_tr < vtxtrack_daughters.size(); i_tr++){
	      vtxtrack_cosmic -> Fill(vtxtrack_daughters.at(i_tr));
	    }
	    for (int i_sh=0; i_sh < vtxshwr_daughters.size(); i_sh++){
	      vtxshwr_cosmic -> Fill(vtxshwr_daughters.at(i_sh));
	    }
         }
         else if(tpcobj_origin==2) {
            len_mixed -> Fill(maxlen);
            tracks_mixed -> Fill(track_daughters);
            showers_mixed -> Fill(shower_daughters);
            mips_mixed -> Fill(mip_daughters);
	    for (int i_tr=0; i_tr < vtxtrack_daughters.size(); i_tr++){
	      vtxtrack_mixed -> Fill(vtxtrack_daughters.at(i_tr));
	    }
	    for (int i_sh=0; i_sh < vtxshwr_daughters.size(); i_sh++){
	      vtxshwr_mixed -> Fill(vtxshwr_daughters.at(i_sh));
	    }
         }
         else if(tpcobj_origin!=0) {
            len_unknown -> Fill(maxlen);
            tracks_unknown -> Fill(track_daughters);
            showers_unknown -> Fill(shower_daughters);
            mips_unknown -> Fill(mip_daughters);
	    for (int i_tr=0; i_tr < vtxtrack_daughters.size(); i_tr++){
	      vtxtrack_unknown -> Fill(vtxtrack_daughters.at(i_tr));
	    }
	    for (int i_sh=0; i_sh < vtxshwr_daughters.size(); i_sh++){
	      vtxshwr_unknown -> Fill(vtxshwr_daughters.at(i_sh));
	    }
         }
         else if(!inFV(nu_vtxx -> at(0), nu_vtxy -> at(0), nu_vtxz -> at(0))) {
            len_outFV -> Fill(maxlen);
            tracks_outFV -> Fill(track_daughters);
            showers_outFV -> Fill(shower_daughters);
            mips_outFV -> Fill(mip_daughters);
	    for (int i_tr=0; i_tr < vtxtrack_daughters.size(); i_tr++){
	      vtxtrack_outFV -> Fill(vtxtrack_daughters.at(i_tr));
	    }
	    for (int i_sh=0; i_sh < vtxshwr_daughters.size(); i_sh++){
	      vtxshwr_outFV -> Fill(vtxshwr_daughters.at(i_sh));
	    }
         }
         else if(tpcobj_origin==0) {
	   if(nu_isCC -> at(0) == 0) {
	     len_NC -> Fill(maxlen);
	     tracks_NC -> Fill(track_daughters);
	     showers_NC -> Fill(shower_daughters);
	     mips_NC -> Fill(mip_daughters);
	     for (int i_tr=0; i_tr < vtxtrack_daughters.size(); i_tr++){
	       vtxtrack_NC -> Fill(vtxtrack_daughters.at(i_tr));
	     }
	     for (int i_sh=0; i_sh < vtxshwr_daughters.size(); i_sh++){
	       vtxshwr_NC -> Fill(vtxshwr_daughters.at(i_sh));
	     }
	   }
	   else if(std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 13) == 1 && std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 211) == 1 && MCP_PDG -> size() == 2 + std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 2112) + std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 2212) + std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 2000000101)) {
               if(nu_PDG -> at(0) == 14) {
                  len_muCC1pip -> Fill(maxlen);
                  tracks_muCC1pip -> Fill(track_daughters);
                  showers_muCC1pip -> Fill(shower_daughters);
                  mips_muCC1pip -> Fill(mip_daughters);
		  for (int i_tr=0; i_tr < vtxtrack_daughters.size(); i_tr++){
		    vtxtrack_muCC1pip -> Fill(vtxtrack_daughters.at(i_tr));
		  }
		  for (int i_sh=0; i_sh < vtxshwr_daughters.size(); i_sh++){
		    vtxshwr_muCC1pip -> Fill(vtxshwr_daughters.at(i_sh));
		  }
               }
               else {
                  len_otherCC1pip -> Fill(maxlen);
                  tracks_otherCC1pip -> Fill(track_daughters);
                  showers_otherCC1pip -> Fill(shower_daughters);
                  mips_otherCC1pip -> Fill(mip_daughters);
		  for (int i_tr=0; i_tr < vtxtrack_daughters.size(); i_tr++){
		    vtxtrack_otherCC1pip -> Fill(vtxtrack_daughters.at(i_tr));
		  }
		  for (int i_sh=0; i_sh < vtxshwr_daughters.size(); i_sh++){
		    vtxshwr_otherCC1pip -> Fill(vtxshwr_daughters.at(i_sh));
		  }
               }
            }
            else {
               len_CCother -> Fill(maxlen);
               tracks_CCother -> Fill(track_daughters);
               showers_CCother -> Fill(shower_daughters);
               mips_CCother -> Fill(mip_daughters);
	       for (int i_tr=0; i_tr < vtxtrack_daughters.size(); i_tr++){
		 vtxtrack_CCother -> Fill(vtxtrack_daughters.at(i_tr));
	       }
	       for (int i_sh=0; i_sh < vtxshwr_daughters.size(); i_sh++){
		 vtxshwr_CCother -> Fill(vtxshwr_daughters.at(i_sh));
	       }
            }
         }

         // Track/shower classification...
         // Subtracting particle masses in order to plot kinetic energy
         for(int j = 0; j < Sel_MCP_PDG -> size(); j++) {
            if (Sel_MCP_PDG -> at(j) == 13) {
               if(Sel_PFP_isTrack -> at(j)) {
                  eff_track_muon -> Fill(true, Sel_MCP_E -> at(j) - 0.105658);
               }
               else if(Sel_PFP_isShower -> at(j)) {
                  eff_track_muon -> Fill(false, Sel_MCP_E -> at(j) - 0.105658);
               }
            }
            else if(Sel_MCP_PDG -> at(j) == 211) {
               if(Sel_PFP_isTrack -> at(j)) {
                  eff_track_pion -> Fill(true, Sel_MCP_E -> at(j) - 0.139570);
               }
               else if(Sel_PFP_isShower -> at(j)) {
                  eff_track_pion -> Fill(false, Sel_MCP_E -> at(j) - 0.139570);
               }
            }
            else if(Sel_MCP_PDG -> at(j) == 2212) {
               if(Sel_PFP_isTrack -> at(j)) {
                  eff_track_proton -> Fill(true, Sel_MCP_E -> at(j) - 0.938272);
               }
               else if(Sel_PFP_isShower -> at(j)) {
                  eff_track_proton -> Fill(false, Sel_MCP_E -> at(j) - 0.938272);
               }
            }

	    // Look at what particles are being associated with showers
	    if (Sel_PFP_isShower -> at(j)){
	      // If true CC1pip
	      if(isSignal) {
		if (Sel_MCP_PDG -> at(j) == 13){ // muminus
		  shwrtrue_muminus -> Fill(shower_daughters,1.0/shower_daughters);
		}
		else if (Sel_MCP_PDG -> at(j) == -13){ // muplus
		  shwrtrue_muplus -> Fill(shower_daughters,1.0/shower_daughters);
		}
		else if (Sel_MCP_PDG -> at(j) == 211){ // piplus
		  shwrtrue_piplus -> Fill(shower_daughters,1.0/shower_daughters);
		}
		else if (Sel_MCP_PDG -> at(j) == -211){ // piminus
		  shwrtrue_piminus -> Fill(shower_daughters,1.0/shower_daughters);
		}
		else if (Sel_MCP_PDG -> at(j) == 111){ // pizero
		  shwrtrue_pizero -> Fill(shower_daughters,1.0/shower_daughters);
		}
		else if (Sel_MCP_PDG -> at(j) == 2212){ // proton
		  shwrtrue_proton -> Fill(shower_daughters,1.0/shower_daughters);
		}
		else if (Sel_MCP_PDG -> at(j) == 11){ // electron
		  shwrtrue_electron -> Fill(shower_daughters,1.0/shower_daughters);
		}
		else if (Sel_MCP_PDG -> at(j) == -11){ // positron
		shwrtrue_positron -> Fill(shower_daughters,1.0/shower_daughters);
		}
		else if (Sel_MCP_PDG -> at(j) == 22){ // photon
		  shwrtrue_photon -> Fill(shower_daughters,1.0/shower_daughters);
		}
		else if (Sel_MCP_PDG -> at(j) == 2112){ // neutron
		  shwrtrue_neutron -> Fill(shower_daughters,1.0/shower_daughters);
		}
		else if (!(Sel_MCP_PDG -> at (j) == -9999)){ // other
		  shwrtrue_other -> Fill(shower_daughters);
		  std::cout << Sel_MCP_PDG -> at (j) << std::endl;
		}	
	      }
	      else{
		if (Sel_MCP_PDG -> at(j) == 13){ // muminus
		  shwrtrue_bg_muminus -> Fill(shower_daughters,1.0/shower_daughters);
		}
		else if (Sel_MCP_PDG -> at(j) == -13){ // muplus
		  shwrtrue_bg_muplus -> Fill(shower_daughters,1.0/shower_daughters);
		}
		else if (Sel_MCP_PDG -> at(j) == 211){ // piplus
		  shwrtrue_bg_piplus -> Fill(shower_daughters,1.0/shower_daughters);
		}
		else if (Sel_MCP_PDG -> at(j) == -211){ // piminus
		  shwrtrue_bg_piminus -> Fill(shower_daughters,1.0/shower_daughters);
		}
		else if (Sel_MCP_PDG -> at(j) == 111){ // pizero
		  shwrtrue_bg_pizero -> Fill(shower_daughters,1.0/shower_daughters);
		}
		else if (Sel_MCP_PDG -> at(j) == 2212){ // proton
		  shwrtrue_bg_proton -> Fill(shower_daughters,1.0/shower_daughters);
		}
		else if (Sel_MCP_PDG -> at(j) == 11){ // electron
		  shwrtrue_bg_electron -> Fill(shower_daughters,1.0/shower_daughters);
		}
		else if (Sel_MCP_PDG -> at(j) == -11){ // positron
		  shwrtrue_bg_positron -> Fill(shower_daughters,1.0/shower_daughters);
		}
		else if (Sel_MCP_PDG -> at(j) == 22){ // photon
		  shwrtrue_bg_photon -> Fill(shower_daughters,1.0/shower_daughters);
		}
		else if (Sel_MCP_PDG -> at(j) == 2112){ // neutron
		  shwrtrue_bg_neutron -> Fill(shower_daughters,1.0/shower_daughters);
		}
		else if (!(Sel_MCP_PDG -> at (j) == -9999)){ // other
		  shwrtrue_bg_other -> Fill(shower_daughters);
		  std::cout << Sel_MCP_PDG -> at (j) << std::endl;
		}
	      }
	    } // end if (Sel_PFP_isShower -> at(j))  
	 } // end loop over Sel_MCP_PDG (j)
      }
   }

   TCanvas *c1 = new TCanvas("c1", "c1");
/*
   mu_pi_length -> Draw("COLZ");
   c1 -> SaveAs("mu_pi_length.eps");
   mu_mom_length -> Draw("COLZ");
   c1 -> SaveAs("mu_mom_length.eps");
   pi_mom_length -> Draw("COLZ");
   c1 -> SaveAs("pi_mom_length.eps");
   std::cout << "Signal events: " << nsignal << std::endl;
   std::cout << "Contained signal: " << nsignal_contained << std::endl;
   std::cout << "Muon longer: " << (1. * longermuon)/(mu_pi_length->GetEntries()) << std::endl;
*/
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

   len_unknown -> Scale(norm);
   len_mixed -> Scale(norm);
   len_cosmic -> Scale(norm);
   len_outFV -> Scale(norm);
   len_NC -> Scale(norm);
   len_CCother -> Scale(norm);
   len_otherCC1pip -> Scale(norm);
   len_muCC1pip -> Scale(norm);
//   len_stack -> Add(len_unknown);
   len_stack -> Add(len_mixed);
   len_stack -> Add(len_cosmic);
   len_stack -> Add(len_outFV);
   len_stack -> Add(len_NC);
   len_stack -> Add(len_CCother);
//   len_stack -> Add(len_otherCC1pip);
   len_stack -> Add(len_muCC1pip);
   len_stack -> Draw();
   auto len_legend = new TLegend(0.6,0.6,0.9,0.9);
   len_legend -> AddEntry(len_muCC1pip, "#nu_{#mu} CC1#pi^{+}", "f");
//   len_legend -> AddEntry(len_otherCC1pip, "#bar{#nu}_{#mu}, #nu_{e}, #bar{#nu}_{e} CC1#pi^{+}", "f");
   len_legend -> AddEntry(len_CCother, "CC-Other", "f");
   len_legend -> AddEntry(len_NC, "NC", "f");
   len_legend -> AddEntry(len_outFV, "Out of FV", "f");
   len_legend -> AddEntry(len_cosmic, "Cosmic", "f");
   len_legend -> AddEntry(len_mixed, "Mixed", "f");
//   len_legend -> AddEntry(len_unknown, "Unknown Origin", "f");
   len_legend -> Draw("SAME");
   c1 -> SaveAs(TString::Format("%s_THStack_len.eps",SaveString.c_str()));

   tracks_unknown -> Scale(norm);
   tracks_mixed -> Scale(norm);
   tracks_cosmic -> Scale(norm);
   tracks_outFV -> Scale(norm);
   tracks_NC -> Scale(norm);
   tracks_CCother -> Scale(norm);
   tracks_otherCC1pip -> Scale(norm);
   tracks_muCC1pip -> Scale(norm);
//   tracks_stack -> Add(tracks_unknown);
   tracks_stack -> Add(tracks_mixed);
   tracks_stack -> Add(tracks_cosmic);
   tracks_stack -> Add(tracks_outFV);
   tracks_stack -> Add(tracks_NC);
   tracks_stack -> Add(tracks_CCother);
//   tracks_stack -> Add(tracks_otherCC1pip);
   tracks_stack -> Add(tracks_muCC1pip);
   tracks_stack -> Draw();
   auto tracks_legend = new TLegend(0.6,0.6,0.9,0.9);
   tracks_legend -> AddEntry(tracks_muCC1pip, "#nu_{#mu} CC1#pi^{+}", "f");
//   tracks_legend -> AddEntry(tracks_otherCC1pip, "#bar{#nu}_{#mu}, #nu_{e}, #bar{#nu}_{e} CC1#pi^{+}", "f");
   tracks_legend -> AddEntry(tracks_CCother, "CC-Other", "f");
   tracks_legend -> AddEntry(tracks_NC, "NC", "f");
   tracks_legend -> AddEntry(tracks_outFV, "Out of FV", "f");
   tracks_legend -> AddEntry(tracks_cosmic, "Cosmic", "f");
   tracks_legend -> AddEntry(tracks_mixed, "Mixed", "f");
//   tracks_legend -> AddEntry(tracks_unknown, "Unknown Origin", "f");
   tracks_legend -> Draw("SAME");
   c1 -> SaveAs(TString::Format("%s_THStack_tracks.eps",SaveString.c_str()));

   showers_unknown -> Scale(norm);
   showers_mixed -> Scale(norm);
   showers_cosmic -> Scale(norm);
   showers_outFV -> Scale(norm);
   showers_NC -> Scale(norm);
   showers_CCother -> Scale(norm);
   showers_otherCC1pip -> Scale(norm);
   showers_muCC1pip -> Scale(norm);
//   showers_stack -> Add(showers_unknown);
   showers_stack -> Add(showers_mixed);
   showers_stack -> Add(showers_cosmic);
   showers_stack -> Add(showers_outFV);
   showers_stack -> Add(showers_NC);
   showers_stack -> Add(showers_CCother);
//   showers_stack -> Add(showers_otherCC1pip);
   showers_stack -> Add(showers_muCC1pip);
   showers_stack -> Draw();
   auto showers_legend = new TLegend(0.6,0.6,0.9,0.9);
   showers_legend -> AddEntry(showers_muCC1pip, "#nu_{#mu} CC1#pi^{+}", "f");
//   showers_legend -> AddEntry(showers_otherCC1pip, "#bar{#nu}_{#mu}, #nu_{e}, #bar{#nu}_{e} CC1#pi^{+}", "f");
   showers_legend -> AddEntry(showers_CCother, "CC-Other", "f");
   showers_legend -> AddEntry(showers_NC, "NC", "f");
   showers_legend -> AddEntry(showers_outFV, "Out of FV", "f");
   showers_legend -> AddEntry(showers_cosmic, "Cosmic", "f");
   showers_legend -> AddEntry(showers_mixed, "Mixed", "f");
//   showers_legend -> AddEntry(showers_unknown, "Unknown Origin", "f");
   showers_legend -> Draw("SAME");
   c1 -> SaveAs(TString::Format("%s_THStack_showers.eps",SaveString.c_str()));

   mips_unknown -> Scale(norm);
   mips_mixed -> Scale(norm);
   mips_cosmic -> Scale(norm);
   mips_outFV -> Scale(norm);
   mips_NC -> Scale(norm);
   mips_CCother -> Scale(norm);
   mips_otherCC1pip -> Scale(norm);
   mips_muCC1pip -> Scale(norm);
//   mips_stack -> Add(mips_unknown);
   mips_stack -> Add(mips_mixed);
   mips_stack -> Add(mips_cosmic);
   mips_stack -> Add(mips_outFV);
   mips_stack -> Add(mips_NC);
   mips_stack -> Add(mips_CCother);
//   mips_stack -> Add(mips_otherCC1pip);
   mips_stack -> Add(mips_muCC1pip);
   mips_stack -> Draw();
   auto mips_legend = new TLegend(0.6,0.6,0.9,0.9);
   mips_legend -> AddEntry(mips_muCC1pip, "#nu_{#mu} CC1#pi^{+}", "f");
//   mips_legend -> AddEntry(mips_otherCC1pip, "#bar{#nu}_{#mu}, #nu_{e}, #bar{#nu}_{e} CC1#pi^{+}", "f");
   mips_legend -> AddEntry(mips_CCother, "CC-Other", "f");
   mips_legend -> AddEntry(mips_NC, "NC", "f");
   mips_legend -> AddEntry(mips_outFV, "Out of FV", "f");
   mips_legend -> AddEntry(mips_cosmic, "Cosmic", "f");
   mips_legend -> AddEntry(mips_mixed, "Mixed", "f");
//   mips_legend -> AddEntry(mips_unknown, "Unknown Origin", "f");
   mips_legend -> Draw("SAME");
   c1 -> SaveAs(TString::Format("%s_THStack_mips.eps",SaveString.c_str()));

   vtxtrack_unknown -> Scale(norm);
   vtxtrack_mixed -> Scale(norm);
   vtxtrack_cosmic -> Scale(norm);
   vtxtrack_outFV -> Scale(norm);
   vtxtrack_NC -> Scale(norm);
   vtxtrack_CCother -> Scale(norm);
   vtxtrack_otherCC1pip -> Scale(norm);
   vtxtrack_muCC1pip -> Scale(norm);
//   vtxtrack_stack -> Add(vtxtrack_unknown);
   vtxtrack_stack -> Add(vtxtrack_mixed);
   vtxtrack_stack -> Add(vtxtrack_cosmic);
   vtxtrack_stack -> Add(vtxtrack_outFV);
   vtxtrack_stack -> Add(vtxtrack_NC);
   vtxtrack_stack -> Add(vtxtrack_CCother);
//   vtxtrack_stack -> Add(vtxtrack_otherCC1pip);
   vtxtrack_stack -> Add(vtxtrack_muCC1pip);
   vtxtrack_stack -> Draw();
   auto vtxtrack_legend = new TLegend(0.6,0.6,0.9,0.9);
   vtxtrack_legend -> AddEntry(vtxtrack_muCC1pip, "#nu_{#mu} CC1#pi^{+}", "f");
//   vtxtrack_legend -> AddEntry(vtxtrack_otherCC1pip, "#bar{#nu}_{#mu}, #nu_{e}, #bar{#nu}_{e} CC1#pi^{+}", "f");
   vtxtrack_legend -> AddEntry(vtxtrack_CCother, "CC-Other", "f");
   vtxtrack_legend -> AddEntry(vtxtrack_NC, "NC", "f");
   vtxtrack_legend -> AddEntry(vtxtrack_outFV, "Out of FV", "f");
   vtxtrack_legend -> AddEntry(vtxtrack_cosmic, "Cosmic", "f");
   vtxtrack_legend -> AddEntry(vtxtrack_mixed, "Mixed", "f");
//   vtxtrack_legend -> AddEntry(vtxtrack_unknown, "Unknown Origin", "f");
   vtxtrack_legend -> Draw("SAME");
   c1 -> SaveAs(TString::Format("%s_THStack_vtxtrack.eps",SaveString.c_str()));

   vtxshwr_unknown -> Scale(norm);
   vtxshwr_mixed -> Scale(norm);
   vtxshwr_cosmic -> Scale(norm);
   vtxshwr_outFV -> Scale(norm);
   vtxshwr_NC -> Scale(norm);
   vtxshwr_CCother -> Scale(norm);
   vtxshwr_otherCC1pip -> Scale(norm);
   vtxshwr_muCC1pip -> Scale(norm);
//   vtxshwr_stack -> Add(vtxshwr_unknown);
   vtxshwr_stack -> Add(vtxshwr_mixed);
   vtxshwr_stack -> Add(vtxshwr_cosmic);
   vtxshwr_stack -> Add(vtxshwr_outFV);
   vtxshwr_stack -> Add(vtxshwr_NC);
   vtxshwr_stack -> Add(vtxshwr_CCother);
//   vtxshwr_stack -> Add(vtxshwr_otherCC1pip);
   vtxshwr_stack -> Add(vtxshwr_muCC1pip);
   vtxshwr_stack -> Draw();
   auto vtxshwr_legend = new TLegend(0.6,0.6,0.9,0.9);
   vtxshwr_legend -> AddEntry(vtxshwr_muCC1pip, "#nu_{#mu} CC1#pi^{+}", "f");
//   vtxshwr_legend -> AddEntry(vtxshwr_otherCC1pip, "#bar{#nu}_{#mu}, #nu_{e}, #bar{#nu}_{e} CC1#pi^{+}", "f");
   vtxshwr_legend -> AddEntry(vtxshwr_CCother, "CC-Other", "f");
   vtxshwr_legend -> AddEntry(vtxshwr_NC, "NC", "f");
   vtxshwr_legend -> AddEntry(vtxshwr_outFV, "Out of FV", "f");
   vtxshwr_legend -> AddEntry(vtxshwr_cosmic, "Cosmic", "f");
   vtxshwr_legend -> AddEntry(vtxshwr_mixed, "Mixed", "f");
//   vtxshwr_legend -> AddEntry(vtxshwr_unknown, "Unknown Origin", "f");
   vtxshwr_legend -> Draw("SAME");
   c1 -> SaveAs(TString::Format("%s_THStack_vtxshwr.eps",SaveString.c_str()));
   
   shwrtrue_muplus -> Scale(norm);
   shwrtrue_muminus -> Scale(norm);
   shwrtrue_piplus -> Scale(norm);
   shwrtrue_piminus -> Scale(norm);
   shwrtrue_pizero -> Scale(norm);
   shwrtrue_proton -> Scale(norm);
   shwrtrue_electron -> Scale(norm);
   shwrtrue_positron -> Scale(norm);
   shwrtrue_neutron -> Scale(norm);
   shwrtrue_photon -> Scale(norm);
   shwrtrue_other -> Scale(norm);
   shwrtrue_stack -> Add(shwrtrue_muplus);
   shwrtrue_stack -> Add(shwrtrue_muminus);
   shwrtrue_stack -> Add(shwrtrue_piplus);
   shwrtrue_stack -> Add(shwrtrue_piminus);
   shwrtrue_stack -> Add(shwrtrue_pizero);
   shwrtrue_stack -> Add(shwrtrue_proton);
   shwrtrue_stack -> Add(shwrtrue_electron);
   shwrtrue_stack -> Add(shwrtrue_positron);
   shwrtrue_stack -> Add(shwrtrue_neutron);
   shwrtrue_stack -> Add(shwrtrue_photon);
   shwrtrue_stack -> Add(shwrtrue_other);
   shwrtrue_stack -> Draw("hist");
   auto shwrtrue_legend = new TLegend(0.6,0.6,0.9,0.9);
   shwrtrue_legend -> AddEntry(shwrtrue_muminus, "#mu^{-}", "f");
   shwrtrue_legend -> AddEntry(shwrtrue_muplus, "#mu^{+}", "f");
   shwrtrue_legend -> AddEntry(shwrtrue_piplus, "#pi^{+}", "f");
   shwrtrue_legend -> AddEntry(shwrtrue_piminus, "#pi^{-}", "f");
   shwrtrue_legend -> AddEntry(shwrtrue_pizero, "#pi^{0}", "f");
   shwrtrue_legend -> AddEntry(shwrtrue_proton, "p", "f");
   shwrtrue_legend -> AddEntry(shwrtrue_electron, "e^{-}", "f");
   shwrtrue_legend -> AddEntry(shwrtrue_positron, "e^{+}", "f");
   shwrtrue_legend -> AddEntry(shwrtrue_neutron, "n", "f");
   shwrtrue_legend -> AddEntry(shwrtrue_photon, "#gamma", "f");
   shwrtrue_legend -> AddEntry(shwrtrue_other, "other", "f");
   shwrtrue_legend -> Draw("SAME");
   c1 -> SaveAs(TString::Format("%s_THStack_signal_shwrtrue.eps",SaveString.c_str()));
   
   shwrtrue_bg_muplus -> Scale(norm);
   shwrtrue_bg_muminus -> Scale(norm);
   shwrtrue_bg_piplus -> Scale(norm);
   shwrtrue_bg_piminus -> Scale(norm);
   shwrtrue_bg_pizero -> Scale(norm);
   shwrtrue_bg_proton -> Scale(norm);
   shwrtrue_bg_electron -> Scale(norm);
   shwrtrue_bg_positron -> Scale(norm);
   shwrtrue_bg_neutron -> Scale(norm);
   shwrtrue_bg_photon -> Scale(norm);
   shwrtrue_bg_other -> Scale(norm);
   shwrtrue_bg_stack -> Add(shwrtrue_bg_muplus);
   shwrtrue_bg_stack -> Add(shwrtrue_bg_muminus);
   shwrtrue_bg_stack -> Add(shwrtrue_bg_piplus);
   shwrtrue_bg_stack -> Add(shwrtrue_bg_piminus);
   shwrtrue_bg_stack -> Add(shwrtrue_bg_pizero);
   shwrtrue_bg_stack -> Add(shwrtrue_bg_proton);
   shwrtrue_bg_stack -> Add(shwrtrue_bg_electron);
   shwrtrue_bg_stack -> Add(shwrtrue_bg_positron);
   shwrtrue_bg_stack -> Add(shwrtrue_bg_neutron);
   shwrtrue_bg_stack -> Add(shwrtrue_bg_photon);
   shwrtrue_bg_stack -> Add(shwrtrue_bg_other);
   shwrtrue_bg_stack -> Draw("hist");
   auto shwrtrue_bg_legend = new TLegend(0.6,0.6,0.9,0.9);
   shwrtrue_bg_legend -> AddEntry(shwrtrue_bg_muminus, "#mu^{-}", "f");
   shwrtrue_bg_legend -> AddEntry(shwrtrue_bg_muplus, "#mu^{+}", "f");
   shwrtrue_bg_legend -> AddEntry(shwrtrue_bg_piplus, "#pi^{+}", "f");
   shwrtrue_bg_legend -> AddEntry(shwrtrue_bg_piminus, "#pi^{-}", "f");
   shwrtrue_bg_legend -> AddEntry(shwrtrue_bg_pizero, "#pi^{0}", "f");
   shwrtrue_bg_legend -> AddEntry(shwrtrue_bg_proton, "p", "f");
   shwrtrue_bg_legend -> AddEntry(shwrtrue_bg_electron, "e^{-}", "f");
   shwrtrue_bg_legend -> AddEntry(shwrtrue_bg_positron, "e^{+}", "f");
   shwrtrue_bg_legend -> AddEntry(shwrtrue_bg_neutron, "n", "f");
   shwrtrue_bg_legend -> AddEntry(shwrtrue_bg_photon, "#gamma", "f");
   shwrtrue_bg_legend -> AddEntry(shwrtrue_bg_other, "other", "f");
   shwrtrue_bg_legend -> Draw("SAME");
   c1 -> SaveAs(TString::Format("%s_THStack_bg_shwrtrue.eps",SaveString.c_str()));

   eff_track_muon -> Draw("AP");
   eff_track_pion -> Draw("P SAME");
   eff_track_proton -> Draw("P SAME");
   gPad->Update();
   eff_track_muon->SetMarkerColor(kRed);
   eff_track_muon->SetMarkerStyle(21);
   eff_track_muon->SetMarkerSize(1);
   graph = eff_track_muon->GetPaintedGraph(); 
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   eff_track_pion->SetMarkerColor(kBlue);
   eff_track_pion->SetMarkerStyle(21);
   eff_track_pion->SetMarkerSize(1);
   graph = eff_track_pion->GetPaintedGraph(); 
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   eff_track_proton->SetMarkerColor(kGreen);
   eff_track_proton->SetMarkerStyle(21);
   eff_track_proton->SetMarkerSize(1);
   graph = eff_track_proton->GetPaintedGraph(); 
   graph->SetMinimum(0);
   graph->SetMaximum(1);
   gPad->Update();
   c1 -> SaveAs(TString::Format("%s_eff_track.eps",SaveString.c_str()));

   double TotalEntries = 1. * (((TH1*)(len_stack -> GetStack() -> Last())) -> GetEntries());
   std::cout << "muCC1pip: " << len_muCC1pip -> GetEntries() / TotalEntries << std::endl;
//   std::cout << "otherCC1pip: " << len_otherCC1pip -> GetEntries() / TotalEntries << std::endl; 
   std::cout << "CCother: " << len_CCother -> GetEntries() / TotalEntries << std::endl;
   std::cout << "NC: " << len_NC -> GetEntries() / TotalEntries << std::endl;
   std::cout << "outFV: " << len_outFV -> GetEntries() / TotalEntries << std::endl;
   std::cout << "cosmic: " << len_cosmic -> GetEntries() / TotalEntries << std::endl;
   std::cout << "mixed: " << len_mixed -> GetEntries() / TotalEntries << std::endl;
//   std::cout << "unknown: " << len_unknown -> GetEntries() / TotalEntries << std::endl;

   std::cout << "Total Efficiency: " << (1. * (eff_nuE -> GetPassedHistogram() -> GetEntries()))/(eff_nuE -> GetTotalHistogram() -> GetEntries()) << std::endl;
   std::cout << "Total Purity: " << (1. * (pur_nuE -> GetPassedHistogram() -> GetEntries()))/(pur_nuE -> GetTotalHistogram() -> GetEntries()) << std::endl;

   std::cout << "Muon Track Efficiency: " << (1. * (eff_track_muon -> GetPassedHistogram() -> GetEntries()))/(eff_track_muon -> GetTotalHistogram() -> GetEntries()) << std::endl;
   std::cout << "Pion Track Efficiency: " << (1. * (eff_track_pion -> GetPassedHistogram() -> GetEntries()))/(eff_track_pion -> GetTotalHistogram() -> GetEntries()) << std::endl;
   std::cout << "Proton Track Efficiency: " << (1. * (eff_track_proton -> GetPassedHistogram() -> GetEntries()))/(eff_track_proton -> GetTotalHistogram() -> GetEntries()) << std::endl;

   gDirectory->GetList()->Delete();

   return;
}
