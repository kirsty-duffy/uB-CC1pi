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

   bool isData;
   unsigned int run_num;
   unsigned int subrun_num;
   unsigned int event_num;

   std::map<std::string,bool> *Marco_cutflow = nullptr;
   bool Marco_selected = 0;

   NuIntTopology Truth_topology;
   std::vector<double> *TPCObj_PFP_track_length = nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_track_start = nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_track_end = nullptr;
   std::vector<double> *TPCObj_PFP_track_dqdx_truncmean = nullptr;
   std::vector<bool> *TPCObj_PFP_isMIP = nullptr;
   std::vector<double> *TPCObj_PFP_shower_length = nullptr;
   std::vector<std::vector<double>> *TPCObj_PFP_shower_start = nullptr;
   int TPCObj_NPFPs = 0;
   int TPCObj_NTracks = 0;
   int TPCObj_NShowers = 0;
   std::vector<bool> *TPCObj_PFP_isTrack = nullptr;
   std::vector<bool> *TPCObj_PFP_isShower = nullptr;
   std::vector<bool> *TPCObj_PFP_isDaughter = nullptr;
   std::vector<int> *TPCObj_PFP_id = nullptr;
   std::vector<int> *TPCObj_PFP_MCPid = nullptr;
   std::vector<int> *TPCObj_PFP_truePDG = nullptr;
   std::vector<double> *TPCObj_PFP_trueE = nullptr;
   std::vector<double> *TPCObj_PFP_trueKE = nullptr;
   int TPCObj_origin = 0;
   int TPCObj_origin_extra = 0;
   std::vector<double> *TPCObj_reco_vtx = nullptr;

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

   std::vector<double> *nu_vtx = nullptr;
   std::vector<double> *nu_vtx_spacecharge = nullptr;
   bool nu_isCC = 0;
   int nu_PDG = 0;
   double nu_E = 0;

   std::map<std::string,bool> *CC1picutflow = nullptr;

   t -> SetBranchStatus("*",0);
   // t -> SetBranchStatus("isData",1);
   // t -> SetBranchAddress("isData", &isData);
   // t -> SetBranchStatus("run_num",1);
   // t -> SetBranchAddress("run_num", &run_num);
   // t -> SetBranchStatus("subrun_num",1);
   // t -> SetBranchAddress("subrun_num", &subrun_num);
   // t -> SetBranchStatus("event_num",1);
   // t -> SetBranchAddress("event_num", &event_num);
   //
   // t -> SetBranchStatus("Marco_cutflow",1);
   // t -> SetBranchAddress("Marco_cutflow", &Marco_cutflow);
   // t -> SetBranchStatus("Marco_selected",1);
   // t -> SetBranchAddress("Marco_selected", &Marco_selected);
   //
   // t -> SetBranchStatus("Truth_topology",1);
   // t -> SetBranchAddress("Truth_topology", &Truth_topology);
   // t -> SetBranchStatus("TPCObj_PFP_track_length",1);
   // t -> SetBranchAddress("TPCObj_PFP_track_length", &TPCObj_PFP_track_length);
   // t -> SetBranchStatus("TPCObj_PFP_track_start",1);
   // t -> SetBranchAddress("TPCObj_PFP_track_start", &TPCObj_PFP_track_start);
   // t -> SetBranchStatus("TPCObj_PFP_track_end",1);
   // t -> SetBranchAddress("TPCObj_PFP_track_end", &TPCObj_PFP_track_end);
   // t -> SetBranchStatus("TPCObj_PFP_track_dqdx_truncmean",1);
   // t -> SetBranchAddress("TPCObj_PFP_track_dqdx_truncmean", &TPCObj_PFP_track_dqdx_truncmean);
   // t -> SetBranchStatus("TPCObj_PFP_isMIP",1);
   // t -> SetBranchAddress("TPCObj_PFP_isMIP", &TPCObj_PFP_isMIP);
   // t -> SetBranchStatus("TPCObj_PFP_shower_length",1);
   // t -> SetBranchAddress("TPCObj_PFP_shower_length", &TPCObj_PFP_shower_length);
   // t -> SetBranchStatus("TPCObj_PFP_shower_start",1);
   // t -> SetBranchAddress("TPCObj_PFP_shower_start", &TPCObj_PFP_shower_start);
   // t -> SetBranchStatus("TPCObj_NPFPs",1);
   // t -> SetBranchAddress("TPCObj_NPFPs", &TPCObj_NPFPs);
   // t -> SetBranchStatus("TPCObj_NTracks",1);
   // t -> SetBranchAddress("TPCObj_NTracks", &TPCObj_NTracks);
   // t -> SetBranchStatus("TPCObj_NShowers",1);
   // t -> SetBranchAddress("TPCObj_NShowers", &TPCObj_NShowers);
   // t -> SetBranchStatus("TPCObj_PFP_isTrack",1);
   // t -> SetBranchAddress("TPCObj_PFP_isTrack", &TPCObj_PFP_isTrack);
   // t -> SetBranchStatus("TPCObj_PFP_isShower",1);
   // t -> SetBranchAddress("TPCObj_PFP_isShower", &TPCObj_PFP_isShower);
   // t -> SetBranchStatus("TPCObj_PFP_isDaughter",1);
   // t -> SetBranchAddress("TPCObj_PFP_isDaughter", &TPCObj_PFP_isDaughter);
   // t -> SetBranchStatus("TPCObj_PFP_id",1);
   // t -> SetBranchAddress("TPCObj_PFP_id", &TPCObj_PFP_id);
   // t -> SetBranchStatus("TPCObj_PFP_MCPid",1);
   // t -> SetBranchAddress("TPCObj_PFP_MCPid", &TPCObj_PFP_MCPid);
   // t -> SetBranchStatus("TPCObj_PFP_truePDG",1);
   // t -> SetBranchAddress("TPCObj_PFP_truePDG", &TPCObj_PFP_truePDG);
   // t -> SetBranchStatus("TPCObj_PFP_trueE",1);
   // t -> SetBranchAddress("TPCObj_PFP_trueE", &TPCObj_PFP_trueE);
   // t -> SetBranchStatus("TPCObj_PFP_trueKE",1);
   // t -> SetBranchAddress("TPCObj_PFP_trueKE", &TPCObj_PFP_trueKE);
   // t -> SetBranchStatus("TPCObj_origin",1);
   // t -> SetBranchAddress("TPCObj_origin", &TPCObj_origin);
   // t -> SetBranchStatus("TPCObj_origin_extra",1);
   // t -> SetBranchAddress("TPCObj_origin_extra", &TPCObj_origin_extra);
   // t -> SetBranchStatus("TPCObj_reco_vtx",1);
   // t -> SetBranchAddress("TPCObj_reco_vtx", &TPCObj_reco_vtx);
   //
   // t -> SetBranchStatus("MCP_PDG",1);
   // t -> SetBranchAddress("MCP_PDG", &MCP_PDG);
   // t -> SetBranchStatus("MCP_length",1);
   // t -> SetBranchAddress("MCP_length", &MCP_length);
   // t -> SetBranchStatus("MCP_process",1);
   // t -> SetBranchAddress("MCP_process", &MCP_process);
   // t -> SetBranchStatus("MCP_endprocess",1);
   // t -> SetBranchAddress("MCP_endprocess", &MCP_endprocess);
   // t -> SetBranchStatus("MCP_numdaughters",1);
   // t -> SetBranchAddress("MCP_numdaughters", &MCP_numdaughters);
   // t -> SetBranchStatus("MCP_P",1);
   // t -> SetBranchAddress("MCP_P", &MCP_P);
   // t -> SetBranchStatus("MCP_Px",1);
   // t -> SetBranchAddress("MCP_Px", &MCP_Px);
   // t -> SetBranchStatus("MCP_Py",1);
   // t -> SetBranchAddress("MCP_Py", &MCP_Py);
   // t -> SetBranchStatus("MCP_Pz",1);
   // t -> SetBranchAddress("MCP_Pz", &MCP_Pz);
   // t -> SetBranchStatus("MCP_E",1);
   // t -> SetBranchAddress("MCP_E", &MCP_E);
   // t -> SetBranchStatus("MCP_KE",1);
   // t -> SetBranchAddress("MCP_KE", &MCP_KE);
   // t -> SetBranchStatus("MCP_isContained",1);
   // t -> SetBranchAddress("MCP_isContained", &MCP_isContained);
   //
   // t -> SetBranchStatus("nu_vtx",1);
   // t -> SetBranchAddress("nu_vtx", &nu_vtx);
   // t -> SetBranchStatus("nu_vtx_spacecharge",1);
   // t -> SetBranchAddress("nu_vtx_spacecharge", &nu_vtx_spacecharge);
   // t -> SetBranchStatus("nu_isCC",1);
   // t -> SetBranchAddress("nu_isCC", &nu_isCC);
   // t -> SetBranchStatus("nu_PDG",1);
   // t -> SetBranchAddress("nu_PDG", &nu_PDG);
   // t -> SetBranchStatus("nu_E",1);
   // t -> SetBranchAddress("nu_E", &nu_E);

   t -> SetBranchStatus("CC1picutflow",1);
   t -> SetBranchAddress("CC1picutflow", &CC1picutflow);


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
   StackedHistTopology* vtxtrack_stack = new StackedHistTopology("vtxtrack_stack", ";Distance between vertex and start of #nu daughter tracks in selected TPCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)", 20, 0, 20);
   StackedHistTopology* pfps_stack = new StackedHistTopology("pfps_stack", ";Number of #nu daughter PFPs in selected TPCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)", 10, 0, 10);

   StackedHistPDGCode* len_stack_byPDG = new StackedHistPDGCode("len_stack_byPDG", ";Longest Track Length [cm];Selected Events (normalised to 3.782 #times 10^{19} POT)", 30, 0, 700);
   StackedHistPDGCode* isMIP_stack_byPDG = new StackedHistPDGCode("isMIP_stack_byPDG", ";#nu daughter track is MIP-like?;Selected Events (normalised to 3.782 #times 10^{19} POT)", 2, 0, 2);
   StackedHistPDGCode* mips_stack_byPDG = new StackedHistPDGCode("mips_stack_byPDG", ";Number of MIP-like #nu daughter tracks in selected TPCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)", 10, 0, 10);
   StackedHistPDGCode* pfps_stack_byPDG = new StackedHistPDGCode("pfps_stack_byPDG", ";Number of #nu daughter PFPs in selected TPCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)", 10, 0, 10);

   TH2D* longest2tracks_byPDG_signal = new TH2D("longest2tracks_byPDG_signal",";True PDG of longest track;True PDG of second longest track",4,0,4,4,0,4);
   TH2D* longest2tracks_byPDG_bg = new TH2D("longest2tracks_byPDG_bg",";True PDG of longest track;True PDG of second longest track",4,0,4,4,0,4);

//   TH1D* vtx_truereco_dist = new TH1D("vtx_truereco_dist",";Distance between true and reco vertices;Selected Events (normalised to 3.782 #times 10^{19} POT)",50,0,50);
   TH1D* vtx_truereco_spacecharge_dist = new TH1D("vtx_truereco_spacecharge_dist",";Distance between true and reco vertices (including spacecharge);Selected Events (normalised to 3.782 #times 10^{19} POT)",50,0,300);

   const int nentries = t -> GetEntries();

   for (int i = 0; i < nentries; i++) {
      std::cout << i << std::endl;
      t -> GetEntry(i);
      std::cout << i << std::endl;
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
      //      if (Truth_topology == kCC1piplus0p || Truth_topology == kCC1piplus1p || Truth_topology == kCC1piplusNp) {
      if(nu_isCC && nu_PDG == 14 && inFV(nu_vtx -> at(0), nu_vtx -> at(1), nu_vtx -> at(2)) && std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 13) == 1 && std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 211) == 1 && MCP_PDG -> size() == 2 + std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 2112) + std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 2212) + std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 2000000101) + std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 1000180400)) {

         isSignal = true;

         // Check neutrino flavour (just in case)
         if (nu_PDG != 14) std::cout << "WARNING: Signal event has neutrino with PDG code " << nu_PDG << " (should be 14)" << std::endl;

         // Efficiency...
         eff_nuE -> Fill(SelectedEvent, nu_E);

         for (int j = 0; j < MCP_PDG -> size(); j++) {
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

         for (int j = 0; j < MCP_PDG -> size(); j++) {
            if (MCP_PDG -> at(j) == 13) {
               pur_muP -> Fill(isSignal, MCP_P -> at(j));
               TVector3 muP(MCP_Px -> at(j), MCP_Py -> at(j), MCP_Pz -> at(j));
               pur_muCosT -> Fill(isSignal, muP.CosTheta());
               pur_muPhi -> Fill(isSignal, muP.Phi());
               break;
            }
         }

         // THStacks...

         //Find longest and second longest tracks
         double maxlen = 0, maxlen2 = 0;
         PDGCode maxPDG, maxPDG2;
         for(int j = 0; j < TPCObj_NPFPs; j++) {
            if (TPCObj_PFP_track_length -> at(j) > maxlen) {
               maxlen2 = maxlen;
               maxlen = TPCObj_PFP_track_length -> at(j);
               maxPDG2 = maxPDG;
               maxPDG = (PDGCode)TPCObj_PFP_truePDG -> at(j);
            }
            else if (TPCObj_PFP_track_length -> at(j) > maxlen2) {
               maxlen2 = TPCObj_PFP_track_length -> at(j);
               maxPDG2 = (PDGCode)TPCObj_PFP_truePDG -> at(j);
            }
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
            }
            if(TPCObj_PFP_isMIP -> at(j) && TPCObj_PFP_isDaughter -> at(j)) {
               mip_daughters++;
            }
         }

         len_stack -> Fill(Truth_topology, maxlen);
         mips_stack -> Fill(Truth_topology, mip_daughters);
         pfps_stack -> Fill(Truth_topology, track_daughters);
         len_stack_byPDG -> Fill(maxPDG, maxlen);

         for (int j = 0; j < track_daughters; j++) {
            vtxtrack_stack -> Fill(Truth_topology, vtxtrack_daughters.at(j));
         }
         for(int j = 0; j < TPCObj_NPFPs; j++) {
            if(TPCObj_PFP_isDaughter -> at(j)) {
               isMIP_stack_byPDG -> Fill((PDGCode)TPCObj_PFP_truePDG->at(j), TPCObj_PFP_isMIP->at(j));
               pfps_stack_byPDG -> Fill((PDGCode)TPCObj_PFP_truePDG->at(j), track_daughters, 1./track_daughters);
               if(mip_daughters>0) mips_stack_byPDG -> Fill((PDGCode)TPCObj_PFP_truePDG->at(j), mip_daughters, 1./mip_daughters);
            }
         }

         //Longest 2 daughter tracks by PDG...
         if(track_daughters >= 2) {
            if(isSignal) {
               if(maxPDG == kMuMinus) {
                  if(maxPDG2 == kMuMinus) longest2tracks_byPDG_signal -> Fill(0., 0);
                  else if(maxPDG2 == kPiPlus) longest2tracks_byPDG_signal -> Fill(0.,1);
                  else if(maxPDG2 == kProton) longest2tracks_byPDG_signal -> Fill(0.,2);
                  else longest2tracks_byPDG_signal -> Fill(0.,3);
               }
               else if(maxPDG == kPiPlus) {
                  if(maxPDG2 == kMuMinus) longest2tracks_byPDG_signal -> Fill(1,0);
                  else if(maxPDG2 == kPiPlus) longest2tracks_byPDG_signal -> Fill(1,1);
                  else if(maxPDG2 == kProton) longest2tracks_byPDG_signal -> Fill(1,2);
                  else longest2tracks_byPDG_signal -> Fill(1,3);
               }
               else if(maxPDG == kProton) {
                  if(maxPDG2 == kMuMinus) longest2tracks_byPDG_signal -> Fill(2,0);
                  else if(maxPDG2 == kPiPlus) longest2tracks_byPDG_signal -> Fill(2,1);
                  else if(maxPDG2 == kProton) longest2tracks_byPDG_signal -> Fill(2,2);
                  else longest2tracks_byPDG_signal -> Fill(2,3);
               }
               else {
                  if(maxPDG2 == kMuMinus) longest2tracks_byPDG_signal -> Fill(3,0);
                  else if(maxPDG2 == kPiPlus) longest2tracks_byPDG_signal -> Fill(3,1);
                  else if(maxPDG2 == kProton) longest2tracks_byPDG_signal -> Fill(3,2);
                  else longest2tracks_byPDG_signal -> Fill(3,3);
               }
            }
            else {
               if(maxPDG == kMuMinus) {
                  if(maxPDG2 == kMuMinus) longest2tracks_byPDG_bg -> Fill(0., 0);
                  else if(maxPDG2 == kPiPlus) longest2tracks_byPDG_bg -> Fill(0.,1);
                  else if(maxPDG2 == kProton) longest2tracks_byPDG_bg -> Fill(0.,2);
                  else longest2tracks_byPDG_bg -> Fill(0.,3);
               }
               else if(maxPDG == kPiPlus) {
                  if(maxPDG2 == kMuMinus) longest2tracks_byPDG_bg -> Fill(1,0);
                  else if(maxPDG2 == kPiPlus) longest2tracks_byPDG_bg -> Fill(1,1);
                  else if(maxPDG2 == kProton) longest2tracks_byPDG_bg -> Fill(1,2);
                  else longest2tracks_byPDG_bg -> Fill(1,3);
               }
               else if(maxPDG == kProton) {
                  if(maxPDG2 == kMuMinus) longest2tracks_byPDG_bg -> Fill(2,0);
                  else if(maxPDG2 == kPiPlus) longest2tracks_byPDG_bg -> Fill(2,1);
                  else if(maxPDG2 == kProton) longest2tracks_byPDG_bg -> Fill(2,2);
                  else longest2tracks_byPDG_bg -> Fill(2,3);
               }
               else {
                  if(maxPDG2 == kMuMinus) longest2tracks_byPDG_bg -> Fill(3,0);
                  else if(maxPDG2 == kPiPlus) longest2tracks_byPDG_bg -> Fill(3,1);
                  else if(maxPDG2 == kProton) longest2tracks_byPDG_bg -> Fill(3,2);
                  else longest2tracks_byPDG_bg -> Fill(3,3);
               }
            }
         } //end longest 2 daughter tracks by PDG

         if (Truth_topology == kOutFV) {
//            vtx_truereco_dist -> Fill(std::hypot(std::hypot(TPCObj_reco_vtx->at(0) - nu_vtx->at(0), TPCObj_reco_vtx->at(1) - nu_vtx->at(1)), TPCObj_reco_vtx->at(2) - nu_vtx->at(2)));
            vtx_truereco_spacecharge_dist -> Fill(std::hypot(std::hypot(TPCObj_reco_vtx->at(0) - nu_vtx_spacecharge->at(0), TPCObj_reco_vtx->at(1) - nu_vtx_spacecharge->at(1)), TPCObj_reco_vtx->at(2) - nu_vtx_spacecharge->at(2)));
         }

      } // end if(SelectedEvent)
   } //end loop over entries

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

   len_stack_byPDG -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THStack_trklen_byPDG.eps",SaveString.c_str()));

   isMIP_stack_byPDG -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THSTack_isMIP_byPDG.eps",SaveString.c_str()));

   mips_stack_byPDG -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THSTack_mips_byPDG.eps",SaveString.c_str()));

   pfps_stack_byPDG -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THSTack_pfps_byPDG.eps",SaveString.c_str()));

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

   std::cout << "Total Efficiency: " << (1. * (eff_nuE -> GetPassedHistogram() -> GetEntries()))/(eff_nuE -> GetTotalHistogram() -> GetEntries()) << std::endl;
   std::cout << "Total Purity: " << (1. * (pur_nuE -> GetPassedHistogram() -> GetEntries()))/(pur_nuE -> GetTotalHistogram() -> GetEntries()) << std::endl;

   gDirectory->GetList()->Delete();

   return;
}
