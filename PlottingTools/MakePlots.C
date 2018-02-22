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

//TPC boundary + FV cut
const double FVxmin = 0 + 12;
const double FVxmax = 256.35 - 12;
const double FVymin = -115.53 + 35;
const double FVymax = 117.47 - 35;
const double FVzmin = 0.1 + 25;
const double FVzmax = 1036.9 - 85;

//Cut out additional section of mostly dead wires in z
const double deadzmin = 675;
const double deadzmax = 775;

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

   std::map<std::string,bool> *Marco_cutflow = NULL;
   bool Marco_selected;

   NuIntTopology TPCObj_beamnu_topology;
   std::vector<double> *TPCObj_PFP_track_length = NULL;
   std::vector<std::vector<double>> *TPCObj_PFP_track_start = NULL;
   std::vector<std::vector<double>> *TPCObj_PFP_track_end = NULL;
   std::vector<double> *TPCObj_PFP_track_dqdx_truncmean = NULL;
   std::vector<bool> *TPCObj_PFP_isMIP = NULL;
   std::vector<double> *TPCObj_PFP_shower_length = NULL;
   std::vector<std::vector<double>> *TPCObj_PFP_shower_start = NULL;
   int TPCObj_NPFPs;
   int TPCObj_NTracks;
   int TPCObj_NShowers;
   std::vector<bool> *TPCObj_PFP_isTrack = NULL;
   std::vector<bool> *TPCObj_PFP_isShower = NULL;
   std::vector<bool> *TPCObj_PFP_isDaughter = NULL;
   std::vector<int> *TPCObj_PFP_id = NULL;
   std::vector<int> *TPCObj_PFP_MCPid = NULL;
   std::vector<int> *TPCObj_PFP_truePDG = NULL;
   std::vector<double> *TPCObj_PFP_trueE = NULL;
   std::vector<double> *TPCObj_PFP_trueKE = NULL;
   int TPCObj_origin;
   int TPCObj_origin_extra;
   std::vector<double> *TPCObj_reco_vtx = NULL;

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
   std::vector<double> *MCP_KE = NULL;
   std::vector<bool> *MCP_isContained = NULL;

   std::vector<double> *nu_vtxx = NULL;
   std::vector<double> *nu_vtxy = NULL;
   std::vector<double> *nu_vtxz = NULL;
   std::vector<bool> *nu_isCC = NULL;
   std::vector<int> *nu_PDG = NULL;
   std::vector<double> *nu_E = NULL;

   std::map<std::string,bool> *CC1picutflow = NULL;

   t -> SetBranchAddress("isData", &isData);
   t -> SetBranchAddress("run_num", &run_num);
   t -> SetBranchAddress("subrun_num", &subrun_num);
   t -> SetBranchAddress("event_num", &event_num);

   t -> SetBranchAddress("Marco_cutflow", &Marco_cutflow);
   t -> SetBranchAddress("Marco_selected", &Marco_selected);

   t -> SetBranchAddress("TPCObj_beamnu_topology", &TPCObj_beamnu_topology);
   t -> SetBranchAddress("TPCObj_PFP_track_length", &TPCObj_PFP_track_length);
   t -> SetBranchAddress("TPCObj_PFP_track_start", &TPCObj_PFP_track_start);
   t -> SetBranchAddress("TPCObj_PFP_track_end", &TPCObj_PFP_track_end);
   t -> SetBranchAddress("TPCObj_PFP_track_dqdx_truncmean", &TPCObj_PFP_track_dqdx_truncmean);
   t -> SetBranchAddress("TPCObj_PFP_isMIP", &TPCObj_PFP_isMIP);
   t -> SetBranchAddress("TPCObj_PFP_shower_length", &TPCObj_PFP_shower_length);
   t -> SetBranchAddress("TPCObj_PFP_shower_start", &TPCObj_PFP_shower_start);
   t -> SetBranchAddress("TPCObj_NPFPs", &TPCObj_NPFPs);
   t -> SetBranchAddress("TPCObj_NTracks", &TPCObj_NTracks);
   t -> SetBranchAddress("TPCObj_NShowers", &TPCObj_NShowers);
   t -> SetBranchAddress("TPCObj_PFP_isTrack", &TPCObj_PFP_isTrack);
   t -> SetBranchAddress("TPCObj_PFP_isShower", &TPCObj_PFP_isShower);
   t -> SetBranchAddress("TPCObj_PFP_isDaughter", &TPCObj_PFP_isDaughter);
   t -> SetBranchAddress("TPCObj_PFP_id", &TPCObj_PFP_id);
   t -> SetBranchAddress("TPCObj_PFP_MCPid", &TPCObj_PFP_MCPid);
   t -> SetBranchAddress("TPCObj_PFP_truePDG", &TPCObj_PFP_truePDG);
   t -> SetBranchAddress("TPCObj_PFP_trueE", &TPCObj_PFP_trueE);
   t -> SetBranchAddress("TPCObj_PFP_trueKE", &TPCObj_PFP_trueKE);
   t -> SetBranchAddress("TPCObj_origin", &TPCObj_origin);
   t -> SetBranchAddress("TPCObj_origin_extra", &TPCObj_origin_extra);
   t -> SetBranchAddress("TPCObj_reco_vtx", &TPCObj_reco_vtx);

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
   t -> SetBranchAddress("MCP_KE", &MCP_KE);
   t -> SetBranchAddress("MCP_isContained", &MCP_isContained);

   t -> SetBranchAddress("nu_vtxx", &nu_vtxx);
   t -> SetBranchAddress("nu_vtxy", &nu_vtxy);
   t -> SetBranchAddress("nu_vtxz", &nu_vtxz);
   t -> SetBranchAddress("nu_isCC", &nu_isCC);
   t -> SetBranchAddress("nu_PDG", &nu_PDG);
   t -> SetBranchAddress("nu_E", &nu_E);

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
   StackedHistTopology* tracks_stack = new StackedHistTopology("tracks_stack", ";Number of #nu daughter tracks in selected TPCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)", 10, 0, 10);
   StackedHistTopology* showers_stack = new StackedHistTopology("showers_stack", ";Number of #nu daughter showers in selected TCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)", 10, 0, 10);
   StackedHistTopology* mips_stack = new StackedHistTopology("mips_stack", ";Number of MIP-like #nu daughter tracks in selected TPCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)", 10, 0, 10);
   StackedHistTopology* vtxtrack_stack = new StackedHistTopology("vtxtrack_stack", ";Distance between vertex and start of #nu daughter tracks in selected TPCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)", 50, 0, 50);
   StackedHistTopology* vtxshwr_stack = new StackedHistTopology("vtxshwr_stack", ";Distance between vertex and start of #nu daughter showers in selected TPCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)", 50, 0, 50);
   StackedHistTopology* pfps_stack = new StackedHistTopology("pfps_stack", ";Number of #nu daughter PFPs in selected TPCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)", 10, 0, 10);

   StackedHistPDGCode* shwrtrue_stack = new StackedHistPDGCode("shwrtrue_sig",";True CC1pi: Number of shower #nu daughter tracks in selected TPCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)",10,0,10);
   StackedHistPDGCode* shwrtrue_bg_stack = new StackedHistPDGCode("shwrtrue_bg_stack", ";Background: Number of shower #nu daughter tracks in selected TPCObject;Selected Events (normalised to 3.782 #times 10^{19} POT)",10,0,10);
   StackedHistPDGCode* len_stack_byPDG = new StackedHistPDGCode("len_stack_byPDG", ";Longest Track Length [cm];Selected Events (normalised to 3.782 #times 10^{19} POT)", 30, 0, 700);
   StackedHistPDGCode* isMIP_stack = new StackedHistPDGCode("isMIP_stack", ";#nu daughter track is MIP-like?;Selected Events (normalised to 3.782 #times 10^{19} POT)",2,0,2);

   TEfficiency* eff_track_muon = new TEfficiency("eff_track_muon",";True Kinetic Energy [GeV];",10,0,2);
   TEfficiency* eff_track_pion = new TEfficiency("eff_track_pion",";True Kinetic Energy [GeV];",10,0,2);
   TEfficiency* eff_track_proton = new TEfficiency("eff_track_proton",";True Kinetic Energy [GeV];",10,0,2);

   TH2D* longest2tracks_byPDG = new TH2D("longest2tracks_byPDG",";True PDG of longest track;True PDG of second longest track",4,0,4,4,0,4);

   const int nentries = t -> GetEntries();

   for (int i = 0; i < nentries; i++) {

      t -> GetEntry(i);

      // Set topology (check for cosmic, mixed, OFV)
      // Move this to module! 
      if (TPCObj_origin == 1) TPCObj_beamnu_topology = kCosmic;
      else if (TPCObj_origin == 2) TPCObj_beamnu_topology = kMixed;
      else if (TPCObj_origin == 0 && !inFV(nu_vtxx -> at(0), nu_vtxy -> at(0), nu_vtxz -> at(0))) TPCObj_beamnu_topology = kOutFV;
      else if (TPCObj_origin == 0 && PDGCode(abs(nu_PDG -> at(0))) == kNuE && TPCObj_beamnu_topology != kNC) TPCObj_beamnu_topology = kCCNue;


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
      //      if (TPCObj_beamnu_topology == kCC1piplus0p || TPCObj_beamnu_topology == kCC1piplus1p || TPCObj_beamnu_topology == kCC1piplusNp) {
      if(nu_isCC -> at(0) && nu_PDG -> at(0) == 14 && inFV(nu_vtxx -> at(0), nu_vtxy -> at(0), nu_vtxz -> at(0)) && std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 13) == 1 && std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 211) == 1 && MCP_PDG -> size() == 2 + std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 2112) + std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 2212) + std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 2000000101) + std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 1000180400)) {

         isSignal = true;

         // Check neutrino flavour (just in case)
         if (nu_PDG -> at(0) != 14) std::cout << "WARNING: Signal event has neutrino with PDG code " << nu_PDG -> at(0) << " (should be 14)" << std::endl;

         // Efficiency...
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

         // Purity...
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

         int shower_daughters = 0;
         int track_daughters = 0;
         int mip_daughters = 0;
         std::vector<double> vtxtrack_daughters;
         std::vector<double> vtxshwr_daughters;

         for(int j = 0; j < TPCObj_NPFPs; j++) {
            if(TPCObj_PFP_isShower -> at(j) && TPCObj_PFP_isDaughter -> at(j)){
               shower_daughters++;
               vtxshwr_daughters.emplace_back(std::hypot(std::hypot(TPCObj_reco_vtx->at(0) - TPCObj_PFP_shower_start->at(j).at(0), TPCObj_reco_vtx->at(1) - TPCObj_PFP_shower_start->at(j).at(1)), TPCObj_reco_vtx->at(2) - TPCObj_PFP_shower_start->at(j).at(2)));
            }
            if(TPCObj_PFP_isTrack -> at(j) && TPCObj_PFP_isDaughter -> at(j)){
               track_daughters++;
               vtxtrack_daughters.emplace_back(std::hypot(std::hypot(TPCObj_reco_vtx->at(0) - TPCObj_PFP_track_start->at(j).at(0), TPCObj_reco_vtx->at(1) - TPCObj_PFP_track_start->at(j).at(1)), TPCObj_reco_vtx->at(2) - TPCObj_PFP_track_start->at(j).at(2)));
            }
            if(TPCObj_PFP_isMIP -> at(j) && TPCObj_PFP_isDaughter -> at(j)) {
               mip_daughters++;
            }
         }

         len_stack -> Fill(TPCObj_beamnu_topology, maxlen);
         tracks_stack -> Fill(TPCObj_beamnu_topology, track_daughters);
         showers_stack -> Fill(TPCObj_beamnu_topology, shower_daughters);
         mips_stack -> Fill(TPCObj_beamnu_topology, mip_daughters);
         pfps_stack -> Fill(TPCObj_beamnu_topology, track_daughters+shower_daughters);
         for (int j = 0; j < track_daughters; j++) {
            vtxtrack_stack -> Fill(TPCObj_beamnu_topology, vtxtrack_daughters.at(j));
         }
         for (int j = 0; j < shower_daughters; j++) {
            vtxshwr_stack -> Fill(TPCObj_beamnu_topology, vtxshwr_daughters.at(j));
         }
         len_stack_byPDG -> Fill(maxPDG, maxlen);
         for(int j = 0; j < TPCObj_NPFPs; j++) {
            if(TPCObj_PFP_isDaughter -> at(j) && TPCObj_PFP_isTrack -> at(j)) isMIP_stack -> Fill((PDGCode)TPCObj_PFP_truePDG->at(j),TPCObj_PFP_isMIP->at(j));
         }

         //Longest 2 tracks by PDG...
         if(TPCObj_NTracks >= 2) {
            if(maxPDG == kMuMinus) {
               if(maxPDG2 == kMuMinus) longest2tracks_byPDG -> Fill(0., 0);
               else if(maxPDG2 == kPiPlus) longest2tracks_byPDG -> Fill(0.,1);
               else if(maxPDG2 == kProton) longest2tracks_byPDG -> Fill(0.,2);
               else longest2tracks_byPDG -> Fill(0.,3);
            }
            else if(maxPDG == kPiPlus) {
               if(maxPDG2 == kMuMinus) longest2tracks_byPDG -> Fill(1,0);
               else if(maxPDG2 == kPiPlus) longest2tracks_byPDG -> Fill(1,1);
               else if(maxPDG2 == kProton) longest2tracks_byPDG -> Fill(1,2);
               else longest2tracks_byPDG -> Fill(1,3);
            }
            else if(maxPDG == kProton) {
               if(maxPDG2 == kMuMinus) longest2tracks_byPDG -> Fill(2,0);
               else if(maxPDG2 == kPiPlus) longest2tracks_byPDG -> Fill(2,1);
               else if(maxPDG2 == kProton) longest2tracks_byPDG -> Fill(2,2);
               else longest2tracks_byPDG -> Fill(2,3);
            }
            else {
               if(maxPDG2 == kMuMinus) longest2tracks_byPDG -> Fill(3,0);
               else if(maxPDG2 == kPiPlus) longest2tracks_byPDG -> Fill(3,1);
               else if(maxPDG2 == kProton) longest2tracks_byPDG -> Fill(3,2);
               else longest2tracks_byPDG -> Fill(3,3);
            }
         }

         // Track/shower classification...
         for(int j = 0; j < TPCObj_NPFPs; j++) {
            if (TPCObj_PFP_truePDG -> at(j) == 13) {
               if(TPCObj_PFP_isTrack -> at(j)) {
                  eff_track_muon -> Fill(true, TPCObj_PFP_trueKE -> at(j));
               }
               else if(TPCObj_PFP_isShower -> at(j)) {
                  eff_track_muon -> Fill(false, TPCObj_PFP_trueKE -> at(j));
               }
            }
            else if(TPCObj_PFP_truePDG -> at(j) == 211) {
               if(TPCObj_PFP_isTrack -> at(j)) {
                  eff_track_pion -> Fill(true, TPCObj_PFP_trueKE -> at(j));
               }
               else if(TPCObj_PFP_isShower -> at(j)) {
                  eff_track_pion -> Fill(false, TPCObj_PFP_trueKE -> at(j));
               }
            }
            else if(TPCObj_PFP_truePDG -> at(j) == 2212) {
               if(TPCObj_PFP_isTrack -> at(j)) {
                  eff_track_proton -> Fill(true, TPCObj_PFP_trueKE -> at(j));
               }
               else if(TPCObj_PFP_isShower -> at(j)) {
                  eff_track_proton -> Fill(false, TPCObj_PFP_trueKE -> at(j));
               }
            }

            // Look at what particles are being associated with showers
            if (TPCObj_PFP_isShower -> at(j) && TPCObj_PFP_isDaughter -> at(j)){
               PDGCode mcp_pdgcode = (PDGCode)TPCObj_PFP_truePDG->at(j);
               if(isSignal) { // !!!Change this to use topology!!! If true CC1pip
                  shwrtrue_stack->Fill(mcp_pdgcode, shower_daughters, 1.0/shower_daughters);
               }
               else{
                  shwrtrue_bg_stack->Fill(mcp_pdgcode, shower_daughters, 1.0/shower_daughters);
               }
            } // end if (TPCObj_PFP_isShower -> at(j))  
         } // end loop over PFPs
      }
   }

   gStyle->SetOptStat(0);

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

   tracks_stack -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THStack_tracks.eps",SaveString.c_str()));

   showers_stack -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THStack_showers.eps",SaveString.c_str()));

   mips_stack -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THStack_mips.eps",SaveString.c_str()));

   pfps_stack -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THStack_pfps.eps",SaveString.c_str()));

   vtxtrack_stack -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THStack_vtxtrack.eps",SaveString.c_str()));

   vtxshwr_stack -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THStack_vtxshwr.eps",SaveString.c_str()));

   shwrtrue_stack -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THStack_signal_shwrtrue.eps",SaveString.c_str()));

   shwrtrue_bg_stack -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THStack_bg_shwrtrue.eps",SaveString.c_str()));

   len_stack_byPDG -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THStack_trklen_byPDG.eps",SaveString.c_str()));

   isMIP_stack -> DrawStack(norm, c1);
   c1 -> SaveAs(TString::Format("%s_THSTack_isMIP.eps",SaveString.c_str()));

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

   longest2tracks_byPDG->GetXaxis()->SetBinLabel(1,PDGenum2str(kMuMinus).c_str());
   longest2tracks_byPDG->GetXaxis()->SetBinLabel(2,PDGenum2str(kPiPlus).c_str());
   longest2tracks_byPDG->GetXaxis()->SetBinLabel(3,PDGenum2str(kProton).c_str());
   longest2tracks_byPDG->GetXaxis()->SetBinLabel(4,"Other");
   longest2tracks_byPDG->GetYaxis()->SetBinLabel(1,PDGenum2str(kMuMinus).c_str());
   longest2tracks_byPDG->GetYaxis()->SetBinLabel(2,PDGenum2str(kPiPlus).c_str());
   longest2tracks_byPDG->GetYaxis()->SetBinLabel(3,PDGenum2str(kProton).c_str());
   longest2tracks_byPDG->GetYaxis()->SetBinLabel(4,"Other");
   longest2tracks_byPDG->Draw("COLZ TEXT");
   c1 -> SaveAs(TString::Format("%s_longest2tracks_byPDG.eps",SaveString.c_str()));

   std::cout << "Total Efficiency: " << (1. * (eff_nuE -> GetPassedHistogram() -> GetEntries()))/(eff_nuE -> GetTotalHistogram() -> GetEntries()) << std::endl;
   std::cout << "Total Purity: " << (1. * (pur_nuE -> GetPassedHistogram() -> GetEntries()))/(pur_nuE -> GetTotalHistogram() -> GetEntries()) << std::endl;

   std::cout << "Muon Track Efficiency: " << (1. * (eff_track_muon -> GetPassedHistogram() -> GetEntries()))/(eff_track_muon -> GetTotalHistogram() -> GetEntries()) << std::endl;
   std::cout << "Pion Track Efficiency: " << (1. * (eff_track_pion -> GetPassedHistogram() -> GetEntries()))/(eff_track_pion -> GetTotalHistogram() -> GetEntries()) << std::endl;
   std::cout << "Proton Track Efficiency: " << (1. * (eff_track_proton -> GetPassedHistogram() -> GetEntries()))/(eff_track_proton -> GetTotalHistogram() -> GetEntries()) << std::endl;

   gDirectory->GetList()->Delete();

   return;
   }
