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

// Can swap to using Marco's FV calculator to make this cleaner
// Can also just put this in the tree once we settle on a particular FV def
// This appears to be what's in v3.0.2, but Marco's latest presentation also had a chunk of dead region
// in z cut out, so double check (especially if we see higher event rate)
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

   if (x > FVxmin && x < FVxmax && y > FVymin && y < FVymax && z > FVzmin && z < FVzmax) return true;
   else return false;

}


void MakePlots(std::string Cut, bool Passes, std::string SaveString, TString FileName) {

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

   std::vector<double> *nu_vtxx = NULL;
   std::vector<double> *nu_vtxy = NULL;
   std::vector<double> *nu_vtxz = NULL;
   std::vector<bool> *nu_isCC = NULL;
   std::vector<int> *nu_PDG = NULL;
   std::vector<double> *nu_E = NULL;

   std::vector<bool> *MIPConsistency = NULL;
   std::vector<double> *dqdx_trunc_uncalib = NULL;

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

   t -> SetBranchAddress("nu_vtxx", &nu_vtxx);
   t -> SetBranchAddress("nu_vtxy", &nu_vtxy);
   t -> SetBranchAddress("nu_vtxz", &nu_vtxz);
   t -> SetBranchAddress("nu_isCC", &nu_isCC);
   t -> SetBranchAddress("nu_PDG", &nu_PDG);
   t -> SetBranchAddress("nu_E", &nu_E);

   t -> SetBranchAddress("MIPConsistency", &MIPConsistency);
   t -> SetBranchAddress("dqdx_trunc_uncalib", &dqdx_trunc_uncalib);

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

   THStack* len_stack = new THStack("len_stack", ";Longest Track Length [cm];Selected Events (normalised to 3.728 #times 10^{19} POT)");
   TH1F* len_muCC1pip = new TH1F("len_muCC1pip", ";;", 30, 0, 700);
   TH1F* len_otherCC1pip = new TH1F("len_otherCC1pip", ";;", 30, 0, 700);
   TH1F* len_CCother = new TH1F("len_CCother", ";;", 30, 0, 700);
   TH1F* len_NC = new TH1F("len_NC", ";;", 30, 0, 700);
   TH1F* len_cosmic = new TH1F("len_cosmic", ";;", 30, 0, 700);
   TH1F* len_mixed = new TH1F("len_mixed", ";;", 30, 0, 700);
   TH1F* len_unknown = new TH1F("len_unknown", ";;", 30, 0, 700);
   TH1F* len_outFV = new TH1F("len_outFV", ";;", 30, 0, 700);

   len_muCC1pip -> SetFillColor(kRed);
   len_otherCC1pip -> SetFillColor(kOrange);
   len_CCother -> SetFillColor(kMagenta);
   len_NC -> SetFillColor(kGray);
   len_cosmic -> SetFillColor(kBlue);
   len_mixed -> SetFillColor(kCyan);
   len_unknown -> SetFillColor(kBlack);
   len_outFV -> SetFillColor(kGreen);

   TEfficiency* eff_track_muon = new TEfficiency("eff_track_muon",";True Energy [GeV];",10,0,2);
   TEfficiency* eff_track_pion = new TEfficiency("eff_track_pion",";True Energy [GeV];",10,0,2);
   TEfficiency* eff_track_proton = new TEfficiency("eff_track_proton",";True Energy [GeV];",10,0,2);

   const int nentries = t -> GetEntries();

   for (int i = 0; i < nentries; i++) {

      t -> GetEntry(i);

      // Replace this with a cutflow map: MakePlots should take a std::map<std::string,bool>
      // And check for each of the given cuts (strings) whether the pass/fail value matches the given bool
      if(Passes) {
         if(Cut.compare("ExactlyTwoMIPCut") == 0) {
            if(CC1piSelecFailureReason -> compare("") != 0) continue;
         }
         else if(Cut.compare("TwoMIPCut") == 0) {
            if((CC1piSelecFailureReason -> compare("") != 0) && (CC1piSelecFailureReason -> compare("ExactlyTwoMIPCut") != 0)) continue;
         }
         else if(Cut.compare("TwoTrackCut") == 0) {
            if((CC1piSelecFailureReason -> compare("") != 0) && (CC1piSelecFailureReason -> compare("ExactlyTwoMIPCut") != 0) && (CC1piSelecFailureReason -> compare("TwoMIPCut") != 0)) continue;
         }
         else if(Cut.compare("MarcosSelec") == 0) {
            if(!isSelected) continue;
         }
      }
      else {
         if(Cut.compare("ExactlyTwoMIPCut") == 0) {
            if(!(CC1piSelecFailureReason -> compare("") != 0)) continue;
         }
         else if(Cut.compare("TwoMIPCut") == 0) {
            if(!((CC1piSelecFailureReason -> compare("") != 0) && (CC1piSelecFailureReason -> compare("ExactlyTwoMIPCut") != 0))) continue;
         }
         else if(Cut.compare("TwoTrackCut") == 0) {
            if(!((CC1piSelecFailureReason -> compare("") != 0) && (CC1piSelecFailureReason -> compare("ExactlyTwoMIPCut") != 0) && (CC1piSelecFailureReason -> compare("TwoMIPCut") != 0))) continue;
         }
         else if(Cut.compare("MarcosSelec") == 0) {
            if(isSelected) continue;
         }
      }


      bool isSignal = false;

      // CC; muon neutrino; vtx in FV; with exactly 1 muon, exactly 1 pi+, any nucleons/Ar-40, and nothing else
      if(nu_isCC -> at(0) && nu_PDG -> at(0) == 14 && inFV(nu_vtxx -> at(0), nu_vtxy -> at(0), nu_vtxz -> at(0)) && std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 13) == 1 && std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 211) == 1 && MCP_PDG -> size() == 2 + std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 2112) + std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 2212) + std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 2000000101)) {

         isSignal = true;

         //Efficiency...
         eff_nuE -> Fill(isSelected, nu_E -> at(0));

         for (int j = 0; j < MCP_PDG -> size(); j++) {
            if (MCP_PDG -> at(j) == 13) {
               eff_muP -> Fill(isSelected, MCP_P -> at(j));
               TVector3 muP(MCP_Px -> at(j), MCP_Py -> at(j), MCP_Pz -> at(j));
               eff_muCosT -> Fill(isSelected, muP.CosTheta());
               eff_muPhi -> Fill(isSelected, muP.Phi());
               break;
            }
         }
      }

      if(isSelected) {

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

         //THStack...
         double maxlen = *std::max_element(track_length -> begin(), track_length -> end());

         if(tpcobj_origin==1) {
            len_cosmic -> Fill(maxlen);
         }
         else if(tpcobj_origin==2) {
            len_mixed -> Fill(maxlen);
         }
         else if(tpcobj_origin!=0) {
            len_unknown -> Fill(maxlen);
         }
         else if(!inFV(nu_vtxx -> at(0), nu_vtxy -> at(0), nu_vtxz -> at(0))) {
            len_outFV -> Fill(maxlen);
         }
         else if(tpcobj_origin==0) {
            if(nu_isCC -> at(0) == 0) {
               len_NC -> Fill(maxlen);
            }
            else if(std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 13) == 1 && std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 211) == 1 && MCP_PDG -> size() == 2 + std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 2112) + std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 2212) + std::count(MCP_PDG -> begin(), MCP_PDG -> end(), 2000000101)) {
               if(nu_PDG -> at(0) == 14) {
                  len_muCC1pip -> Fill(maxlen);
               }
               else {
                  len_otherCC1pip -> Fill(maxlen);
               }
            }
            else {
               len_CCother -> Fill(maxlen);
            }
         }

         //Track/shower classification...
         for(int j = 0; j < Sel_MCP_PDG -> size(); j++) {
            if (Sel_MCP_PDG -> at(j) == 13) {
               if(Sel_PFP_isTrack -> at(j)) {
                  eff_track_muon -> Fill(true, Sel_MCP_E -> at(j));
               }
               else if(Sel_PFP_isShower -> at(j)) {
                  eff_track_muon -> Fill(false, Sel_MCP_E -> at(j));
               }
            }
            else if(Sel_MCP_PDG -> at(j) == 211) {
               if(Sel_PFP_isTrack -> at(j)) {
                  eff_track_pion -> Fill(true, Sel_MCP_E -> at(j));
               }
               else if(Sel_PFP_isShower -> at(j)) {
                  eff_track_pion -> Fill(false, Sel_MCP_E -> at(j));
               }
            }
            else if(Sel_MCP_PDG -> at(j) == 2212) {
               if(Sel_PFP_isTrack -> at(j)) {
                  eff_track_proton -> Fill(true, Sel_MCP_E -> at(j));
               }
               else if(Sel_PFP_isShower -> at(j)) {
                  eff_track_proton -> Fill(false, Sel_MCP_E -> at(j));
               }
            }
         }
      }
   }


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

   double pot_marco = 3.728e+19;
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
   len_stack -> Add(len_unknown);
   len_stack -> Add(len_mixed);
   len_stack -> Add(len_cosmic);
   len_stack -> Add(len_outFV);
   len_stack -> Add(len_NC);
   len_stack -> Add(len_CCother);
   len_stack -> Add(len_otherCC1pip);
   len_stack -> Add(len_muCC1pip);
   len_stack -> Draw();
   auto len_legend = new TLegend(0.6,0.6,0.9,0.9);
   len_legend -> AddEntry(len_muCC1pip, "#nu_{#mu} CC1#pi^{+}", "f");
   len_legend -> AddEntry(len_otherCC1pip, "#bar{#nu}_{#mu}, #nu_{e}, #bar{#nu}_{e} CC1#pi^{+}", "f");
   len_legend -> AddEntry(len_CCother, "CC-Other", "f");
   len_legend -> AddEntry(len_NC, "NC", "f");
   len_legend -> AddEntry(len_outFV, "Out of FV", "f");
   len_legend -> AddEntry(len_cosmic, "Cosmic", "f");
   len_legend -> AddEntry(len_mixed, "Mixed", "f");
   len_legend -> AddEntry(len_unknown, "Unknown Origin", "f");
   len_legend -> Draw("SAME");
   c1 -> SaveAs(TString::Format("%s_THStack_len.eps",SaveString.c_str()));

   std::cout << "muCC1pip: " << len_muCC1pip -> GetEntries() << std::endl;
   std::cout << "otherCC1pip: " << len_otherCC1pip -> GetEntries() << std::endl; 
   std::cout << "CCother: " << len_CCother -> GetEntries() << std::endl;
   std::cout << "NC: " << len_NC -> GetEntries() << std::endl;
   std::cout << "outFV: " << len_outFV -> GetEntries() << std::endl;
   std::cout << "cosmic: " << len_cosmic -> GetEntries() << std::endl;
   std::cout << "mixed: " << len_mixed -> GetEntries() << std::endl;
   std::cout << "unknown: " << len_unknown -> GetEntries() << std::endl;

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

   std::cout << "Total Efficiency: " << (1. * (eff_nuE -> GetPassedHistogram() -> GetEntries()))/(eff_nuE -> GetTotalHistogram() -> GetEntries()) << std::endl;
   std::cout << "Total Purity: " << (1. * (pur_nuE -> GetPassedHistogram() -> GetEntries()))/(pur_nuE -> GetTotalHistogram() -> GetEntries()) << std::endl;

   return;
}
