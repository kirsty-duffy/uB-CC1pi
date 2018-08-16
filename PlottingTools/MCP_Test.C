#include "../Algorithms/PDGEnums.h"

void MCP_Test(TString FileName) {

   //Setup tree...
   TChain *t = new TChain("cc1piselec/outtree");
   t -> Add(FileName);

   std::vector<int> *MCP_PDG=nullptr;
   std::vector<double> *MCP_length=nullptr;
   std::vector<std::string> *MCP_process=nullptr;
   std::vector<std::string> *MCP_endprocess=nullptr;
   std::vector<int> *MCP_numdaughters=nullptr;
   std::vector<double> *MCP_P=nullptr;
   std::vector<double> *MCP_Px=nullptr;
   std::vector<double> *MCP_Py=nullptr;
   std::vector<double> *MCP_Pz=nullptr;
   std::vector<double> *MCP_E=nullptr;
   std::vector<double> *MCP_KE=nullptr;
   std::vector<bool> *MCP_isContained=nullptr;
   std::vector<int> *MCP_ID=nullptr;
   std::vector<std::vector<int>> *MCP_DaughterIDs=nullptr;
   std::vector<int> *MCP_StatusCode=nullptr;
   std::vector<int> *MCP_MotherID=nullptr;
   std::vector<std::vector<double>> *MCP_StartPosition=nullptr;
   std::vector<std::vector<double>> *MCP_EndPosition=nullptr;

   t -> SetBranchStatus("*",0);

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
   t -> SetBranchStatus("MCP_ID",1);
   t -> SetBranchAddress("MCP_ID", &MCP_ID);
   t -> SetBranchStatus("MCP_DaughterIDs",1);
   t -> SetBranchAddress("MCP_DaughterIDs", &MCP_DaughterIDs);
   t -> SetBranchStatus("MCP_StatusCode",1);
   t -> SetBranchAddress("MCP_StatusCode", &MCP_StatusCode);
   t -> SetBranchStatus("MCP_MotherID",1);
   t -> SetBranchAddress("MCP_MotherID", &MCP_MotherID);
   t -> SetBranchStatus("MCP_StartPosition",1);
   t -> SetBranchAddress("MCP_StartPosition", &MCP_StartPosition);
   t -> SetBranchStatus("MCP_EndPosition",1);
   t -> SetBranchAddress("MCP_EndPosition", &MCP_EndPosition);

   const int nentries = t -> GetEntries();

   for (int i_entry = 0; i_entry < nentries; i_entry++) {
      t -> GetEntry(i_entry);

      for (int i_MCP = 0; i_MCP < MCP_ID->size(); i_MCP++) {

         if(MCP_MotherID->at(i_MCP)==0 && (MCP_PDG->at(i_MCP)==kMuMinus || MCP_PDG->at(i_MCP)==kPiPlus)) {
            std::cout << "Found neutrino daughter with PDG " << MCP_PDG->at(i_MCP) << " and EndProcess " << MCP_endprocess->at(i_MCP) << std::endl;
            std::cout << "Daughters: " << std::endl;
            for (int i_Daughter = 0; i_Daughter < MCP_DaughterIDs->at(i_MCP).size(); i_Daughter++) {
               std::cout << "ID: " << MCP_DaughterIDs->at(i_MCP).at(i_Daughter);
               int Daughter_index = std::distance(MCP_ID->begin(),std::find(MCP_ID->begin(), MCP_ID->end(), MCP_DaughterIDs->at(i_MCP).at(i_Daughter)));
               if(Daughter_index >= MCP_ID->size()) {
                  std::cout << " - Daughter not found?" << std::endl;
                  continue;
               }
               std::cout << " – PDG: " << MCP_PDG->at(Daughter_index);
               std::cout << " – Start Process: " << MCP_process->at(Daughter_index);
               std::cout << " – End Process: " << MCP_endprocess->at(Daughter_index);
               double distance = std::hypot(std::hypot(MCP_EndPosition->at(i_MCP).at(0)-MCP_StartPosition->at(Daughter_index).at(0),MCP_EndPosition->at(i_MCP).at(1)-MCP_StartPosition->at(Daughter_index).at(1)),MCP_EndPosition->at(i_MCP).at(2)-MCP_StartPosition->at(Daughter_index).at(2));
               std::cout << " - Distance from end: " << distance;
               std::cout << std::endl;
            }
         }
      } // End loop over MCPs
      std::cout << std::endl;
   } // End loop over entries

/*
// Just print out basic information: All the MCP IDs and all their daughter IDs
   for (int i_entry = 0; i_entry < 5; i_entry++) {
      t -> GetEntry(i_entry);

      std::cout << std::endl << "Event " << i_entry << std::endl;
      for (int i_MCP = 0; i_MCP < MCP_ID->size(); i_MCP++) {
         std::cout << "ID: " << MCP_ID -> at(i_MCP) << std::endl;
         for (int i_Daughter = 0; i_Daughter < MCP_DaughterIDs->at(i_MCP).size(); i_Daughter++) {
            std::cout << "          ID: " << MCP_DaughterIDs->at(i_MCP).at(i_Daughter) << std::endl;
         }
      }
   }
*/
}
