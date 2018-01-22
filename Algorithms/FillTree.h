// Add includes!

// Struct to hold tree variables for CC1pi analysis

// Note: if you want to add a variable to the tree, you have to do three things:
//   1) add the data object in cc1pianavars
//   2) add a default value for the data object in Clear() (inside cc1pianavars)
//   3) add a line to set the branch in MakeAnaBranches

struct cc1pianavars{

   art::RunNumber_t run_num;
   art::SubRunNumber_t subrun_num;
   art::EventNumber_t event_num;

   std::map<std::string,bool> cutflow;
   bool isSelected;
   std::vector<float> track_length;
   std::vector<float> shower_length;
   int NPFPs;
   int NTracks;
   int NShowers;
   std::vector<bool> Sel_PFP_isTrack;
   std::vector<bool> Sel_PFP_isShower;
   std::vector<int> Sel_PFP_ID;
   std::vector<int> Sel_MCP_ID;
   std::vector<int> Sel_MCP_PDG;
   std::vector<float> Sel_MCP_E;
   int tpcobj_origin;
   int tpcobj_origin_extra;

   std::vector<int> MCP_PDG;
   std::vector<float> MCP_length;
   std::vector<std::string> MCP_process;
   std::vector<std::string> MCP_endprocess;
   std::vector<int> MCP_numdaughters;
   std::vector<float> MCP_P;
   std::vector<float> MCP_E;
   
   std::vector<float> nu_vtxx;
   std::vector<float> nu_vtxy;
   std::vector<float> nu_vtxz;
   std::vector<bool> nu_isCC;
   std::vector<int> nu_PDG;
   std::vector<float> nu_E;

   std::vector<bool> MIPConsistency;

   void Clear();

   void SetReco2Vars(art::Event &evt);

};

void MakeAnaBranches(TTree *t, cc1pianavars *vars);
