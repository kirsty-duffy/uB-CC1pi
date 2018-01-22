// Add includes!

// Struct to hold tree variables for CC1pi analysis

// Note: if you want to add a variable to the tree, you have to do three things:
//   1) add the data object in cc1pianavars
//   2) add a default value for the data object in Clear() (inside cc1pianavars)
//   3) add a line to set the branch in MakeAnaBranches

struct cc1pianavars{

  // Change these to whatever you like
  Int_t evtnum;
  bool MIPConsistency;

  void Clear();

  void SetReco2Vars(art::Event &evt);

};

void MakeAnaBranches(TTree *t, cc1pianavars *vars);
