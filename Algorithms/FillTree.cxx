#include "FillTree.h"

// probably need some art includes here too
// and some root ones

void cc1pianavars::Clear(){
  
  // Set default values (that will be used if the event isn't filled)
  evtnum = -999;
  MIPConsistency = false;
  
}



void cc1pianavars::SetReco2Vars(art::Event &evt){
  
  // Set all the values that are in the reco2 file
  // This just cleans the module up - move all these lines to here instead of in the module
  evtnum = e.EvtNum(); // almost definitely wrong, I don't remember the syntax, but you know what I mean...
  
}



void MakeAnaBranches(TTree *t, cc1pianavars *vars){

  // Note: it's very important that we use the syntax &(vars->evtnum) to make sure
  // we get the value from reference, otherwise it won't work!
  t->Branch("evtnum", &(vars->evtnum), "evtnum/I");
  t->Branch("MIPConsistency", &(vars->MIPConsistency), "MIPConsistency/O");

}
