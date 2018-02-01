#include "Wrapper.h"

//run by doing: root 'Wrapper.C("MyFileName.root")'

void Wrapper(TString FileName) {

   // MakePlots takes: 
   // A map of cut names and whether you want to look at events that passed or failed them
   // A string to append to plot save names
   // And the file name

   std::map<std::string,bool> CC1picutflow;

   CC1picutflow["MarcosSelec"] = true;
   MakePlots(CC1picutflow, "PassMarcosSelec", FileName);

   CC1picutflow["TwoTrackCut"] = true;
   MakePlots(CC1picutflow, "PassTwoTrackCut", FileName);

   CC1picutflow["TwoMIPCut"] = true;
   MakePlots(CC1picutflow, "PassTwoMIPCut", FileName);

   CC1picutflow["ExactlyTwoMIPCut"] = true;
   MakePlots(CC1picutflow, "PassExactlyTwoMIPCut", FileName);

}
