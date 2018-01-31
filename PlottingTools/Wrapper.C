#include "Wrapper.h"

//run by doing: root 'Wrapper.C("MyFileName.root")'

void Wrapper(TString FileName) {

   // MakePlots takes: the name of a cut (currently just these hardcoded four)
   // A bool for whether you want to look at events that pass or fail that cut
   // A string to append to plot save names
   // And the file name
   // 
   // Eventually the former two should be replaced with passing it a cutflow map

   MakePlots("MarcosSelec", true, "PassMarcosSelec", FileName);
   //MakePlots("TwoTrackCut", true, "PassTwoTrackCut", FileName);
   //MakePlots("TwoMIPCut", true, "PassTwoMIPCut", FileName);
   //MakePlots("ExactlyTwoMIPCut", true, "PassExactlyTwoMIPCut", FileName);

}
