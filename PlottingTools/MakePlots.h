#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TEfficiency.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TString.h"

#include <vector>
#include <iostream>


double GetPOT(TString FileName);

bool inFV(double x, double y, double z);

void MakePlots(std::string Cut, bool Passes, std::string SaveString, TString FileName);
