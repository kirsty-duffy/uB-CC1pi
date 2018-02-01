#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TEfficiency.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TString.h"

#include <vector>
#include <map>
#include <iostream>


double GetPOT(TString FileName);

bool inFV(double x, double y, double z);

void MakePlots(std::map<std::string,bool> SelectionCutflow, std::string SaveString, TString FileName);
