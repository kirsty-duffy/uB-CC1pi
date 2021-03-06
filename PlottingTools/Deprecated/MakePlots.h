#include "TChain.h"
#include "TTree.h"
#include "TFile.h"
#include "TEfficiency.h"
#include "TVector3.h"
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TGraphAsymmErrors.h"

#include <vector>
#include <map>
#include <iostream>
#include <algorithm>

#include "../Algorithms/TopologyEnums.h"
#include "../Algorithms/PDGEnums.h"
#include "StackedHistPDGCode.h"
#include "StackedHistTopology.h"


double GetPOT(TString FileName);

void MakePlots(std::map<std::string,bool> SelectionCutflow, std::string SaveString, TString FileName);
