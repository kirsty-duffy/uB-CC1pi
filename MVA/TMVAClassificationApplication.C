/**********************************************************************************
 * Project   : TMVA - a Root-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Exectuable: TMVAClassificationApplication                                      *
 *                                                                                *
 * This macro provides a simple example on how to use the trained classifiers     *
 * within an analysis module                                                      *
 **********************************************************************************/

#include <cstdlib>
#include <vector>
#include <iostream>
#include <map>
#include <string>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TSystem.h"
#include "TROOT.h"
#include "TStopwatch.h"
#include "TGraph.h"
#include "TCanvas.h"
#include "TMarker.h"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Reader.h"
#include "TMVA/MethodCuts.h"
#endif

using namespace TMVA;
using namespace std;

void TMVAClassificationApplication( TString myMethodList = "" ) 
{   
#ifdef __CINT__
   gROOT->ProcessLine( ".O0" ); // turn off optimization in CINT
#endif

   //---------------------------------------------------------------

   // This loads the library
   TMVA::Tools::Instance();

   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

   ofstream outfile("pass.txt");

   // --- Cut optimisation
   Use["Cuts"]            = 0;
   Use["CutsD"]           = 0;
   Use["CutsPCA"]         = 0;
   Use["CutsGA"]          = 0;
   Use["CutsSA"]          = 0;
   // 
   // --- 1-dimensional likelihood ("naive Bayes estimator")
   Use["Likelihood"]      = 0;
   Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
   Use["LikelihoodPCA"]   = 0; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
   Use["LikelihoodKDE"]   = 0;
   Use["LikelihoodMIX"]   = 0;
   //
   // --- Mutidimensional likelihood and Nearest-Neighbour methods
   Use["PDERS"]           = 0;
   Use["PDERSD"]          = 0;
   Use["PDERSPCA"]        = 0;
   Use["PDEFoam"]         = 0;
   Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
   Use["KNN"]             = 0; // k-nearest neighbour method
   //
   // --- Linear Discriminant Analysis
   Use["LD"]              = 0; // Linear Discriminant identical to Fisher
   Use["Fisher"]          = 0;
   Use["FisherG"]         = 0;
   Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
   Use["HMatrix"]         = 0;
   //
   // --- Function Discriminant analysis
   Use["FDA_GA"]          = 0; // minimisation of user-defined function using Genetics Algorithm
   Use["FDA_SA"]          = 0;
   Use["FDA_MC"]          = 0;
   Use["FDA_MT"]          = 0;
   Use["FDA_GAMT"]        = 0;
   Use["FDA_MCMT"]        = 0;
   //
   // --- Neural Networks (all are feed-forward Multilayer Perceptrons)
   Use["MLP"]             = 0; // Recommended ANN
   Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
   Use["MLPBNN"]          = 0; // Recommended ANN with BFGS training method and bayesian regulator
   Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
   Use["TMlpANN"]         = 0; // ROOT's own ANN
   //
   // --- Support Vector Machine 
   Use["SVM"]             = 0;
   // 
   // --- Boosted Decision Trees
   Use["BDT"]             = 0; // uses Adaptive Boost
   Use["BDTG"]            = 1; // uses Gradient Boost
   Use["BDTB"]            = 0; // uses Bagging
   Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
   // 
   // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 0;
   // ---------------------------------------------------------------
   Use["Plugin"]          = 0;
   Use["Category"]        = 0;
   Use["SVM_Gauss"]       = 0;
   Use["SVM_Poly"]        = 0;
   Use["SVM_Lin"]         = 0;

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassificationApplication" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod 
                      << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
               std::cout << it->first << " ";
            }
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // --- Create the Reader object

   TMVA::Reader *reader = new TMVA::Reader( "!Color:!Silent" );    

   // Create a set of variables and declare them to the reader
   // - the variable names MUST corresponds in name and type to those given in the weight file(s) used
   float evtcharge;
   float ntrack;
   float maxtrklength, avgtrklength;
   float trkdedx;
   float trkrch;
   float trkrt;
   float trkfr;
   float trkpida;
   float nshower;
   float showerdedx;
   float eshower;
   float frshower;
   float nhitspershw;
   float shwlength;
   float shwmax;
   float fract_5_wires;
   float fract_10_wires;
   float fract_50_wires;
   float fract_100_wires;
   float shwdisx;
   float shwdisy;
   float shwdisz;
   float shwcosx;
   float shwcosy;
   float shwcosz;
   float trkcosx;
   float trkcosy;
   float trkcosz, ET;
   Float_t var1, var2;
   Float_t var3, var4;

   bool SelectNuE = true;
   bool SelectNuMu = false;
   //bool SelectNuE = false;
   //bool SelectNuMu = true;

   reader->AddVariable("evtcharge", &evtcharge);
   reader->AddVariable("ntrack", &ntrack);
   reader->AddVariable("maxtrklength", &maxtrklength);
   reader->AddVariable("avgtrklength", &avgtrklength);
   reader->AddVariable("trkdedx", &trkdedx);
   reader->AddVariable("trkrch", &trkrch);
   reader->AddVariable("trkrt", &trkrt);
   reader->AddVariable("trkfr", &trkfr);
   reader->AddVariable("trkpida", &trkpida);
   reader->AddVariable("fract_5_wires", &fract_5_wires);
   reader->AddVariable("fract_10_wires", &fract_10_wires);
   reader->AddVariable("fract_50_wires", &fract_50_wires);
   reader->AddVariable("fract_100_wires", &fract_100_wires);
   reader->AddVariable("trkcosx", &trkcosx);
   reader->AddVariable("trkcosy", &trkcosy);
   reader->AddVariable("trkcosz", &trkcosz);
   reader->AddVariable("ET", &ET);

   reader->AddVariable("nshower", &nshower);
   reader->AddVariable("showerdedx", &showerdedx);
   reader->AddVariable("eshower", &eshower);
   reader->AddVariable("frshower", &frshower);
   reader->AddVariable("nhitspershw", &nhitspershw);
   reader->AddVariable("shwlength", &shwlength);
   reader->AddVariable("shwmax", &shwmax);
   reader->AddVariable("shwdisx", &shwdisx);
   reader->AddVariable("shwdisy", &shwdisy);
   reader->AddVariable("shwdisz", &shwdisz);
   reader->AddVariable("shwcosx", &shwcosx);
   reader->AddVariable("shwcosy", &shwcosy);
   reader->AddVariable("shwcosz", &shwcosz);

   
   // Spectator variables declared in the training have to be added to the reader, too
//   Float_t spec1,spec2;
//   reader->AddSpectator( "spec1 := var1*2",   &spec1 );
//   reader->AddSpectator( "spec2 := var1*3",   &spec2 );

   Float_t Category_cat1, Category_cat2, Category_cat3;
   if (Use["Category"]){
      // Add artificial spectators for distinguishing categories
      reader->AddSpectator( "Category_cat1 := var3<=0",             &Category_cat1 );
      reader->AddSpectator( "Category_cat2 := (var3>0)&&(var4<0)",  &Category_cat2 );
      reader->AddSpectator( "Category_cat3 := (var3>0)&&(var4>=0)", &Category_cat3 );
   }

   // --- Book the MVA methods

   TString dir    = "weights/";
   TString prefix = "TMVAClassification";

   // Book method(s)
   for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) {
      if (it->second) {
         TString methodName = TString(it->first) + TString(" method");
         TString weightfile = dir + prefix + TString("_") + TString(it->first) + TString(".weights.xml");
         reader->BookMVA( methodName, weightfile ); 
      }
   }
   
   // Book output histograms
   UInt_t nbin = 100;
   TH1F   *histLk(0), *histLkD(0), *histLkPCA(0), *histLkKDE(0), *histLkMIX(0), *histPD(0), *histPDD(0);
   TH1F   *histPDPCA(0), *histPDEFoam(0), *histPDEFoamErr(0), *histPDEFoamSig(0), *histKNN(0), *histHm(0);
   TH1F   *histFi(0), *histFiG(0), *histFiB(0), *histLD(0), *histNn(0),*histNnbfgs(0),*histNnbnn(0);
   TH1F   *histNnC(0), *histNnT(0), *histBdt(0), *histBdtG(0), *histBdtD(0), *histRf(0), *histSVMG(0);
   TH1F   *histSVMP(0), *histSVML(0), *histFDAMT(0), *histFDAGA(0), *histCat(0), *histPBdt(0);

   if (Use["Likelihood"])    histLk      = new TH1F( "MVA_Likelihood",    "MVA_Likelihood",    nbin, -1, 1 );
   if (Use["LikelihoodD"])   histLkD     = new TH1F( "MVA_LikelihoodD",   "MVA_LikelihoodD",   nbin, -1, 0.9999 );
   if (Use["LikelihoodPCA"]) histLkPCA   = new TH1F( "MVA_LikelihoodPCA", "MVA_LikelihoodPCA", nbin, -1, 1 );
   if (Use["LikelihoodKDE"]) histLkKDE   = new TH1F( "MVA_LikelihoodKDE", "MVA_LikelihoodKDE", nbin,  -0.00001, 0.99999 );
   if (Use["LikelihoodMIX"]) histLkMIX   = new TH1F( "MVA_LikelihoodMIX", "MVA_LikelihoodMIX", nbin,  0, 1 );
   if (Use["PDERS"])         histPD      = new TH1F( "MVA_PDERS",         "MVA_PDERS",         nbin,  0, 1 );
   if (Use["PDERSD"])        histPDD     = new TH1F( "MVA_PDERSD",        "MVA_PDERSD",        nbin,  0, 1 );
   if (Use["PDERSPCA"])      histPDPCA   = new TH1F( "MVA_PDERSPCA",      "MVA_PDERSPCA",      nbin,  0, 1 );
   if (Use["KNN"])           histKNN     = new TH1F( "MVA_KNN",           "MVA_KNN",           nbin,  0, 1 );
   if (Use["HMatrix"])       histHm      = new TH1F( "MVA_HMatrix",       "MVA_HMatrix",       nbin, -0.95, 1.55 );
   if (Use["Fisher"])        histFi      = new TH1F( "MVA_Fisher",        "MVA_Fisher",        nbin, -4, 4 );
   if (Use["FisherG"])       histFiG     = new TH1F( "MVA_FisherG",       "MVA_FisherG",       nbin, -1, 1 );
   if (Use["BoostedFisher"]) histFiB     = new TH1F( "MVA_BoostedFisher", "MVA_BoostedFisher", nbin, -2, 2 );
   if (Use["LD"])            histLD      = new TH1F( "MVA_LD",            "MVA_LD",            nbin, -2, 2 );
   if (Use["MLP"])           histNn      = new TH1F( "MVA_MLP",           "MVA_MLP",           nbin, -1.25, 1.5 );
   if (Use["MLPBFGS"])       histNnbfgs  = new TH1F( "MVA_MLPBFGS",       "MVA_MLPBFGS",       nbin, -1.25, 1.5 );
   if (Use["MLPBNN"])        histNnbnn   = new TH1F( "MVA_MLPBNN",        "MVA_MLPBNN",        nbin, -1.25, 1.5 );
   if (Use["CFMlpANN"])      histNnC     = new TH1F( "MVA_CFMlpANN",      "MVA_CFMlpANN",      nbin,  0, 1 );
   if (Use["TMlpANN"])       histNnT     = new TH1F( "MVA_TMlpANN",       "MVA_TMlpANN",       nbin, -1.3, 1.3 );
   if (Use["BDT"])           histBdt     = new TH1F( "MVA_BDT",           "MVA_BDT",           nbin, -0.8, 0.8 );
   if (Use["BDTD"])          histBdtD    = new TH1F( "MVA_BDTD",          "MVA_BDTD",          nbin, -0.8, 0.8 );
   if (Use["BDTG"])          histBdtG    = new TH1F( "MVA_BDTG",          "MVA_BDTG",          nbin, -1.0, 1.0 );
   if (Use["RuleFit"])       histRf      = new TH1F( "MVA_RuleFit",       "MVA_RuleFit",       nbin, -2.0, 2.0 );
   if (Use["SVM_Gauss"])     histSVMG    = new TH1F( "MVA_SVM_Gauss",     "MVA_SVM_Gauss",     nbin,  0.0, 1.0 );
   if (Use["SVM_Poly"])      histSVMP    = new TH1F( "MVA_SVM_Poly",      "MVA_SVM_Poly",      nbin,  0.0, 1.0 );
   if (Use["SVM_Lin"])       histSVML    = new TH1F( "MVA_SVM_Lin",       "MVA_SVM_Lin",       nbin,  0.0, 1.0 );
   if (Use["FDA_MT"])        histFDAMT   = new TH1F( "MVA_FDA_MT",        "MVA_FDA_MT",        nbin, -2.0, 3.0 );
   if (Use["FDA_GA"])        histFDAGA   = new TH1F( "MVA_FDA_GA",        "MVA_FDA_GA",        nbin, -2.0, 3.0 );
   if (Use["Category"])      histCat     = new TH1F( "MVA_Category",      "MVA_Category",      nbin, -2., 2. );
   if (Use["Plugin"])        histPBdt    = new TH1F( "MVA_PBDT",          "MVA_BDT",           nbin, -0.8, 0.8 );

   // PDEFoam also returns per-event error, fill in histogram, and also fill significance
   if (Use["PDEFoam"]) {
      histPDEFoam    = new TH1F( "MVA_PDEFoam",       "MVA_PDEFoam",              nbin,  0, 1 );
      histPDEFoamErr = new TH1F( "MVA_PDEFoamErr",    "MVA_PDEFoam error",        nbin,  0, 1 );
      histPDEFoamSig = new TH1F( "MVA_PDEFoamSig",    "MVA_PDEFoam significance", nbin,  0, 10 );
   }

   // Book example histogram for probability (the other methods are done similarly)
   TH1F *probHistFi(0), *rarityHistFi(0);
   if (Use["Fisher"]) {
      probHistFi   = new TH1F( "MVA_Fisher_Proba",  "MVA_Fisher_Proba",  nbin, 0, 1 );
      rarityHistFi = new TH1F( "MVA_Fisher_Rarity", "MVA_Fisher_Rarity", nbin, 0, 1 );
   }

   // Prepare input tree (this must be replaced by your data source)
   // in this example, there is a toy tree with signal and one with background events
   // we'll later on use only the "signal" events for the test in this example.
   //   
   TFile *input(0);
   //TString fname = "./histtrain_nue_1.root";   
   TString fname = "numu.train.root";   

   if (!gSystem->AccessPathName( fname )) 
      input = TFile::Open( fname ); // check if file in local directory exists
   else    
      input = TFile::Open( "http://root.cern.ch/files/tmva_class_example.root" ); // if not: download from ROOT server
   
   if (!input) {
      std::cout << "ERROR: could not open data file" << std::endl;
      exit(1);
   }
   std::cout << "--- TMVAClassificationApp    : Using input file: " << input->GetName() << std::endl;
   
   // --- Event loop

   // Prepare the event tree
   // - here the variable names have to corresponds to your tree
   // - you can use the same variables as above which is slightly faster,
   //   but of course you can use different ones and copy the values inside the event loop
   //
   std::cout << "--- Select signal sample" << std::endl;
   TTree* theTree = (TTree*)input->Get("mvaselect/all");
   //Float_t userVar1, userVar2;
   //theTree->SetBranchAddress( "var1", &userVar1 );
   float ntrk;
   float nshw;
   int itype;
   float weight;
   int run, subrun, event;
   float showerlength;
   theTree->SetBranchAddress("itype",&itype);
   theTree->SetBranchAddress("evtcharge", &evtcharge);
   theTree->SetBranchAddress("ntrack", &ntrk);
   theTree->SetBranchAddress("maxtrklength", &maxtrklength);
   theTree->SetBranchAddress("avgtrklength", &avgtrklength);
   theTree->SetBranchAddress("trkdedx", &trkdedx);
   theTree->SetBranchAddress("trkrch", &trkrch);
   theTree->SetBranchAddress("trkrt", &trkrt);
   theTree->SetBranchAddress("trkfr", &trkfr);
   theTree->SetBranchAddress("trkpida", &trkpida);
   theTree->SetBranchAddress("fract_5_wires",&fract_5_wires);
   theTree->SetBranchAddress("fract_10_wires",&fract_10_wires);
   theTree->SetBranchAddress("fract_50_wires",&fract_50_wires);
   theTree->SetBranchAddress("fract_100_wires",&fract_100_wires);
   theTree->SetBranchAddress("trkcosx",&trkcosx);
   theTree->SetBranchAddress("trkcosy",&trkcosy);
   theTree->SetBranchAddress("trkcosz",&trkcosz);
   theTree->SetBranchAddress("ET",&ET);
   theTree->SetBranchAddress("nshower", &nshw);
   theTree->SetBranchAddress("showerdedx", &showerdedx);
   theTree->SetBranchAddress("eshower", &eshower);
   theTree->SetBranchAddress("frshower", &frshower);
   theTree->SetBranchAddress("nhitspershw",&nhitspershw);
   theTree->SetBranchAddress("shwlength",&showerlength);
   theTree->SetBranchAddress("shwmax",&shwmax);
   theTree->SetBranchAddress("shwdisx",&shwdisx);
   theTree->SetBranchAddress("shwdisy",&shwdisy);
   theTree->SetBranchAddress("shwdisz",&shwdisz);
   theTree->SetBranchAddress("shwcosx",&shwcosx);
   theTree->SetBranchAddress("shwcosy",&shwcosy);
   theTree->SetBranchAddress("shwcosz",&shwcosz);
   theTree->SetBranchAddress("weight",&weight);
   theTree->SetBranchAddress("run",&run);
   theTree->SetBranchAddress("subrun",&subrun);
   theTree->SetBranchAddress("event",&event);

   TH1F *bdt_sig = new TH1F("bdt_sig","bdt_sig",100,-1.25,1.25);
   TH1F *bdt_nc = new TH1F("bdt_nc","bdt_nc",100,-1.25,1.25);
   TH1F *bdt_cc = new TH1F("bdt_cc","bdt_cc",100,-1.25,1.25);
   TH1F *bdt_bnue = new TH1F("bdt_bnue","bdt_bnue",100,-1.25,1.25);
   TH1F *bdt_nutau = new TH1F("bdt_nutau","bdt_nutau",100,-1.25,1.25);

   // Efficiency calculator for cut method
   Int_t    nSelCutsGA = 0;
   Double_t effS       = 0.7;
   const int ncuts = 100;
   float cut[ncuts];
   float sigall = 0, ncall = 0, ccall = 0, bnueall = 0, nutauall = 0;
   float sig[ncuts];
   float nc[ncuts];
   float cc[ncuts];
   float bnue[ncuts];
   float nutau[ncuts];
   for (int i = 0; i<ncuts; ++i){
     cut[i] = -1+i*0.02;
     sig[i] = 0;
     nc[i] = 0;
     cc[i] = 0;
     bnue[i] = 0;
     nutau[i] = 0;
   }
   std::vector<Float_t> vecVar(4); // vector for EvaluateMVA tests

   std::cout << "--- Processing: " << theTree->GetEntries() << " events" << std::endl;
   TStopwatch sw;
   sw.Start();
   Long64_t nent = theTree->GetEntries();
   if( nent > 50000 ) nent = 50000;
   for (Long64_t ievt=0; ievt<nent;ievt++) {

      if (ievt%1000 == 0) std::cout << "--- ... Processing event: " << ievt << std::endl;

      theTree->GetEntry(ievt);
      ntrack = ntrk;
      nshower = nshw;
      shwlength = showerlength;
      //var1 = userVar1 + userVar2;
      //var2 = userVar1 - userVar2;

      // --- Return the MVA outputs and fill into histograms

      if (Use["CutsGA"]) {
         // Cuts is a special case: give the desired signal efficienciy
         Bool_t passed = reader->EvaluateMVA( "CutsGA method", effS );
         if (passed) nSelCutsGA++;
      }

      if (Use["Likelihood"   ])   histLk     ->Fill( reader->EvaluateMVA( "Likelihood method"    ) );
      if (Use["LikelihoodD"  ])   histLkD    ->Fill( reader->EvaluateMVA( "LikelihoodD method"   ) );
      if (Use["LikelihoodPCA"])   histLkPCA  ->Fill( reader->EvaluateMVA( "LikelihoodPCA method" ) );
      if (Use["LikelihoodKDE"])   histLkKDE  ->Fill( reader->EvaluateMVA( "LikelihoodKDE method" ) );
      if (Use["LikelihoodMIX"])   histLkMIX  ->Fill( reader->EvaluateMVA( "LikelihoodMIX method" ) );
      if (Use["PDERS"        ])   histPD     ->Fill( reader->EvaluateMVA( "PDERS method"         ) );
      if (Use["PDERSD"       ])   histPDD    ->Fill( reader->EvaluateMVA( "PDERSD method"        ) );
      if (Use["PDERSPCA"     ])   histPDPCA  ->Fill( reader->EvaluateMVA( "PDERSPCA method"      ) );
      if (Use["KNN"          ])   histKNN    ->Fill( reader->EvaluateMVA( "KNN method"           ) );
      if (Use["HMatrix"      ])   histHm     ->Fill( reader->EvaluateMVA( "HMatrix method"       ) );
      if (Use["Fisher"       ])   histFi     ->Fill( reader->EvaluateMVA( "Fisher method"        ) );
      if (Use["FisherG"      ])   histFiG    ->Fill( reader->EvaluateMVA( "FisherG method"       ) );
      if (Use["BoostedFisher"])   histFiB    ->Fill( reader->EvaluateMVA( "BoostedFisher method" ) );
      if (Use["LD"           ])   histLD     ->Fill( reader->EvaluateMVA( "LD method"            ) );
      if (Use["MLP"          ])   histNn     ->Fill( reader->EvaluateMVA( "MLP method"           ) );
      if (Use["MLPBFGS"      ])   histNnbfgs ->Fill( reader->EvaluateMVA( "MLPBFGS method"       ) );
      if (Use["MLPBNN"       ])   histNnbnn  ->Fill( reader->EvaluateMVA( "MLPBNN method"        ) );
      if (Use["CFMlpANN"     ])   histNnC    ->Fill( reader->EvaluateMVA( "CFMlpANN method"      ) );
      if (Use["TMlpANN"      ])   histNnT    ->Fill( reader->EvaluateMVA( "TMlpANN method"       ) );
      if (Use["BDT"          ])   histBdt    ->Fill( reader->EvaluateMVA( "BDT method"           ) );
      if (Use["BDTD"         ])   histBdtD   ->Fill( reader->EvaluateMVA( "BDTD method"          ) );
      if (Use["BDTG"         ])   histBdtG   ->Fill( reader->EvaluateMVA( "BDTG method"          ) );
      if (Use["RuleFit"      ])   histRf     ->Fill( reader->EvaluateMVA( "RuleFit method"       ) );
      if (Use["SVM_Gauss"    ])   histSVMG   ->Fill( reader->EvaluateMVA( "SVM_Gauss method"     ) );
      if (Use["SVM_Poly"     ])   histSVMP   ->Fill( reader->EvaluateMVA( "SVM_Poly method"      ) );
      if (Use["SVM_Lin"      ])   histSVML   ->Fill( reader->EvaluateMVA( "SVM_Lin method"       ) );
      if (Use["FDA_MT"       ])   histFDAMT  ->Fill( reader->EvaluateMVA( "FDA_MT method"        ) );
      if (Use["FDA_GA"       ])   histFDAGA  ->Fill( reader->EvaluateMVA( "FDA_GA method"        ) );
      if (Use["Category"     ])   histCat    ->Fill( reader->EvaluateMVA( "Category method"      ) );
      if (Use["Plugin"       ])   histPBdt   ->Fill( reader->EvaluateMVA( "P_BDT method"         ) );

      if (itype==0) bdt_sig->Fill(reader->EvaluateMVA( "BDTG method"           ), weight);
      if (itype==1) bdt_nc->Fill(reader->EvaluateMVA( "BDTG method"           ), weight);
      if (itype==2) bdt_cc->Fill(reader->EvaluateMVA( "BDTG method"           ), weight);
      if (itype==3) bdt_bnue->Fill(reader->EvaluateMVA( "BDTG method"           ), weight);
      if (itype==4) bdt_nutau->Fill(reader->EvaluateMVA( "BDTG method"           ), weight);

      //std::cout<<run<<" "<<event<<" "<<reader->EvaluateMVA( "BDTG method"           )<<std::endl;

      if(SelectNuE){
	if (itype==0) sigall+= weight;
	if (itype==2) ccall+= weight;
      } else {
	if (itype==2) sigall+= weight;
	if (itype==0) ccall+= weight;
      }
      if (itype==1) ncall+= weight;
      if (itype==3) bnueall+= weight;
      if (itype==4) nutauall+= weight;
      //std::cout<<evtcharge<<" "<<ntrack<<" "<<maxtrklength<<" "<<avgtrklength<<" "<<trkdedx<<" "<<trkrch<<" "<<trkrt<<" "<<trkfr<<" "<<trkpida<<" "<<fract_5_wires<<" "<<fract_10_wires<<" "<<fract_50_wires<<" "<<fract_100_wires<<" "<<trkcosx<<" "<<trkcosy<<" "<<trkcosz<<" "<<ET<<" "<<nshower<<" "<<showerdedx<<" "<<eshower<<" "<<frshower<<" "<<nhitspershw<<" "<<shwlength<<" "<<shwmax<<" "<<shwdisx<<" "<<shwdisy<<" "<<shwdisz<<" "<<shwcosx<<" "<<shwcosy<<" "<<shwcosz<<std::endl;
      //cout<<ntrack<<" "<<trkdedx<<" "<<trkrch<<" "<<reader->EvaluateMVA( "BDTG method"           )<<endl;
      //break;
      for (int icut = 0; icut<ncuts; ++icut){
	if (//maxtrklength<300&&
	    //evtcharge>1e5&&
	    //evtcharge<5e6&&
	    //eshower>0.5&&
	    reader->EvaluateMVA( "BDTG method"           )>cut[icut]){
      if(SelectNuE){
	  if (itype==0) sig[icut]+= weight;
	  if (itype==2) cc[icut]+= weight;
      } else {
	  if (itype==2) sig[icut]+= weight;
	  if (itype==0) cc[icut]+= weight;
      }
	  if (itype==1) nc[icut]+= weight;
	  if (itype==3) bnue[icut]+= weight;
	  if (itype==4) nutau[icut]+= weight;
	  if (icut==ncuts/2) outfile<<run<<" "<<subrun<<" "<<event<<endl;
	}
      }	
      //if (itype==0&&nshower==0) outfile<<run<<" "<<subrun<<" "<<event<<endl;
      // Retrieve also per-event error
      if (Use["PDEFoam"]) {
         Double_t val = reader->EvaluateMVA( "PDEFoam method" );
         Double_t err = reader->GetMVAError();
         histPDEFoam   ->Fill( val );
         histPDEFoamErr->Fill( err );         
         if (err>1.e-50) histPDEFoamSig->Fill( val/err );
      }         

      // Retrieve probability instead of MVA output
      if (Use["Fisher"])   {
         probHistFi  ->Fill( reader->GetProba ( "Fisher method" ) );
         rarityHistFi->Fill( reader->GetRarity( "Fisher method" ) );
      }
   }

   // Get elapsed time
   sw.Stop();
   std::cout << "--- End of event loop: "; sw.Print();

   // Get efficiency for cuts classifier
   if (Use["CutsGA"]) std::cout << "--- Efficiency for CutsGA method: " << double(nSelCutsGA)/theTree->GetEntries()
                                << " (for a required signal efficiency of " << effS << ")" << std::endl;

   if (Use["CutsGA"]) {

      // test: retrieve cuts for particular signal efficiency
      // CINT ignores dynamic_casts so we have to use a cuts-secific Reader function to acces the pointer  
      TMVA::MethodCuts* mcuts = reader->FindCutsMVA( "CutsGA method" ) ;

      if (mcuts) {      
         std::vector<Double_t> cutsMin;
         std::vector<Double_t> cutsMax;
         mcuts->GetCuts( 0.7, cutsMin, cutsMax );
         std::cout << "--- -------------------------------------------------------------" << std::endl;
         std::cout << "--- Retrieve cut values for signal efficiency of 0.7 from Reader" << std::endl;
         for (UInt_t ivar=0; ivar<cutsMin.size(); ivar++) {
            std::cout << "... Cut: " 
                      << cutsMin[ivar] 
                      << " < \"" 
                      << mcuts->GetInputVar(ivar)
                      << "\" <= " 
                      << cutsMax[ivar] << std::endl;
         }
         std::cout << "--- -------------------------------------------------------------" << std::endl;
      }
   }

   // --- Write histograms

   TFile *target  = new TFile( "TMVApp.root","RECREATE" );
   if (Use["Likelihood"   ])   histLk     ->Write();
   if (Use["LikelihoodD"  ])   histLkD    ->Write();
   if (Use["LikelihoodPCA"])   histLkPCA  ->Write();
   if (Use["LikelihoodKDE"])   histLkKDE  ->Write();
   if (Use["LikelihoodMIX"])   histLkMIX  ->Write();
   if (Use["PDERS"        ])   histPD     ->Write();
   if (Use["PDERSD"       ])   histPDD    ->Write();
   if (Use["PDERSPCA"     ])   histPDPCA  ->Write();
   if (Use["KNN"          ])   histKNN    ->Write();
   if (Use["HMatrix"      ])   histHm     ->Write();
   if (Use["Fisher"       ])   histFi     ->Write();
   if (Use["FisherG"      ])   histFiG    ->Write();
   if (Use["BoostedFisher"])   histFiB    ->Write();
   if (Use["LD"           ])   histLD     ->Write();
   if (Use["MLP"          ])   histNn     ->Write();
   if (Use["MLPBFGS"      ])   histNnbfgs ->Write();
   if (Use["MLPBNN"       ])   histNnbnn  ->Write();
   if (Use["CFMlpANN"     ])   histNnC    ->Write();
   if (Use["TMlpANN"      ])   histNnT    ->Write();
   if (Use["BDT"          ])   histBdt    ->Write();
   if (Use["BDTD"         ])   histBdtD   ->Write();
   if (Use["BDTG"         ])   histBdtG   ->Write(); 
   if (Use["RuleFit"      ])   histRf     ->Write();
   if (Use["SVM_Gauss"    ])   histSVMG   ->Write();
   if (Use["SVM_Poly"     ])   histSVMP   ->Write();
   if (Use["SVM_Lin"      ])   histSVML   ->Write();
   if (Use["FDA_MT"       ])   histFDAMT  ->Write();
   if (Use["FDA_GA"       ])   histFDAGA  ->Write();
   if (Use["Category"     ])   histCat    ->Write();
   if (Use["Plugin"       ])   histPBdt   ->Write();
   bdt_sig->Write();
   bdt_nc->Write();
   bdt_cc->Write();
   bdt_bnue->Write();
   bdt_nutau->Write();
   // Write also error and significance histos
   if (Use["PDEFoam"]) { histPDEFoam->Write(); histPDEFoamErr->Write(); histPDEFoamSig->Write(); }

   // Write also probability hists
   if (Use["Fisher"]) { if (probHistFi != 0) probHistFi->Write(); if (rarityHistFi != 0) rarityHistFi->Write(); }
   target->Close();

   std::cout << "--- Created root file: \"TMVApp.root\" containing the MVA output histograms" << std::endl;
  
   delete reader;
    
   std::cout << "==> TMVAClassificationApplication is done!" << endl << std::endl;
//   cout<<sig<<"/"<<sigall<<"="<<sig/sigall<<endl;
//   cout<<nc<<"/"<<ncall<<"="<<nc/ncall<<endl;
//   cout<<cc<<"/"<<ccall<<"="<<cc/ccall<<endl;
//   cout<<bnue<<"/"<<bnueall<<"="<<bnue/bnueall<<endl;
//   cout<<nutau<<"/"<<nutauall<<"="<<nutau/nutauall<<endl;
   outfile.close();
   for (int i = 0; i<ncuts; ++i){
     sig[i]/=sigall;
     nc[i]/=ncall; nc[i] = 1-nc[i];
     cc[i]/=ccall; cc[i] = 1-cc[i];
     bnue[i]/=bnueall; bnue[i] = 1-bnue[i];
     nutau[i]/=nutauall; nutau[i] = 1-nutau[i];
   }
   TCanvas *c1 = new TCanvas("c1","c1");
   TGraph *rocnc = new TGraph(ncuts,sig,nc);
   rocnc->SetLineWidth(3);
   rocnc->SetTitle("");
   rocnc->GetXaxis()->SetTitle("Signal Efficiency");
   rocnc->GetYaxis()->SetTitle("NC Background Rejection");
   rocnc->GetXaxis()->SetRangeUser(0,1.1);
   rocnc->GetYaxis()->SetRangeUser(0,1.1);
   rocnc->Draw("ac");
   cout<<rocnc->Eval(0.8)<<endl;
   TMarker *m1 = new TMarker(0.8,0.99,29);
   m1->SetMarkerSize(2.5);
   m1->SetMarkerColor(2);
   //m1->Draw();

   TCanvas *c2 = new TCanvas("c2","c2");
   TGraph *roccc = new TGraph(ncuts,sig,cc);
   roccc->SetLineWidth(3);
   roccc->SetTitle("");
   roccc->GetXaxis()->SetTitle("Signal Efficiency");
   roccc->GetYaxis()->SetTitle("CC Background Rejection");
   roccc->GetXaxis()->SetRangeUser(0,1.1);
   roccc->GetYaxis()->SetRangeUser(0,1.1);
   roccc->Draw("ac");
   cout<<roccc->Eval(0.8)<<endl;
   TMarker *m2 = new TMarker(0.8,0.999,29);
   m2->SetMarkerSize(2.5);
   m2->SetMarkerColor(2);
   //m2->Draw();

   TCanvas *c3 = new TCanvas("c3","c3");
   TGraph *rocbnue = new TGraph(ncuts,sig,bnue);
   rocbnue->SetLineWidth(3);
   rocbnue->GetXaxis()->SetTitle("Signal Efficiency");
   rocbnue->GetYaxis()->SetTitle("Beam nue Background Rejection");
   rocbnue->Draw("ac");
   cout<<rocbnue->Eval(0.8)<<endl;

   TCanvas *c4 = new TCanvas("c4","c4");
   TGraph *rocnutau = new TGraph(ncuts,sig,nutau);
   rocnutau->SetLineWidth(3);
   rocnutau->GetXaxis()->SetTitle("Signal Efficiency");
   rocnutau->GetYaxis()->SetTitle("Nutau Background Rejection");
   rocnutau->Draw("ac");
   cout<<rocnutau->Eval(0.8)<<endl;

   c1->Print("rocnc.pdf");
   c2->Print("roccc.pdf");
   c3->Print("rocbnue.pdf");
   c4->Print("rocnutau.pdf");
} 
