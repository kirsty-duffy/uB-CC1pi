// @(#)root/tmva $Id$
/**********************************************************************************
 * Project   : TMVA - a ROOT-integrated toolkit for multivariate data analysis    *
 * Package   : TMVA                                                               *
 * Root Macro: TMVAClassification                                                 *
 *                                                                                *
 * This macro provides examples for the training and testing of the               *
 * TMVA classifiers.                                                              *
 *                                                                                *
 * As input data is used a toy-MC sample consisting of four Gaussian-distributed  *
 * and linearly correlated input variables.                                       *
 *                                                                                *
 * The methods to be used can be switched on and off by means of booleans, or     *
 * via the prompt command, for example:                                           *
 *                                                                                *
 *    root -l ./TMVAClassification.C\(\"Fisher,Likelihood\"\)                     *
 *                                                                                *
 * (note that the backslashes are mandatory)                                      *
 * If no method given, a default set of classifiers is used.                      *
 *                                                                                *
 * The output file "TMVA.root" can be analysed with the use of dedicated          *
 * macros (simply say: root -l <macro.C>), which can be conveniently              *
 * invoked through a GUI that will appear at the end of the run of this macro.    *
 * Launch the GUI via the command:                                                *
 *                                                                                *
 *    root -l ./TMVAGui.C                                                         *
 *                                                                                *
 **********************************************************************************/

#include <cstdlib>
#include <iostream>
#include <map>
#include <string>
#include <math.h>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"


#if not defined(__CINT__) || defined(__MAKECINT__)
// needs to be included when makecint runs (ACLIC)
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#endif

void TMVAClassification( TString myMethodList = "" )
{
   // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
   // if you use your private .rootrc, or run from a different directory, please copy the
   // corresponding lines from .rootrc

   // methods to be processed can be given as an argument; use format:
   //
   // mylinux~> root -l TMVAClassification.C\(\"myMethod1,myMethod2,myMethod3\"\)
   //
   // if you like to use a method via the plugin mechanism, we recommend using
   //
   // mylinux~> root -l TMVAClassification.C\(\"P_myMethod\"\)
   // (an example is given for using the BDT as plugin (see below),
   // but of course the real application is when you write your own
   // method based)

   //---------------------------------------------------------------
   // This loads the library
   TMVA::Tools::Instance();

   // to get access to the GUI and all tmva macros
    // TString tmva_dir(TString(gRootDir) + "/tmva");
    // if(gSystem->Getenv("TMVASYS"))
    //    tmva_dir = TString(gSystem->Getenv("TMVASYS"));
    TString tmva_dir("/uboone/app/users/kduffy/CC1pi/CC1pi_uboonecode/srcs/uboonecode/uboone/CC1pi/MVA/:");
    // gROOT->SetMacroPath(tmva_dir + gROOT->GetMacroPath() );
    gROOT->SetMacroPath(tmva_dir);
    gROOT->ProcessLine(".L TMVAGui.C");



   // Default MVA methods to be trained + tested
   std::map<std::string,int> Use;

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
   Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting
   //
   // --- Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
   Use["RuleFit"]         = 0;
   // ---------------------------------------------------------------

   std::cout << std::endl;
   std::cout << "==> Start TMVAClassification" << std::endl;

   // Select methods (don't look at this code - not of interest)
   if (myMethodList != "") {
      for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

      std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
      for (UInt_t i=0; i<mlist.size(); i++) {
         std::string regMethod(mlist[i]);

         if (Use.find(regMethod) == Use.end()) {
            std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
            for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
            std::cout << std::endl;
            return;
         }
         Use[regMethod] = 1;
      }
   }

   // --------------------------------------------------------------------------------------------------

   // --- Here the preparation phase begins

   // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
   TString outfileName( "TMVA_mupi_output.root" );
   TFile* outputFile = TFile::Open( outfileName, "RECREATE" );

   // Create the factory object. Later you can choose the methods
   // whose performance you'd like to investigate. The factory is
   // the only TMVA object you have to interact with
   //
   // The first argument is the base of the name of all the
   // weightfiles in the directory weight/
   //
   // The second argument is the output file for the training results
   // All TMVA output can be suppressed by removing the "!" (not) in
   // front of the "Silent" argument in the option string
   TMVA::Factory *factory_cont = new TMVA::Factory( "TMVAClassification_cont", outputFile, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Auto" );
   TMVA::DataLoader *dataloader_cont = new TMVA::DataLoader("dataset_cont");
   TMVA::Factory *factory_exiting = new TMVA::Factory( "TMVAClassification_exiting", outputFile, "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Auto" );
   TMVA::DataLoader *dataloader_exiting = new TMVA::DataLoader("dataset_exiting");

   // If you wish to modify default settings
   // (please check "src/Config.h" to see all available global options)
   //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
   //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

   // Define the input variables that shall be used for the MVA training
   // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
   // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]

   // branch, name, unit, type
   dataloader_cont->AddVariable("length_over_startend", "length_over_startend", "", 'D');
   dataloader_cont->AddVariable("track_length_over_longest","track_length_over_longest","",'D');
   dataloader_cont->AddVariable("ndaughters", "ndaughters", "", 'D');
   dataloader_cont->AddVariable("lnLmipovermu", "lnLmipovermu", "", 'D');
   dataloader_cont->AddVariable("perc_used_hits", "perc_used_hits", "", 'D');
   // dataloader_cont->AddVariable("residuals_mean", "residuals_mean", "", 'D');
   // dataloader_cont->AddVariable("residuals_stddev", "residuals_stddev", "", 'D');
   dataloader_cont->AddVariable("MCS_pi_maxScatter", "MCS_pi_maxScatter", "", 'D');
   dataloader_cont->AddVariable("MCS_pi_meanScatter", "MCS_pi_meanScatter", "", 'D');
   dataloader_cont->AddVariable("n_unused_hits_nearend", "n_unused_hits_nearend", "", 'D');
   dataloader_cont->AddVariable("unmatched_charge_nearend_plane2", "unmatched_charge_nearend_plane2", "", 'D');

   // You can add so-called "Spectator variables", which are not used in the MVA training,
   // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
   // input variables, the response values of all trained MVAs, and the spectator variables
   // factory->AddSpectator( "spec1 := var1*2",  "Spectator 1", "units", 'F' );
   // factory->AddSpectator( "spec2 := var1*3",  "Spectator 2", "units", 'F' );
   dataloader_cont->AddSpectator("track_length", "track_length", "", 'F');
   dataloader_cont->AddSpectator("pandoraclassedastrack", "pandoraclassedastrack", "", 'F');
   dataloader_cont->AddSpectator("lnLmipoverpi", "lnLmipoverpi", "", 'F');

   // Now do the same for exiting sample
   // branch, name, unit, type
   dataloader_exiting->AddVariable("length_over_startend", "length_over_startend", "", 'D');
   dataloader_exiting->AddVariable("track_length_over_longest","track_length_over_longest","",'D');
   dataloader_exiting->AddVariable("perc_used_hits", "perc_used_hits", "", 'D');
   // dataloader_exiting->AddVariable("residuals_mean", "residuals_mean", "", 'D');
   // dataloader_exiting->AddVariable("residuals_stddev", "residuals_stddev", "", 'D');
   dataloader_exiting->AddVariable("MCS_pi_maxScatter", "MCS_pi_maxScatter", "", 'D');
   dataloader_exiting->AddVariable("MCS_pi_meanScatter", "MCS_pi_meanScatter", "", 'D');

   dataloader_exiting->AddSpectator("track_length", "track_length", "", 'F'); dataloader_exiting->AddSpectator("pandoraclassedastrack", "pandoraclassedastrack", "", 'F');
   dataloader_exiting->AddSpectator("lnLmipovermu", "lnLmipovermu", "", 'F');
   dataloader_exiting->AddSpectator("lnLmipoverpi", "lnLmipoverpi", "", 'F');
   dataloader_exiting->AddSpectator("ndaughters", "ndaughters", "", 'F');
   dataloader_exiting->AddSpectator("n_unused_hits_nearend", "n_unused_hits_nearend", "", 'F');
   dataloader_exiting->AddSpectator("unmatched_charge_nearend_plane2", "unmatched_charge_nearend_plane2", "", 'F');

   // Read training and test data
   // (it is also possible to use ASCII format as input -> see TMVA Users Guide)

   //TString fname = "./histtrain_nue.root";
   //TString fname = "./train.anue.root";
   //TString fname = "./numu.train.root";
   //TString fname = "/dune/data/users/dbrailsf/oscillations/nu_mu/v06_11_00_B_prodgenie_nuall_dune10kt_1x2x6_hadd.root";
   //TString fname = "/dune/data/users/dbrailsf/oscillations/nu_mu/production/v06_11_00/B/mc/fdsenseopt/v06_11_00_B_prodgenie_nuall_dune10kt_1x2x6_hadd.root";
   //TString fname = "/dune/data/users/dbrailsf/oscillations/nu_mu/production/v06_11_00/B/mc/fdsenseopt/v06_11_00_B_prodgenie_nuall_dune10kt_1x2x6_hadd.root";
   //TString fname = "/dune/data/users/dbrailsf/oscillations/nu_mu/production/v06_18_00/A/mc/fdsenseopt/v06_18_00_A_prodgenie_nuall_dune10kt_1x2x6_hadd.root";
   //TString fname = "/dune/data/users/dbrailsf/oscillations/nu_mu/production/v06_34_00/B/mc/fdsenseopt/v06_34_00_B_prodgenie_nuall_dune10kt_1x2x6_hadd.root";
   //TString fname = "/dune/data/users/dbrailsf/oscillations/nu_mu/production/v06_34_00/D/mc/fdsenseopt/v06_34_00_D_prodgenie_nuall_dune10kt_1x2x6_hadd.root";
   //TString fname_mu = "/sbnd/app/users/dbrailsf/warwick_pid/select/flat_tree_mu.root";
   //TString fname_proton = "/sbnd/app/users/dbrailsf/warwick_pid/select/flat_tree_proton.root";
   //TString fname_electron = "/sbnd/app/users/dbrailsf/warwick_pid/select/flat_tree_electron.root";
   //TString fname_pi0 = "/sbnd/app/users/dbrailsf/warwick_pid/select/flat_tree_pi0.root";

   //TString fname = "/dune/data/users/dbrailsf/oscillations/nu_mu/production/v06_18_00/A/mc/fdsenseopt/v06_18_00_A_prodgenie_anuall_dune10kt_1x2x6_hadd.root";
   //TString fname = "/dune/data/users/dbrailsf/oscillations/nu_mu/production/v06_18_00/E/mc/fdsenseopt/v06_18_00_E_prodgenie_nuall_dune10kt_1x2x6_hadd.root";

   //TString fname = "./anumu.train.root";
   //TString fname = "./nuetrain.def.root";
   //TString fname = "./nuetrain.3mm.root";
   //TString fname = "./nuetrain.45deg.root";

   //TString fname = "/pnfs/dune/persistent/TaskForce_AnaTree/far/train_nue_disambig/anahist.root";
   //TString fname = "/pnfs/dune/persistent/TaskForce_AnaTree/far/train_nue/anahist.root";
   //TString fname = "/pnfs/dune/persistent/TaskForce_AnaTree/far/train_numu/anahist.root";
   //TString fname = "/pnfs/dune/persistent/TaskForce_AnaTree/far/train_anue/anahist.root";
   //TString fname = "/pnfs/dune/persistent/TaskForce_AnaTree/far/train_anumu/anahist.root";

   //TString fname = "/pnfs/dune/persistent/TaskForce_AnaTree/far/train_numu/anahist.root";

//   if (gSystem->AccessPathName( fname ))  // file does not exist in local directory
//      gSystem->Exec("wget http://root.cern.ch/files/tmva_class_example.root");

   TFile *input = TFile::Open("/uboone/app/users/kduffy/CC1pi/CC1pi_uboonecode/srcs/uboonecode/uboone/CC1pi/ignore/2019-01-28/MVA_Trees_mupi.root");



//   std::cout << "--- TMVAClassification       : Using input file: " << input->GetName() << std::endl;

   // --- Register the training and test trees
/*
   TTree *muon_contained = (TTree*) input->Get("muon_contained");
   TTree *muon_uncontained = (TTree*) input->Get("muon_uncontained");
   TTree *pion_contained = (TTree*) input->Get("pion_contained");
   TTree *pion_uncontained = (TTree*) input->Get("pion_uncontained");
*/
   TTree *muon_cont = (TTree*) input->Get("muon_cont");
   TTree *pion_cont = (TTree*) input->Get("pion_cont");
   TTree *muon_exit = (TTree*) input->Get("muon_exit");
   TTree *pion_exit = (TTree*) input->Get("pion_exit");


   // global event weights per tree (see below for setting event-wise weights)
   Double_t signalWeight     = 1.0;
   Double_t backgroundWeight = 1.0;

   // You can add an arbitrary number of signal or background trees
/*
   dataloader->AddSignalTree(muon_contained, 1./(muon_contained->GetEntries()));
   dataloader->AddSignalTree(muon_uncontained, 1./(muon_uncontained->GetEntries()));
   dataloader->AddSignalTree(pion_contained, 1./(pion_contained->GetEntries()));
   dataloader->AddSignalTree(pion_uncontained, 1./(pion_uncontained->GetEntries()));
*/
   dataloader_cont->AddSignalTree(pion_cont, 1./(pion_cont->GetEntries()));
   dataloader_cont->AddBackgroundTree(muon_cont, 1./(muon_cont->GetEntries()));

   dataloader_exiting->AddSignalTree(pion_exit, 1./(pion_exit->GetEntries()));
   dataloader_exiting->AddBackgroundTree(muon_exit, 1./(muon_exit->GetEntries()));


   // To give different trees for training and testing, do as follows:
   //    factory->AddSignalTree( signalTrainingTree, signalTrainWeight, "Training" );
   //    factory->AddSignalTree( signalTestTree,     signalTestWeight,  "Test" );

   // Use the following code instead of the above two or four lines to add signal and background
   // training and test events "by hand"
   // NOTE that in this case one should not give expressions (such as "var1+var2") in the input
   //      variable definition, but simply compute the expression before adding the event
   //
   //     // --- begin ----------------------------------------------------------
   //     std::vector<Double_t> vars( 4 ); // vector has size of number of input variables
   //     Float_t  treevars[4], weight;
   //
   //     // Signal
   //     for (UInt_t ivar=0; ivar<4; ivar++) signal->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
   //     for (UInt_t i=0; i<signal->GetEntries(); i++) {
   //        signal->GetEntry(i);
   //        for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
   //        // add training and test events; here: first half is training, second is testing
   //        // note that the weight can also be event-wise
   //        if (i < signal->GetEntries()/2.0) factory->AddSignalTrainingEvent( vars, signalWeight );
   //        else                              factory->AddSignalTestEvent    ( vars, signalWeight );
   //     }
   //
   //     // Background (has event weights)
   //     background->SetBranchAddress( "weight", &weight );
   //     for (UInt_t ivar=0; ivar<4; ivar++) background->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
   //     for (UInt_t i=0; i<background->GetEntries(); i++) {
   //        background->GetEntry(i);
   //        for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
   //        // add training and test events; here: first half is training, second is testing
   //        // note that the weight can also be event-wise
   //        if (i < background->GetEntries()/2) factory->AddBackgroundTrainingEvent( vars, backgroundWeight*weight );
   //        else                                factory->AddBackgroundTestEvent    ( vars, backgroundWeight*weight );
   //     }
         // --- end ------------------------------------------------------------
   //
   // --- end of tree registration

   // Set individual event weights (the variables must exist in the original TTree)
   //    for signal    : factory->SetSignalWeightExpression    ("weight1*weight2");
   //    for background: factory->SetBackgroundWeightExpression("weight1*weight2");
   //factory->SetSignalWeightExpression( "weight" );
   //factory->SetBackgroundWeightExpression( "weight" );
   //factory->SetSignalWeightExpression( "projected_weight" );
   //factory->SetBackgroundWeightExpression( "projected_weight" );

   //factory->SetSignalWeightExpression( "projected_pot*energy_weight" );
   //factory->SetBackgroundWeightExpression( "projected_pot*energy_weight" );



   // Apply additional cuts on the signal and background samples (can be different)
   // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
   // for example: TCut mycutb = "abs(var1)<0.5";
   TCut mycuts = "!isnan(residuals_mean) && !isnan(residuals_stddev) && !isnan(length_over_startend) && !isnan(track_length) && !isnan(ndaughters) && !isnan(lnLmipovermu) && !isnan(lnLmipoverpi) && !isnan(pandoraclassedastrack) && !isnan(perc_used_hits) && !isnan(MCS_pi_maxScatter) && !isnan(MCS_pi_meanScatter) && track_length_over_longest<=1";
   TCut mycutb = mycuts;


   // Tell the factory how to use the training and testing events
   //
   // Use 1000 events for training and the rest for testing:
   dataloader_cont->PrepareTrainingAndTestTree( mycuts, mycutb,  "nTrain_Signal=1000:nTrain_Background=1000:nTest_Signal=0:nTest_Background=0:SplitMode=Random:NormMode=None:!V");
   dataloader_exiting->PrepareTrainingAndTestTree( mycuts, mycutb,  "nTrain_Signal=1000:nTrain_Background=1000:nTest_Signal=0:nTest_Background=0:SplitMode=Random:NormMode=None:!V");


   // To also specify the number of testing events, use:
   //factory->PrepareTrainingAndTestTree( mycuts,mycutb,
   //					"nTrain_Signal=1500:nTrain_Background=1500:nTest_Signal=498:nTest_Background=498:SplitMode=Random:!V" );
   //factory->PrepareTrainingAndTestTree( mycuts, mycutb,
   //"nTrain_Signal=0:nTrain_Background=0:SplitMode=Random:NormMode=NumEvents:!V" );

   // ---- Book MVA methods
   //
   // Please lookup the various method configuration options in the corresponding cxx files, eg:
   // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
   // it is possible to preset ranges in the option string in which the cut optimisation should be done:
   // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

   // Cut optimisation
   if (Use["Cuts"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kCuts, "Cuts",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kCuts, "Cuts",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );
      }

   if (Use["CutsD"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kCuts, "CutsD",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kCuts, "CutsD",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );
      }

   if (Use["CutsPCA"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kCuts, "CutsPCA",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kCuts, "CutsPCA",
                           "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );
      }

   if (Use["CutsGA"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kCuts, "CutsGA",
                           "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kCuts, "CutsGA",
                           "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

      }

   if (Use["CutsSA"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kCuts, "CutsSA",
                           "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kCuts, "CutsSA",
                           "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );
      }

   // Likelihood ("naive Bayes estimator")
   if (Use["Likelihood"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kLikelihood, "Likelihood",
                           "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kLikelihood, "Likelihood",
                           "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );
      }

   // Decorrelated likelihood
   if (Use["LikelihoodD"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kLikelihood, "LikelihoodD",
                           "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kLikelihood, "LikelihoodD",
                           "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" );
      }

   // PCA-transformed likelihood
   if (Use["LikelihoodPCA"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kLikelihood, "LikelihoodPCA",
                           "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kLikelihood, "LikelihoodPCA",
                           "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" );
      }

   // Use a kernel density estimator to approximate the PDFs
   if (Use["LikelihoodKDE"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kLikelihood, "LikelihoodKDE",
                           "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kLikelihood, "LikelihoodKDE",
                           "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" );
      }

   // Use a variable-dependent mix of splines and kernel density estimator
   if (Use["LikelihoodMIX"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kLikelihood, "LikelihoodMIX",
                           "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kLikelihood, "LikelihoodMIX",
                           "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" );
      }

   // Test the multi-dimensional probability density estimator
   // here are the options strings for the MinMax and RMS methods, respectively:
   //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
   //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
   if (Use["PDERS"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kPDERS, "PDERS",
                           "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kPDERS, "PDERS",
                           "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );
      }

   if (Use["PDERSD"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kPDERS, "PDERSD",
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kPDERS, "PDERSD",
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );
      }

   if (Use["PDERSPCA"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kPDERS, "PDERSPCA",
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kPDERS, "PDERSPCA",
                           "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );
      }

   // Multi-dimensional likelihood estimator using self-adapting phase-space binning
   if (Use["PDEFoam"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kPDEFoam, "PDEFoam",
                           "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kPDEFoam, "PDEFoam",
                           "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );
      }

   if (Use["PDEFoamBoost"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kPDEFoam, "PDEFoamBoost",
                           "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kPDEFoam, "PDEFoamBoost",
                           "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );
      }

   // K-Nearest Neighbour classifier (KNN)
   if (Use["KNN"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kKNN, "KNN",
                           "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kKNN, "KNN",
                           "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );
      }

   // H-Matrix (chi2-squared) method
   if (Use["HMatrix"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kHMatrix, "HMatrix", "!H:!V:VarTransform=None" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kHMatrix, "HMatrix", "!H:!V:VarTransform=None" );
   }

   // Linear discriminant (same as Fisher discriminant)
   if (Use["LD"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );
   }

   // Fisher discriminant (same as LD)
   if (Use["Fisher"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );
   }

   // Fisher with Gauss-transformed input variables
   if (Use["FisherG"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );
   }

   // Composite classifier: ensemble (tree) of boosted Fisher classifiers
   if (Use["BoostedFisher"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kFisher, "BoostedFisher",
                           "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kFisher, "BoostedFisher",
                           "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring" );
      }

   // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
   if (Use["FDA_MC"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kFDA, "FDA_MC",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kFDA, "FDA_MC",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );
      }

   if (Use["FDA_GA"]){ // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kFDA, "FDA_GA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kFDA, "FDA_GA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=300:Cycles=3:Steps=20:Trim=True:SaveBestGen=1" );
      }

   if (Use["FDA_SA"]){ // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kFDA, "FDA_SA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kFDA, "FDA_SA",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );
      }

   if (Use["FDA_MT"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kFDA, "FDA_MT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kFDA, "FDA_MT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );
      }

   if (Use["FDA_GAMT"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kFDA, "FDA_GAMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kFDA, "FDA_GAMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );
      }

   if (Use["FDA_MCMT"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kFDA, "FDA_MCMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kFDA, "FDA_MCMT",
                           "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );
      }

   // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
   if (Use["MLP"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );
   }

   if (Use["MLPBFGS"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );
   }

   if (Use["MLPBNN"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators
   }

   // CF(Clermont-Ferrand)ANN
   if (Use["CFMlpANN"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=2000:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...
   }

   // Tmlp(Root)ANN
   if (Use["TMlpANN"]){
      //factory->BookMethod(dataloader, TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=1000:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.01"  ); // n_cycles:#nodes:#nodes:...
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=1000:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.01"  ); // n_cycles:#nodes:#nodes:...
   }

   // Support Vector Machine
   if (Use["SVM"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );
   }

   // Boosted Decision Trees
   if (Use["BDTG"]){ // Gradient Boost
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kBDT, "BDTG",
                           //"!H:!V:NTrees=2000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );
			   "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2:CreateMVAPdfs" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kBDT, "BDTG",
                           //"!H:!V:NTrees=2000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );
			   "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2:CreateMVAPdfs" );
   }

   if (Use["BDT"]){  // Adaptive Boost
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20:CreateMVAPdfs" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kBDT, "BDT",
                           "!H:!V:NTrees=850:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20:CreateMVAPdfs" );
   }

   if (Use["BDTB"]){ // Bagging
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kBDT, "BDTB",
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20:CreateMVAPdfs" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kBDT, "BDTB",
                           "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20:CreateMVAPdfs" );
   }

   if (Use["BDTD"]){ // Decorrelation + Adaptive Boost
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kBDT, "BDTD",
                           "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate:CreateMVAPdfs" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kBDT, "BDTD",
                           "!H:!V:NTrees=400:MinNodeSize=5%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate:CreateMVAPdfs" );
   }

   if (Use["BDTF"]){  // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kBDT, "BDTMitFisher",
                           "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:CreateMVAPdfs" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kBDT, "BDTMitFisher",
                           "!H:!V:NTrees=50:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20:CreateMVAPdfs" );
   }

   // RuleFit -- TMVA implementation of Friedman's method
   if (Use["RuleFit"]){
      factory_cont->BookMethod(dataloader_cont, TMVA::Types::kRuleFit, "RuleFit",
                           "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );
      factory_exiting->BookMethod(dataloader_exiting, TMVA::Types::kRuleFit, "RuleFit",
                           "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );
   }

   // For an example of the category classifier usage, see: TMVAClassificationCategory

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can optimize the setting (configuration) of the MVAs using the set of training events

   // ---- STILL EXPERIMENTAL and only implemented for BDT's !
   // factory->OptimizeAllMethods("SigEffAt001","Scan");
   // factory->OptimizeAllMethods("ROCIntegral","FitGA");

   // --------------------------------------------------------------------------------------------------

   // ---- Now you can tell the factory to train, test, and evaluate the MVAs

   // Train MVAs using the set of training events
   factory_cont->TrainAllMethods();
   factory_exiting->TrainAllMethods();

   // ---- Evaluate all MVAs using the set of test events
   factory_cont->TestAllMethods();
   factory_exiting->TestAllMethods();

   // ----- Evaluate and compare performance of all configured MVAs
   factory_cont->EvaluateAllMethods();
   factory_exiting->EvaluateAllMethods();

   // --------------------------------------------------------------

   // Save the output
   outputFile->Close();

   std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
   std::cout << "==> TMVAClassification is done!" << std::endl;

   delete factory_cont;
   delete factory_exiting;

   // Launch the GUI for the root macros
   if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );
}
