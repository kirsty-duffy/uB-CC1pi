// Macro to make the BDT output plots I'm interested in without having to use the GUI
// Adapted largely from TMVAGui.C to call the root macros that are already supplied

#include <iostream>
#include <vector>

#include "TList.h"
#include "TROOT.h"
#include "TKey.h"
#include "TString.h"

#include "/uboone/app/users/kduffy/CC1pi/CC1pi_uboonecode/srcs/uboonecode/uboone/CC1pi/MVA/tmvaglob.C"

// some global lists
static TList*               TMVAGui_keyContent;
static std::vector<TString> TMVAGui_inactiveButtons;

TList* GetKeyList( const TString& pattern )
{
   TList* list = new TList();

   TIter next( TMVAGui_keyContent );
   TKey* key(0);
   while ((key = (TKey*)next())) {
      if (TString(key->GetName()).Contains( pattern )) {
       list->Add( new TObjString( key->GetName() ) ); }
   }
   return list;
}

// main function
void TMVAPlots_noGui(const char* fName = "TMVA.root")
{

   TString curMacroPath(gROOT->GetMacroPath());
   // uncomment next line for macros submitted to next root version
   gROOT->SetMacroPath(curMacroPath+":./:$ROOTSYS/tmva/test/:");
   // This line is specific to Kirsty's setup
   gROOT->SetMacroPath(curMacroPath+":/uboone/app/users/kduffy/CC1pi/CC1pi_uboonecode/srcs/uboonecode/uboone/CC1pi/MVA/");

   // check if file exist
   TFile* file = TFile::Open( fName );
   if (!file) {
      cout << "==> Abort TMVAGui, please verify filename" << endl;
      return;
   }

   // find all references
   TMVAGui_keyContent = (TList*)file->GetListOfKeys()->Clone();

   // If there is only one thing at the top level and it is a TDirectoryFile, go inside that TDirectory and look for keys in there
   // I have no idea how TMVAGui doesn't have to do this
   TString firstleveldir = "";
   if (TMVAGui_keyContent->GetEntries()<10){
      TDirectory *dir;
      TIter next(TMVAGui_keyContent);
      TKey *key(0);
      while ((key = (TKey*)next())){
         if (!strcmp(key->GetClassName(),"TDirectoryFile")){
            firstleveldir = (TString)(key->GetName());
            gDirectory->cd(key->GetName());
            TMVAGui_keyContent = (TList*)gDirectory->GetListOfKeys()->Clone();

            // Plots for all input variables
            // Call variables.C
            TList* keylist = GetKeyList( "InputVariables" );
            TListIter it( keylist );
            TObjString* str = 0;
            char ch = 'a';
            while ((str = (TObjString*)it())) {
               TString tmp   = str->GetString();
               TString title = Form( "Input variables '%s'-transformed (training sample)", tmp.ReplaceAll("InputVariables_","").Data() );
               if (tmp.Contains( "Id" )) title = "Input variables (training sample)";
               TString macro = Form(".x variables.C(\"%s\",\"%s\",\"%s\",%s,%s,\"%s\")",fName, str->GetString().Data(), title.Data(), "kFALSE", "kTRUE", firstleveldir.Data());
               gROOT->ProcessLine(macro);
            }

            // Correlations plots
            // Call correlations.C
            TString macro = Form(".x correlations.C(\"%s\",%s,%s,%s,\"%s\")",fName,"kFALSE","kFALSE", "kTRUE", firstleveldir.Data() );
            gROOT->ProcessLine(macro);

            // MVA-specific background and signal overlay plots
            // Call mvas.C
            macro = Form(".x mvas.C(\"%s\",MVAType,%s,\"%s\")", fName, "kTRUE", firstleveldir.Data() );
            gROOT->ProcessLine(macro);
            macro = Form(".x mvas.C(\"%s\",ProbaType,%s,\"%s\")", fName, "kTRUE", firstleveldir.Data() );
            gROOT->ProcessLine(macro);
            macro = Form(".x mvas.C(\"%s\",RarityType,%s,\"%s\")", fName, "kTRUE", firstleveldir.Data() );
            gROOT->ProcessLine(macro);
            macro = Form(".x mvas.C(\"%s\",CompareType,%s,\"%s\")", fName, "kTRUE", firstleveldir.Data() );
            gROOT->ProcessLine(macro);

            // Signal efficiency and purity as a function of cut value on MVA response
            // Call mvaeffs.C
            // macro = Form(".x mvaeffs.C++(\"%s\",%s,\"%s\",\"%s\")", fName, "kTRUE", "S/sqrt(S+B)", firstleveldir.Data());
            // gROOT->ProcessLine(macro);

            // Background rejection vs signal efficiency (ROC curve)
            // Call efficiencies.C
            macro = Form(".x efficiencies.C(\"%s\",1,kTRUE,\"%s\")", fName, firstleveldir.Data());
            gROOT->ProcessLine(macro);
            macro = Form(".x efficiencies.C(\"%s\",2,kTRUE,\"%s\")", fName, firstleveldir.Data());
            gROOT->ProcessLine(macro);
            macro = Form(".x efficiencies.C(\"%s\",3,kTRUE,\"%s\")", fName, firstleveldir.Data());
            gROOT->ProcessLine(macro);

            // PDFs of classifiers (will only produce plots if "CreateMVAPdfs" option is set)
            // Call probas.C
            macro = Form(".x probas.C(\"%s\",kTRUE,\"%s\")", fName, firstleveldir.Data());
            gROOT->ProcessLine(macro);

            // BDT control plots
            // Call BDTControlPlots.C
            macro = Form(".x BDTControlPlots.C(\"%s\",kTRUE,\"%s\")", fName, firstleveldir.Data());
            gROOT->ProcessLine(macro);
         }
      } // end loop over top-level directories
   } // end if nkeys<10 (bad standin for if you're looking at the directories or if you have the output already)


   // close file
   file->Close();
}
