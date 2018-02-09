#include "uboone/CC1pi/Algorithms/PDGEnums.h"

class StackedHistPDGCode{

 public:
  StackedHistPDGCode(char* histname, char* title, int nbins, double lowlimit, double highlimit);// Constructor

 protected:
  THStack *stack;
  TH1F *hist_muminus;

}

// Define functions here instead of a .C file so we don't need to many includes

StackedHistPDGCode::StackedHistPDGCode(char* histname, char* title, int nbins, double lowlimit, double highlimit)
{
  stack = new THStack(histname,title);
  std::string histname_pdg = std::string(histname)+std::string("_muminus");
  hist_muminus = new TH1F(histname_pdg,"",nbins,lowlimit,highlimit);
}
