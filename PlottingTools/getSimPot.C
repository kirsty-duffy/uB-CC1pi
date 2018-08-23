// Script to get POT from an MC file
// NOTE: only works for MC! For data, see dumpRunSubRun.C

void getSimPot(std::string inputfile){

  double sr_pot;

  TFile *_file0 = new TFile(inputfile.c_str(),"read");
  TTree* tree = (TTree*)_file0->Get("cc1piselec/pottree");

  tree->SetBranchAddress("pot", &sr_pot);

  double totPot = 0;
  for (int i = 0; i < tree->GetEntries(); i++){

    tree->GetEntry(i);

    totPot += sr_pot;
  }

  std::cout << "Total MC POT: " <<  totPot << std::endl;

}
