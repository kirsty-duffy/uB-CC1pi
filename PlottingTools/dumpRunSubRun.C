// Script to produce a text file containing run/subrun numbers: run_subrun.txt (produced in current directory)
// This can then be used with Zarko's script to get the total POT: /uboone/app/users/zarko/getDataInfo.py -v2 --run-subrun-list run_subrun.txt
// NOTE: this only works for data!
// For MC, use getSimPot.C.

void dumpRunSubRun(std::string inputfile){

  UInt_t run;
  UInt_t sub_run;

  TFile *_file0 = new TFile(inputfile.c_str(),"read");
  TTree* tree = (TTree*)_file0->Get("cc1piselec/pottree");
  // tree->Print();

  tree->SetBranchStatus("*",0);
  tree->SetBranchStatus("run_num",1);
  tree->SetBranchAddress("run_num", &run);
  tree->SetBranchStatus("subrun_num",1);
  tree->SetBranchAddress("subrun_num", &sub_run);

  std::vector<std::pair<int, int>> rsr_pairs;

  for (int i = 0; i < tree->GetEntries(); i++){
    if (i%1000==0) std::cout << i << "/" << tree->GetEntries() << std::endl;

    tree->GetEntry(i);

    std::pair<int, int> rsr_pair;
    rsr_pair.first = run;
    rsr_pair.second = sub_run;
    rsr_pairs.push_back(rsr_pair);

  }

  // open output txt file
  std::ofstream ofs("run_subrun.txt", std::ofstream::out);

  std::vector<std::pair<int, int>> outvec;
  for (int i = 0; i < rsr_pairs.size(); i++){
    if (i%1000==0) std::cout << "-- " << i << "/" << rsr_pairs.size() << std::endl;

    // get run subrun info
    int runt1 = rsr_pairs.at(i).first;
    int subrunt1 = rsr_pairs.at(i).second;

    bool isWritten = false;
    for (int j = 0; j < outvec.size(); j++){

      int runt2 = outvec.at(j).first;
      int subrunt2 = outvec.at(j).second;

      if (runt1 == runt2 && subrunt1 == subrunt2){
        isWritten = true;
      }

    }
    if (isWritten == false) outvec.push_back(rsr_pairs.at(i));

  }

  for (int i = 0; i < outvec.size(); i++){
    if (i%1000==0) std::cout << i << "/" << outvec.size() << std::endl;
    ofs << outvec.at(i).first << " " << outvec.at(i).second << std::endl;

  }

  ofs.close();

}
