#include <iostream>
#include <assert.h>

#include <root/TFile.h>
#include <root/TTree.h>


#include <Detectors/MuonScatter.h>
#include "IBMuonTTreeReader.h"


#include "testing-prototype.h"


using namespace uLib;

int main()
{

  BEGIN_TESTING(IBMuon);

  TFile* hbfile = new TFile("/var/local/data/root/run_1363/muRadio_1363.root");
  if (hbfile->IsZombie()) {
      exit(-1);
  }
  TTree* tree = (TTree*) hbfile->Get("RADMU");
  assert(tree);

  IBMuonReader reader(tree,Demo_LNL);

  std::cout << "read->num of events = " << reader.GetNumberOfEvents() << "\n";

  MuonScatter *muon = reader.GetNext();





  END_TESTING;
}




