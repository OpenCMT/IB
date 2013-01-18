#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include "IBMuonEventTTreeReader.h"
#include "IBMuonEventTTreeR3DmcReader.h"

using namespace uLib;

int main() {

    // errors //
    IBMuonError sigma(11.93,2.03, 18.53,2.05);
    // reader //
    TFile* f = new TFile
            ("/var/local/data/root/run_20130124/muSteel_PDfit_20130124_v12.root");
    TTree* t = (TTree*)f->Get("n");
    IBMuonEventTTreeReader* reader = IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::R3D_MC);
    reader->setAcquisitionTime(40);
    reader->setTTree(t);
    reader->setError(sigma);
    reader->setMomentum(0.7);
    reader->selectionCode(IBMuonEventTTreeR3DmcReader::All);
    std::cout << "Number of Events: " << reader->getNumberOfEvents() << std::endl;

    exit(0);

}
