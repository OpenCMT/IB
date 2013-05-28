#include "TTree.h"
#include "TFile.h"

#include "IBMuonEventTTreeReader.h"

#include "IBMuonEventTTreeR3DmcReader.h"
#include "IBMuonEventTTreeR2DmcReader.h"
#include "IBMuonEventTTreeLNLmcReader.h"
#include "IBMuonEventTTreeLNLdataReader.h"

using namespace uLib;


IBMuonEventTTreeReader *IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::IBMuonEventTTreeReaderSources S)
{
    switch (S) {
    case R3D_MC:
        return new IBMuonEventTTreeR3DmcReader;
    case R2D_MC:
        return new IBMuonEventTTreeR2DmcReader;
    case LNL_MC:
        return new IBMuonEventTTreeLNLmcReader;
    case LNL_DATA:
        return new IBMuonEventTTreeLNLdataReader;
    }
}

IBMuonEventTTreeReader *IBMuonEventTTreeReader::New(TFile* f)
{
    IBMuonEventTTreeReader* reader;
    if (f->Get("RADMU")) {
        reader = IBMuonEventTTreeReader::New(LNL_DATA);
        reader->setTFile(f);
        return reader;
    }
    if (f->Get("n")) {
        TTree* t = (TTree*)f->Get("n");
        if (t->FindBranch("flagMC")) {
            reader = IBMuonEventTTreeReader::New(R3D_MC);
            reader->setTTree(t);
            return reader;
        }
        if (t->GetLeaf("nptx1")) {
            reader = IBMuonEventTTreeReader::New(LNL_MC);
            reader->setTTree(t);
            return reader;
        }
    }
    printf("File or TTree not found! Aborting...\n");
    exit(1);
}
