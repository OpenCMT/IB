/*//////////////////////////////////////////////////////////////////////////////
// CMT Cosmic Muon Tomography project //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

  Copyright (c) 2014, Universita' degli Studi di Padova, INFN sez. di Padova

  Coordinators: Prof. Gianni Zumerle < gianni.zumerle@pd.infn.it >
                Paolo Checchia       < paolo.checchia@pd.infn.it >

  Authors: Andrea Rigoni Garola < andrea.rigoni@pd.infn.it >
           Matteo Furlan        < nuright@gmail.com >
           Sara Vanini          < sara.vanini@pd.infn.it >

  All rights reserved
  ------------------------------------------------------------------

  This file can not be copied and/or distributed without the express
  permission of  Prof. Gianni Zumerle  < gianni.zumerle@pd.infn.it >

//////////////////////////////////////////////////////////////////////////////*/



#include "TTree.h"
#include "TFile.h"

#include "IBMuonEventTTreeReader.h"

#include "IBMuonEventTTreeLNLdataReader.h"

using namespace uLib;


IBMuonEventTTreeReader *IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::IBMuonEventTTreeReaderSources S)
{
    switch (S) {
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
