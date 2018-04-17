/*////////////////////////////////////////////////////////////////////////////
 Copyright 2018 Istituto Nazionale di Fisica Nucleare

 Licensed under the EUPL, Version 1.2 or - as soon they will be approved by
 the European Commission - subsequent versions of the EUPL (the "Licence").
 You may not use this work except in compliance with the Licence.

 You may obtain a copy of the Licence at:

 https://joinup.ec.europa.eu/software/page/eupl

 Unless required by applicable law or agreed to in writing, software
 distributed under the Licence is distributed on an "AS IS" basis, WITHOUT
 WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
 Licence for the specific language governing permissions and limitations under
 the Licence.
////////////////////////////////////////////////////////////////////////////*/


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
