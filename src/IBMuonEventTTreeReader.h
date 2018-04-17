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



#ifndef IBMUONEVENTTTREEREADER_H
#define IBMUONEVENTTTREEREADER_H

#include "Detectors/MuonScatter.h"
#include "IBMuonError.h"

using namespace uLib;

class IBMuonEventTTreeReader {

public:
    enum IBMuonEventTTreeReaderSources {
        R3D_MC,
        R2D_MC,
        LNL_MC,
        LNL_DATA
    };

    static IBMuonEventTTreeReader* New(enum IBMuonEventTTreeReaderSources S);
    static IBMuonEventTTreeReader* New (class TFile* f);
    virtual void setTTree(class TTree* tree)        = 0;
    virtual void setTFile(class TFile* file)        = 0;
    virtual void setHitCuts(int nx_cut, int nz_cut) = 0;
    virtual void setMomentum(Scalarf p)             = 0;
    virtual void setError(IBMuonError &e)           = 0;
    virtual void selectionCode(short code)          = 0;
    virtual void setAcquisitionTime(float min)       {}
    virtual void setStartTime(float min)             {}
    virtual void setReadPCut(float pcut)             {}
    virtual void readPguess(bool yn=true)            {}
    virtual unsigned long getNumberOfEvents()  = 0;
    virtual unsigned long getCurrentPosition() = 0;

    virtual bool readNext(uLib::MuonScatter* event) = 0;

    virtual ~IBMuonEventTTreeReader() {}
protected:
    IBMuonEventTTreeReader() {}

};

#endif // IBMUONEVENTTTREEREADER_H
