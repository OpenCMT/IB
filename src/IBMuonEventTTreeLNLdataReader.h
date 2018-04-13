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


#ifndef IBMUONEVENTTTREELNLDATAREADER_H
#define IBMUONEVENTTTREELNLDATAREADER_H

#include "IBMuonEventTTreeReader.h"

class IBMuonEventTTreeLNLdataReader : public IBMuonEventTTreeReader
{
public:
    IBMuonEventTTreeLNLdataReader();
    ~IBMuonEventTTreeLNLdataReader();
    void setTTree(class TTree* tree);
    void setTFile(class TFile* file);
    void setHitCuts(int nx_cut, int nz_cut);
    void setMomentum(Scalarf p);
    void selectionCode(short code);
    void setError(IBMuonError &e);
    void setAcquisitionTime(float min);
    void setStartTime(float min);

    void setAlignmentFromData(float min = 0.0);
    void setAlignment(Matrix4f align);
    Matrix4f getAlignment();


    unsigned long getNumberOfEvents();
    unsigned long getCurrentPosition();

    bool readNext(uLib::MuonScatter *event);
    
private:
    friend class IBMuonEventTTreeLNLdataReaderPimpl;
    class IBMuonEventTTreeLNLdataReaderPimpl *d;
};

#endif // IBMUONEVENTTTREELNLDATAREADER_H
