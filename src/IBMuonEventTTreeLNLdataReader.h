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

    unsigned long getNumberOfEvents();
    unsigned long getCurrentPosition();

    bool readNext(uLib::MuonScatter *event);

private:
    friend class IBMuonEventTTreeLNLdataReaderPimpl;
    class IBMuonEventTTreeLNLdataReaderPimpl *d;
};

#endif // IBMUONEVENTTTREELNLDATAREADER_H
