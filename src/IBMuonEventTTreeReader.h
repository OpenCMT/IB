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

    virtual void setTTree(class TTree* tree)        = 0;
    virtual void setHitCuts(int nx_cut, int nz_cut) = 0;
    virtual void setMomentum(Scalarf p)             = 0;
    virtual void setError(IBMuonError &e)           = 0;
    virtual void selectionCode(short code)          = 0;
    virtual void setAcquisitionTime(float min)       {}
    virtual void setStartTime(float min)             {}
    virtual void setReadPCut(float pcut)             {}
    virtual unsigned long getNumberOfEvents()  = 0;
    virtual unsigned long getCurrentPosition() = 0;

    virtual bool readNext(uLib::MuonScatter* event) = 0;

protected:
    IBMuonEventTTreeReader() {}
    virtual ~IBMuonEventTTreeReader() {}

};

#endif // IBMUONEVENTTTREEREADER_H
