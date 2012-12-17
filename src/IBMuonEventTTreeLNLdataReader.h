#ifndef IBMUONEVENTTTREELNLDATAREADER_H
#define IBMUONEVENTTTREELNLDATAREADER_H

#include "IBMuonEventTTreeReader.h"

class IBMuonEventTTreeLNLdataReader : public IBMuonEventTTreeReader
{
public:
    IBMuonEventTTreeLNLdataReader();
    ~IBMuonEventTTreeLNLdataReader();
    void setTTree(class TTree* tree);
    void setHitCuts(int nx_cut, int nz_cut);
    void setMomentum(Scalarf p);
    void selectionCode(short code);
    void setError(IBMuonError e);

    unsigned long getNumberOfEvents();
    unsigned long getCurrentPosition();

    bool readNext(uLib::MuonScatter *event);

private:
    friend class IBMuonEventTTreeLNLdataReaderPimpl;
    class IBMuonEventTTreeLNLdataReaderPimpl *d;
};

#endif // IBMUONEVENTTTREELNLDATAREADER_H
