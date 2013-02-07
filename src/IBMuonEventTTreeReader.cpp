
#include "TTree.h"

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
