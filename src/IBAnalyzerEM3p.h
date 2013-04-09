#ifndef IBANALYZEREM3P_H
#define IBANALYZEREM3P_H

#include "IBAnalyzerEM.h"

class IBAnalyzerEM3p : public IBAnalyzerEM {
    typedef IBAnalyzerEM BaseClass;

public:
    IBAnalyzerEM3p(IBVoxCollection &voxels) : BaseClass(voxels) {}

    void SetPocaAlgorithm(IBPocaEvaluator *PocaAlgorithm);

    virtual void AddMuon(const MuonScatterData &muon);

};



#endif // IBANALYZEREM3P_H


