#ifndef IBANALYZERPOCA_H
#define IBANALYZERPOCA_H

#include "IBAnalyzer.h"

using namespace uLib;

class IBPocaEvaluator;

class IBAnalyzerPoca : public IBAnalyzer {

public:
    IBAnalyzerPoca();
    ~IBAnalyzerPoca();

    void AddMuon(const MuonScatterData &event);

    void Run(unsigned int iterations = 1, float muons_ratio = 1);

    void SetPocaAlgorithm(IBPocaEvaluator *poca);

private:
    class IBAnalyzerPocaPimpl *d;
};




#endif // IBANALYZERPOCA_H
