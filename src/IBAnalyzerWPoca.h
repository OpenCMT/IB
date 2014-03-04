#ifndef IBANALYZERWPOCA_H
#define IBANALYZERWPOCA_H

#include "IBAnalyzer.h"

using namespace uLib;

class IBPocaEvaluator;
class IBMinimizationVariablesEvaluator;

class IBAnalyzerWPoca : public IBAnalyzer
{
    uLibTypeMacro(IBAnalyzerWPoca,IBAnalyzer)
public:
    IBAnalyzerWPoca();
    ~IBAnalyzerWPoca();

    bool AddMuon(const MuonScatterData &event);

    void Run(unsigned int iteration = 1, float muons_ratio = 1);

    void SetPocaAlgorithm(IBPocaEvaluator *poca);

    void SetVarAlgorithm(IBMinimizationVariablesEvaluator *evaluator);

private:
    friend class IBAnalyzerWPocaPimpl;
    class IBAnalyzerWPocaPimpl *d;
};

#endif // IBANALYZERWPOCA_H
