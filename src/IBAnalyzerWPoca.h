#ifndef IBANALYZERWPOCA_H
#define IBANALYZERWPOCA_H

#include "IBAnalyzer.h"

using namespace uLib;

class IBPocaEvaluator;
class IBMinimizationVariablesEvaluator;

class IBAnalyzerWPoca : public IBAnalyzer
{
public:
    IBAnalyzerWPoca();
    ~IBAnalyzerWPoca();

    void AddMuon(const MuonScatterData &event);

    void Run(unsigned int iteration = 1, float muons_ratio = 1);

    void SetPocaAlgorithm(IBPocaEvaluator *poca);

    void SetVariablesAlgorithm(IBMinimizationVariablesEvaluator *evaluator);

private:
    friend class IBAnalyzerWPocaPimpl;
    class IBAnalyzerWPocaPimpl *d;
};

#endif // IBANALYZERWPOCA_H
