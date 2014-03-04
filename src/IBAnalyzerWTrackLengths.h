#ifndef IBANALYZERWTRACKLENGTHS_H
#define IBANALYZERWTRACKLENGTHS_H

#include "IBAnalyzer.h"

#include "IBVoxRaytracer.h"
#include "IBPocaEvaluator.h"

#include "IBMinimizationVariablesEvaluator.h"

using namespace uLib;


class IBAnalyzerWTrackLengths : public IBAnalyzer {
    uLibTypeMacro(IBAnalyzerWTrackLengths,IBAnalyzer)
public:
    IBAnalyzerWTrackLengths();
    ~IBAnalyzerWTrackLengths();

    bool AddMuon(const MuonScatterData &muon);

    void Run(unsigned int iterations = 1, float muons_ratio = 1);

    void SetRayAlgorithm(IBVoxRaytracer *raytracer);

    void SetPocaAlgorithm(IBPocaEvaluator *evaluator);

    void SetVarAlgorithm(IBMinimizationVariablesEvaluator *algorithm);

    void SetPocaProximity(float sigma = 0); // TODO: //

private:
    class IBAnalyzerWTrackLengthsPimpl *d;
};



#endif // IBANALYZERWTRACKLENGTHS_H
