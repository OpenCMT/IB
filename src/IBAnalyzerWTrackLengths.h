#ifndef IBANALYZERWTRACKLENGTHS_H
#define IBANALYZERWTRACKLENGTHS_H

#include "IBAnalyzer.h"

#include "IBVoxRaytracer.h"
#include "IBPocaEvaluator.h"

#include "IBMinimizationVariablesEvaluator.h"

using namespace uLib;


class IBAnalyzerWTrackLengths : public IBAnalyzer {

public:
    IBAnalyzerWTrackLengths();
    ~IBAnalyzerWTrackLengths();

    bool AddMuon(const MuonScatterData &muon);

    void Run(unsigned int iterations = 1, float muons_ratio = 1);

    void SetRaytracer(IBVoxRaytracer *raytracer);

    void SetPocaAlgorithm(IBPocaEvaluator *evaluator);

    void SetVaraiblesAlgorithm(IBMinimizationVariablesEvaluator *algorithm);

private:
    class IBAnalyzerWTrackLengthsPimpl *d;
};



#endif // IBANALYZERWTRACKLENGTHS_H
