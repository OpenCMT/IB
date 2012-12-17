#ifndef IBANALYZERTRACKCOUNT_H
#define IBANALYZERTRACKCOUNT_H


#include "IBAnalyzer.h"
#include "IBVoxRaytracer.h"
#include "IBPocaEvaluator.h"

using namespace uLib;


class IBAnalyzerTrackCount : public IBAnalyzer {

public:
    IBAnalyzerTrackCount();
    ~IBAnalyzerTrackCount();

    void AddMuon(MuonScatterData &muon);

    void Run(unsigned int iterations = 1, float muons_ratio = 1);

    void SetRaytracer(IBVoxRaytracer *raytracer);

    void SetPocaAlgorithm(IBPocaEvaluator *evaluator);

private:
    class IBAnalyzerTrackCountPimpl *d;
};


#endif // IBANALYZERTRACKCOUNT_H
