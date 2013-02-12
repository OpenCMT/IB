#ifndef IBANALYZERTRACKLENGTHS_H
#define IBANALYZERTRACKLENGTHS_H


#include "IBAnalyzer.h"

#include "IBVoxRaytracer.h"
#include "IBPocaEvaluator.h"

using namespace uLib;


class IBAnalyzerTrackLengths : public IBAnalyzer {

public:
    IBAnalyzerTrackLengths();
    ~IBAnalyzerTrackLengths();

    void AddMuon(const MuonScatterData &muon);

    void Run(unsigned int iterations = 1, float muons_ratio = 1);

    void SetRaytracer(IBVoxRaytracer *raytracer);

    void SetPocaAlgorithm(IBPocaEvaluator *evaluator);

private:
    class IBAnalyzerTrackLengthsPimpl *d;
};



#endif // IBANALYZERTRACKLENGTHS_H
