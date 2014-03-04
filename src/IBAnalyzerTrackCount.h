#ifndef IBANALYZERTRACKCOUNT_H
#define IBANALYZERTRACKCOUNT_H


#include "IBAnalyzer.h"
#include "IBVoxRaytracer.h"
#include "IBPocaEvaluator.h"

using namespace uLib;


class IBAnalyzerTrackCount : public IBAnalyzer {    
    uLibTypeMacro(IBAnalyzerTrackCount,IBAnalyzer)
public:

    IBAnalyzerTrackCount();
    ~IBAnalyzerTrackCount();

    bool AddMuon(const MuonScatterData &muon);

    void Run(unsigned int iterations = 1, float muons_ratio = 1);

    void SetRayAlgorithm(IBVoxRaytracer *raytracer);

    void SetPocaAlgorithm(IBPocaEvaluator *evaluator);

    void Clear();

    unsigned int Size() const;

    void SetMuonCollection(IBMuonCollection *muons);
private:
    class IBAnalyzerTrackCountPimpl *d;
};



#endif // IBANALYZERTRACKCOUNT_H
