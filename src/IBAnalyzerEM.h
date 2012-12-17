#ifndef IBANALYZEREM_H
#define IBANALYZEREM_H

#include "IBAnalyzer.h"
#include "IBVoxRaytracer.h"

class IBPocaEvaluator;
class IBMinimizationVariablesEvaluator;

class IBAnalyzerEM : public IBAnalyzer {
    typedef IBAnalyzer BaseClass;
public:
    enum Algorithm {
        PXTZ,
        PX,
        TZ,
        PT,
        XZ,
        P,
        T,
        X,
        Z
    };

    ULIB_OBJECT_PARAMETERS(BaseClass)
    {
        enum Algorithm algorithm;
    };

public:
    IBAnalyzerEM();
    ~IBAnalyzerEM();

    void AddMuon(MuonScatterData &muon);

    unsigned int Size();

    void Run(unsigned int iterations, float muons_ratio);

    void SetPocaAlgorithm(IBPocaEvaluator *evaluator);

    void SetVariablesAlgorithm(IBMinimizationVariablesEvaluator *evaluator);

    void SetRaytracer(IBVoxRaytracer *raytracer);

    void SijCut(float threshold);

private:
    friend class IBAnalyzerEMPimpl;
    class IBAnalyzerEMPimpl *d;
};


inline void IBAnalyzerEM::init_parameters() {
    Parameters &p = this->parameters();
    p.algorithm = PXTZ;
}


#endif // IBANALYZEREM_H
