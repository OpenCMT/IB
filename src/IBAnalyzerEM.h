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
        PXTZ = 0,
        PX,
        TZ,
        PT,
        XZ,
        P,
        T,
        X,
        Z
    };

    enum PWeigthAlgorithm {
        PWeigth_disabled = 0,
        PWeigth_pw,
        PWeigth_sw,
        PWeigth_cw
    };

    ULIB_OBJECT_PARAMETERS(BaseClass)
    {
        Scalarf nominal_momentum;
        enum Algorithm algorithm;
        enum PWeigthAlgorithm pweigth;
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

    void UpdatePW(enum PWeigthAlgorithm algorithm);

private:
    friend class IBAnalyzerEMPimpl;
    class IBAnalyzerEMPimpl *d;
};


inline void IBAnalyzerEM::init_parameters() {
    Parameters &p = this->parameters();
    p.nominal_momentum = 3;
    p.algorithm = PXTZ;
    p.pweigth = PWeigth_pw;
}


#endif // IBANALYZEREM_H
