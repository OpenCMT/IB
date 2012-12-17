#ifndef IBANALYZER_H
#define IBANALYZER_H

#include <Core/Object.h>

#include <Detectors/MuonScatter.h>
#include <Detectors/MuonError.h>

#include "IBExperiment.h"

using namespace uLib;

class IBAnalyzerPoca;
class IBAnalyzerWPoca;
class IBAnalyzerEM;

class IBPocaEvaluator;
class IBMinimizationVariablesEvaluator;

class IBVoxCollectionCap;


class IBAnalyzer : public Object {

public:
    enum IBAnalyzerFactoryType {
        Poca  = 0,
        WPoca = 1,
        EM    = 2
    };

public:

    static IBAnalyzer * New(enum IBAnalyzerFactoryType id);

    uLibGetSetMacro(Experiment,IBExperiment *)
    uLibGetSetMacro(VoxCollection,IBVoxCollectionCap *)

    virtual void AddMuon(MuonScatterData &event) = 0;

    virtual unsigned int Size() {}

    virtual void Run(unsigned int iterations, float muons_ratio) = 0;

    virtual void SetPocaAlgorithm(IBPocaEvaluator *poca) {}
    virtual void SetVariablesAlgorithm(IBMinimizationVariablesEvaluator *evaluator) {}

protected:
    IBAnalyzer() {}
    virtual ~IBAnalyzer() {}

private:    
    IBExperiment       *m_Experiment;
    IBVoxCollectionCap *m_VoxCollection;
};







#endif // IBANALYZER_H



















