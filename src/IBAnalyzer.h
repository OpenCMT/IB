#ifndef IBANALYZER_H
#define IBANALYZER_H

#include <Core/Object.h>

#include <Detectors/MuonScatter.h>
#include <Detectors/MuonError.h>

#include "IBMuonCollection.h"
#include "IBExperiment.h"

using namespace uLib;

class IBAnalyzerPoca;
class IBAnalyzerWPoca;
class IBAnalyzerEM;

class IBPocaEvaluator;
class IBMinimizationVariablesEvaluator;

class IBVoxCollection;


class IBAnalyzer : public Object {

public:

    uLibGetSetMacro(Experiment,IBExperiment *)

    uLibGetSetMacro(VoxCollection,IBVoxCollection *)

    uLibGetMacro(MuonCollection,IBMuonCollection *)

    virtual uLibSetMacro(MuonCollection,IBMuonCollection *)

    virtual bool AddMuon(const MuonScatterData &event) = 0;

    virtual unsigned int Size() {}

    virtual void Run(unsigned int iterations, float muons_ratio) = 0;

protected:
    IBAnalyzer() :
        m_Experiment(NULL),
        m_VoxCollection(NULL),
        m_MuonCollection(NULL)
    {}

    virtual ~IBAnalyzer() {}

private:    
    IBExperiment       *m_Experiment;
    IBVoxCollection    *m_VoxCollection;
    IBMuonCollection   *m_MuonCollection;
};







#endif // IBANALYZER_H



















