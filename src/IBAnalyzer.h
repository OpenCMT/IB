/*//////////////////////////////////////////////////////////////////////////////
// CMT Cosmic Muon Tomography project //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

  Copyright (c) 2014, Universita' degli Studi di Padova, INFN sez. di Padova

  Coordinators: Prof. Gianni Zumerle < gianni.zumerle@pd.infn.it >
                Paolo Checchia       < paolo.checchia@pd.infn.it >

  Authors: Andrea Rigoni Garola < andrea.rigoni@pd.infn.it >
           Matteo Furlan        < nuright@gmail.com >
           Sara Vanini          < sara.vanini@pd.infn.it >

  All rights reserved
  ------------------------------------------------------------------

  This file can not be copied and/or distributed without the express
  permission of  Prof. Gianni Zumerle  < gianni.zumerle@pd.infn.it >

//////////////////////////////////////////////////////////////////////////////*/



#ifndef IBANALYZER_H
#define IBANALYZER_H

#include <Core/Object.h>

#include <Detectors/MuonScatter.h>
#include <Detectors/MuonError.h>

#include "IBMuonCollection.h"

using namespace uLib;

class IBExperiment;
class IBAnalyzerPoca;
class IBAnalyzerWPoca;
class IBAnalyzerEM;

class IBPocaEvaluator;
class IBMinimizationVariablesEvaluator;

class IBVoxCollection;


class IBAnalyzer : public Object {
    uLibTypeMacro(IBAnalyzer,Object)

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



















