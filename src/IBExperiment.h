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



#ifndef IBEXPERIMENT_H
#define IBEXPERIMENT_H

//#include <iostream>
//#include <fstream>

//#include <TFile.h>
//#include <TTree.h>

//#include "Core/Object.h"
//#include "Core/Options.h"
//#include "Detectors/MuonScatter.h"

//#include "IB.h"

//#include "IBMuonError.h"
#include "IBVoxCollection.h"
#include "IBMuonCollection.h"

#include "IBVoxRaytracer.h"
#include "IBPocaEvaluator.h"
#include "IBMinimizationVariablesEvaluator.h"
#include "IBVoxFilters.h"

#include "IBMuonEventTTreeReader.h"

#include "IBAnalyzer.h"


////////////////////////////////////////////////////////////////////////////////
// EXPERIMENT PARAMETERS //

namespace IB {

class Experiment : public uLib::Object {
public:
    Experiment() :
        m_Reader(NULL),
        m_Poca(NULL),
        m_Tracer(NULL),
        m_Analyzer(NULL),
        m_Voxels(NULL)
    {}

    ~Experiment() {}

    uLibRefMacro(m_Muons,IBMuonCollection)
    uLibConstRefMacro(m_Muons,IBMuonCollection)

    uLibGetSetMacro(m_Reader,IBMuonEventTTreeReader*)
    uLibGetSetMacro(m_Poca,IBPocaEvaluator*)
    uLibGetSetMacro(m_Tracer,IBVoxRaytracer*)
    uLibGetSetMacro(m_Analyzer,IBAnalyzer*)
    uLibGetSetMacro(m_Voxels,IBVoxCollection*)

protected:
    // members //
    IBMuonCollection       m_Muons;
    IBMuonEventTTreeReader *m_Reader;
    IBPocaEvaluator        *m_Poca;
    IBVoxRaytracer         *m_Tracer;
    IBAnalyzer             *m_Analyzer;
    IBVoxCollection        *m_Voxels;
};

} // IB




#endif // IBEXPERIMENT_H
