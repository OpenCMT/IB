/*////////////////////////////////////////////////////////////////////////////
 Copyright 2018 Istituto Nazionale di Fisica Nucleare

 Licensed under the EUPL, Version 1.2 or - as soon they will be approved by
 the European Commission - subsequent versions of the EUPL (the "Licence").
 You may not use this work except in compliance with the Licence.

 You may obtain a copy of the Licence at:

 https://joinup.ec.europa.eu/software/page/eupl

 Unless required by applicable law or agreed to in writing, software
 distributed under the Licence is distributed on an "AS IS" basis, WITHOUT
 WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
 Licence for the specific language governing permissions and limitations under
 the Licence.
////////////////////////////////////////////////////////////////////////////*/



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
