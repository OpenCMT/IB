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



#ifndef IBPOCAEVALUATOR_H
#define IBPOCAEVALUATOR_H

#include "Core/Macros.h"
#include "Math/Dense.h"
#include "Math/Utils.h"
#include "Detectors/MuonScatter.h"

using namespace uLib;

class IBTiltedAxisPocaEvaluator;
class IBLineDistancePocaEvaluator;

class IBPocaEvaluator
{
public:
    virtual ~IBPocaEvaluator() {}

    enum IBPocaEvaluationAlgorithms {
        TiltedAxis,
        LineDistance
    };

    static IBPocaEvaluator* New(enum IBPocaEvaluationAlgorithms S);

    virtual bool evaluate(MuonScatterData muon) = 0;
    virtual HPoint3f getPoca() = 0;

protected:
    IBPocaEvaluator() {}
};

#endif // IBPOCAEVALUATOR_H
