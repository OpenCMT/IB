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



#include "IBPocaEvaluator.h"
#include "IBTiltedAxisPocaEvaluator.h"
#include "IBLineDistancePocaEvaluator.h"

using namespace uLib;

IBPocaEvaluator *IBPocaEvaluator::New(IBPocaEvaluator::IBPocaEvaluationAlgorithms S)
{
    switch (S) {
    case TiltedAxis:
        return new IBTiltedAxisPocaEvaluator;
    case LineDistance:
        return new IBLineDistancePocaEvaluator;
    }
}
