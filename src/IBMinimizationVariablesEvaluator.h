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
#ifndef IBMINIMIZATIONVARIABLESEVALUATOR_H
#define IBMINIMIZATIONVARIABLESEVALUATOR_H

#include "Core/Macros.h"
#include "Math/Dense.h"
#include "Detectors/MuonScatter.h"
#include "IBVoxRaytracer.h"
#include "IBMuonError.h"

using namespace uLib;

class IBNormalPlaneMinimizationVariablesEvaluator;
class IBSimpleTwoViewsMinimizationVariablesEvaluator;

class IBMinimizationVariablesEvaluator : public Object
{
    uLibTypeMacro(IBMinimizationVariablesEvaluator,Object)
public:
    enum IBMinVarEvaluatorAlgorithm {
        NormalPlane,
        SimpleTwoViews
    };

    static IBMinimizationVariablesEvaluator* New(IBMinVarEvaluatorAlgorithm S);

    virtual bool evaluate(MuonScatterData muon)         = 0;
    virtual Vector4f getDataVector()                    = 0;
    virtual Scalarf  getDataVector(int i)               = 0;
    virtual Matrix4f getCovarianceMatrix()              = 0;
    virtual Scalarf  getCovarianceMatrix(int i, int j)  = 0;

    virtual void setRaytracer(IBVoxRaytracer* tracer)   = 0;
    virtual void setDisplacementScatterOnly(bool,bool,bool) = 0;
    
    // virtual void setConfiguration();

    virtual ~IBMinimizationVariablesEvaluator() {} // make this protected again!!

protected:
    IBMinimizationVariablesEvaluator() {}
};

#endif // IBMINIMIZATIONVARIABLESEVALUATOR_H
