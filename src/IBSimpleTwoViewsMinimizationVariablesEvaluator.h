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



#ifndef IBSIMPLETWOVIEWSMINIMIZATIONVARIABLESEVALUATOR_H
#define IBSIMPLETWOVIEWMINIMIZATIONVARIABLESEVALUATOR_H

#include "IBMinimizationVariablesEvaluator.h"

using namespace uLib;

class IBSimpleTwoViewsMinimizationVariablesEvaluator : public IBMinimizationVariablesEvaluator
{
public:
    IBSimpleTwoViewsMinimizationVariablesEvaluator();
    ~IBSimpleTwoViewsMinimizationVariablesEvaluator();

    bool evaluate(MuonScatterData muon);

    Vector4f getDataVector();
    Scalarf  getDataVector(int i);
    Matrix4f getCovarianceMatrix();
    Scalarf  getCovarianceMatrix(int i, int j);
    void setRaytracer(IBVoxRaytracer *tracer);

    // virtual void setConfiguration();
private:
    friend class IBSimpleTwoViewsMinimizationVariablesEvaluatorPimpl;
    class IBSimpleTwoViewsMinimizationVariablesEvaluatorPimpl *d;

};

#endif // IBSimpleTwoViewsMINIMIZATIONVARIABLESEVALUATOR_H
