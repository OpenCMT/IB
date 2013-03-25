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
