#ifndef IBNORMALPLANEMINIMIZATIONVARIABLESEVALUATOR_H
#define IBNORMALPLANEMINIMIZATIONVARIABLESEVALUATOR_H

#include "IBMinimizationVariablesEvaluator.h"

using namespace uLib;

class IBNormalPlaneMinimizationVariablesEvaluator : public IBMinimizationVariablesEvaluator
{
public:
    IBNormalPlaneMinimizationVariablesEvaluator();
    ~IBNormalPlaneMinimizationVariablesEvaluator();

    bool evaluate(MuonScatterData muon);

    Vector4f getDataVector();
    Scalarf  getDataVector(int i);
    Matrix4f getCovarianceMatrix();
    Scalarf  getCovarianceMatrix(int i, int j);
    void setRaytracer(IBVoxRaytracer *tracer);

    // virtual void setConfiguration();
private:
    friend class IBNormalPlaneMinimizationVariablesEvaluatorPimpl;
    class IBNormalPlaneMinimizationVariablesEvaluatorPimpl *d;

};

#endif // IBNORMALPLANEMINIMIZATIONVARIABLESEVALUATOR_H
