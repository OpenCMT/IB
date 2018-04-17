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
