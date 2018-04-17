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



#ifndef IBNORMALPLANEMINIMIZATIONVARIABLESEVALUATOR_H
#define IBNORMALPLANEMINIMIZATIONVARIABLESEVALUATOR_H

#include "Core/Object.h"
#include "IBMinimizationVariablesEvaluator.h"

using namespace uLib;

class IBNormalPlaneMinimizationVariablesEvaluator : public IBMinimizationVariablesEvaluator
{
    uLibTypeMacro(IBNormalPlaneMinimizationVariablesEvaluator,IBMinimizationVariablesEvaluator)
public:
    properties() {
        bool    use_free_rotation;
        Scalarf alphaXZ;
    };

    bool m_scatterOnly, m_displacementOnly, m_oneD;

    
    IBNormalPlaneMinimizationVariablesEvaluator();
    ~IBNormalPlaneMinimizationVariablesEvaluator();

    bool evaluate(MuonScatterData muon);

    Vector4f getDataVector();
    Scalarf  getDataVector(int i);
    Matrix4f getCovarianceMatrix();
    Scalarf  getCovarianceMatrix(int i, int j);
    void setRaytracer(IBVoxRaytracer *tracer);
    void setDisplacementScatterOnly(bool,bool,bool);
    // virtual void setConfiguration();
private:
    friend class IBNormalPlaneMinimizationVariablesEvaluatorPimpl;
    class IBNormalPlaneMinimizationVariablesEvaluatorPimpl *d;
};

inline void IBNormalPlaneMinimizationVariablesEvaluator::init_properties() {
    $_init();
    $$.use_free_rotation = 0;
    $$.alphaXZ = 1; // 0 = X ---> 1 = Z
}


#endif // IBNORMALPLANEMINIMIZATIONVARIABLESEVALUATOR_H
