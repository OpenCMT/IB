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
    virtual HPoint3f getInTrackPoca() = 0;
    virtual HPoint3f getOutTrackPoca() = 0;
    virtual Scalarf getDistance() = 0;
protected:
    IBPocaEvaluator() {}
};

#endif // IBPOCAEVALUATOR_H
