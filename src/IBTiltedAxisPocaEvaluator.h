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



#ifndef IBTILTEDAXISPOCAEVALUATOR_H
#define IBTILTEDAXISPOCAEVALUATOR_H

#include "IBPocaEvaluator.h"
#include "IBMuonError.h"

class IBTiltedAxisPocaEvaluator : public IBPocaEvaluator
{
public:

    IBTiltedAxisPocaEvaluator();
    ~IBTiltedAxisPocaEvaluator();

    bool     evaluate(MuonScatterData muon);
    HPoint3f getPoca();
    HPoint3f getInTrackPoca(){ return HPoint3f();}
    HPoint3f getOutTrackPoca(){ return HPoint3f();}

    // dummy functions for the moment... to be implemented
    inline Scalarf getDistance() {return 0;}

private:

    Matrix4f getRotationMatrix(const HVector3f &track_direction);
    Matrix4f compileYRotation(Scalarf angle);
    Matrix4f compileZRotation(Scalarf angle);
    HVector3f getDirectorCosines(const HVector3f &track_direction);

    MuonScatterData m_muon;
    HPoint3f        m_poca;
    bool            m_integrity;

};

#endif // IBTILTEDAXISPOCAEVALUATOR_H
