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



#ifndef IBLINEDISTANCEPOCAEVALUATOR_H
#define IBLINEDISTANCEPOCAEVALUATOR_H

#include "IBPocaEvaluator.h"

class IBLineDistancePocaEvaluator : public IBPocaEvaluator
{
public:
    IBLineDistancePocaEvaluator();
    ~IBLineDistancePocaEvaluator();

    bool evaluate(MuonScatterData muon);
    HPoint3f getPoca();
    HPoint3f getOutTrackPoca();
    HPoint3f getInTrackPoca();
    void setDistanceCut(Scalarf length);
    Scalarf getDistance();

private:

    HVector3f getDirectorCosines(const HVector3f &track_direction);

    //TODO probably useless
    bool m_integrity;

    Scalarf m_cutlength;
    Scalarf m_dist;
    MuonScatterData m_muon;
    HPoint3f m_poca;
    HPoint3f m_inPoca;
    HPoint3f m_outPoca;
};

#endif // IBLINEDISTANCEPOCAEVALUATOR_H
