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


#include "IBLineDistancePocaEvaluator.h"
#include <iostream>


IBLineDistancePocaEvaluator::IBLineDistancePocaEvaluator() {
        m_integrity = true;
        m_cutlength = 0.f;
        m_dist = NAN;
        m_poca << 0,0,0,1;
}

IBLineDistancePocaEvaluator::~IBLineDistancePocaEvaluator()
{}

HVector3f IBLineDistancePocaEvaluator::getDirectorCosines(const HVector3f &track_direction)
{
    return track_direction / track_direction.head(3).norm();
}


bool IBLineDistancePocaEvaluator::evaluate(MuonScatterData muon)
{
    m_muon = muon;
    m_integrity = true;

    m_dist = NAN;

    HPoint3f p    = m_muon.LineIn().origin;
    HPoint3f q    = m_muon.LineOut().origin;
    Vector4f v    = getDirectorCosines(m_muon.LineIn().direction);
    Vector4f w    = getDirectorCosines(m_muon.LineOut().direction);
    Vector4f diff = q - p;
    Scalarf  prod = v.transpose() * w;
    Scalarf  den  = 1./(1. - prod*prod);
    if (unlikely(!isFinite(den))) {
            float pr = diff.dot(v);
            Vector3f proj = (v * pr).head(3);
            Vector3f dist = diff.head(3);
            dist = proj.cross(dist);
            dist /= pr;
            m_dist = dist.head(3).norm();
            if (unlikely(m_cutlength != 0.f && dist.head(3).norm()>m_cutlength)) {
                m_integrity = false;
            }
            m_integrity = false;
    }
    else {
        Scalarf lambda = (diff.transpose()*(v*den-(w*prod*den)));
        Scalarf mu     = (diff.transpose()*((v*prod)*den-w*den));
        m_inPoca       = p + v*lambda;
        m_outPoca      = q + w*mu;
        Vector4f mdseg = (m_outPoca-m_inPoca);
        if (unlikely(m_cutlength != 0.f && mdseg.head(3).norm() > m_cutlength)) {
            m_integrity = false;
        }
        m_dist = mdseg.head(3).norm();
        m_poca = m_inPoca + (mdseg*0.5);
        m_integrity   = true;
    }

    return m_integrity;
}

HPoint3f IBLineDistancePocaEvaluator::getPoca()
{
    return m_poca;
}

HPoint3f IBLineDistancePocaEvaluator::getInTrackPoca()
{
    return m_inPoca;
}

void IBLineDistancePocaEvaluator::setDistanceCut(Scalarf length)
{
    m_cutlength = length;
}

Scalarf IBLineDistancePocaEvaluator::getDistance()
{
    return m_dist;
}

HPoint3f IBLineDistancePocaEvaluator::getOutTrackPoca()
{
    return m_outPoca;
}
