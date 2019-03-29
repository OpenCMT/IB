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


#include "IBTiltedAxisPocaEvaluator.h"

using namespace uLib;

IBTiltedAxisPocaEvaluator::IBTiltedAxisPocaEvaluator()
{
    m_integrity = true;
}

IBTiltedAxisPocaEvaluator::~IBTiltedAxisPocaEvaluator()
{}

bool IBTiltedAxisPocaEvaluator::evaluate(MuonScatterData muon)
{
    m_muon = muon;
    m_integrity = true;

    Vector4f  diff = m_muon.LineIn().origin - m_muon.LineOut().origin;
    HVector3f new_axis_direction(diff(0)/diff(1), 1, diff(2)/diff(1));
    Matrix4f  rotation = getRotationMatrix(new_axis_direction);

    Scalarf   dtphi, dttheta, distance, yp, yt, np, nt, s2yp, s2yt;
    Scalarf   e2pi = m_muon.ErrorIn().direction_error(0),
              e2po = m_muon.ErrorOut().direction_error(0),
              e2ti = (m_muon.ErrorIn().direction_error(2) ==0.f) ? m_muon.ErrorIn().direction_error(1)
                                                                 : m_muon.ErrorIn().direction_error(2),
              e2to = (m_muon.ErrorOut().direction_error(2)==0.f) ? m_muon.ErrorOut().direction_error(1)
                                                                 : m_muon.ErrorOut().direction_error(2);
    e2pi *= e2pi;
    e2po *= e2po;
    e2ti *= e2ti;
    e2to *= e2to;

    Scalarf   poca_x=0, poca_y=0, poca_z=0;

    HVector3f new_in_direction  = rotation * this->getDirectorCosines(m_muon.LineIn().direction);
    HVector3f new_out_direction = rotation * this->getDirectorCosines(m_muon.LineOut().direction);
    new_in_direction(0)  /= new_in_direction(1);
    new_in_direction(2)  /= new_in_direction(1);
    new_out_direction(0) /= new_out_direction(1);
    new_out_direction(2) /= new_out_direction(1);

    dtphi    = new_in_direction(0) - new_out_direction(0);
    dttheta  = new_in_direction(2) - new_out_direction(2);
    distance = diff.head(3).norm();
    yp = (new_in_direction(0) *distance) / dtphi;
    np = (new_out_direction(0)*distance) / dtphi;
    yt = (new_in_direction(2) *distance) / dttheta;
    nt = (new_out_direction(2)*distance) / dttheta;


    if (unlikely(!isFinite(yp))) {
        poca_y = yt;
        if (unlikely(!isFinite(yt))) {
            m_integrity = false;
        }
    } else  if (unlikely(!isFinite(yt))) {
            poca_y = yp;
    } else {
        s2yp = (np*np*e2pi)/(dtphi*dtphi) + (yp*yp*e2po)/(dtphi*dtphi);
        s2yt = (nt*nt*e2ti)/(dttheta*dttheta) + (yt*yt*e2to)/(dttheta*dttheta);
        poca_y = (yp/s2yp + yt/s2yt)/(1./s2yp + 1./s2yt);
    }
    if (unlikely(!isFinite(poca_y))){
        m_integrity = false;
    }
    poca_x = (pow(cos(atan(new_in_direction(0))),2)  * (new_in_direction(0)*(poca_y-distance)) +
              pow(cos(atan(new_out_direction(0))),2) * (new_out_direction(0)* (poca_y)))       /
             (pow(cos(atan(new_in_direction(0))),2)  + pow(cos(atan(new_out_direction(0))),2));

    poca_z = (pow(cos(atan(new_in_direction(2))),2)  * (new_in_direction(2)*(poca_y-distance)) +
              pow(cos(atan(new_out_direction(2))),2) * (new_out_direction(2)*(poca_y)))        /
             (pow(cos(atan(new_in_direction(2))),2)  + pow(cos(atan(new_out_direction(2))),2));

    HVector3f poca(poca_x, poca_y, poca_z);
    Matrix4f inverse;
    bool invertible = false;
    rotation.computeInverseWithCheck(inverse, invertible); // let it transpose it!
    if (unlikely(!invertible)){
        m_integrity = false;
    }
    poca = inverse * poca;
    poca += m_muon.LineOut().origin;
    m_poca(0) = poca(0);
    m_poca(1) = poca(1);
    m_poca(2) = poca(2);
    if (unlikely(!isFinite(poca(0))||!isFinite(poca(1))||!isFinite(poca(2)))) {
        m_integrity = false;
    }
    return m_integrity;
}

HPoint3f IBTiltedAxisPocaEvaluator::getPoca()
{
    return m_poca;
}

Matrix4f IBTiltedAxisPocaEvaluator::getRotationMatrix(const HVector3f &track_direction)
{
    HVector3f dc = this->getDirectorCosines(track_direction);

    Scalarf theta = acos(dc(1));
    Scalarf phi = atan2(dc(2),dc(0));

    Matrix4f y_rotation  = compileYRotation(phi);
    Matrix4f z_rotation  = compileZRotation(theta);
    Matrix4f out = z_rotation * y_rotation;
    return out;
}

Matrix4f IBTiltedAxisPocaEvaluator::compileYRotation(Scalarf angle)
{
    Matrix4f out;
    out << cos(angle), 0, sin(angle), 0,
           0,          1, 0,          0,
          -sin(angle), 0, cos(angle), 0,
           0,          0, 0,          1;
    return out;
}

Matrix4f IBTiltedAxisPocaEvaluator::compileZRotation(Scalarf angle)
{
    Matrix4f out;
    out << cos(angle), -sin(angle), 0, 0,
           sin(angle),  cos(angle), 0, 0,
           0,           0,          1, 0,
           0,           0,          0, 1;
    return out;
}

HVector3f IBTiltedAxisPocaEvaluator::getDirectorCosines(const HVector3f &track_direction)
{
    return track_direction / track_direction.head(3).norm();
}

