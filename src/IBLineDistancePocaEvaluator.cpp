#include "IBLineDistancePocaEvaluator.h"
#include <iostream>

class IBLineDistancePocaEvaluatorPimpl
{
public:
    IBLineDistancePocaEvaluatorPimpl()
    {
        m_integrity = true;
        m_poca << 0,0,0;
    }

    void evaluatePoCA()
    {
        HPoint3f p    = m_muon.LineIn().origin;
        HPoint3f q    = m_muon.LineOut().origin;
        Vector4f v    = getDirectorCosines(m_muon.LineIn().direction);
        Vector4f w    = getDirectorCosines(m_muon.LineOut().direction);
        Vector4f diff = q - p;
        Scalarf  prod = v.transpose() * w;
        Scalarf  den  = 1./(1 - prod*prod);
        m_inPoca      = p + v*(diff.transpose()*((v-(w*prod))*den));
        m_outPoca     = q + w*(diff.transpose()*(((v*prod)-w)*den));
        m_poca        = m_inPoca + (m_outPoca-m_inPoca)*0.5;
        m_integrity   = true;
    }

    HVector3f getDirectorCosines(const HVector3f &track_direction)
    {
        return track_direction / track_direction.head(3).norm();
    }

public:
    bool m_integrity;
    MuonScatterData m_muon;
    HPoint3f m_poca;
    HPoint3f m_inPoca;
    HPoint3f m_outPoca;
};


IBLineDistancePocaEvaluator::IBLineDistancePocaEvaluator() :
    d(new IBLineDistancePocaEvaluatorPimpl)
{
}

IBLineDistancePocaEvaluator::~IBLineDistancePocaEvaluator()
{
    delete d;
}

bool IBLineDistancePocaEvaluator::evaluate(MuonScatterData muon)
{
    d->m_muon = muon;
    d->m_integrity = true;
    d->evaluatePoCA();
    return d->m_integrity;
}

HPoint3f IBLineDistancePocaEvaluator::getPoca()
{
    return d->m_poca;
}

HPoint3f IBLineDistancePocaEvaluator::getInTrackPoca()
{
    return d->m_inPoca;
}

HPoint3f IBLineDistancePocaEvaluator::getOutTrackPoca()
{
    return d->m_outPoca;
}
