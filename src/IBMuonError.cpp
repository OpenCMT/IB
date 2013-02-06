#include <IBMuonError.h>

using namespace uLib;

IBMuonError::IBMuonError(Scalarf xA, Scalarf zA, Scalarf ratio)
{
    m_Axi = xA;
    m_Azi = zA;
    m_Axo = ratio * xA;
    m_Azo = ratio * zA;
    IBMESimpler sp(this);
    m_simpler = &sp;
    m_shader  = NULL;
}

bool IBMuonError::evaluate(MuonScatter &event, int i, int j)
{
    if (m_shader == NULL) {
        return m_simpler->evaluate(event, i, j);
    } else {
        return m_shader->evaluate(event, i, j);
    }
    return false;
}

void IBMuonError::setScrapsImage(IBLightCollection &image, bool evPM)
{
    m_Axo = m_Axi;
    m_Azo = m_Axi;
    IBPocaEvaluator * pproc = IBPocaEvaluator::New(IBPocaEvaluator::LineDistance);
    IBVoxRaytracer trace(image);
    IBMEShader sh(this);
    m_shader = &sh;
    m_shader->m_image  = &image;
    m_shader->m_tracer = &trace;
    m_shader->m_pproc  = pproc;
    m_shader->m_evPM   = evPM;
}


///////////////////////////////////////////////////////////////////
//// SIMPLER //////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

bool IBMuonError::IBMESimpler::evaluate(MuonScatter &event, int i, int j)
{
    event.ErrorIn().direction_error(0)  = d->mpdEval(d->m_Axi, event.GetMomentum(), event.LineIn().direction(0));
    event.ErrorIn().direction_error(i)  = d->mpdEval(d->m_Azi, event.GetMomentum(), event.LineIn().direction(i));
    event.ErrorOut().direction_error(0) = d->mpdEval(d->m_Axo, event.GetMomentum(), event.LineOut().direction(0));
    event.ErrorOut().direction_error(j) = d->mpdEval(d->m_Azo, event.GetMomentum(), event.LineOut().direction(j));
    return true;
}


///////////////////////////////////////////////////////////////////
//// SHADER ///////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

bool IBMuonError::IBMEShader::evaluate(MuonScatter &event, int i, int j)
{
    IBVoxRaytracer::RayData ray;
    {
        HPoint3f entry_pt, poca, exit_pt;
        if (!m_tracer->GetEntryPoint(event.LineIn(),entry_pt) ||
                !m_tracer->GetExitPoint(event.LineOut(),exit_pt))
            return false;
        bool test = m_pproc->evaluate(event);
        poca = m_pproc->getPoca();
        if (test && m_image->IsInsideBounds(poca)) {
            ray = m_tracer->TraceBetweenPoints(entry_pt, poca);
            ray.AppendRay(m_tracer->TraceBetweenPoints(poca, exit_pt));
        } else {
            ray = m_tracer->TraceBetweenPoints(entry_pt, exit_pt);
        }
    }
    float cL = 0;
    for (int ii=0; ii<ray.Data().size(); ++ii) {
        const IBVoxRaytracer::RayData::Element *el = &ray.Data().at(ii);
        float L   = el->L;
        float val = (m_image->operator [](el->vox_id)).Value;
        cL += val*L;
    }
    float pi = event.GetMomentum();
    float po = sqrt((pi*pi)/(1+pi*pi*cL));
    event.ErrorIn().direction_error(0)  = d->mpdEval(d->m_Axi, pi, event.LineIn().direction(0));
    event.ErrorIn().direction_error(i)  = d->mpdEval(d->m_Azi, pi, event.LineIn().direction(i));
    event.ErrorOut().direction_error(0) = d->mpdEval(d->m_Axo, po, event.LineOut().direction(0));
    event.ErrorOut().direction_error(j) = d->mpdEval(d->m_Azo, po, event.LineOut().direction(j));
    event.SetMomentum((pi+po)/2.);

}

Scalarf IBMuonError::mpdEval(Scalarf a, Scalarf p, Scalarf d)
{
    return 1E-3 * sqrt(a*a/(p*p)) * (1+d*d);
}

