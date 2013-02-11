
#include <IBMuonError.h>

#include <IBLineDistancePocaEvaluator.h>


using namespace uLib;

IBMuonError::IBMuonError(Scalarf xA, Scalarf zA, Scalarf ratio) :
    m_simpler(new IBMESimpler(this)),
    m_shader(NULL)
{
    m_Ax = xA;
    m_Az = zA;
    m_pratio = ratio;
    m_averPcorr = false;
    m_azimPcorr = false;
}

IBMuonError::~IBMuonError()
{
    delete m_simpler;
    if(m_shader) {
        delete m_shader->m_pproc;
        delete m_shader->m_tracer;
        delete m_shader;
    }
}

bool IBMuonError::evaluate(MuonScatter &event, int i, int j)
{
    if (m_shader) {
        return m_shader->evaluate(event, i, j);
    } else {
        return m_simpler->evaluate(event, i, j);
    }
    return false;
}

void IBMuonError::azimuthalMomentumCorrection(bool enable)
{
    m_azimPcorr = enable;
}

void IBMuonError::averageMomentumCorrection(bool enable)
{
    m_averPcorr = enable;
}

void IBMuonError::setScrapsImage(IBLightCollection &image, bool evPM)
{
    if(m_shader) {
        delete m_shader->m_pproc;
        delete m_shader->m_tracer;
        delete m_shader;
    }

    m_shader = new IBMEShader(this);
    m_shader->m_pproc = new IBLineDistancePocaEvaluator();
    m_shader->m_tracer = new IBVoxRaytracer(image);
    m_shader->m_image = &image;
    m_shader->m_evPM = evPM;
}


///////////////////////////////////////////////////////////////////
//// SIMPLER //////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////

bool IBMuonError::IBMESimpler::evaluate(MuonScatter &event, int i, int j)
{
    if (unlikely(event.GetMomentum()!=0.f)) {
        event.SetMomentumPrime(event.GetMomentum()/d->m_pratio);
    } else {
        if (d->m_azimPcorr) {
            float azAngl = atan(sqrt((event.LineIn().direction(0)*event.LineIn().direction(0)+
                                      event.LineIn().direction(i)*event.LineIn().direction(i))));
            azAngl = cos(azAngl);
            float pSqInv = -0.7022*(azAngl*azAngl)+2.0807*azAngl+0.6215;
            event.SetMomentumPrime(sqrt(1./pSqInv));
        }
        event.SetMomentum(event.GetMomentumPrime()*d->m_pratio);
    }
    event.ErrorIn().direction_error(0)  = d->mpdEval(d->m_Ax, event.GetMomentum(), event.LineIn().direction(0));
    event.ErrorIn().direction_error(i)  = d->mpdEval(d->m_Az, event.GetMomentum(), event.LineIn().direction(i));
    event.ErrorOut().direction_error(0) = d->mpdEval(d->m_Ax, event.GetMomentumPrime(), event.LineOut().direction(0));
    event.ErrorOut().direction_error(j) = d->mpdEval(d->m_Az, event.GetMomentumPrime(), event.LineOut().direction(j));
    if (d->m_averPcorr) event.SetMomentum((event.GetMomentum()+event.GetMomentumPrime())/2.);
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
    event.ErrorIn().direction_error(0)  = d->mpdEval(d->m_Ax, pi, event.LineIn().direction(0));
    event.ErrorIn().direction_error(i)  = d->mpdEval(d->m_Az, pi, event.LineIn().direction(i));
    event.ErrorOut().direction_error(0) = d->mpdEval(d->m_Ax, po, event.LineOut().direction(0));
    event.ErrorOut().direction_error(j) = d->mpdEval(d->m_Az, po, event.LineOut().direction(j));
    if(m_evPM) event.SetMomentum((pi+po)/2.);

}

Scalarf IBMuonError::mpdEval(Scalarf a, Scalarf p, Scalarf d)
{
    return 1E-3 * sqrt(a*a/(p*p)) * (1+d*d);
}

