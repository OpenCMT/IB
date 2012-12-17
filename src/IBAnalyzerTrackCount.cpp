#include <Core/Vector.h>
#include "IBVoxCollectionCap.h"
#include "IBAnalyzerTrackCount.h"
#include "IBVoxRaytracer.h"

using namespace uLib;



///// PIMPL ////////////////////////////////////////////////////////////////////

class IBAnalyzerTrackCountPimpl {
public:
    struct Event {
        struct Element {
            IBVoxel *voxel;
        };
        Vector<Element> elements;
    };
public:
    IBAnalyzerTrackCountPimpl() :
        m_RayAlgorithm(NULL),
        m_PocaAlgorithm(NULL)
    {}


    void Project(Event *evc) {
        IBVoxel *vox;
        for (unsigned int j = 0; j < evc->elements.size(); ++j) {
            vox = evc->elements[j].voxel;
            vox->Value += 1; // track count in density value
        }
    }

    // members //
    Vector<Event>    m_Events;
    VoxRaytracer    *m_RayAlgorithm;
    IBPocaEvaluator *m_PocaAlgorithm;
};





IBAnalyzerTrackCount::IBAnalyzerTrackCount() :
    d(new IBAnalyzerTrackCountPimpl)
{}

IBAnalyzerTrackCount::~IBAnalyzerTrackCount()
{
    delete d;
}



void IBAnalyzerTrackCount::AddMuon(MuonScatterData &muon)
{
    if(!d->m_RayAlgorithm || !d->m_PocaAlgorithm) return;
    IBAnalyzerTrackCountPimpl::Event evc;

    IBVoxRaytracer::RayData ray;
    { // Get RayTrace RayData //
        HPoint3f entry_pt,poca,exit_pt;
        if( !d->m_RayAlgorithm->GetEntryPoint(muon.LineIn(),entry_pt) ||
                !d->m_RayAlgorithm->GetExitPoint(muon.LineOut(),exit_pt) )
            return;
        bool test = d->m_PocaAlgorithm->evaluate(muon);
        poca = d->m_PocaAlgorithm->getPoca();
        if(test && this->GetVoxCollection()->IsInsideBounds(poca)) {
            poca = d->m_PocaAlgorithm->getPoca();
            ray = d->m_RayAlgorithm->TraceBetweenPoints(entry_pt,poca);
            ray.AppendRay( d->m_RayAlgorithm->TraceBetweenPoints(poca,exit_pt) );
        }
        else {
            ray = d->m_RayAlgorithm->TraceBetweenPoints(entry_pt,exit_pt);
        }
    }

    IBAnalyzerTrackCountPimpl::Event::Element elc;
    Scalarf T = ray.TotalLength();
    for(int i=0; i<ray.Data().size(); ++i)
    {
        const IBVoxRaytracer::RayData::Element *el = &ray.Data().at(i);
        elc.voxel = &this->GetVoxCollection()->operator [](el->vox_id);
        evc.elements.push_back(elc);
    }
    d->m_Events.push_back(evc);
}

void IBAnalyzerTrackCount::Run(unsigned int iterations, float muons_ratio)
{
    for(int i=0; i<d->m_Events.size(); ++i)
        d->Project(&d->m_Events[i]);
}

void IBAnalyzerTrackCount::SetRaytracer(IBVoxRaytracer *raytracer)
{
    d->m_RayAlgorithm = raytracer;
}

void IBAnalyzerTrackCount::SetPocaAlgorithm(IBPocaEvaluator *evaluator)
{
    d->m_PocaAlgorithm = evaluator;
}





