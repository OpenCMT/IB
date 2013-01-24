
#include <Core/Vector.h>
#include "IBVoxCollectionCap.h"
#include "IBAnalyzerWTrackLengths.h"
#include "IBVoxRaytracer.h"
#include "IBPocaEvaluator.h"

using namespace uLib;


///// PIMPL ////////////////////////////////////////////////////////////////////


class IBAnalyzerWTrackLengthsPimpl {
public:
    struct Event {
        struct Element {
            Matrix4f Wij;
            IBVoxel *voxel;
        };
        Vector4f Variables;
        Vector<Element> elements;
    };

    IBAnalyzerWTrackLengthsPimpl() :
        m_RayAlgorithm(NULL),
        m_PocaAlgorithm(NULL),
        m_VarAlgorithm(NULL)
    {}

    void Project(Event *evc) {
        IBVoxel *vox;
        for (unsigned int j = 0; j < evc->elements.size(); ++j) {
            vox = evc->elements[j].voxel;
            vox->Value += ( pow(fabs(evc->Variables(0)),2) + pow(fabs(evc->Variables(2)),2)  ) *
                    evc->elements[j].Wij(0,0);
            vox->Count++;
        }
    }

    // members //
    Vector<Event> m_Events;
    VoxRaytracer *m_RayAlgorithm;
    IBPocaEvaluator *m_PocaAlgorithm;
    IBMinimizationVariablesEvaluator *m_VarAlgorithm;
};




IBAnalyzerWTrackLengths::IBAnalyzerWTrackLengths() :
    d(new IBAnalyzerWTrackLengthsPimpl)
{}

IBAnalyzerWTrackLengths::~IBAnalyzerWTrackLengths()
{
    delete d;
}

void IBAnalyzerWTrackLengths::AddMuon(MuonScatterData &muon)
{
    if(!d->m_RayAlgorithm || !d->m_PocaAlgorithm || !d->m_VarAlgorithm) {
        std::cerr << "not all parameters setted\n";
        return;
    }
    IBAnalyzerWTrackLengthsPimpl::Event evc;

    // VARIABLES //
    if(likely(d->m_VarAlgorithm->evaluate(muon))) {
        evc.Variables = d->m_VarAlgorithm->getDataVector();

    }
    else return;

    // RAY //
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

    // LENGTHS //
    IBAnalyzerWTrackLengthsPimpl::Event::Element elc;
    Scalarf T = ray.TotalLength();
    for(int i=0; i<ray.Data().size(); ++i)
    {
        const IBVoxRaytracer::RayData::Element *el = &ray.Data().at(i);
        elc.voxel = &this->GetVoxCollection()->operator [](el->vox_id);
        Scalarf L = el->L;
        T -= L;

        Matrix2f wij_block;
        wij_block << L ,          L*L/2 + L*T,
                     L*L/2 + L*T, L*L*L/3 + L*L*T + L*T*T;
        elc.Wij = Matrix4f::Zero();
        elc.Wij.block<2,2>(0,0) = wij_block;
        elc.Wij.block<2,2>(2,2) = wij_block;

        evc.elements.push_back(elc);
    }
    d->m_Events.push_back(evc);
}

void IBAnalyzerWTrackLengths::Run(unsigned int iterations, float muons_ratio)
{
    for(int i=0; i<d->m_Events.size(); ++i)
        d->Project(&d->m_Events[i]);
}

void IBAnalyzerWTrackLengths::SetRaytracer(IBVoxRaytracer *raytracer)
{
    d->m_RayAlgorithm = raytracer;
}

void IBAnalyzerWTrackLengths::SetPocaAlgorithm(IBPocaEvaluator *evaluator)
{
    d->m_PocaAlgorithm = evaluator;
}

void IBAnalyzerWTrackLengths::SetVaraiblesAlgorithm(IBMinimizationVariablesEvaluator *algorithm)
{
    d->m_VarAlgorithm = algorithm;
}
