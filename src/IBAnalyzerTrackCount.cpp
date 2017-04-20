/*//////////////////////////////////////////////////////////////////////////////
// CMT Cosmic Muon Tomography project //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

  Copyright (c) 2014, Universita' degli Studi di Padova, INFN sez. di Padova

  Coordinators: Prof. Gianni Zumerle < gianni.zumerle@pd.infn.it >
                Paolo Checchia       < paolo.checchia@pd.infn.it >

  Authors: Andrea Rigoni Garola < andrea.rigoni@pd.infn.it >
           Matteo Furlan        < nuright@gmail.com >
           Sara Vanini          < sara.vanini@pd.infn.it >

  All rights reserved
  ------------------------------------------------------------------

  This file can not be copied and/or distributed without the express
  permission of  Prof. Gianni Zumerle  < gianni.zumerle@pd.infn.it >

//////////////////////////////////////////////////////////////////////////////*/
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
        m_PocaAlgorithm(NULL),
        m_detSgnZ(0)
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
    //20170420 select detector based on Z coordinate: -1, 1, 0 if not used
    int m_detSgnZ;
};


IBAnalyzerTrackCount::IBAnalyzerTrackCount() :
    d(new IBAnalyzerTrackCountPimpl)
{}

IBAnalyzerTrackCount::~IBAnalyzerTrackCount()
{
    delete d;
}

bool IBAnalyzerTrackCount::AddMuon(const MuonScatterData &muon)
{
    if(!d->m_RayAlgorithm || !d->m_PocaAlgorithm) return false;
    IBAnalyzerTrackCountPimpl::Event evc;

    IBVoxRaytracer::RayData ray;
    // ENTRY and EXIT point present
    if( !std::isnan(muon.LineOut().origin.prod()) )
    { // Get RayTrace RayData //
        HPoint3f entry_pt,
                poca,
                exit_pt;
        if( !d->m_RayAlgorithm->GetEntryPoint(muon.LineIn(),entry_pt) ||
                !d->m_RayAlgorithm->GetExitPoint(muon.LineOut(),exit_pt) )
            return false;

        // 20170420 SV - select only one detector based on Z coordinate
        if(d->m_detSgnZ  && entry_pt[2]*d->m_detSgnZ < 0)
            return false;

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
    } else // Get RayTrace Data for stopping muon //
        ray = d->m_RayAlgorithm->TraceLine(muon.LineIn());

    IBAnalyzerTrackCountPimpl::Event::Element elc;
    Scalarf T = ray.TotalLength();
    for(int i=0; i<ray.Data().size(); ++i)
    {
        const IBVoxRaytracer::RayData::Element *el = &ray.Data().at(i);
        elc.voxel = &this->GetVoxCollection()->operator [](el->vox_id);
        evc.elements.push_back(elc);
    }
    d->m_Events.push_back(evc);
    return true;
}

void IBAnalyzerTrackCount::SetMuonCollection(IBMuonCollection *muons)
{
    uLibAssert(muons);
    d->m_Events.clear();
    for(int i=0; i<muons->size(); ++i)
    {
        this->AddMuon(muons->At(i));
    }
    BaseClass::SetMuonCollection(muons);
}



void IBAnalyzerTrackCount::Run(unsigned int iterations, float muons_ratio)
{
    for(int i=0; i<d->m_Events.size(); ++i)
        d->Project(&d->m_Events[i]);
}

void IBAnalyzerTrackCount::SetRayAlgorithm(IBVoxRaytracer *raytracer)
{
    d->m_RayAlgorithm = raytracer;
}

void IBAnalyzerTrackCount::SetPocaAlgorithm(IBPocaEvaluator *evaluator)
{
    d->m_PocaAlgorithm = evaluator;
}

void IBAnalyzerTrackCount::SetDetectorZSelection(int selectZ)
{
    d->m_detSgnZ = selectZ;
}

void IBAnalyzerTrackCount::Clear()
{
    d->m_Events.clear();
}

unsigned int IBAnalyzerTrackCount::Size() const
{
    return d->m_Events.size();
}





