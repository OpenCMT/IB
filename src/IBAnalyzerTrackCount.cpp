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
#include "IBVoxCollectionCap.h"
#include "IBAnalyzerTrackCount.h"
#include "IBVoxRaytracer.h"

using namespace uLib;

IBAnalyzerTrackCount::IBAnalyzerTrackCount() :
    m_RayAlgorithm(NULL),
    m_PocaAlgorithm(NULL),
    m_detSgnZ(0)
{}

IBAnalyzerTrackCount::~IBAnalyzerTrackCount()
{}

bool IBAnalyzerTrackCount::AddMuon(const MuonScatterData &muon)
{
    if(!m_RayAlgorithm || !m_PocaAlgorithm) return false;
    IBAnalyzerTrackCount::Event evc;

    IBVoxRaytracer::RayData ray;
    // ENTRY and EXIT point present
    if( !std::isnan(muon.LineOut().origin.prod()) )
    { // Get RayTrace RayData //
        HPoint3f entry_pt,
                poca,
                exit_pt;
        if( !m_RayAlgorithm->GetEntryPoint(muon.LineIn(),entry_pt) ||
                !m_RayAlgorithm->GetExitPoint(muon.LineOut(),exit_pt) )
            return false;

        // 20170420 SV - select only one detector based on Z coordinate
        if(m_detSgnZ  && entry_pt[2]*m_detSgnZ < 0)
            return false;

        bool test = m_PocaAlgorithm->evaluate(muon);
        poca = m_PocaAlgorithm->getPoca();
        if(test && this->GetVoxCollection()->IsInsideBounds(poca)) {
            poca = m_PocaAlgorithm->getPoca();
            ray = m_RayAlgorithm->TraceBetweenPoints(entry_pt,poca);
            ray.AppendRay( m_RayAlgorithm->TraceBetweenPoints(poca,exit_pt) );
        }
        else {
            ray = m_RayAlgorithm->TraceBetweenPoints(entry_pt,exit_pt);
        }
    } else // Get RayTrace Data for stopping muon //
        ray = m_RayAlgorithm->TraceLine(muon.LineIn());

    IBAnalyzerTrackCount::Event::Element elc;
    Scalarf T = ray.TotalLength();
    for(int i=0; i<ray.Data().size(); ++i)
    {
        const IBVoxRaytracer::RayData::Element *el = &ray.Data().at(i);
        elc.voxel = &this->GetVoxCollection()->operator [](el->vox_id);
        evc.elements.push_back(elc);
    }
    m_Events.push_back(evc);
    return true;
}

void IBAnalyzerTrackCount::SetMuonCollection(IBMuonCollection *muons)
{
    uLibAssert(muons);
    m_Events.clear();
    for(int i=0; i<muons->size(); ++i)
    {
        this->AddMuon(muons->At(i));
    }
    BaseClass::SetMuonCollection(muons);
}



void IBAnalyzerTrackCount::Run(unsigned int iterations, float muons_ratio)
{
    for(int i=0; i<m_Events.size(); ++i) {

        IBAnalyzerTrackCount::Event *evc = &m_Events[i];
        IBVoxel *vox;

        for (unsigned int j = 0; j < evc->elements.size(); ++j) {
            vox = evc->elements[j].voxel;
            vox->Value += 1; // track count in density value
        }
    }
}

void IBAnalyzerTrackCount::SetRayAlgorithm(IBVoxRaytracer *raytracer)
{
    m_RayAlgorithm = raytracer;
}

void IBAnalyzerTrackCount::SetPocaAlgorithm(IBPocaEvaluator *evaluator)
{
    m_PocaAlgorithm = evaluator;
}

void IBAnalyzerTrackCount::SetDetectorZSelection(int selectZ)
{
    m_detSgnZ = selectZ;
}

void IBAnalyzerTrackCount::Clear()
{
    m_Events.clear();
}

unsigned int IBAnalyzerTrackCount::Size() const
{
    return m_Events.size();
}





