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



#include <Core/Vector.h>
#include "IBVoxCollectionCap.h"
#include "IBAnalyzerWTrackLengths.h"
#include "IBVoxRaytracer.h"
#include "IBPocaEvaluator.h"

using namespace uLib;


IBAnalyzerWTrackLengths::IBAnalyzerWTrackLengths()
{}

IBAnalyzerWTrackLengths::~IBAnalyzerWTrackLengths()
{}

bool IBAnalyzerWTrackLengths::AddMuon(const MuonScatterData &muon)
{
    if(!m_RayAlgorithm || !m_PocaAlgorithm /*|| !m_VarAlgorithm*/) {
        std::cerr << "not all parameters set\n";
        return false;
    }

    IBAnalyzerWTrackLengths::Event evc;

    evc.Momentum = muon.GetMomentum();

    // VARIABLES //
    if(likely(m_VarAlgorithm && m_VarAlgorithm->evaluate(muon))) {
        evc.Variables = m_VarAlgorithm->getDataVector();
    }
    else {
        Vector3f in, out;
        in  = muon.LineIn().direction.head(3);
        out = muon.LineOut().direction.head(3);
        float a = in.transpose() * out;
        a = acos(a / (in.norm() * out.norm()) );
        if(uLib::isFinite(a)) evc.Variables(0) = pow(a,2);
        else evc.Variables(0) = 0;
    }

    // RAY //
    IBVoxRaytracer::RayData ray;
    { // Get RayTrace RayData //
        HPoint3f entry_pt,poca,exit_pt;
        if( !m_RayAlgorithm->GetEntryPoint(muon.LineIn(),entry_pt) ||
                !m_RayAlgorithm->GetExitPoint(muon.LineOut(),exit_pt) )
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
    }

    // LENGTHS //
    IBAnalyzerWTrackLengths::Event::Element elc;
    Scalarf T = ray.TotalLength();
    evc.Length = T;
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
//    m_Events.push_back(evc);
    Project(&evc);
    return true;
}

void IBAnalyzerWTrackLengths::Run(unsigned int iterations, float muons_ratio)
{
    std::cerr << "WARNING: Run function does nothing ... as the projection is "
                 "made inside AddMuon funcion (without accounting events).";
    //    for(int i=0; i<m_Events.size(); ++i)
    //        Project(&m_Events[i]);
}

void IBAnalyzerWTrackLengths::SetRayAlgorithm(IBVoxRaytracer *raytracer)
{
    m_RayAlgorithm = raytracer;
}

void IBAnalyzerWTrackLengths::SetPocaAlgorithm(IBPocaEvaluator *evaluator)
{
    m_PocaAlgorithm = evaluator;
}

void IBAnalyzerWTrackLengths::SetVarAlgorithm(IBMinimizationVariablesEvaluator *algorithm)
{
    m_VarAlgorithm = algorithm;
}

void IBAnalyzerWTrackLengths::SetPocaProximity(float sigma)
{
    // FINIRE //
    //    m_poca_proximity_rms = p;
}

void IBAnalyzerWTrackLengths::Project(IBAnalyzerWTrackLengths::Event *evc) {
    IBVoxel *vox;
    for (unsigned int j = 0; j < evc->elements.size(); ++j) {
        vox = evc->elements[j].voxel;
        float a;
        if(m_VarAlgorithm)
            a = ( pow(fabs(evc->Variables(0)),2) + pow(fabs(evc->Variables(2)),2)  ) *
                    evc->elements[j].Wij(0,0) * pow(evc->Momentum,2) * 1.5E-6 / evc->Length;
        else
            a = evc->Variables(0) * pow(evc->Momentum,2) * 1.5E-6 * evc->elements[j].Wij(0,0) / evc->Length;
        // FINIRE //
        //  if(m_poca_proximity_rms > 0) {
        //    float d = 1/sqrt(2*M_PI)/fabs(m_poca_proximity_rms) ;// ...
        //  }
        if(a<1E-6) {
            vox->Value += a;
            vox->Count++;
        }
    }
}

