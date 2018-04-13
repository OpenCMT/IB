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



#include <assert.h>

#include <Core/Vector.h>
//#include <Math/Utils.h>
#include "IBAnalyzerPoca.h"
#include "IBPocaEvaluator.h"
#include "IBVoxCollectionCap.h"

using namespace uLib;

class IBAnalyzerPocaPimpl {
public:
    IBAnalyzerPocaPimpl(IBAnalyzerPoca *pt) :
        m_pt(pt),
        m_PocaAlgorithm(NULL)
    {}

    bool CollectMuon(const MuonScatterData &muon)
    {
        assert(m_PocaAlgorithm);
        if(m_PocaAlgorithm->evaluate(muon)) {

            Vector4f poca = m_PocaAlgorithm->getPoca();

            //---- Check that the POCA is valid
            Vector4f in, out;
            in  = poca - muon.LineIn().origin();
            out = muon.LineOut().origin() - poca;
            float poca_prj = in.transpose() * out;
            bool validPoca = poca_prj > 0;

            if(validPoca)
                m_Data.push_back(m_PocaAlgorithm->getPoca());
            return true;
        } else return false;
    }

    void SetVoxels(IBVoxCollection *voxels)
    {
        for(int i=0; i<m_Data.size(); ++i) {
            Vector3i id = voxels->Find(m_Data[i]);
            if(voxels->IsInsideGrid(id))
                voxels->operator [](id).Value += 1;
        }
    }

    // members //
    IBAnalyzerPoca  *m_pt;
    IBPocaEvaluator *m_PocaAlgorithm;
    Vector<Vector4f> m_Data;
};




IBAnalyzerPoca::IBAnalyzerPoca() :
    d(new IBAnalyzerPocaPimpl(this))
{}

IBAnalyzerPoca::~IBAnalyzerPoca()
{
    delete d;
}


bool IBAnalyzerPoca::AddMuon(const MuonScatterData &event) {
    return d->CollectMuon(event);
}

void IBAnalyzerPoca::Run(unsigned int iterations, float muons_ratio) {
    d->SetVoxels((IBVoxCollection *)this->GetVoxCollection());
}

void IBAnalyzerPoca::SetPocaAlgorithm(IBPocaEvaluator *poca)
{   d->m_PocaAlgorithm = poca;  }

void IBAnalyzerPoca::SetMuonCollection(IBMuonCollection *muons)
{
    uLibAssert(muons);
    for(int i=0; i<muons->size(); ++i)
    {
        this->AddMuon(muons->At(i));
    }
    BaseClass::SetMuonCollection(muons);
}

