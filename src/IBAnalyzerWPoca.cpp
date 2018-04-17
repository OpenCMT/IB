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


#include "TTree.h"
#include "TFile.h"

#include "IBAnalyzerWPoca.h"
#include "IBMinimizationVariablesEvaluator.h"
#include "IBPocaEvaluator.h"
#include "IBVoxCollectionCap.h"

class IBAnalyzerWPocaPimpl {
    struct Data {
        HPoint3f poca;
        Scalarf  weight;
    };

public:
    IBAnalyzerWPocaPimpl()
    {
        m_PocaAlgorithm = NULL;
        m_Minimizator = NULL;

#ifndef NDEBUG
        m_out   = new TFile("dump_poca_weighting.root", "RECREATE");
        m_tree = new TTree("w", "weighting");
        m_tree->Branch("w",   &tmp.weight, "w");
#endif
    }

    bool CollectMuon(const MuonScatterData &muon)
    {
        assert(m_PocaAlgorithm);
        if (m_PocaAlgorithm->evaluate(muon)) {
            tmp.poca   = m_PocaAlgorithm->getPoca();
            if(m_Minimizator && m_Minimizator->evaluate(muon)) {
                // weight with two views scattering angles
                Scalarf t_w_1 = pow(tan((m_Minimizator->getDataVector(0))), 2);
                Scalarf t_w_2 = pow(tan((m_Minimizator->getDataVector(2))), 2);
                tmp.weight = (t_w_1 + t_w_2) * pow(muon.GetMomentum(),2) * 1.5E-6;
                //tmp.weight = (t_w_1 + t_w_2);
            }
            else {
                // weight with angle in 3D space
                Vector3f in, out;
                in  = muon.LineIn().direction.head(3);
                out = muon.LineOut().direction.head(3);
                float a = in.transpose() * out;
                a = fabs( acos(a / (in.norm() * out.norm())) );
                if(uLib::isFinite(a))
                    tmp.weight = pow(a * muon.GetMomentum(),2) * 1.5E-6;
                    //tmp.weight = pow(a,2);
                else tmp.weight = 0;
            }
            m_Data.push_back(tmp);
#           ifndef NDEBUG
            m_tree->Fill();
#           endif
            return true;
        }
        return false;
    }

    void SetVoxels(IBVoxCollection *voxels) {
        for (int i=0; i<m_Data.size(); ++i) {
            Vector3i id = voxels->Find(m_Data[i].poca);
            if (voxels->IsInsideGrid(id)) {
                IBVoxel &vox = voxels->operator [](id);
                vox.Value += m_Data[i].weight;
                vox.Count++;
            }
        }
    }

    IBPocaEvaluator                  *m_PocaAlgorithm;
    IBMinimizationVariablesEvaluator *m_Minimizator;
    Vector<Data>                      m_Data;

    Data    tmp;
    Scalarf t_w_1,
            t_w_2;

#ifndef NDEBUG
    TFile* m_out;
    TTree* m_tree;
#endif

};

IBAnalyzerWPoca::IBAnalyzerWPoca() :
    d(new IBAnalyzerWPocaPimpl)
{}

IBAnalyzerWPoca::~IBAnalyzerWPoca()
{
#ifndef NDEBUG
    d->m_out->cd();
    d->m_tree->Write();
    d->m_out->Close();
#endif
    delete d;
}

bool IBAnalyzerWPoca::AddMuon(const MuonScatterData &event)
{
    return d->CollectMuon(event);
}

void IBAnalyzerWPoca::Run(unsigned int iteration, float muons_ratio)
{
    d->SetVoxels((IBVoxCollection*)this->GetVoxCollection());
}

void IBAnalyzerWPoca::SetPocaAlgorithm(IBPocaEvaluator *poca)
{
    d->m_PocaAlgorithm = poca;
}

void IBAnalyzerWPoca::SetVarAlgorithm(IBMinimizationVariablesEvaluator *evaluator)
{
    d->m_Minimizator = evaluator;
}
