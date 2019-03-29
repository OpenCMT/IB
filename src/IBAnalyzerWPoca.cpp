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


IBAnalyzerWPoca::IBAnalyzerWPoca() {
    m_PocaAlgorithm = NULL;
    m_Minimizator = NULL;

#ifndef NDEBUG
    m_out   = new TFile("dump_poca_weighting.root", "RECREATE");
    m_tree = new TTree("w", "weighting");
    m_tree->Branch("w",   &tmp.weight, "w");
#endif
}

IBAnalyzerWPoca::~IBAnalyzerWPoca()
{
#ifndef NDEBUG
    m_out->cd();
    m_tree->Write();
    m_out->Close();
#endif
}

bool IBAnalyzerWPoca::AddMuon(const MuonScatterData &event)
{
    assert(m_PocaAlgorithm);
    if (m_PocaAlgorithm->evaluate(event)) {
        tmp.poca   = m_PocaAlgorithm->getPoca();
        if(m_Minimizator && m_Minimizator->evaluate(event)) {
            // weight with two views scattering angles
            Scalarf t_w_1 = pow(tan((m_Minimizator->getDataVector(0))), 2);
            Scalarf t_w_2 = pow(tan((m_Minimizator->getDataVector(2))), 2);
            tmp.weight = (t_w_1 + t_w_2) * pow(event.GetMomentum(),2) * 1.5E-6;
            //tmp.weight = (t_w_1 + t_w_2);
        }
        else {
            // weight with angle in 3D space
            Vector3f in, out;
            in  = event.LineIn().direction.head(3);
            out = event.LineOut().direction.head(3);
            float a = in.transpose() * out;
            a = fabs( acos(a / (in.norm() * out.norm())) );
            if(uLib::isFinite(a))
                tmp.weight = pow(a * event.GetMomentum(),2) * 1.5E-6;
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

void IBAnalyzerWPoca::Run(unsigned int iteration, float muons_ratio) {
    IBVoxCollection *voxels = (IBVoxCollection*)this->GetVoxCollection();

    for (int i=0; i<m_Data.size(); ++i) {
        Vector3i id = voxels->Find(m_Data[i].poca);
        if (voxels->IsInsideGrid(id)) {
            IBVoxel &vox = voxels->operator [](id);
            vox.Value += m_Data[i].weight;
            vox.Count++;
        }
    }

}

void IBAnalyzerWPoca::SetPocaAlgorithm(IBPocaEvaluator *poca)
{
    m_PocaAlgorithm = poca;
}

void IBAnalyzerWPoca::SetVarAlgorithm(IBMinimizationVariablesEvaluator *evaluator)
{
    m_Minimizator = evaluator;
}
