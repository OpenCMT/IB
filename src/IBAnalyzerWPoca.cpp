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
                Scalarf t_w_1 = pow(tan((m_Minimizator->getDataVector(0))), 2);
                Scalarf t_w_2 = pow(tan((m_Minimizator->getDataVector(2))), 2);
                tmp.weight = (t_w_1 + t_w_2) * pow(muon.GetMomentum(),2) * 1.5E-6;
            }
            else {
                Vector3f in, out;
                in  = muon.LineIn().direction.head(3);
                out = muon.LineOut().direction.head(3);
                float a = in.transpose() * out;
                a = fabs( acos(a / (in.norm() * out.norm())) );
                if(uLib::isFinite(a)) tmp.weight = pow(a * muon.GetMomentum(),2) * 1.5E-6;
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
