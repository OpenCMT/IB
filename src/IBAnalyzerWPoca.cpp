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
        m_tree->Branch("w_x", &t_w_1,      "w_x");
        m_tree->Branch("w_z", &t_w_2,      "w_z");
        m_tree->Branch("w",   &tmp.weight, "w");
#endif
    }

    bool CollectMuon(const MuonScatterData &muon)
    {
        assert(m_PocaAlgorithm);
        assert(m_Minimizator);
        if (m_PocaAlgorithm->evaluate(muon) && m_Minimizator->evaluate(muon)) {
            tmp.poca   = m_PocaAlgorithm->getPoca();
            t_w_1 = tan((m_Minimizator->getDataVector(0)));
            t_w_2 = tan((m_Minimizator->getDataVector(2)));
            t_w_1 *= t_w_1;
            t_w_2 *= t_w_2;
//            t_w_1 = (t_w_1>0.5) ? 0.5 : t_w_1; // <<< HARDCODED!!
//            t_w_2 = (t_w_2>0.5) ? 0.5 : t_w_2; // Think about delta angle!

            tmp.weight = 1E-6*(sqrt(t_w_1 + t_w_2));

            m_Data.push_back(tmp);

#ifndef NDEBUG
            m_tree->Fill();
#endif
            return true;
        }
        return false;
    }

    void SetVoxels(IBVoxCollection *voxels) {
        for (int i=0; i<m_Data.size(); ++i) {
            Vector3i id = voxels->Find(m_Data[i].poca);
            if (voxels->IsInsideGrid(id))
                voxels->operator [](id).Value += m_Data[i].weight;
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

void IBAnalyzerWPoca::SetVariablesAlgorithm(IBMinimizationVariablesEvaluator *evaluator)
{
    d->m_Minimizator = evaluator;
}
