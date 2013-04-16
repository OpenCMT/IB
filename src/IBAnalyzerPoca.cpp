
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
    Vector<HPoint3f> m_Data;
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



