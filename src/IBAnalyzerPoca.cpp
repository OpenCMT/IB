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

