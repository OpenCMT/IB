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

#include "IBAnalyzerTF.h"

#include "tensorflow/cc/client/client_session.h"
#include "tensorflow/cc/ops/standard_ops.h"
#include "tensorflow/core/framework/tensor.h"

using tensorflow::DT_FLOAT;


// ****************************************************************************
// Graph variables
// ****************************************************************************

TFVariableSet::TFVariableSet(const Scope &scope, const Vector3i &size) :
    tf_vars(size.prod()),
    StructuredData(size)
{
    for(int k; k<size.prod(); k++)
    {
        tf_vars[k] = new Variable(scope, {1}, DT_FLOAT);
    }
}

TFVariableSet::~TFVariableSet()
{
    for(int k; k<tf_vars.size(); k++)
    {
        delete tf_vars[k];
    }
}

Variable& TFVariableSet::At(int i) const
{
    return *tf_vars.at(i);
}

Variable& TFVariableSet::At(const Vector3i &id) const
{
    return *tf_vars.at(Map(id));
}




// ****************************************************************************
// Plugin
// ****************************************************************************

IBAnalyzerTF::IBAnalyzerTF(IBVoxCollection &voxels,
                           const Scope &scope,
                           IBPocaEvaluator *poca_algo,
                           IBMinimizationVariablesEvaluator *var_algo,
                           IBVoxRaytracer *ray_algo) :
    tf_scope(scope),
    m_Events(),
    m_PocaAlgorithm(poca_algo),
    m_VarAlgorithm(var_algo),
    m_RayAlgorithm(ray_algo),
    nominal_momentum(3),
    tf_Variables(TFVariableSet(tf_scope, voxels.GetDims()))
{
}

IBAnalyzerTF::~IBAnalyzerTF()
{}

void IBAnalyzerTF::SetMuonCollection(IBMuonCollection *muons)
{
    m_Events.clear();

    std::vector<MuonScatterData>::iterator itr = muons->Data().begin();
    std::vector<std::vector<HPoint3f> >::iterator path_itr = muons->FullPath().begin();

    while(itr != muons->Data().end())
    {
        if(!AddMuonFullPath(*itr, *path_itr))
        {
            std::swap(*itr, muons->Data().back());
            muons->Data().pop_back();

            if(muons->FullPath().size() > 0)
            {
                std::swap(*path_itr, muons->FullPath().back());
                muons->FullPath().pop_back();
            }
        }
        else
        {
            itr++;
            path_itr++;
        }
    }

    BaseClass::SetMuonCollection(muons);
}

bool IBAnalyzerTF::AddMuonFullPath(const MuonScatterData &muon,
                                   std::vector<HPoint3f>& muonPath)
{
    if(!m_VarAlgorithm->evaluate(muon))
    {
        return false;
    }
    if(muon.GetMomentum() == 0)
    {
        return false;
    }

    Event evn;
    evn.Di = m_VarAlgorithm->getDataVector();
    evn.E = m_VarAlgorithm->getCovarianceMatrix();
    evn.InitialSqrP = pow(nominal_momentum/muon.GetMomentum(), 2);
    evn.elements = std::vector<Event::Element>(0);

    if(evn.Di[0] != NAN)
    {
        return false;
    }
    std::vector<HPoint3f> pts;

    HPoint3f entry_pt, exit_pt;
    if(!m_RayAlgorithm->GetEntryPoint(muon.LineIn(), entry_pt))
    {
        return false;
    }
    if(!m_RayAlgorithm->GetExitPoint(muon.LineOut(), exit_pt))
    {
        return false;
    }

    pts.push_back(entry_pt);

    if(m_PocaAlgorithm->evaluate(muon))
    {
        HPoint3f poca  = m_PocaAlgorithm->getPoca();
        HVector3f in   = poca - muon.LineIn().origin;
        HVector3f out  = muon.LineOut().origin - poca;
        float poca_prj = in.transpose() * out;

        if(poca_prj > 0 && GetVoxCollection()->IsInsideBounds(poca))
        {
            pts.push_back(poca);
        }

        pts.push_back(exit_pt);
    }

    std::map<int,std::vector<HPoint3f> > voxelMap;
    std::vector<int> voxelOrder;
    Scalarf totalLength = 0.;

    for(int i = 0; i < pts.size() - 1; ++i)
    {
        HPoint3f& pt1 = pts[i];
        HPoint3f& pt2 = pts[i+1];

        Scalarf  rayLength = (pt2-pt1).norm();
        HPoint3f rayDir    = (pt2-pt1)/rayLength;

        IBVoxRaytracer::RayData ray = m_RayAlgorithm->TraceBetweenPoints(pt1,pt2);

        float cumulativeLength = 0.;
        foreach(const IBVoxRaytracer::RayData::Element &item, ray.Data())
        {
            HPoint3f pti = pt1 + cumulativeLength * rayDir;
            cumulativeLength += item.L;
            HPoint3f ptj = pt1 + cumulativeLength * rayDir;

            if(voxelMap.find(item.vox_id) == voxelMap.end())
            {
                voxelMap[item.vox_id] = std::vector<HPoint3f>();
                voxelMap[item.vox_id].push_back(pti);
                voxelMap[item.vox_id].push_back(ptj);
                voxelOrder.push_back(item.vox_id);
            }
            voxelMap[item.vox_id][1] = ptj;
        }
        totalLength += rayLength;
    }


    if(voxelOrder.size() < 2)
    {
        return false;
    }

    Scalarf normIn = muon.LineIn().direction.norm();
    Scalarf H = muon.LineIn().direction.transpose() * (exit_pt - entry_pt);
    H = H/normIn;

    for(std::vector<int>::const_iterator it = voxelOrder.begin(); it != voxelOrder.end(); it++)
    {
        const HPoint3f& pt1 = voxelMap[*it][0];
        const HPoint3f& pt2 = voxelMap[*it][1];

        Event::Element elc;
        elc.v_idx = *it;

        Scalarf L = (pt2 - pt1).norm();
        Scalarf h = muon.LineIn().direction.transpose() * (pt2 - entry_pt);
        Scalarf T = H - h/normIn;
        if(T < 0) T = 0.;

        elc.Wij << L,
                   L*L/2. + L*T,
                   L*L/2. + L*T,
                   L*L*L/3. + L*L*T + L*T*T;
        elc.pw = evn.InitialSqrP;

        evn.elements.push_back(elc);
    }

    if(evn.elements.size() > 1)
    {
        m_Events.push_back(evn);
        return true;
    }
    return false;
}

void IBAnalyzerTF::Run(unsigned int iterations, float muons_ratio)
{}

unsigned int IBAnalyzerTF::Size()
{
    return m_Events.size();
}
