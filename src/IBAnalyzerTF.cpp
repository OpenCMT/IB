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
#include "tensorflow/cc/framework/gradients.h"
#include "tensorflow/core/framework/tensor.h"

using namespace tensorflow;
using namespace tensorflow::ops;


// ****************************************************************************
// Plugin
// ****************************************************************************

IBAnalyzerTF::IBAnalyzerTF(IBVoxCollection &voxels,
                           const Scope &scope,
                           IBPocaEvaluator *poca_algo,
                           IBMinimizationVariablesEvaluator *var_algo,
                           IBVoxRaytracer *ray_algo,
                           float learn_rate) :
    tf_scope(scope),
    m_Events(),
    m_PocaAlgorithm(poca_algo),
    m_VarAlgorithm(var_algo),
    m_RayAlgorithm(ray_algo),
    learning_rate(learn_rate),
    nominal_momentum(3)
{
    SetVoxCollection(&voxels);
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

    if(evn.Di[0] == NAN)
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
{
    auto init_density = Const(tf_scope, { 7.f * 1.E-8f }); // air density

    int m2a[4][4] =
    {
        { 0, 1, 2, 3 },
        { 1, 4, 5, 6 },
        { 2, 5, 7, 8 },
        { 3, 6, 8, 9 }
    };

    auto ident_mat = OneHot(tf_scope, { 0, 1, 2, 3 }, 4, 1.f, 0.f);

    Output oh_sym[10];
    oh_sym[m2a[0][0]] = OneHot(tf_scope, { 0, -1, -1, -1 }, 4, 1.f, 0.f);
    oh_sym[m2a[1][1]] = OneHot(tf_scope, { -1, 1, -1, -1 }, 4, 1.f, 0.f);
    oh_sym[m2a[2][2]] = OneHot(tf_scope, { -1, -1, 2, -1 }, 4, 1.f, 0.f);
    oh_sym[m2a[3][3]] = OneHot(tf_scope, { -1, -1, -1, 3 }, 4, 1.f, 0.f);
    oh_sym[m2a[0][1]] = OneHot(tf_scope, { 1, 0, -1, -1 }, 4, 1.f, 0.f);
    oh_sym[m2a[0][2]] = OneHot(tf_scope, { 2, -1, 0, -1 }, 4, 1.f, 0.f);
    oh_sym[m2a[0][3]] = OneHot(tf_scope, { 3, -1, -1, 0 }, 4, 1.f, 0.f);
    oh_sym[m2a[1][2]] = OneHot(tf_scope, { -1, 2, 1, -1 }, 4, 1.f, 0.f);
    oh_sym[m2a[1][3]] = OneHot(tf_scope, { -1, 3, -1, 1 }, 4, 1.f, 0.f);
    oh_sym[m2a[2][3]] = OneHot(tf_scope, { -1, -1, 3, 2 }, 4, 1.f, 0.f);

    Output oh_det[24] =
    {
        // positive terms
        OneHot(tf_scope,{ 0, 1, 2, 3 }, 4, 1.f, 0.f),
        OneHot(tf_scope,{ 0, 2, 3, 1 }, 4, 1.f, 0.f),
        OneHot(tf_scope,{ 0, 3, 1, 2 }, 4, 1.f, 0.f),
        OneHot(tf_scope,{ 1, 0, 3, 2 }, 4, 1.f, 0.f),
        OneHot(tf_scope,{ 1, 2, 0, 3 }, 4, 1.f, 0.f),
        OneHot(tf_scope,{ 1, 3, 2, 0 }, 4, 1.f, 0.f),
        OneHot(tf_scope,{ 2, 0, 1, 3 }, 4, 1.f, 0.f),
        OneHot(tf_scope,{ 2, 1, 3, 0 }, 4, 1.f, 0.f),
        OneHot(tf_scope,{ 2, 3, 0, 1 }, 4, 1.f, 0.f),
        OneHot(tf_scope,{ 3, 0, 2, 1 }, 4, 1.f, 0.f),
        OneHot(tf_scope,{ 3, 1, 0, 2 }, 4, 1.f, 0.f),
        OneHot(tf_scope,{ 3, 2, 1, 0 }, 4, 1.f, 0.f),
        //negative terms
        OneHot(tf_scope,{ 0, 1, 3, 2 }, 4, 1.f, 0.f),
        OneHot(tf_scope,{ 0, 2, 1, 3 }, 4, 1.f, 0.f),
        OneHot(tf_scope,{ 0, 3, 2, 1 }, 4, 1.f, 0.f),
        OneHot(tf_scope,{ 1, 0, 2, 3 }, 4, 1.f, 0.f),
        OneHot(tf_scope,{ 1, 2, 3, 0 }, 4, 1.f, 0.f),
        OneHot(tf_scope,{ 1, 3, 0, 2 }, 4, 1.f, 0.f),
        OneHot(tf_scope,{ 2, 0, 3, 1 }, 4, 1.f, 0.f),
        OneHot(tf_scope,{ 2, 1, 0, 3 }, 4, 1.f, 0.f),
        OneHot(tf_scope,{ 2, 3, 1, 0 }, 4, 1.f, 0.f),
        OneHot(tf_scope,{ 3, 0, 1, 2 }, 4, 1.f, 0.f),
        OneHot(tf_scope,{ 3, 1, 2, 0 }, 4, 1.f, 0.f),
        OneHot(tf_scope,{ 3, 2, 0, 1 }, 4, 1.f, 0.f)
    };

    Output oh_wij[3] =
    {
        OneHot(tf_scope,{ 0, -1, 2, -1 }, 4, 1.f, 0.f),
        OneHot(tf_scope,{ 1, 0, 3, 2 }, 4, 1.f, 0.f),
        OneHot(tf_scope,{ -1, 1, -1, 3 }, 4, 1.f, 0.f)
    };

    // ************************************************************************
    // Likelihood function composition
    // ************************************************************************

    std::vector<Output> likeli_parts;

    for(IBAnalyzerTF::Event& evn_item : m_Events)
    {
        auto d_input = Const(tf_scope,
        {
            { evn_item.Di(0) },
            { evn_item.Di(1) },
            { evn_item.Di(2) },
            { evn_item.Di(3) }
        });

        std::vector<Output> sigma_parts;
        for(IBAnalyzerTF::Event::Element& el_item : evn_item.elements)
        {

            std::vector<Output> wij_parts;
            wij_parts.push_back(Multiply(tf_scope, oh_wij[0], el_item.Wij(0,0)));
            wij_parts.push_back(Multiply(tf_scope, oh_wij[1], el_item.Wij(0,1)));
            wij_parts.push_back(Multiply(tf_scope, oh_wij[2], el_item.Wij(1,1)));
            auto all_wij = AddN(tf_scope, wij_parts);

            auto tmp_term = Multiply(tf_scope, all_wij, evn_item.InitialSqrP);

            Variable* tmpvar = nullptr;
            auto tmpitem = tf_Variables.find(el_item.v_idx);
            if(tmpitem == tf_Variables.end())
            {
                tmpvar = new Variable(tf_scope, {1}, DT_FLOAT);
                tf_Variables[el_item.v_idx] = tmpvar;
            }
            else
            {
                tmpvar = tmpitem->second;
            }

            auto s_term = Multiply(tf_scope, tmp_term, *tmpvar);

            sigma_parts.push_back(s_term);
        }
        auto tmp_sigma = AddN(tf_scope, sigma_parts);

        std::vector<Output> err_parts;
        for(int k = 0; k < 4; k++)
        {
            for(int j = 0; j < k; j++)
            {
                err_parts.push_back(Multiply(tf_scope,
                                            oh_sym[m2a[k][j]],
                                            evn_item.E(k,j)));
            }
        }
        auto err_mat = AddN(tf_scope, err_parts);

        auto s_tensor = Add(tf_scope, tmp_sigma, err_mat);

        // ************************************************************************
        // Determinant
        // ************************************************************************

        std::vector<Output> det_parts;
        for(int k = 0; k < 24; k++)
        {
            auto tmp_mul = Multiply(tf_scope, oh_det[k], s_tensor);
            auto tmp_term = Prod(tf_scope, Sum(tf_scope, tmp_mul, 0), 0);
            if(k < 12)
            {
                det_parts.push_back(tmp_term);
            }
            else
            {
                det_parts.push_back(Negate(tf_scope, tmp_term));
            }

        }
        auto s_deter = AddN(tf_scope, det_parts);

        // ************************************************************************
        // matrix inverse
        // https://en.wikipedia.org/wiki/Invertible_matrix#Inversion_of_4_%C3%97_4_matrices
        // ************************************************************************

        auto s_tensor2 = MatMul(tf_scope, s_tensor, s_tensor);
        auto s_tensor3 = MatMul(tf_scope, s_tensor2, s_tensor);

        auto s_trace = Sum(tf_scope, Sum(tf_scope, Multiply(tf_scope, s_tensor, ident_mat), 0), 0);
        auto s_trace2 = Multiply(tf_scope, s_trace, s_trace);
        auto s_trace3 = Multiply(tf_scope, s_trace, s_trace2);
        auto s2_trace = Sum(tf_scope, Sum(tf_scope, Multiply(tf_scope, s_tensor2, ident_mat), 0), 0);
        auto s3_trace = Sum(tf_scope, Sum(tf_scope, Multiply(tf_scope, s_tensor3, ident_mat), 0), 0);

        std::vector<Output> c0_parts;
        c0_parts.push_back(s_trace3);
        c0_parts.push_back(Multiply(tf_scope,
                            Const(tf_scope, -3.f), Multiply(tf_scope, s_trace, s2_trace)));
        c0_parts.push_back(Multiply(tf_scope, Const(tf_scope, 2.f), s3_trace));

        std::vector<Output> c1_parts;
        c1_parts.push_back(Multiply(tf_scope, Const(tf_scope, 3.f), s_trace2));
        c1_parts.push_back(Multiply(tf_scope, Const(tf_scope, -3.f), s2_trace));

        auto c2 = Multiply(tf_scope, Const(tf_scope, 6.f), s_trace);
        auto c3 = Const(tf_scope, -6.f);

        std::vector<Output> inv_parts;
        inv_parts.push_back(Multiply(tf_scope, ident_mat, AddN(tf_scope, c0_parts)));
        inv_parts.push_back(Multiply(tf_scope, s_tensor, AddN(tf_scope, c1_parts)));
        inv_parts.push_back(Multiply(tf_scope, s_tensor2, c2));
        inv_parts.push_back(Multiply(tf_scope, s_tensor3, c3));

        auto s_inv = Multiply(tf_scope, AddN(tf_scope, inv_parts), Inv(tf_scope, s_deter));
        auto m_comp = MatMul(tf_scope,
                             MatMul(tf_scope, d_input, s_inv, MatMul::TransposeA(true)),
                             d_input);

        likeli_parts.push_back(Add(tf_scope, Log(tf_scope, s_deter), m_comp));
        std::cout << "Built graph module " << likeli_parts.size() << std::endl;
    }

    if(likeli_parts.size() == 0) return;

    auto likeli_fnct = AddN(tf_scope, likeli_parts);

    // ************************************************************************
    // Gradients setup
    // ************************************************************************

    std::vector<Output> v_assigns;
    std::vector<Output> algo_apply;
    int m_cnt = 0;
    for(auto& item : tf_Variables)
    {
        v_assigns.push_back(Assign(tf_scope, *(item.second), init_density));

        std::vector<Output> grad_outputs;
        TF_CHECK_OK(AddSymbolicGradients(tf_scope, { likeli_fnct }, { *(item.second) }, &grad_outputs));
        algo_apply.push_back(ApplyGradientDescent(tf_scope, *(item.second),
                                                  learning_rate, { grad_outputs[0] }));
        std::cout << "Gradient " << m_cnt << "/" << tf_Variables.size() << std::endl;
        m_cnt++;
    }

    // ************************************************************************
    // Minimization step
    // ************************************************************************

    ClientSession session(tf_scope);
    TF_CHECK_OK(session.Run(v_assigns, nullptr));

    std::vector<Tensor> likeli_val;
    for (int i = 0; i < iterations; ++i)
    {
        TF_CHECK_OK(session.Run({ likeli_fnct }, &likeli_val));
        std::cout << "Loss after " << i << " steps " << likeli_val[0].scalar<float>() << std::endl;
        TF_CHECK_OK(session.Run(algo_apply, nullptr));
    }

}

unsigned int IBAnalyzerTF::Size()
{
    return m_Events.size();
}

