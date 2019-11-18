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

#ifndef IBANALYZERTF_H
#define IBANALYZERTF_H

#include <Math/Dense.h>
#include <Math/StructuredData.h>
#include "IBAnalyzer.h"
#include "IBMinimizationVariablesEvaluator.h"
#include "IBMuonCollection.h"
#include "IBPocaEvaluator.h"
#include "IBVoxCollection.h"
#include "IBVoxRaytracer.h"

#include <tensorflow/cc/ops/state_ops.h>

using uLib::Scalarf;
using uLib::StructuredGrid;
using tensorflow::Scope;
using tensorflow::ops::Variable;


class IBAnalyzerTF : public IBAnalyzer
{
public:
    typedef IBAnalyzer BaseClass;

    struct Event {
        struct Element {
            Matrix2f Wij;
            int v_idx;
            Scalarf pw;
        };

        Vector4f Di;
        Matrix4f E;
        Scalarf  InitialSqrP;
        std::vector<Element> elements;
    };

    IBAnalyzerTF(IBVoxCollection &voxels,
                 const Scope &scope,
                 IBPocaEvaluator *poca_algo,
                 IBMinimizationVariablesEvaluator *var_algo,
                 IBVoxRaytracer *ray_algo,
                 float learn_rate = 0.01f);

    virtual ~IBAnalyzerTF();

    inline virtual const char *type_name() const { return "IBAnalyzerTF"; }

    // deprecated
    bool AddMuon(const MuonScatterData &muon) { return false; }

    void SetMuonCollection(IBMuonCollection *muons);

    bool AddMuonFullPath(const MuonScatterData &muon, std::vector<HPoint3f>& muonPath);

    void Run(unsigned int iterations, float muons_ratio);

    unsigned int Size();

private:

    Scalarf                            nominal_momentum;
    IBPocaEvaluator                   *m_PocaAlgorithm;
    IBMinimizationVariablesEvaluator  *m_VarAlgorithm;
    IBVoxRaytracer                    *m_RayAlgorithm;

    Scope                              tf_scope;
    std::unordered_map<int, Variable*> tf_Variables;
    float                              learning_rate;
    std::vector<IBAnalyzerTF::Event>   m_Events;
};


#endif // IBANALYZERTF_H
