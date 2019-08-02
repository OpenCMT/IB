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



#ifndef IBANALYZERWPOCA_H
#define IBANALYZERWPOCA_H

#include "IBAnalyzer.h"

using namespace uLib;

class IBPocaEvaluator;
class IBMinimizationVariablesEvaluator;

class IBAnalyzerWPoca : public IBAnalyzer
{
public:
    IBAnalyzerWPoca();
    ~IBAnalyzerWPoca();

    inline virtual const char *type_name() const { return "IBAnalyzerWPoca"; }

    bool AddMuon(const MuonScatterData &event);

    void Run(unsigned int iteration = 1, float muons_ratio = 1);

    void SetPocaAlgorithm(IBPocaEvaluator *poca);

    void SetVarAlgorithm(IBMinimizationVariablesEvaluator *evaluator);

private:

    struct Data {
        HPoint3f poca;
        Scalarf  weight;
    };

    IBPocaEvaluator                  *m_PocaAlgorithm;
    IBMinimizationVariablesEvaluator *m_Minimizator;
    std::vector<Data>                 m_Data;

    Data    tmp;
    Scalarf t_w_1;
    Scalarf t_w_2;

#ifndef NDEBUG
    TFile* m_out;
    TTree* m_tree;
#endif
};

#endif // IBANALYZERWPOCA_H
