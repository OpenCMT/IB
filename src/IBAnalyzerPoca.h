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



#ifndef IBANALYZERPOCA_H
#define IBANALYZERPOCA_H

#include "IBAnalyzer.h"

using namespace uLib;

class IBPocaEvaluator;

class IBAnalyzerPoca : public IBAnalyzer {
    uLibTypeMacro(IBAnalyzerPoca,IBAnalyzer)
public:
    IBAnalyzerPoca();
    ~IBAnalyzerPoca();

    bool AddMuon(const MuonScatterData &event);

    void Run(unsigned int iterations = 1, float muons_ratio = 1);

    void SetPocaAlgorithm(IBPocaEvaluator *poca);

    void SetMuonCollection(IBMuonCollection *muons);

private:
    IBPocaEvaluator *m_PocaAlgorithm;
    Vector<HPoint3f> m_Data;
};




#endif // IBANALYZERPOCA_H
