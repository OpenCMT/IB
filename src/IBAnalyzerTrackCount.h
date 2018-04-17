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



#ifndef IBANALYZERTRACKCOUNT_H
#define IBANALYZERTRACKCOUNT_H


#include "IBAnalyzer.h"
#include "IBVoxRaytracer.h"
#include "IBPocaEvaluator.h"

using namespace uLib;


class IBAnalyzerTrackCount : public IBAnalyzer {    
    uLibTypeMacro(IBAnalyzerTrackCount,IBAnalyzer)
public:

    IBAnalyzerTrackCount();
    ~IBAnalyzerTrackCount();

    bool AddMuon(const MuonScatterData &muon);

    void Run(unsigned int iterations = 1, float muons_ratio = 1);

    void SetRayAlgorithm(IBVoxRaytracer *raytracer);

    void SetPocaAlgorithm(IBPocaEvaluator *evaluator);

    void SetDetectorZSelection(int selectZ);

    void Clear();

    unsigned int Size() const;

    void SetMuonCollection(IBMuonCollection *muons);


private:
    class IBAnalyzerTrackCountPimpl *d;
};



#endif // IBANALYZERTRACKCOUNT_H
