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


#ifndef IBANALYZERWTRACKLENGTHS_H
#define IBANALYZERWTRACKLENGTHS_H

#include "IBAnalyzer.h"

#include "IBVoxRaytracer.h"
#include "IBPocaEvaluator.h"

#include "IBMinimizationVariablesEvaluator.h"

using namespace uLib;


class IBAnalyzerWTrackLengths : public IBAnalyzer {
    uLibTypeMacro(IBAnalyzerWTrackLengths,IBAnalyzer)
public:
    IBAnalyzerWTrackLengths();
    ~IBAnalyzerWTrackLengths();

    bool AddMuon(const MuonScatterData &muon);

    void Run(unsigned int iterations = 1, float muons_ratio = 1);

    void SetRayAlgorithm(IBVoxRaytracer *raytracer);

    void SetPocaAlgorithm(IBPocaEvaluator *evaluator);

    void SetVarAlgorithm(IBMinimizationVariablesEvaluator *algorithm);

    void SetPocaProximity(float sigma = 0); // TODO: //

private:
    class IBAnalyzerWTrackLengthsPimpl *d;
};



#endif // IBANALYZERWTRACKLENGTHS_H
