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


#ifndef IBANALYZEREMALGORITHM_H
#define IBANALYZEREMALGORITHM_H

#include "Core/Object.h"
#include "IBAnalyzerEM.h"

using namespace uLib;



class IBAnalyzerEMAlgorithm : public Object {
protected:
    typedef struct IBAnalyzerEM::Event Event;
    uLibTypeMacro(IBAnalyzerEMAlgorithm, uLib::Object)
public:
    ULIB_props() {
        Scalarf inertia;
    };

    IBAnalyzerEMAlgorithm() { init_properties(); }

    virtual void evaluate(Matrix4f &Sigma, Event *evc) = 0;

    virtual bool ComputeSigma(Matrix4f &Sigma, Event *evc);

protected:
    virtual ~IBAnalyzerEMAlgorithm() {}

};


inline void IBAnalyzerEMAlgorithm::init_properties() {
    $_init();
    $$.inertia = 1;
}

#endif // IBANALYZEREMALGORITHM_H
