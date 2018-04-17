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



#ifndef IBANALYZEREMALGORITHMSGA_H
#define IBANALYZEREMALGORITHMSGA_H

#include "IBAnalyzerEMAlgorithm.h"


typedef IBAnalyzerEMAlgorithm IBAnalyzerEMAlgorithmSGA;


class IBAnalyzerEMAlgorithmSGA_PXTZ : public IBAnalyzerEMAlgorithmSGA {
public:
    void evaluate(Matrix4f &Sigma, Event *evc);
};

class IBAnalyzerEMAlgorithmSGA_PXTZ2 : public IBAnalyzerEMAlgorithmSGA {
public:
    void evaluate(Matrix4f &Sigma, Event *evc);
};

class IBAnalyzerEMAlgorithmSGA_PXTZ3 : public IBAnalyzerEMAlgorithmSGA_PXTZ {
    typedef IBAnalyzerEMAlgorithmSGA_PXTZ BaseClass;
public:
    IBAnalyzerEMAlgorithmSGA_PXTZ3() : m_Factor(0) {}

    uLibGetSetMacro(Factor,Scalarf)

    bool ComputeSigma(Matrix4f &Sigma, Event *evc);
private:
    Scalarf m_Factor;
};

class IBAnalyzerEMAlgorithmSGA_PXTZ4 : public IBAnalyzerEMAlgorithmSGA {
    typedef IBAnalyzerEMAlgorithmSGA_PXTZ BaseClass;
public:
    IBAnalyzerEMAlgorithmSGA_PXTZ4() : m_AR(0,0,0), m_MA(0,0,0) {}

    uLibGetSetMacro(AR,Vector3f)
    uLibGetSetMacro(MA,Vector3f)

    void evaluate(Matrix4f &Sigma, Event *evc);

private:
    Vector3f m_AR;
    Vector3f m_MA;
};


class IBAnalyzerEMAlgorithmSGA_PX : public IBAnalyzerEMAlgorithmSGA {
public:
    void evaluate(Matrix4f &Sigma, Event *evc);
};

class IBAnalyzerEMAlgorithmSGA_PXT : public IBAnalyzerEMAlgorithmSGA {
public:
    void evaluate(Matrix4f &Sigma, Event *evc);
};

class IBAnalyzerEMAlgorithmSGA_TZ : public IBAnalyzerEMAlgorithmSGA {
public:
    void evaluate(Matrix4f &Sigma, Event *evc);
};

class IBAnalyzerEMAlgorithmSGA_PT : public IBAnalyzerEMAlgorithmSGA {
public:
    void evaluate(Matrix4f &Sigma, Event *evc);
};

class IBAnalyzerEMAlgorithmSGA_XZ : public IBAnalyzerEMAlgorithmSGA {
public:
    void evaluate(Matrix4f &Sigma, Event *evc);
};

class IBAnalyzerEMAlgorithmSGA_P : public IBAnalyzerEMAlgorithmSGA {
public:
    void evaluate(Matrix4f &Sigma, Event *evc);
};

class IBAnalyzerEMAlgorithmSGA_T : public IBAnalyzerEMAlgorithmSGA {
public:
    void evaluate(Matrix4f &Sigma, Event *evc);
};

class IBAnalyzerEMAlgorithmSGA_X : public IBAnalyzerEMAlgorithmSGA {
public:
    void evaluate(Matrix4f &Sigma, Event *evc);
};

class IBAnalyzerEMAlgorithmSGA_Z : public IBAnalyzerEMAlgorithmSGA {
public:
    void evaluate(Matrix4f &Sigma, Event *evc);
};




class IBAnalyzerEMAlgorithmSGA_M : public IBAnalyzerEMAlgorithmSGA {
public:
    void evaluate(Matrix4f &Sigma, Event *evc);
};

class IBAnalyzerEMAlgorithmSGA_M_PX : public IBAnalyzerEMAlgorithmSGA {
public:
    void evaluate(Matrix4f &Sigma, Event *evc);
};









#endif // IBANALYZEREMALGORITHMSGA_H
