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
