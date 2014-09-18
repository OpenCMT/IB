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
    properties() {
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
