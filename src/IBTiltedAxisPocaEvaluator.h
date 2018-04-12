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



#ifndef IBTILTEDAXISPOCAEVALUATOR_H
#define IBTILTEDAXISPOCAEVALUATOR_H

#include "IBPocaEvaluator.h"
#include "IBMuonError.h"

class IBTiltedAxisPocaEvaluator : public IBPocaEvaluator
{
public:
    IBTiltedAxisPocaEvaluator();
    ~IBTiltedAxisPocaEvaluator();

    bool     evaluate(MuonScatterData muon);
    Vector4f getPoca();
    Vector4f getInTrackPoca(){ return HPoint3f();}
    Vector4f getOutTrackPoca(){ return HPoint3f();}

    // dummy functions for the moment... to be implemented
    inline Scalarf getDistance() {return 0;}

private:
    friend class IBTiltedAxisPocaEvaluatorPimpl;
    class IBTiltedAxisPocaEvaluatorPimpl *d;
};

#endif // IBTILTEDAXISPOCAEVALUATOR_H
