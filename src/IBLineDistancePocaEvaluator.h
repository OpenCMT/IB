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



#ifndef IBLINEDISTANCEPOCAEVALUATOR_H
#define IBLINEDISTANCEPOCAEVALUATOR_H

#include "IBPocaEvaluator.h"

class IBLineDistancePocaEvaluator : public IBPocaEvaluator
{
public:
    IBLineDistancePocaEvaluator();
    ~IBLineDistancePocaEvaluator();

    bool evaluate(MuonScatterData muon);
    HPoint3f getPoca();
    HPoint3f getOutTrackPoca();
    HPoint3f getInTrackPoca();
    void setDistanceCut(Scalarf length);
private:
    friend class IBLineDistancePocaEvaluatorPimpl;
    class IBLineDistancePocaEvaluatorPimpl *d;
};

#endif // IBLINEDISTANCEPOCAEVALUATOR_H
