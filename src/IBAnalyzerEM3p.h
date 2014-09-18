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



#ifndef IBANALYZEREM3P_H
#define IBANALYZEREM3P_H

#include "IBAnalyzerEM.h"

class IBAnalyzerEM3p : public IBAnalyzerEM {
    typedef IBAnalyzerEM BaseClass;

public:
    IBAnalyzerEM3p(IBVoxCollection &voxels) : BaseClass(voxels) {}

    void SetPocaAlgorithm(IBPocaEvaluator *PocaAlgorithm);

    virtual bool AddMuon(const MuonScatterData &muon);

};



#endif // IBANALYZEREM3P_H


