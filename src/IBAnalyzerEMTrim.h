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



#ifndef IBANALYZEREMTRIM_H
#define IBANALYZEREMTRIM_H

#include "IBAnalyzerEM.h"

using namespace uLib;

class IBAnalyzerEMTrim : public IBAnalyzerEM {

    typedef IBAnalyzerEM BaseClass;
public:

    IBAnalyzerEMTrim(IBVoxCollection &voxels);
    ~IBAnalyzerEMTrim();

    void Run(unsigned int iterations, float muons_ratio, float a=0, float b=0);

    void SetMLAlgorithm(IBAnalyzerEMAlgorithm *MLAlgorithm);

private:
    friend class IBAnalyzerEMTrimPimpl;
    class IBAnalyzerEMTrimPimpl *d;
};



#endif // IBANALYZEREMTRIM_H
