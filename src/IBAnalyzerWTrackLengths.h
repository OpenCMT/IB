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
