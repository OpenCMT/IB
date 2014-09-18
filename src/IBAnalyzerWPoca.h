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



#ifndef IBANALYZERWPOCA_H
#define IBANALYZERWPOCA_H

#include "IBAnalyzer.h"

using namespace uLib;

class IBPocaEvaluator;
class IBMinimizationVariablesEvaluator;

class IBAnalyzerWPoca : public IBAnalyzer
{
    uLibTypeMacro(IBAnalyzerWPoca,IBAnalyzer)
public:
    IBAnalyzerWPoca();
    ~IBAnalyzerWPoca();

    bool AddMuon(const MuonScatterData &event);

    void Run(unsigned int iteration = 1, float muons_ratio = 1);

    void SetPocaAlgorithm(IBPocaEvaluator *poca);

    void SetVarAlgorithm(IBMinimizationVariablesEvaluator *evaluator);

private:
    friend class IBAnalyzerWPocaPimpl;
    class IBAnalyzerWPocaPimpl *d;
};

#endif // IBANALYZERWPOCA_H
