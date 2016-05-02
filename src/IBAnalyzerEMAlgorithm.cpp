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



#include "IBAnalyzerEMAlgorithm.h"

// Algorithms Headers //

#include "IBAnalyzerEMAlgorithmSGA.h"
#include "IBAnalyzerEMAlgorithmMGA.h"



bool IBAnalyzerEMAlgorithm::ComputeSigma(Matrix4f &Sigma,
                                         IBAnalyzerEMAlgorithm::Event *evc)
{
  //  std::cout << "-------------------------" << std::endl;
  Matrix2f _Sigma = Matrix2f::Zero();
    for (unsigned int j = 0; j < evc->elements.size(); ++j) {
	evc->elements[j].lambda = fabs (evc->elements[j].voxel->Value ); // fabs needed to cope with negative (fixed) lambdas //
        //_Sigma += evc->elements[j].Wij * evc->elements[j].lambda * evc->elements[j].pw;
	_Sigma += evc->elements[j].Wij * evc->elements[j].lambda;
    }

    //    if(isnan(evc->header.InitialSqrP)) std::cout << "sto calcolando Sigma e ho trovato 1/p2 a nan \n" << std::flush;
    _Sigma *= evc->header.InitialSqrP;

    Sigma.block<2,2>(0,0) = _Sigma;
    Sigma.block<2,2>(2,2) = _Sigma;
    Sigma += evc->header.E;
    return true;
}
