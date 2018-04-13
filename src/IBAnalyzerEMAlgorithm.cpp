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


#include "IBAnalyzerEMAlgorithm.h"

// Algorithms Headers //

#include "IBAnalyzerEMAlgorithmSGA.h"
#include "IBAnalyzerEMAlgorithmMGA.h"



bool IBAnalyzerEMAlgorithm::ComputeSigma(Matrix4f &Sigma,
                                         IBAnalyzerEMAlgorithm::Event *evc)
{
  Matrix2f _Sigma = Matrix2f::Zero();
    for (unsigned int j = 0; j < evc->elements.size(); ++j) {
	evc->elements[j].lambda = fabs (evc->elements[j].voxel->Value ); // fabs needed to cope with negative (fixed) lambdas //
        _Sigma += evc->elements[j].Wij * evc->elements[j].lambda * evc->elements[j].pw;	
	//_Sigma += evc->elements[j].Wij * evc->elements[j].lambda;
    }

    //    if(isnan(evc->header.InitialSqrP)) std::cout << "sto calcolando Sigma e ho trovato 1/p2 a nan \n" << std::flush;
    //    _Sigma *= evc->header.InitialSqrP;

    Sigma.block<2,2>(0,0) = _Sigma;
    Sigma.block<2,2>(2,2) = _Sigma;
    Sigma += evc->header.E;
    return true;
}
