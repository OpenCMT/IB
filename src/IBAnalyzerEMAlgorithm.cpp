#include "IBAnalyzerEMAlgorithm.h"

// Algorithms Headers //

#include "IBAnalyzerEMAlgorithmSGA.h"
#include "IBAnalyzerEMAlgorithmMGA.h"



bool IBAnalyzerEMAlgorithm::ComputeSigma(Matrix4f &Sigma,
                                         IBAnalyzerEMAlgorithm::Event *evc)
{
    Matrix2f _Sigma = Matrix2f::Zero();
//    float RLen = 0;
    for (unsigned int j = 0; j < evc->elements.size(); ++j) {
        evc->elements[j].lambda = evc->elements[j].voxel->Value;
//        float L = evc->elements[j].Wij(0,0);
//        float T = (Wij(0,1) - L*L/2)/L;
        _Sigma += evc->elements[j].Wij * evc->elements[j].lambda;

    }

    _Sigma *= evc->header.InitialSqrP;

    Sigma.block<2,2>(0,0) = _Sigma;
    Sigma.block<2,2>(2,2) = _Sigma;
    Sigma += evc->header.E;

    return true;
}
