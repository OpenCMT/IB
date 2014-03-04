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
        _Sigma += evc->elements[j].Wij * evc->elements[j].lambda;
    }

    _Sigma *= evc->header.InitialSqrP;

    Sigma.block<2,2>(0,0) = _Sigma;
    Sigma.block<2,2>(2,2) = _Sigma;
    Sigma += evc->header.E;

    return true;
}
