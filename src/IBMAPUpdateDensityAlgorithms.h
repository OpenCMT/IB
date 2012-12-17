#ifndef IBMAPUPDATEDENSITYALGORITHMS_H
#define IBMAPUPDATEDENSITYALGORITHMS_H

#include "IBVoxCollection.h"

using namespace uLib;

////////////////////////////////////////////////////////////////////////////////
/////// GAUSSIAN MAP UPDATE  ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class IBMAPPriorGaussianUpdateAlgorithm :
        public IBAbstract::IBVoxCollectionMAPAlgorithm {
public:
    IBMAPPriorGaussianUpdateAlgorithm(Scalarf beta) :
        m_beta(beta), m_DensityPrior(NULL) {}

    uLibSetMacro(DensityPrior, IBVoxCollection *)

    void UpdateDensity(IBVoxCollection *voxels, unsigned int threshold);

    float GetDensity(IBVoxel &voxel);

    float GetDensity(IBVoxel &voxel, IBVoxel &prior_voxel);

private:
    int CheckVoxelsCorrectness(IBVoxCollection *voxels);

    Scalarf m_beta;
    IBVoxCollection *m_DensityPrior;
};


inline void
IBMAPPriorGaussianUpdateAlgorithm::UpdateDensity(IBVoxCollection *voxels,
                                                 unsigned int threshold)
{
    // TODO : use of threshold !
    if(CheckVoxelsCorrectness(voxels)) {
        for(unsigned int i=0; i< voxels->Data().size(); ++i) {
            IBVoxel &voxel = voxels->Data()[i];
            IBVoxel &prior = m_DensityPrior->Data()[i]; // warning !
            float tau = static_cast<float>(voxel.Count) / m_beta;
            if( tau > prior.Value * prior.Value / 3 )
                voxel.Value = GetDensity(voxel,prior);
            else
                voxel.Value = GetDensity(voxel);
        }
    }
    else {
        for(unsigned int i=0; i< voxels->Data().size(); ++i) {
            IBVoxel &voxel = voxels->Data()[i];
            voxel.Value = GetDensity(voxel);
        }
    }
}

#include <Math/Utils.h>

inline float
IBMAPPriorGaussianUpdateAlgorithm::GetDensity(IBVoxel &voxel)
{
    // tau = 1/betaj = sigma0 (lambda distribution) //
    float tau = voxel.Count / m_beta;
    // t assumes that Analyzer has already produced his value //
    float t = voxel.Value;
    float a = tau * t / 2;
    float b = tau * sqrt( (t * t / 4) + ( tau / 27 ) );
    float a_b = a-b;
    float ret = (pow(a+b, 1./3) + fastSign(a_b) * pow(fabs(a_b), 1./3));
    return ret;
}

inline float
IBMAPPriorGaussianUpdateAlgorithm::GetDensity(IBVoxel &voxel, IBVoxel &prior)
{
    // tau = 1/betaj = sigma0 (lambda distribution) //
    float tau = static_cast<float>(voxel.Count) / m_beta;
    float lambda0 = prior.Value;
    float p = tau - lambda0 * lambda0 / 3;
    float q = lambda0 * tau * tau / 3 - 2 * lambda0 * lambda0 * lambda0 / 27 -
            tau * tau * voxel.Value;
    float D = sqrt(q*q/4 + p*p*p/27);

    return lambda0/3 + fastSign(-q/2+D) * pow(fabs(-q/2+D), 1./3) + fastSign(-q/2-D) * pow(fabs(-q/2-D), 1./3);
}

inline int
IBMAPPriorGaussianUpdateAlgorithm::CheckVoxelsCorrectness(IBVoxCollection *voxels)
{
    // WARNING: should also check that prior has the same axes order //
    if(m_DensityPrior)
        return voxels->GetDims() == m_DensityPrior->GetDims();
    else
        return false;
}










////////////////////////////////////////////////////////////////////////////////
/////// LAPLACIAN MAP UPDATE  //////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class IBMAPPriorLaplacianUpdateAlgorithm :
        public IBAbstract::IBVoxCollectionMAPAlgorithm {
public:
    IBMAPPriorLaplacianUpdateAlgorithm(Scalarf beta) :
        m_beta(beta) {}

    void UpdateDensity(IBVoxCollection *voxels, unsigned int threshold);

private:
    Scalarf m_beta;
};

inline void
IBMAPPriorLaplacianUpdateAlgorithm::UpdateDensity(IBVoxCollection *voxels,
                                                unsigned int threshold)
{
    for(unsigned int i=0; i< voxels->Data().size(); ++i) {
        IBVoxel &voxel = voxels->Data()[i];
        float tau = static_cast<float>(voxel.Count) / m_beta; //  1/bj
        float t = voxel.Value;
        voxel.Value = (sqrt(1+4*t*tau)-1)/(2*tau);
    }
}





#endif // IBMAPUPDATEDENSITYALGORITHMS_H
