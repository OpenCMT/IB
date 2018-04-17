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



#ifndef IBMAPUPDATEDENSITYALGORITHMS_H
#define IBMAPUPDATEDENSITYALGORITHMS_H

#include <math.h>

#include "IBVoxCollection.h"
#include "IBVoxFilters.h"
#include <Math/Utils.h>

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




////////////////////////////////////////////////////////////////////////////////
////// GIANNI TOTAL WEIGHT   ///////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class IBMAPPriorTotalWeigth :
        public IBAbstract::IBVoxCollectionMAPAlgorithm {

public:
    IBMAPPriorTotalWeigth(Scalarf mean_weight, Scalarf var_weight) :
        m_weigth(mean_weight), m_variance(var_weight) {}

    void UpdateDensity(IBVoxCollection *voxels, unsigned int threshold);

protected:

    Scalarf ComputeFactorA(IBVoxCollection *voxels);

private:
    Scalarf m_weigth, m_variance;
};


inline void
IBMAPPriorTotalWeigth::UpdateDensity(IBVoxCollection *voxels,
                                     unsigned int threshold)
{
    Scalarf a = this->ComputeFactorA(voxels);
    std::cout << "MAP GW .. [a= "<< a << "] \n";
    for (int i=0; i< voxels->Data().size(); ++i)
    {
        IBVoxel &vox = voxels->Data()[i];
        vox.Value -= a * vox.Value * vox.Value / vox.Count;
    }
}

inline Scalarf
IBMAPPriorTotalWeigth::ComputeFactorA(IBVoxCollection *voxels)
{
    Scalarf sum = 0;
    for(int i=0; i< voxels->Data().size(); ++i)
        sum += voxels->At(i).Value;
    std::cout << " [lmean=" << sum/voxels->Data().size() << "] ";
    return (sum - m_weigth)/m_variance;
}





////////////////////////////////////////////////////////////////////////////////
//////  GIANNI NEIGHBOUR DENSITY GAUSSIAN PRIOR  ///////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class IBMAPPriorNeighbourDensity :
        public IBAbstract::IBVoxCollectionMAPAlgorithm
{
public:
    IBMAPPriorNeighbourDensity(const IBVoxCollection &image, float sigma) :
        m_NeighMap(image), m_Sigma(sigma) {}

    uLibSetMacro(Filter,Abstract::VoxImageFilter *)

    void UpdateDensity(IBVoxCollection *voxels, unsigned int threshold);


protected:
    Abstract::VoxImageFilter *m_Filter;
    IBVoxCollection           m_NeighMap;
    Scalarf                   m_Sigma;
};


void IBMAPPriorNeighbourDensity::UpdateDensity(IBVoxCollection *voxels,
                                               unsigned int threshold)
{
    //    std::cout << "\n == UPDATE DENSITY == \n";
    int rcount = 0;
    for (int i=0; i < voxels->Data().size(); ++i)
    {
        IBVoxel &vox = voxels->Data()[i];
        const IBVoxel &vox_int = this->m_NeighMap.At(i);
        float Mj = vox.Count;

        // 1) check sigma
        if(m_Sigma > vox_int.Value / 5)
        {
            // 2) update lambda
            float a = vox_int.Value;
            float b = Mj * m_Sigma * m_Sigma;
            float c = - b * vox.Value;

            float p = -a*a/3 + b;
            float q = -2*a*a*a/27 - a*b/3 + c;
            float D = sqrt(q*q/4 + p*p*p/27);

            vox.Value = a/3 + cbrtf( -q/2+D ) + cbrtf( -q/2-D );
            rcount++;
        }
    }

    // debug //
    std::cout << " MAP active: " << rcount * 100 / voxels->Data().size() << "%\n";


    // 3) save to local image map
    this->m_NeighMap = *voxels;
    this->m_Filter->SetImage(&m_NeighMap);
    this->m_Filter->Run();

}




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Da "Formalismo per ricostruzione tomografica pg. 23,24 (10 febbraio 2013)  //


class IBMAPPriorNeighbourDensity2 :
        public IBAbstract::IBVoxCollectionMAPAlgorithm
{
public:
    IBMAPPriorNeighbourDensity2(const IBVoxCollection &image, Scalarf a) :
        m_NeighMap(image) , m_a(a) {}

    uLibSetMacro(Filter,Abstract::VoxImageFilter *)

    void UpdateDensity(IBVoxCollection *voxels, unsigned int threshold);


protected:
    Abstract::VoxImageFilter *m_Filter;
    IBVoxCollection           m_NeighMap;
    Scalarf                   m_a;
};


void IBMAPPriorNeighbourDensity2::UpdateDensity(IBVoxCollection *voxels,
                                               unsigned int threshold)
{
    for (int i=0; i < voxels->Data().size(); ++i)
    {
        IBVoxel       &vox = voxels->Data()[i];
        const IBVoxel &vox_int = this->m_NeighMap.At(i);
        float Mj = vox.Count;

        // 1) do computation //
        float p = Mj * vox.Value - vox_int.Value / (m_a * m_a);
        float q = 2*(Mj + 1);

        vox.Value = 1/q * (p + sqrt( p*p + 2*q * pow(vox_int.Value/m_a,2) ));
    }


    // 2) save to local image map
    this->m_NeighMap = *voxels;
    this->m_Filter->SetImage(&m_NeighMap);
    this->m_Filter->Run();
}






#endif // IBMAPUPDATEDENSITYALGORITHMS_H




