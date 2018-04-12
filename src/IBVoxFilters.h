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



#ifndef IB_VOXFILTERS_H
#define IB_VOXFILTERS_H

#include <Core/Object.h>
#include <Core/StaticInterface.h>
#include <Math/VoxImageFilter.h>

#include <IB.h>
#include "IBVoxel.h"

#include <IBVoxImageFilterPlasmon.hpp>
#include <IBVoxImageFilterGradient.hpp>




class IBFilterGaussShape : public uLib::Interface::VoxImageFilterShape {
public:
    IBFilterGaussShape(float sigma) : m_sigma2(sigma*sigma) {}

    float operator()(float d) {
        return exp( -d/m_sigma2);
    }

    float operator()(const Vector3f &pos) {
        return compute_gaussian(pos(0)-0.5, pos(0)+0.5, m_sigma2) *
                compute_gaussian(pos(1)-0.5, pos(1)+0.5, m_sigma2) *
                compute_gaussian(pos(2)-0.5, pos(2)+0.5, m_sigma2);
    }

private:        

    static float compute_gaussian(float x0, float x1, float sigma2, int divs =5)
    {
        float sum = 0;
        for (int i=0; i<divs; ++i) {
            float x = x0 + (x1-x0)/divs * (i+0.5);
            sum += 1/(sqrt(2*M_PI * sigma2)) * exp(-0.5*x*x/sigma2);
        }
        return sum;
    }

    float m_sigma2;
};





typedef uLib::VoxFilterAlgorithmLinear<IBVoxel> IBVoxFilter_Linear;

typedef uLib::VoxFilterAlgorithmMedian<IBVoxel> IBVoxFilter_Median;

typedef uLib::VoxFilterAlgorithmAbtrim<IBVoxel> IBVoxFilter_Abtrim;

typedef uLib::VoxFilterAlgorithmBilateral<IBVoxel> IBVoxFilter_Bilateral;

typedef uLib::VoxFilterAlgorithmBilateralTrim<IBVoxel> IBVoxFilter_BilateralTrim;

typedef uLib::VoxFilterAlgorithmSPR<IBVoxel> IBVoxFilter_SPR;

typedef VoxImageFilterAlgorithmPlasmon<IBVoxel> IBVoxFilter_Plasmon;

typedef VoxImageFilterAlgorithmGradient<IBVoxel> IBVoxFilter_Gradient;

typedef uLib::VoxFilterAlgorithm2ndStat<IBVoxel> IBVoxFilter_2ndStat;


#endif // IBVOXFILTERS_H
