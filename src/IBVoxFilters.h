#ifndef IB_VOXFILTERS_H
#define IB_VOXFILTERS_H

#include <Core/Object.h>
#include <Core/StaticInterface.h>
#include <Math/VoxImageFilter.h>
#include "IBVoxel.h"

#include <IBVoxImageFilterPlasmon.hpp>
#include <IBVoxImageFilterGradient.hpp>




class IBFilterGaussShape : public uLib::Interface::VoxImageFilterShape {
public:
    IBFilterGaussShape(float sigma) : m_sigma(sigma) {}

    float operator()(float d) {
        // no need to normalize if autonormalized filter applied //
        return exp( -d/m_sigma);
    }

private:
    float m_sigma;
};


typedef uLib::VoxFilterAlgorithmLinear<IBVoxel> IBVoxFilter_Linear;

typedef uLib::VoxFilterAlgorithmMedian<IBVoxel> IBVoxFilter_Median;

typedef uLib::VoxFilterAlgorithmAbtrim<IBVoxel> IBVoxFilter_Abtrim;

typedef VoxImageFilterAlgorithmPlasmon<IBVoxel> IBVoxFilter_Plasmon;

typedef VoxImageFilterAlgorithmGradient<IBVoxel> IBVoxFilter_Gradient;

typedef uLib::VoxFilterAlgorithm2ndStat<IBVoxel> IBVoxFilter_2ndStat;


#endif // IBVOXFILTERS_H
