#ifndef IB_VOXFILTERS_H
#define IB_VOXFILTERS_H

#include <Core/Object.h>
#include <Core/StaticInterface.h>
#include <Math/VoxImageFilter.h>

#include "IBVoxel.h"


// Direct implementations  .. waiting for parameters //

//class IBVoxFilter : public uLib::VoxImageFilter
//{
//    typedef uLib::VoxImageFilter BaseClass;
//public:
//    enum AlgorithmTypeEnum {
//        Linear,
//        Median,
//        Abtrim
//    };

//    ULIB_OBJECT_PARAMETERS(BaseClass)
//    {
//      enum AlgorithmTypeEnum type;
//    };

//public:
//    static IBVoxFilter * New(enum AlgorithmTypeEnum type);


//protected:
//    IBVoxFilter() {}
//    virtual ~IBVoxFilter() {}

//};

//inline IBVoxFilter *IBVoxFilter::New(IBVoxFilter::AlgorithmTypeEnum type)
//{
//    switch(type) {
//    case IBVoxFilter::Linear:
//        return new uLib::VoxFilterAlgorithmLinear<IBVoxel>();
//        break;
//    case IBVoxFilter::Median:
//        break;
//    case IBVoxFilter::Abtrim:
//        break;
//    }
//}

class IBFilterGaussShape : public uLib::Interface::VoxImageFilterShape {
public:
    IBFilterGaussShape(float sigma) :
        m_sigma(sigma) {}

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






#endif // IBVOXFILTERS_H
