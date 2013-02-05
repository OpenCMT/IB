#ifndef IBVOXEL_H
#define IBVOXEL_H

#include <Math/Dense.h>
#include <Math/VoxImage.h>
#include <Core/Macros.h>



struct IBVoxel {    
    uLib::Scalarf Value;
    uLib::Scalarf SijCap;
    unsigned int  Count;
};

struct IBPMap {
    uLib::Scalarf Value;
    uLib::Scalarf Cj;
};


typedef uLib::VoxImage<uLib::Voxel> IBLightCollection;
typedef uLib::VoxImage<IBPMap> IBShadeCollection;
#endif // IBVOXEL_H
