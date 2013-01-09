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


typedef uLib::VoxImage<uLib::Voxel> IBLightCollection;

#endif // IBVOXEL_H
