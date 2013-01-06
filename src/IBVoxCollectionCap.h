#ifndef IBVOXCOLLECTIONCAP_H
#define IBVOXCOLLECTIONCAP_H

#include <Math/Utils.h>
#include "IBVoxCollection.h"
#include "Math/Utils.h"

using namespace uLib;

class UpdateDensitySijCapAlgorithm :
        public IBInterface::IBVoxCollectionStaticUpdateAlgorithm
{
public:
    static void UpdateDensity(IBVoxCollection *voxels, unsigned int threshold)
    {
        for(unsigned int i=0; i< voxels->Data().size(); ++i) {
            IBVoxel& voxel = voxels->Data()[i];
            unsigned int tcount = voxel.Count;
            if (tcount > 0 && (threshold == 0 || tcount >= threshold)) {
                voxel.Value += voxel.SijCap / static_cast<float>(tcount);
                if(unlikely(!isFinite(voxel.Value) || voxel.Value > 100.E-6)) {  // HARDCODED!!!
                    voxel.Value = 100.E-6;
                } else if (unlikely(voxel.Value < 0.)) voxel.Value = 0.1E-6;
            }
            else
                voxel.Value = 0;
            voxel.SijCap = 0;
        }
    }
};



class IBVoxCollectionCap : public IBVoxCollection
{
    typedef IBVoxCollection BaseClass;
public:
    typedef IBVoxCollection VoxImageType;
    IBVoxCollectionCap(const Vector3i size) :
        BaseClass(size) {}

    IBVoxCollectionCap(const IBVoxCollection &copy) :
        BaseClass(copy) {}

    void UpdateDensity(unsigned int threshold)
    {
        BaseClass::UpdateDensity<UpdateDensitySijCapAlgorithm>(threshold);
    }
};











#endif // IBVOXCOLLECTIONCAP_H
