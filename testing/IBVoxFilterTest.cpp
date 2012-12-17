#include "testing-prototype.h"
#include "IBVoxCollectionCap.h"
#include "../src/IBVoxFilters.h"

int main()
{
    BEGIN_TESTING(IB Vox Filters);

    IBVoxCollectionCap voxels(Vector3i(10,5,10));
    voxels.SetPosition(Vector3f(-5,-3,-5));


    IBVoxel zero = {0.1E-6,0,0};

    voxels.InitLambda(zero);

    voxels[Vector3i(3,3,3)].Value = 5.552368E-6;

    IBVoxFilter_Linear filter(Vector3i(5,5,5));

    IBFilterGaussShape shape(0.2);
    filter.SetKernelSpherical(shape);

    filter.SetImage(&voxels);

    filter.Run();

    voxels.ExportToVtk("filter_tes.vtk",0);


    END_TESTING;
}
