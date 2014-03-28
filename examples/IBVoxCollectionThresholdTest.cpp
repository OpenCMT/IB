#include "testing-prototype.h"
#include "IBVoxCollectionCap.h"

int main()
{
    BEGIN_TESTING(IB Vox Collection Threshold);

    IBVoxCollectionCap voxels(Vector3i(100,100,100));
    voxels.SetSpacing(Vector3f(1,1,1));

    IBVoxel zero = {0.1E-6,0,0};

    voxels.InitLambda(zero);

    voxels[voxels.Find(HPoint3f(2,2,2))].Value = 0.1;
    voxels[voxels.Find(HPoint3f(4,4,4))].Value = 0.2;
    voxels[voxels.Find(HPoint3f(10,10,10))].Value = 1;
    voxels[voxels.Find(HPoint3f(50,50,50))].Value = 10;

    float threshold = 0.1;
    int count = voxels.CountLambdaOverThreshold(threshold);
    std::cout << "found " << count << " voxels over " << threshold << "\n";
    TEST1( count == 4 );

    threshold = 0.2;
    count = voxels.CountLambdaOverThreshold(threshold);
    std::cout << "found " << count << " voxels over " << threshold << "\n";
    TEST1( count == 3 );

    count = voxels.CountLambdaOverThreshold(threshold,HPoint3f(10,10,10), HVector3f(1,1,1));
    std::cout << "found " << count << " voxels over " << threshold << " around point 3 \n";
    TEST1( count == 1 );


    END_TESTING;
}

