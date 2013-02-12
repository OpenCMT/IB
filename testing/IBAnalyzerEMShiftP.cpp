
#include <Math/Dense.h>
#include <Math/VoxImage.h>

#include "IBVoxCollection.h"

using namespace uLib;



int main()
{

    std::vector<std::string> files;

    files.push_back(std::string("./0203_PXTZ_p2_SHIFT_0.vtk"));
    files.push_back(std::string("./0203_PXTZ_p2_SHIFT_1.vtk"));
    files.push_back(std::string("./0203_PXTZ_p2_SHIFT_2.vtk"));
    files.push_back(std::string("./0203_PXTZ_p2_SHIFT_3.vtk"));
    files.push_back(std::string("./0203_PXTZ_p2_SHIFT_4.vtk"));
    files.push_back(std::string("./0203_PXTZ_p2_SHIFT_5.vtk"));

    Vector<IBVoxCollection *> images;

    for(int i=0; i<6; ++i)
    {
        // image //

        IBVoxCollection *image = new IBVoxCollection(Vector3i(0,0,0));
        image->ImportFromVtk((files[i]).c_str());
        images.push_back(image);
        std::cout << "image size: " << image->GetDims().transpose() << "\n";
    }


    IBVoxCollection result = *images[1];

    for(int id=0; id<images[0]->Data().size(); ++id )
    {
        // plane mean //
        float max = 0;
        for(int i=0; i<6; ++i)
        {
            IBVoxCollection &image = *images[i];
            if(image[id].Value > max) max = image[id].Value;
            result[id].Value += image[id].Value;
        }
        result[id].Value -= max;
        result[id].Value /= 5;
    }

    result.ExportToVtk("shifted_result.vtk");

}
