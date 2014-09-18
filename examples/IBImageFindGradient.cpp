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




#include <Math/Dense.h>
#include <Math/VoxImage.h>
#include "IBVoxCollection.h"
#include "IBSubImageGrabber.h"

using namespace uLib;

typedef VoxImage<Voxel> Image;


float get_subimage_mean(IBVoxCollection &image,Vector3i &p0, Vector3i &p1) {
    IBSubImageGrabber<IBVoxCollection> grabber(image);
    Image subimg = grabber.GrabRegion<Image>(p0,p1);

    float mean =0;
    for (int i=0; i<subimg.Data().size(); ++i)
        mean += subimg.Data()[i].Value;
    return mean/subimg.Data().size();

}



int main()
{   
    Vector3i vstart(13,13,9); // starting voxel //
    Vector3i nblocks(10,5,4);   // # of cubes //
    Vector3i bsize(4,4,4);     // in voxels //
    Vector3i bstep(12,10,12);     // in voxels //

    std::string path("./");
    std::vector<std::string> files;

//    files.push_back(std::string("muSteel_PDfit_20130129_v13_500.vtk"));
//    files.push_back(std::string("muSteel_PDfit_20130130_v13_500.vtk"));
//    files.push_back(std::string("muSteel_PDfit_20130131_v13_500.vtk"));
//    files.push_back(std::string("muSteel_PDfit_20130201_v13_500.vtk"));
//    files.push_back(std::string("muSteel_PDfit_20130202_v13_500.vtk"));

    files.push_back(std::string("muSteel_PDfit_20130204_v13_500.vtk"));
    files.push_back(std::string("muSteel_PDfit_20130205_v13_500.vtk"));
    files.push_back(std::string("muSteel_PDfit_20130206_v13_500.vtk"));
    files.push_back(std::string("muSteel_PDfit_20130207_v13_500.vtk"));
    files.push_back(std::string("muSteel_PDfit_20130208_v13_500.vtk"));


    for(int level=0; level<nblocks(1); level++)
    {
        std::cout << path+files[level] << "\n";
        // image //
        IBVoxCollection image(Vector3i(0,0,0));
        image.ImportFromVtk((path+files[level]).c_str());

        std::cout << "image size: " << image.GetDims().transpose() << "\n";

        // plane mean //
        float sum = 0;
        for(int z=0; z < nblocks(2); ++z) {
            for(int x = 0; x < nblocks(0); ++x) {
                Vector3i p0 = vstart + (bstep.cwiseProduct(Vector3i(x,level,z)));
                Vector3i p1 = p0 + bsize;
                std::cout << "getting mean for block [ " <<
                             p0.transpose() << " - " << p1.transpose() << " ] ... ";
                float mean = get_subimage_mean(image,p0,p1);
                std::cout << mean * 1e6 << "\n";
                sum += mean;
            }
        }

        float mean = sum / (nblocks(0) * nblocks(2));
        std::cout << "LEVEL MEAN : " << mean * 1e6 << "\n\n";

    }

}

