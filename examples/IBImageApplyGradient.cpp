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
#include "IBVoxFilters.h"

using namespace uLib;

std::string FileNameRemoveExtension(const std::string& FileName)
{
    std::string file;

    // remove path //
    if(FileName.find_last_of("/") != std::string::npos)
        file = FileName.substr(FileName.find_last_of("/")+1);
    // remove extension //
    if(file.find_last_of(".") != std::string::npos)
        return file.substr(0,file.find_last_of("."));
    return "";
}

int main(int argc, char **argv)
{
    char *filename = argv[1];
    float m = atof(argv[2]);

    IBVoxCollection image(Vector3i(0,0,0));
    image.ImportFromVtk(filename);

    std::cout << "image size: " << image.GetDims().transpose() << "\n";

    IBVoxFilter_Gradient filter;

    filter.SetImage(&image);
    filter.SetCoeff(m);
    filter.Run();


    image.ExportToVtk((std::string(filename)+"_grad.vtk").c_str());
    return 0;

}

