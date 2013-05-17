#include <string>

#include <TFile.h>
#include <TH1F.h>
#include <IBVoxFilters.h>
#include "IBVoxCollection.h"
#include <Vtk/vtkVoxImage.h>

using namespace uLib;



std::string GetFileExtension(const std::string& FileName)
{
    if(FileName.find_last_of(".") != std::string::npos)
        return FileName.substr(FileName.find_last_of(".")+1);
    return "";
}

std::string FileNameRemoveExtension(const std::string& FileName)
{
    if(FileName.find_last_of(".") != std::string::npos)
        return FileName.substr(0,FileName.find_last_of("."));
    return "";
}


int main(int argc, char *argv[])
{

    if(argc == 1)
    {
        std::cerr << "No input filename given ..\n" <<
                     "use: abtrim filename.vtk voxel_size a_trim b_trim scale\n";
        exit(1);
    }

    struct Params {
        char *file;
        int size;
        int atrim;
        int btrim;
        Scalarf scale;
    } parameters = {
        argv[1],
                5, // default size
                0, // default a
                2, // default b
                1  // default scaling factor
    };

    if(argc == 3)
    {
        parameters.scale = atof(argv[2]);
    }

    if(argc == 6)
    {
        parameters.size = atoi(argv[2]);
        parameters.atrim = atoi(argv[3]);
        parameters.btrim = atoi(argv[4]);
        parameters.scale = atof(argv[5]);
    }


    std::cout << "// --------- [abtrim spr] ------------- //\n"
              << "size [ " << parameters.size << "] "
              << " a = " << parameters.atrim
              << " b = " << parameters.btrim << "\n"
              << " scale factor = " << parameters.scale << "\n"
              << "// ------------------------------------ //\n";


    IBVoxCollection image(Vector3i(0,0,0));

    char filename[150];
    sprintf(filename,"%s",parameters.file);

    if(!image.ImportFromVtk(filename))
    {
        std::cerr << "Error: Could not open file\n";
        exit(1);
    }

    int s = parameters.size;
    IBVoxFilter_SPR trim(Vector3i(s,s,s));
    IBFilterGaussShape shape(0.2);
    trim.SetKernelWeightFunction(shape);
    trim.SetABTrim(parameters.atrim, parameters.btrim);
    trim.SetImage(&image);
    trim.Run();

    image *= parameters.scale;

    sprintf(filename,"%s_trim%d%d%dspr_scale%.2f.vtk",
            FileNameRemoveExtension(parameters.file).c_str(),
            parameters.size, parameters.atrim, parameters.btrim,
            parameters.scale);
    image.ExportToVtk(filename,0);


    return 0;
}

