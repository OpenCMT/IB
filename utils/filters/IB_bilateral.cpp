
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
        float ssigma;
        float isigma;
        Scalarf scale;
    } parameters = {
        argv[1],
                5, // default size
                0.5, // default ssigma
                1, // default isigma
                1  // default scaling factor
    };

    if(argc == 3)
    {
        parameters.size = atof(argv[2]);
    }

    if(argc == 6)
    {
        parameters.size = atoi(argv[2]);
        parameters.ssigma = atof(argv[3]);
        parameters.isigma = atof(argv[4]);
        parameters.scale = atof(argv[5]);
    }


    std::cout << "// --------- [abtrim] ----------------- //\n"
              << " size [ " << parameters.size << "] "
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
    IBVoxFilter_Bilateral filter(Vector3i(s,s,s));
    IBFilterGaussShape shape(sqrt(parameters.ssigma));
    filter.SetKernelWeightFunction(shape);
    filter.SetIntensitySigma(parameters.isigma);
    filter.SetImage(&image);
    filter.Run();

    image *= parameters.scale;

    sprintf(filename,"%s_bilateral%d_s%.2f_i%.2f_scale%.2f.vtk",
            FileNameRemoveExtension(parameters.file).c_str(),
            parameters.size, parameters.ssigma, parameters.isigma,
            parameters.scale);
    image.ExportToVtk(filename,0);


    return 0;
}