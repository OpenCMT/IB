/*////////////////////////////////////////////////////////////////////////////
 Copyright 2018 Istituto Nazionale di Fisica Nucleare

 Licensed under the EUPL, Version 1.2 or - as soon they will be approved by
 the European Commission - subsequent versions of the EUPL (the "Licence").
 You may not use this work except in compliance with the Licence.

 You may obtain a copy of the Licence at:

 https://joinup.ec.europa.eu/software/page/eupl

 Unless required by applicable law or agreed to in writing, software
 distributed under the Licence is distributed on an "AS IS" basis, WITHOUT
 WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
 Licence for the specific language governing permissions and limitations under
 the Licence.
////////////////////////////////////////////////////////////////////////////*/



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
                     "use: abtrim filename.vtk voxel_size s_sigma i_sigma a_trim b_trim scale\n";
        exit(1);
    }

    struct Params {
        char *file;
        int size;
        float ssigma;
        float isigma;
        int a,b;
        Scalarf scale;
    } parameters = {
        argv[1],
                5,   // default size
                0.7, // default ssigma
                1,   // default isigma
                0,
                0,
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

    if(argc == 8)
    {
        parameters.size = atoi(argv[2]);
        parameters.ssigma = atof(argv[3]);
        parameters.isigma = atof(argv[4]);
        parameters.a      = atoi(argv[5]);
        parameters.b      = atoi(argv[6]);
        parameters.scale  = atof(argv[7]);
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
    IBVoxFilter_BilateralTrim filter(Vector3i(s,s,s));
    IBFilterGaussShape shape(parameters.ssigma);
    filter.SetKernelWeightFunction(shape);
    filter.SetIntensitySigma(parameters.isigma);
    filter.SetABTrim(parameters.a, parameters.b);
    filter.SetImage(&image);
    filter.Run();

    image *= parameters.scale;

    sprintf(filename,"%s_bilateral_%d%d%d_s%.2f_i%.2f_scale%.2f.vtk",
            FileNameRemoveExtension(parameters.file).c_str(),
            parameters.size,parameters.a,parameters.b, parameters.ssigma,
            parameters.isigma, parameters.scale);
    image.ExportToVtk(filename,0);


    return 0;
}
