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


    std::cout << "// --------- [abtrim] ----------------- //\n"
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
    IBVoxFilter_Abtrim trim(Vector3i(s,s,s));
    std::vector <float> values;
    for (int i=0; i<trim.GetKernelData().GetDims().prod(); ++i) {
        values.push_back(1.);
    }
    trim.SetKernelNumericXZY(values);
    trim.SetABTrim(parameters.atrim, parameters.btrim);
    trim.SetImage(&image);
    trim.Run();

    image *= parameters.scale;

    sprintf(filename,"%s_trim%d%d%du_scale%.2f.vtk",
            FileNameRemoveExtension(parameters.file).c_str(),
            parameters.size, parameters.atrim, parameters.btrim,
            parameters.scale);
    image.ExportToVtk(filename,0);


    return 0;
}
