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



#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <stdio.h>
#include "IBVoxCollectionCap.h"
#include "IBSubImageGrabber.h"
#include "IBVoxImageScanner.h"
#include <IBVoxFilters.h>
//ridurre a scan singolo

void PrintRData(RangeThresholdScan::ScanData d);
void PrintData(SimpleThresholdScan::ScanData d);

int main()
{
    RangeThresholdScan::ScanOption opt;
    for(int i=0; i<10; ++i) {
        SimpleThresholdScan::ScanOption o;
        o.Threshold = (0.1e-6)*i;
        opt.push_back(o);
    }
    char fname[50];
    sprintf(fname, "./roc_set/lead_image_1.vtk");
    IBVoxCollectionCap image(Vector3i(0,0,0));

    image.ImportFromVtk(fname);
    // filter
    IBVoxFilter_Abtrim trim(Vector3i(5,5,5));
    IBFilterGaussShape shape(0.2);
    trim.SetKernelSpherical(shape);
    trim.SetABTrim(0,2);
    trim.SetImage(&image);
    trim.Run();
    // grabbing
    IBSubImageGrabber<IBVoxCollectionCap> grabber(image);
    IBLightCollection box =  grabber.GrabRegion<IBLightCollection>(
                Vector3i(21,14,9),Vector3i(28,22,16));
    IBVoxImageScanner<IBLightCollection> scanner;
    scanner.SetImage(&box);
    RangeThresholdScan::ScanData res = scanner.ScanImage<RangeThresholdScan>(opt);
    PrintRData(res);
    return 0;
}

void PrintData(SimpleThresholdScan::ScanData d)
{
    std::cout << "Percent = " << d.Percent << "\% | Intensity = " << d.Intensity << std::endl;
}

void PrintRData(RangeThresholdScan::ScanData d)
{
    for(int i=0; i<d.size(); ++i) PrintData(d.at(i));
}
