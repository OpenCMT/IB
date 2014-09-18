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



#include <stdio.h>
#include "IBVoxCollectionCap.h"
#include "IBSubImageGrabber.h"
#include "IBVoxImageScanner.h"

int main()
{
    IBVoxCollectionCap image(Vector3i(10,10,10));
    image.SetSpacing(Vector3f(1,1,1));
    image.SetPosition(Vector3f(-5,-5,-5));
    image.SetValue(Vector3i(2,3,4), 5.0e-6);
    image.SetValue(Vector3i(2,2,2), 4.5e-6);

    IBSubImageGrabber<IBVoxCollectionCap> g;
    g.SetImage(image);
    IBVoxCollection c = g.GrabRegion<IBVoxCollection>(
                Vector3i(1,1,1), Vector3i(5,5,5));

    IBVoxImageScanner<IBVoxCollection> scanner(c);
    SimpleThresholdScan::ScanOption    opt;
    SimpleThresholdScan::ScanData      output;
    opt.Threshold = 4.6E-6;
    output = scanner.ScanImage<SimpleThresholdScan>(opt);

    PrintData(output);

    RangeThresholdScan::ScanOption r_opt;
    for(int i=0; i<7; ++i) {
        SimpleThresholdScan::ScanOption o;
        o.Threshold = i*1E-6;
        r_opt.push_back(o);
    }
    RangeThresholdScan::ScanData r_output;

    r_output = scanner.ScanImage<RangeThresholdScan>(r_opt);

    PrintRData(r_output);

    return 0;
}

void PrintData(SimpleThresholdScan::ScanData d)
{
    std::cout << "Percent   = " << d.Percent << "\%"<<std::endl;
    std::cout << "Intensity = " << d.Intensity << std::endl;
    std::cout << "------------" << std::endl;
}
void PrintRData(RangeThresholdScan::ScanData d)
{
    for(int i=0; i<d.size(); ++i) PrintData(d.at(i));
}
