#include <TFile.h>
#include <TTree.h>
#include <stdio.h>
#include "IBVoxCollectionCap.h"
#include "IBSubImageGrabber.h"
#include "IBVoxImageScanner.h"
#include <IBVoxFilters.h>

void PrintRData(RangeThresholdScan::ScanData d);
void PrintData(SimpleThresholdScan::ScanData d);

int main()
{

    RangeThresholdScan::ScanOption opt;
    for(int i=0; i<20; ++i) {
        SimpleThresholdScan::ScanOption o;
        o.Threshold = (0.5e-6)*i;
        opt.push_back(o);
    }
    // load
    IBVoxCollectionCap image(Vector3i(0,0,0));
    image.ImportFromVtk("./roc_set/lead_image_1.vtk");
    // filter
    IBVoxFilter_Abtrim trim(Vector3i(5,5,5));
    IBFilterGaussShape shape(0.2);
    trim.SetKernelSpherical(shape);
    trim.SetABTrim(0,2);
    trim.SetImage(&image);
    trim.Run();
    // grabbing
    Vector3i z(0,0,0);
    IBSubImageGrabber<IBVoxCollectionCap> grabber(image);
    IBLightCollection imgLR_B(z),imgLR_M(z), imgLR_F(z);
    imgLR_B = grabber.GrabRegion<IBLightCollection>(Vector3i(21,14,9),Vector3i(28,22,16));
    imgLR_M = grabber.GrabRegion<IBLightCollection>(Vector3i(46,14,9),Vector3i(53,22,16));
    imgLR_F = grabber.GrabRegion<IBLightCollection>(Vector3i(76,14,9),Vector3i(83,22,16));
    IBLightCollection imgMR_B(z), imgMR_M(z), imgMR_F(z);
    imgMR_B = grabber.GrabRegion<IBLightCollection>(Vector3i(36,31,26),Vector3i( 43,39,33));
    imgMR_M = grabber.GrabRegion<IBLightCollection>(Vector3i(66,31,26),Vector3i( 73,39,33));
    imgMR_F = grabber.GrabRegion<IBLightCollection>(Vector3i(96,31,26),Vector3i(103,39,33));
    IBLightCollection imgHR_B(z), imgHR_M(z), imgHR_F(z);
    imgHR_B = grabber.GrabRegion<IBLightCollection>(Vector3i( 56,49,43),Vector3i( 63,57,50));
    imgHR_M = grabber.GrabRegion<IBLightCollection>(Vector3i( 86,49,43),Vector3i( 93,57,50));
    imgHR_F = grabber.GrabRegion<IBLightCollection>(Vector3i(111,49,43),Vector3i(118,57,50));
    std::vector<IBLightCollection> boxC;
    boxC.push_back(imgLR_B);
    boxC.push_back(imgLR_M);
    boxC.push_back(imgLR_F);
    boxC.push_back(imgMR_B);
    boxC.push_back(imgMR_M);
    boxC.push_back(imgMR_F);
    boxC.push_back(imgHR_B);
    boxC.push_back(imgHR_M);
    boxC.push_back(imgHR_F);
    // scanning
    Vector<RangeThresholdScan::ScanData> sdBox;
    IBVoxImageScanner<IBLightCollection> scanner;
    for(int i=0; i<boxC.size(); ++i) {
        scanner.SetImage(&boxC.at(i));
        RangeThresholdScan::ScanData res = scanner.ScanImage<RangeThresholdScan>(opt);
        sdBox.push_back(res);
    }


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
