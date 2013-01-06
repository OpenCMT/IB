#include <TFile.h>
#include <TH1F.h>
#include <stdio.h>
#include "IBVoxCollectionCap.h"
#include "IBSubImageGrabber.h"
#include "IBVoxImageScanner.h"
#include <IBVoxFilters.h>
#include <Vtk/vtkVoxImage.h>

int main()
{
    // load
    IBVoxCollectionCap image(Vector3i(0,0,0));
    vtkVoxImage loadingMask(image);
    loadingMask.ReadFromVKTFile("original.vtk");
    // filter
    IBVoxFilter_Abtrim trim(Vector3i(5,5,5));
    IBFilterGaussShape shape(0.2);
    trim.SetKernelSpherical(shape);
    trim.SetABTrim(0,2);
    trim.SetImage(&image);
    trim.Run();
    // boxing
    Box xlylzl;
//TO COMPLETE!


}
