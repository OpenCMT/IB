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
#include <TH1F.h>
#include <IBVoxFilters.h>
#include "IBVoxCollectionCap.h"
#include <Vtk/vtkVoxImage.h>

using namespace uLib;

int main()
{
    IBVoxCollectionCap map_img(Vector3i(0,0,0));
    IBVoxCollectionCap image(Vector3i(0,0,0));

    vtkVoxImage map_vtk(map_img);
    vtkVoxImage image_vtk(image);
    map_vtk.ReadFromVKTFile("20121224_v495_2PXTZ_200.vtk");
    image_vtk.ReadFromVKTFile("20121218_v490_2PXTZ_200.vtk");

    IBVoxFilter_Abtrim trim(Vector3i(5,5,5));
    IBFilterGaussShape shape(0.2);
    trim.SetKernelSpherical(shape);
    trim.SetABTrim(0,2);
    trim.SetImage(&image);
    trim.Run();
   // image.ExportToVtk("20121218_trimmed.vtk",0);
    trim.SetABTrim(0,2);
    trim.SetImage(&map_img);
    trim.Run();
   // map_img.ExportToVtk("20121224_trimmed.vtk",0);


    IBVoxFilter_Linear gauss(Vector3i(20,20,20));
    gauss.SetKernelSpherical(shape);
    gauss.SetImage(&map_img);
    gauss.Run();

    IBVoxFilter_Plasmon plas;
    plas.SetImage(&image);
    plas.SetMappingImage(&map_img);
    plas.Run();

    {
        int hi_y_start = 14,  hi_y_stop = 22;
        int hi_z_start = 8,   hi_z_stop = 16;
        int me_y_start = 31,  me_y_stop = 39;
        int me_z_start = 26,  me_z_stop = 33;
        int lo_y_start = 49,  lo_y_stop = 57;
        int lo_z_start = 43,  lo_z_stop = 50;
        TFile f("carroting_20121218_omogenized_smoothBG.root","RECREATE");
        TH1F * hi = new TH1F("lower","lower", image.GetDims()(0), 0, image.GetDims()(0)-1);
        TH1F * me = new TH1F("mider","mider", image.GetDims()(0), 0, image.GetDims()(0)-1);
        TH1F * lo = new TH1F("upper","upper", image.GetDims()(0), 0, image.GetDims()(0)-1);
        for (int i=0; i<image.GetDims()(0); ++i) {
            float tmp=0;
            for (int j=hi_y_start; j<=hi_y_stop; ++j) {
                for (int k=hi_z_start; k<=hi_z_stop; ++k) {
                    tmp+=image.At(Vector3i(i,j,k)).Value;
                }
            }
            hi->Fill(i,tmp);
        }
        for (int i=0; i<image.GetDims()(0); ++i) {
            float tmp=0;
            for (int j=me_y_start; j<=me_y_stop; ++j) {
                for (int k=me_z_start; k<=me_z_stop; ++k) {
                    tmp+=image.At(Vector3i(i,j,k)).Value;
                }
            }
            me->Fill(i,tmp);
        }
        for (int i=0; i<image.GetDims()(0); ++i) {
            float tmp=0;
            for (int j=lo_y_start; j<=lo_y_stop; ++j) {
                for (int k=lo_z_start; k<=lo_z_stop; ++k) {
                    tmp+=image.At(Vector3i(i,j,k)).Value;
                }
            }
            lo->Fill(i,tmp);
        }
        hi->Write();
        me->Write();
        lo->Write();
        f.Close();
    }


    image.ExportToVtk("20121218_PlastFiltered_20121224_smooth.vtk",0);
    return 0;
}


