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




#include "testing-prototype.h"
#include "IBVoxCollectionCap.h"

int main()
{
    BEGIN_TESTING(IB Vox Collection);

    IBVoxCollectionCap voxels(Vector3i(10,5,10));
    voxels.SetPosition(Vector3f(-5,-3,-5));


    IBVoxel zero = {0.1E-6,0,0};

    voxels.InitLambda(zero);

    voxels[voxels.Find(HPoint3f(2,2,2))].Value = 5.552368E-6;

    voxels.ExportToVtkXml("test_vox_collection.vti",0);

    voxels.ExportToVtk("test_vox_collection.vtk",0);


    END_TESTING;
}
