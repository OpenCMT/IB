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
#include "../src/IBVoxFilters.h"

int main()
{
    BEGIN_TESTING(IB Vox Filters);

    IBVoxCollectionCap voxels(Vector3i(3,3,3));

    IBVoxel zero = {0,0,0};
    voxels.InitLambda(zero);

    voxels[Vector3i(1,1,1)].Value = 1E-6;

    IBVoxFilter_Linear filter(Vector3i(3,3,3));

    IBFilterGaussShape shape(2);
    filter.SetKernelSpherical(shape);
//    filter.SetKernelWeightFunction(shape);

    filter.SetImage(&voxels);

    filter.Run();

    voxels.ExportToVtk("filter_tes.vtk",0);


    END_TESTING;
}
