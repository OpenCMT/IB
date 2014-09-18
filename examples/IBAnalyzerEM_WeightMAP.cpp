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




#include <iostream>

#include <TFile.h>
#include <TTree.h>

#include "IBPocaEvaluator.h"
#include "IBAnalyzerEM.h"
#include "IBAnalyzerEMAlgorithm.h"
#include "IBVoxRaytracer.h"
#include "IBMinimizationVariablesEvaluator.h"
#include "IBVoxCollection.h"
#include "IBMuonError.h"
#include "Detectors/MuonScatter.h"
#include "IBMuonEventTTreeReader.h"
#include "IBMuonEventTTreeR3DmcReader.h"
#include "IBVoxFilters.h"
#include "IBAnalyzerWPoca.h"
#include "IBMAPUpdateDensityAlgorithms.h"

#include <Vtk/uLibVtkViewer.h>
#include <Vtk/vtkMuonEvent.h>
#include <Vtk/vtkStructuredGrid.h>

#include "IBAnalyzerEMAlgorithmSGA.h"
#include "IBAnalyzerEMAlgorithmMGA.h"

using namespace uLib;



int main(void) {

    // errors //
    IBMuonError sigma(12.24, 0.0,
                      18.85, 0.0,
                      1.4 );
    // reader //

    TFile* f = new TFile ("/var/local/data/root/ROC_sets/201212/20121223/muSteel_PDfit_20121223_1_v11.root");
    TTree* t = (TTree*)f->Get("n");
    IBMuonEventTTreeReader* reader = IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::R3D_MC);
    reader->setTTree(t);
    reader->setError(sigma);
    reader->setMomentum(0.7);
    reader->selectionCode(IBMuonEventTTreeR3DmcReader::All);

    // voxels //
    IBVoxel zero = {0.22E-6,0,0};
    IBVoxCollection voxels(Vector3i(140,72,60));
    voxels.SetSpacing (Vector3f(5,5,5));
    voxels.SetPosition(Vector3f(-350,-180,-150));
    voxels.InitLambda(zero);

    // MAP Algorithm //
//    IBMAPPriorTotalWeigth weight_MAP(0.19, 300E-6 * 300E-6);
//    voxels.SetMAPAlgorithm(&weight_MAP);

    // poca //
    IBPocaEvaluator* processor = IBPocaEvaluator::New(IBPocaEvaluator::LineDistance);

    // tracer //
    IBVoxRaytracer* tracer = new IBVoxRaytracer(voxels);

    // variables //
    IBMinimizationVariablesEvaluator* minimizator =
            IBMinimizationVariablesEvaluator::New(IBMinimizationVariablesEvaluator::NormalPlane);
    minimizator->setRaytracer(tracer);


    // analyzer //
    IBAnalyzerEM* aem = new IBAnalyzerEM(voxels);
    //aem->SetMLAlgorithm(ml_algorithm);
    aem->SetPocaAlgorithm(processor);
    aem->SetRayAlgorithm(tracer);
    aem->SetVarAlgorithm(minimizator);

    IBAnalyzerEMAlgorithmSGA_PXTZ ml_algorithm;
    aem->SetMLAlgorithm(&ml_algorithm);


    // filter //
    IBVoxFilter_Abtrim trim(Vector3i(5,5,5));
    IBFilterGaussShape shape(0.2);
    trim.SetKernelSpherical(shape);
    trim.SetABTrim(0,2);
    trim.SetImage(&voxels);



    //reader->setAcquisitionTime(5);

    std::cout << "There are " << reader->getNumberOfEvents() << " events!\n";
    int tot=0;

    int ev = 1375250;
    for (int i=0; i<ev; i++) {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            aem->AddMuon(mu);
            tot++;
        }
    }


    char file[100];

    int it   = 50;
    int drop = 10;

    aem->SijCut(60);
    std::cout << "Survived: [" << aem->Size() << "]\n";


    // SGA //

    {
        voxels.InitLambda(zero);

        for (int i=1; i<=it; ++i) {
            aem->Run(drop,1);
            sprintf(file, "20121223_SGA_PXTZ_GW_%i.vtk", i*drop);
            voxels.ExportToVtk(file,0);

            IBVoxCollection filtered = voxels;
            trim.SetImage(&filtered);
            trim.Run();
            sprintf(file, "20121223_SGA_PXTZ_GW_trim_%i.vtk", i*drop);
            filtered.ExportToVtk(file,0);
        }
    }



    //MGA //

//    // 1m scraps
//    float w[5] = { 0.104 , 0.286,  0.340,  0.202,  0.063  };
//    float s[5] = { 0.0006, 0.0027, 0.0082, 0.0259, 0.0718 };

//    // 3m scraps
//    float w[5] = { 0.093 , 0.257,  0.337,  0.234,  0.076  };
//    float s[5] = { 0.0006, 0.0024, 0.0070, 0.0194, 0.0527 };

//    {
//        IBAnalyzerEMAlgorithmMGA_PXTZ<5> ml_algorithm_mga;
//        ml_algorithm_mga.SetGaussians(w,s);
//        aem->SetMLAlgorithm(&ml_algorithm_mga);
//        voxels.InitLambda(zero);

//        std::cout << "MGA PXTZ_5_\n";
//        for (int i=1; i<=it; ++i) {
//            aem->Run(drop,1);
//            sprintf(file, "20121223_MGA_5_EnoCorr_B0AZ24_PXTZ_%i.vtk", i*drop);
//            voxels.ExportToVtk(file,0);
//        }
//    }




    delete aem;
    delete minimizator;
    return 0;

}

