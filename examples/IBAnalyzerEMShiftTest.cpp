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

#include "IBMuonCollection.h"

using namespace uLib;

int main() {



    // errors //
    IBMuonError sigma(12.24,0.0,
                      18.85,0.0,
                      1.4);
    // reader //

//    TFile* f = new TFile ("/var/local/data/root/ROC_sets/201301/20130122.root");
    TFile* f = new TFile ("/var/local/data/root/muSteel_PDfit_20130203_v13.root");
    TTree* t = (TTree*)f->Get("n");
    IBMuonEventTTreeReader* reader = IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::R3D_MC);
    reader->setTTree(t);
    reader->setError(sigma);
    reader->setMomentum(0.7);
    reader->selectionCode(IBMuonEventTTreeR3DmcReader::All);

    // voxels //
    IBVoxel zero = {0.1E-6,0,0};
    IBVoxCollection voxels(Vector3i(140,72,60));
    voxels.SetSpacing (Vector3f(5,5,5));
    voxels.SetPosition(Vector3f(-350,-180,-150));
    voxels.InitLambda(zero);


    // poca //
    IBPocaEvaluator* processor = IBPocaEvaluator::New(IBPocaEvaluator::LineDistance);

    // tracer //
    IBVoxRaytracer* tracer = new IBVoxRaytracer(voxels);

    // variables //
    IBMinimizationVariablesEvaluator* minimizator =
            IBMinimizationVariablesEvaluator::New(IBMinimizationVariablesEvaluator::NormalPlane);
    minimizator->setRaytracer(tracer);

    // ML Algorithm //
    IBAnalyzerEMAlgorithmSGA_PXTZ4 *ml_algorithm =
            new IBAnalyzerEMAlgorithmSGA_PXTZ4;  //////// 4 <<<<<<<<< !!!

    // analyzer //
    IBAnalyzerEM* aem = new IBAnalyzerEM(voxels);
    aem->SetMLAlgorithm(ml_algorithm);
    aem->SetPocaAlgorithm(processor);
    aem->SetRayAlgorithm(tracer);
    aem->SetVarAlgorithm(minimizator);

    // filter //
    IBVoxFilter_Abtrim trim(Vector3i(3,3,3));
    IBFilterGaussShape shape(0.2);
    trim.SetKernelSpherical(shape);
    trim.SetABTrim(0,1);
    trim.SetImage(&voxels);


    int tot=0;


    IBMuonCollection muons;
    //    int ev = 1375250;
    reader->setAcquisitionTime(5);
    std::cout << "There are " << reader->getNumberOfEvents() << " events!\n";
    int ev = reader->getNumberOfEvents();
    for (int i=0; i<ev; i++) {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            muons.AddMuon(mu);
            tot++;
        }
    }
    aem->SetMuonCollection(&muons);


    char file[100];

    int it   = 1;
    int drop = 200;

    aem->SijCut(60);
    std::cout << "Spared: [" << aem->Size() << "]\n";


    voxels.InitLambda(zero);
    IBVoxCollection shift_accum = voxels;
    shift_accum -= zero.Value;

    std::cout << "NO SHIFT: " << voxels.GetPosition().transpose() << "\n";
    aem->Run(drop,1);
    trim.SetImage(&voxels);
    trim.Run();
    shift_accum += voxels;

    // SHIFTING //
    Vector3f pos = voxels.GetPosition();
        for(int s=0; s<6; ++s)
    {
        voxels.InitLambda(zero);
        Vector3f position = pos;
        position(s%3) += (1-2*(s%2)) * voxels.GetSpacing()(s%3) / 5;
        std::cout << "SHIFT: " << position.transpose() << "\n";
        voxels.SetPosition(position);
        aem->SetVoxCollection(&voxels);
        aem->Run(drop,1);
//        trim.Run();

        sprintf(file, "0203_PXTZ_p2_SHIFT%i.vtk", s);
        voxels.ExportToVtk(file,0);
        shift_accum += voxels;
    }
    shift_accum /= 7;

    sprintf(file, "0203_PXTZ_p2_SHIFT.vtk");
    shift_accum.ExportToVtk(file,0);



    delete aem;
    delete minimizator;
    return 0;




}

