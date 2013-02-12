
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
    TFile* f = new TFile ("/var/local/data/root/muSteel_PDfit_20130203_v14.root");
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
    IBAnalyzerEMAlgorithmSGA *ml_algorithm = new IBAnalyzerEMAlgorithmSGA_PXTZ;

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
    trim.SetABTrim(0,0);

    // remove center of filter kernel //
//    int center = trim.GetKernelData().GetCenterData();
//    trim.GetKernelData()[center].Value = 0;


    // MAP Algorithm //
    IBMAPPriorNeighbourDensity MAP( voxels, 5E-6);
    MAP.SetFilter(&trim);
        voxels.SetMAPAlgorithm(&MAP);


    reader->setAcquisitionTime(5);
    std::cout << "There are " << reader->getNumberOfEvents() << " events!\n";
    int tot=0;

    IBMuonCollection muons;
    int ev = reader->getNumberOfEvents();
    for (int i=0; i<ev; i++) {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            muons.AddMuon(mu);
            tot++;
        }
    }

    std::cout << "Reader events: " << tot << "\n";

    muons.PrintSelf(std::cout);
    aem->SetMuonCollection(&muons);

    char file[100];

    int it   = 200;
    int drop = 1;

    aem->SijCut(60);
    std::cout << "Spared: [" << aem->Size() << "]\n";
    voxels.InitLambda(zero);

    // SGA //
    std::cout << "SGA PXTZ\n";
    for (int i=1; i<=it; ++i) {
        aem->Run(drop,1);
        sprintf(file, "20130203_PXTZ_p14_%i.vtk", i*drop);
        voxels.ExportToVtk(file,0);
//        trim.SetImage(&voxels);
//        trim.Run();
    }


//    for (int i=11; i<=it; ++i) {
//        aem->Run(drop,1);
//        sprintf(file, "20130203_PXTZ_p14_%i.vtk", i*drop);
//        voxels.ExportToVtk(file,0);
//        trim.SetImage(&voxels);
//        trim.Run();
//    }





    delete aem;
    delete minimizator;
    return 0;




}

