#include <iostream>

#include <TFile.h>
#include <TTree.h>


#include "IBPocaEvaluator.h"
#include "IBAnalyzerEM.h"
#include "IBVoxRaytracer.h"
#include "IBMinimizationVariablesEvaluator.h"
#include "IBVoxCollectionCap.h"
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

using namespace uLib;

int main() {

    // errors //
    Voxel zero = {0};
    IBLightCollection scraps(Vector3i(140,72,60));
    scraps.SetSpacing (Vector3f(5,5,5));
    scraps.SetPosition(Vector3f(-350,-180,-150));
    scraps.InitVoxels(zero);

    for(int x=10; x < 130; ++x) {
        for (int y=10; y < 62; ++y) {
            for (int z=4; z<56; ++z) {
                Vector3i id(x,y,z);
                scraps[id].Value = 0.00427;
            }
        }
    }
    IBMuonError sigma(12.24,
                      18.85,
                      1.4);
    sigma.setScrapsImage(scraps,1);


    // reader //
    TFile* f = new TFile ("/var/local/data/root/muSteel_PDfit_20130203_v14.root");
//    TFile* f = new TFile ("/var/local/data/root/ROC_sets/201301/20130111/muSteel_PDfit_20130111_10_v12.root");
//    TFile* f = new TFile ("/var/local/data/root/ROC_sets/201212/20121223/muSteel_PDfit_20121223_10_v11.root");

    if (f->IsZombie()) {
        std::cerr << "file not found!\n";
        exit(1);
    }

    TTree* t = (TTree*)f->Get("n");
//    IBMuonEventTTreeReader* reader = IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::R3D_MC);
    IBMuonEventTTreeR3DmcReader *reader = new IBMuonEventTTreeR3DmcReader();
    reader->setTTree(t);
    reader->setError(sigma);
    reader->setMomentum(0.7);
    reader->selectionCode(IBMuonEventTTreeR3DmcReader::All);

    // voxels //

    IBVoxel air = {0.1E-6,0,0};
    IBVoxCollectionCap voxels(Vector3i(140,72,60));
    voxels.SetSpacing (Vector3f(5,5,5));
    voxels.SetPosition(Vector3f(-350,-180,-150));
    voxels.InitLambda(air);



    // poca //
    IBPocaEvaluator* processor = IBPocaEvaluator::New(IBPocaEvaluator::LineDistance);

    // tracer //
    IBVoxRaytracer* tracer = new IBVoxRaytracer(voxels);

    // variables //
    IBMinimizationVariablesEvaluator* minimizator =
            IBMinimizationVariablesEvaluator::New(IBMinimizationVariablesEvaluator::NormalPlane);
    minimizator->setRaytracer(tracer);

    // analyzer //
    IBAnalyzerEM* aem = new IBAnalyzerEM;
    aem->SetVoxCollection(&voxels);
    aem->SetPocaAlgorithm(processor);
    aem->SetRaytracer(tracer);
    aem->SetVariablesAlgorithm(minimizator);


    reader->setAcquisitionTime(5);
    std::cout << "There are " << reader->getNumberOfEvents() << " events!\n";

    int tot=0;    
    int ev = reader->getNumberOfEvents();
    for (int i=0; i<ev; i++) {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            aem->AddMuon(mu);
            tot++;
        }
    }



    char file[100];

    int it   = 100;
    int drop = 5;

    aem->SijCut(60);
    std::cout << "Survived: [" << aem->Size() << "]\n";

    // SGA //
    std::cout << "SGA PXTZ\n";
    for (int i=1; i<=it; ++i) {
        aem->Run(drop,1);
        sprintf(file, "20130203_PXTZ_esauto_%i.vtk", i*drop);
        voxels.ExportToVtk(file,0);
    }


    delete aem;
    delete minimizator;
    return 0;

}

