#include <iostream>

#include <TFile.h>
#include <TTree.h>

#include "IBPocaEvaluator.h"
#include "IBAnalyzerTrackCount.h"
#include "IBVoxRaytracer.h"
#include "IBMinimizationVariablesEvaluator.h"
#include "IBVoxCollectionCap.h"
#include "IBMuonError.h"
#include "Detectors/MuonScatter.h"
#include "IBMuonEventTTreeReader.h"
#include "IBMuonEventTTreeR3DmcReader.h"

#include <Vtk/uLibVtkViewer.h>
#include <Vtk/vtkMuonEvent.h>
#include <Vtk/vtkStructuredGrid.h>


using namespace uLib;

int main() {

    // errors //
    IBMuonError sigma(11.93,2.03, 18.53,2.05);

    // reader //
    TFile* f = new TFile("/var/local/data/root/run_PDfit_201210/muSteel_PDfit_20121016_v10.root");
    TTree* t = (TTree*)f->Get("n");
    IBMuonEventTTreeReader* reader = IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::R3D_MC);
    reader->setTTree(t);
    reader->setError(sigma);
    reader->setMomentum(0.7);
    reader->selectionCode(0xff);

    // voxels //
    IBVoxel zero = {0.1E-6,0,0};
    IBVoxCollectionCap voxels(Vector3i(160,96,60));
    voxels.SetSpacing(Vector3f(5,5,5));
    voxels.SetPosition(Vector3f(-400,-300,-150));
    voxels.InitLambda(zero);


    // pocal //
    IBPocaEvaluator* processor = IBPocaEvaluator::New(IBPocaEvaluator::TiltedAxis);

    // tracer //
    IBVoxRaytracer* tracer = new IBVoxRaytracer(voxels);


    // analyzer //
    IBAnalyzerTrackCount atc;
    atc.SetVoxCollection(&voxels);
    atc.SetPocaAlgorithm(processor);
    atc.SetRaytracer(tracer);


    std::cout << "There are " << reader->getNumberOfEvents() << " events!\n";

    for (int i=0; i<reader->getNumberOfEvents(); i++) {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            atc.AddMuon(mu);
        }
    }

    char file[20];

    atc.Run(1,1);
    sprintf(file, "track_count.vtk");
    voxels.ExportToVtk(file,0);

    return 0;

}

