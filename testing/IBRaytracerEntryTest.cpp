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

#include <Vtk/uLibVtkViewer.h>
#include <Vtk/vtkMuonEvent.h>
#include <Vtk/vtkStructuredGrid.h>

using namespace uLib;

int main() {

    // errors //
    IBMuonError sigma(12.24,18.85);

    // reader //
    TFile* f = new TFile("/var/local/data/root/ROC_sets/201212/lead/muSteel_PDfit_2012122300_v11.root");
    TTree* t = (TTree*)f->Get("n");
    IBMuonEventTTreeReader* reader = IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::R3D_MC);
    reader->setTTree(t);
    reader->setError(sigma);
    reader->setMomentum(0.7);
    reader->selectionCode(IBMuonEventTTreeR3DmcReader::Top2Bottom);

    // voxels //
    IBVoxel zero = {0.1E-6,0,0};
    IBVoxCollectionCap voxels(Vector3i(80,48,30));
    voxels.SetSpacing(Vector3f(10,10,10));
    voxels.SetPosition(Vector3f(-400,-300,-150));
    voxels.InitLambda(zero);


    // pocal //
    IBPocaEvaluator* processor = IBPocaEvaluator::New(IBPocaEvaluator::TiltedAxis);

    // tracer //
    IBVoxRaytracer* tracer = new IBVoxRaytracer(voxels);




    std::cout << "There are " << reader->getNumberOfEvents() << " events!\n";

    std::cout << "Evaluating entry and exit point of container for each track\n";
    HPoint3f entry_pt,exit_pt;
    int out_trace = 0;
    for (int i=0; i<1000000/*reader->getNumberOfEvents()*/; i++) {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            if( !tracer->GetEntryPoint(mu.LineIn(),entry_pt) ||
                    !tracer->GetExitPoint(mu.LineOut(),exit_pt) )
                out_trace ++;
        }
    }

    std::cout << out_trace << " traces are out of voxels container\n";
    return 0;
}

