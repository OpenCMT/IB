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
#include "IBMinimizationVariablesEvaluator.h"

#include "IBAnalyzerTrackLengths.h"
#include "IBAnalyzerWTrackLengths.h"

#include <Vtk/uLibVtkViewer.h>
#include <Vtk/vtkMuonEvent.h>
#include <Vtk/vtkStructuredGrid.h>




using namespace uLib;

int main() {

    // errors //
//    IBMuonError sigma(11.93,2.03, 18.53,2.05);
    IBMuonError sigma(12.24,0, 18.85,0);

    // reader //
    TFile* f = new TFile ("/home/rigoni/muSteel_PDfit_20130127_v13.root");
    TTree* t = (TTree*)f->Get("n");
    IBMuonEventTTreeReader* reader = IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::R3D_MC);
    reader->setTTree(t);
    reader->setError(sigma);
    reader->setMomentum(0.7);
    reader->selectionCode(0xff);

    // voxels //
    IBVoxel zero = {0.1E-6,0,0};
    IBVoxCollectionCap voxels(Vector3i(140,72,60));
    voxels.SetSpacing (Vector3f(5,5,5));
    voxels.SetPosition(Vector3f(-350,-180,-150));
    voxels.InitLambda(zero);


    // pocal //
    IBPocaEvaluator* processor = IBPocaEvaluator::New(IBPocaEvaluator::TiltedAxis);

    // tracer //
    IBVoxRaytracer* tracer = new IBVoxRaytracer(voxels);

    // variables //
    IBMinimizationVariablesEvaluator* minimizator =
            IBMinimizationVariablesEvaluator::New(IBMinimizationVariablesEvaluator::NormalPlane);
    minimizator->setRaytracer(tracer);


    // analyzer //
    IBAnalyzerWTrackLengths atc;
    atc.SetVoxCollection(&voxels);
    atc.SetPocaAlgorithm(processor);
    atc.SetRaytracer(tracer);
    atc.SetVaraiblesAlgorithm(minimizator);

    std::cout << "There are " << reader->getNumberOfEvents() << " events!\n";


    int tot  = 0;
    int tot2 = 0;
    do {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            atc.AddMuon(mu);
            tot++;
            if(tot%10000 == 0 ) std::cout<<tot<<"\n"<<std::flush;
        }
        tot2++;
    } while (tot2<2000000);


    char file[20];

    atc.Run(1,1);
    sprintf(file, "track_count_W2L.vtk");
    voxels.ExportToVtk(file,0);

    return 0;

}

