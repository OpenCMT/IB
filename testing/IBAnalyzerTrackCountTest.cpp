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
    IBMuonError sigma(12.24,18.85);

    // reader //
    TFile* f = new TFile ("/var/local/data/root/muSteel_PDfit_20130203_v14.root");
    //    TFile* f = new TFile ("/var/local/data/root/ROC_sets/201301/20130122/muSteel_PDfit_20130122_1_v12.root");
    //    TFile* f = new TFile ("/var/local/data/root/muSteel_PDfit_20130129_v13.root");
    //    TFile* f = new TFile ("/var/local/data/root/muSteel_PDfit_20130123_v13.root");


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
    IBPocaEvaluator* processor = IBPocaEvaluator::New(IBPocaEvaluator::LineDistance);

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


    reader->setAcquisitionTime(5);
    std::cout << "There are " << reader->getNumberOfEvents() << " events!\n";
    for(int i=0; i<reader->getNumberOfEvents(); ++i) {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            atc.AddMuon(mu);
            i++;
            if(i%10000 == 0 ) std::cout<<i<<"\n"<<std::flush;
        }
    }

    char file[20];

    atc.Run(1,1);
    sprintf(file, "wtrack_count_tb_20130203.vtk");
    voxels.ExportToVtk(file,0);

    return 0;

}

