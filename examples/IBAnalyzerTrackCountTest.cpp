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




using namespace uLib;


int main(int argc, char *argv[]) {

    // errors //
//    IBMuonError sigma(11.93,2.03, 18.53,2.05);
    IBMuonError sigma(13,13);

    // reader //
    TFile* f = new TFile ("/mnt/mu-tom1/geant4/fitPD/20131122_production_v18/muSteel_PDfit_20131203_0_v18.root");
    //    TFile* f = new TFile ("/var/local/data/root/muSteel_PDfit_20130210_0_v15.root");

    //    TFile* f = new TFile ("/var/local/data/root/ROC_sets/201301/20130122/muSteel_PDfit_20130122_1_v12.root");
    //    TFile* f = new TFile ("/var/local/data/root/muSteel_PDfit_20130129_v13.root");
    //    TFile* f = new TFile ("/var/local/data/root/muSteel_PDfit_20130123_v13.root");


    TTree* t = (TTree*)f->Get("n");
    IBMuonEventTTreeR3DmcReader *reader = new IBMuonEventTTreeR3DmcReader();


    //    TTree* t = (TTree*)f->Get("n");
    //    IBMuonEventTTreeReader* reader = IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::R3D_MC);
    reader->setTTree(t);
    reader->setError(sigma);
    reader->setMomentum(0.7);
    reader->selectionCode(0xff);

    // voxels //
    IBVoxel zero = {0,0,0};
    IBVoxel air = {0.1E-6,0,0};

    float vox_size = 5;
    Vector3f vox_bounding(700,360,300); // centered bounding size //

    IBVoxCollection voxels(Vector3i(vox_bounding(0)/vox_size,
                                    vox_bounding(1)/vox_size,
                                    vox_bounding(2)/vox_size));
    voxels.SetSpacing (Vector3f(vox_size,
                                vox_size,
                                vox_size));
    voxels.SetPosition(Vector3f( - 350,
                                 - 180,
                                 - 150 ));
    voxels.InitLambda(air);

    std::cout << "Voxels count " << voxels.GetDims().prod() << "\n";

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
    IBAnalyzerTrackCount atc;
    atc.SetVoxCollection(&voxels);
    atc.SetPocaAlgorithm(processor);
    atc.SetRayAlgorithm(tracer);
    //    atc.SetVaraiblesAlgorithm(minimizator);


    reader->setAcquisitionTime(5.);
    std::cout << "There are " << reader->getNumberOfEvents() << " events!\n";
    for(int i=0; i<reader->getNumberOfEvents(); ++i) {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            atc.AddMuon(mu);
            if(i++%(reader->getNumberOfEvents()/80) == 0) std::cout << "o" << std::flush;
        }
    }
    std::cout << "\n";

    char file[20];

    atc.Run(1,1);
    sprintf(file, "track_count.vtk");
    voxels /= 1E6;
    voxels.ExportToVtk(file,0);

    return 0;

}

