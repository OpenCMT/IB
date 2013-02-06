#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include "IBPocaEvaluator.h"
#include "IBAnalyzerEM.h"
#include "IBVoxRaytracer.h"
#include "IBMinimizationVariablesEvaluator.h"
#include "IBVoxCollectionCap.h"
#include "IBMuonError.h"
#include "IBMuonEventTTreeR3DmcReader.h"
#include "IBVoxFilters.h"
#include "IBAnalyzerWPoca.h"

using namespace uLib;

int main(int argc, char* argv[]) {

    if (argc!=3) {
        Printf("Usage: IBBulkBuild infile.root outfile.vtk\n");
        exit(0);
    }
    char* infile  = argv[1];
    char* outfile = argv[2];
    // errors //
    IBMuonError sigma(12.24,18.85);

    // reader //
    TFile* f = new TFile(infile);
    TTree* t = (TTree*)f->Get("n");
    IBMuonEventTTreeReader* reader = IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::R3D_MC);
    reader->setTTree(t);
    reader->setError(sigma);
    reader->setMomentum(0.7);
    reader->selectionCode(IBMuonEventTTreeR3DmcReader::All);

    // voxels //
    IBVoxel zero = {0.1E-6,0,0};
    IBVoxCollectionCap voxels(Vector3i(140,72,60));
    //3cm: 233 120 100
    //5cm:140 72 60;
    // 18: 70;36;30
    //small3= 60,40,60;
    //BIG = 30 45 30
    //Portal= 80, 45, 30;
    //Fat = 190 58 90
    voxels.SetSpacing(Vector3f(5,5,5));
    voxels.SetPosition(Vector3f(-350,-180,-150));
    //18: -350;-180;-150
    //small= -90 -60 -90;
    //BIG= -150 -270 -150
    //Size= -400, -270, -150

    voxels.InitLambda(zero);

    // poca //
    IBPocaEvaluator* processor = IBPocaEvaluator::New(IBPocaEvaluator::TiltedAxis);

    // tracer //
    IBVoxRaytracer* tracer = new IBVoxRaytracer(voxels);

    // variables //
    IBMinimizationVariablesEvaluator* minimizator =
            IBMinimizationVariablesEvaluator::New(IBMinimizationVariablesEvaluator::NormalPlane);
    minimizator->setRaytracer(tracer);


    // poca analyzer //
    IBAnalyzerWPoca* ap = new IBAnalyzerWPoca();
    ap->SetPocaAlgorithm(processor);
    ap->SetVariablesAlgorithm(minimizator);
    ap->SetVoxCollection(&voxels);

    // analyzer //
    IBAnalyzerEM* aem = new IBAnalyzerEM;
    aem->SetVoxCollection(&voxels);
    aem->SetPocaAlgorithm(processor);
    aem->SetRaytracer(tracer);
    aem->SetVariablesAlgorithm(minimizator);

    int toBread = 2754500;// 1377250; // 2.7545M = 10 min
    for (int i=0; i<toBread; i++) {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            aem->AddMuon(mu);
        }
    }

    aem->SijCut(60);

    aem->Run(200,1);

    voxels.ExportToVtk(outfile,0);



    delete aem;
    delete minimizator;
    return 0;
}
