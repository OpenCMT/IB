#include <iostream>
#include <TFile.h>
#include <TTree.h>
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

int main() {

    // errors //
    IBMuonError sigma(12.24,18.85);

    // reader //
    TFile* f = new TFile
            ("/var/local/data/root/ROC_sets/201212/lead/muSteel_PDfit_2012122300_v11.root");
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
    ap->SetVarAlgorithm(minimizator);
    ap->SetVoxCollection(&voxels);

    // analyzer //
    IBAnalyzerEM* aem = new IBAnalyzerEM;
    aem->SetVoxCollection(&voxels);
    aem->SetPocaAlgorithm(processor);
    aem->SetRaytracer(tracer);
    aem->SetVariablesAlgorithm(minimizator);

    int toBread = 2500000;
    for (int i=0; i<toBread; i++) {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            aem->AddMuon(mu);
        }
        if (i%100==0) printf("\rReading Events...%.2f\%", i*100/toBread);
    }

    aem->SijCut(60);

    int it   = 1;
    int drop = 100;

    char file[100];

    std::cout << "PXTZ\n";
    for (int i=1; i<=it; ++i) {
        aem->Run(drop,1);
 //       sprintf(file, "20121224_v495_2PXTZ_%i.vtk", i*drop);
 //       voxels.ExportToVtk(file,0);
    }

    delete aem;
    delete minimizator;
    return 0;
}
