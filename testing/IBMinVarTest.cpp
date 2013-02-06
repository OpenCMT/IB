#include <omp.h>
#include "TFile.h"
#include "TTree.h"
#include "IBMinimizationVariablesEvaluator.h"
#include "IBMuonEventTTreeReader.h"
#include "IBMuonEventTTreeR3DmcReader.h"
#include "testing-prototype.h"
#include "IBVoxRaytracer.h"
#include "IBVoxCollectionCap.h"

using namespace uLib;

int main() {
    BEGIN_TESTING(IBMinVar);

    // errors //
    IBMuonError sigma(12.24,18.85); // parameters relative to scattering angles NOT measured angles!!

    // reader //
    TFile* f = new TFile
            ("/var/local/data/root/run_PDfit_201210/muSteel_PDfit_20121030_v10.root");
    TTree* t = (TTree*)f->Get("n");
    IBMuonEventTTreeReader* reader = IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::R3D_MC);
    reader->setTTree(t);
    reader->setError(sigma);
    reader->setMomentum(0.7);
    reader->selectionCode(IBMuonEventTTreeR3DmcReader::Top2Bottom);

    // voxels //
    IBVoxel zero = {0.00082E-6,0,0};
    IBVoxCollectionCap voxels(Vector3i(80,45,30));
    voxels.SetSpacing(Vector3f(10,10,10));
    voxels.SetPosition(Vector3f(-400,-270,-150));
    voxels.InitLambda(zero);

    // tracer //
    IBVoxRaytracer* tracer = new IBVoxRaytracer(voxels);

    // variables //
    IBMinimizationVariablesEvaluator* minimizator =
            IBMinimizationVariablesEvaluator::New(IBMinimizationVariablesEvaluator::NormalPlane);
    minimizator->setRaytracer(tracer);
    int tot  = 0;
    int tot2 = 0;
    for (int i = 0; i<500000/*reader->getNumberOfEvents()*/; ++i) {
        MuonScatter event;
        if (reader->readNext(&event)) {
            if(minimizator->evaluate(event)) tot2++;
            tot++;
        }
    }
    delete minimizator;
    std::cout << "Reader has processed " << tot  << " events" << std::endl;
    std::cout << "MinVar has processed " << tot2 << " events" << std::endl;

    END_TESTING
}
