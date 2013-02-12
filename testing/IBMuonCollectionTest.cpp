
#include "IBMuonCollection.h"
#include "IBMuonEventTTreeReader.h"
#include "IBMuonEventTTreeR3DmcReader.h"
#include "IBVoxRaytracer.h"
#include "IBMinimizationVariablesEvaluator.h"
#include "testing-prototype.h"

#include "root/TTree.h"

int main()
{

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


    // tracer //
    IBVoxRaytracer* tracer = new IBVoxRaytracer(voxels);

    // variables //
    IBMinimizationVariablesEvaluator* minimizator =
            IBMinimizationVariablesEvaluator::New(IBMinimizationVariablesEvaluator::NormalPlane);
    minimizator->setRaytracer(tracer);



    reader->setAcquisitionTime(5);
    std::cout << "There are " << reader->getNumberOfEvents() << " events!\n";
    int tot=0;

    IBMuonCollection muons;
    int ev = reader->getNumberOfEvents();
    for (int i=0; i<ev; i++) {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            muons.AddMuon(mu);
//            minimizator->evaluate(mu);
            tot++;
        }
    }

    muons.SetHiPassAngle(0.01);
    std::cout << "After Hipass Cut events: " << muons.size() << "\n";

    muons.SetLowPassAngle(0.01);
    std::cout << "Afetr Lowpass Cut events: " << muons.size() << "\n";






    return 0;
}

