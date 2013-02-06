
#include <iostream>

#include <TFile.h>
#include <TTree.h>

#include "IBMuonError.h"
#include "Detectors/MuonScatter.h"
#include "IBMuonEventTTreeReader.h"
#include "IBMuonEventTTreeR3DmcReader.h"
#include "IBPocaEvaluator.h"
#include "IBMinimizationVariablesEvaluator.h"
#include "IBVoxRaytracer.h"
#include "IBVoxCollectionCap.h"
#include "IBAnalyzerWPoca.h"
#include "IBAnalyzerPoca.h"
#include "IBVoxFilters.h"
#include "IBLineDistancePocaEvaluator.h"


using namespace uLib;

int main() {

    // errors //
//    IBMuonError sigma(11.93,2.03,18.53,2.05);
    IBMuonError sigma(12.24,18.85);

    // reader //
    TFile* f = new TFile ("/var/local/data/root/ROC_sets/201212/20121223/muSteel_PDfit_20121223_1_v11.root");
    TTree* t = (TTree*)f->Get("n");
    IBMuonEventTTreeReader* reader = IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::R3D_MC);
    reader->setTTree(t);
    reader->setError(sigma);
    reader->setMomentum(1.);
    reader->selectionCode(IBMuonEventTTreeR3DmcReader::All);

    // voxels //
    IBVoxel zero = {0,0,0};
    IBVoxCollectionCap voxels(Vector3i(140,72,60));
    voxels.SetSpacing (Vector3f(5,5,5));
    voxels.SetPosition(Vector3f(-350,-180,-150));
    voxels.InitLambda(zero);

    IBVoxCollectionCap voxels2(voxels);
    // filter
    IBVoxFilter_Abtrim abfilt(Vector3i(5,5,5));
    IBFilterGaussShape shape(0.2);
    abfilt.SetKernelSpherical(shape);
    abfilt.SetImage(&voxels);
    abfilt.SetABTrim(0,0);

    // tracer //
    IBVoxRaytracer* tracer = new IBVoxRaytracer(voxels);

    // variables //
    IBMinimizationVariablesEvaluator* minimizator =
            IBMinimizationVariablesEvaluator::New(IBMinimizationVariablesEvaluator::NormalPlane);
    minimizator->setRaytracer(tracer);

    //poca
    IBPocaEvaluator* processor = IBPocaEvaluator::New(IBPocaEvaluator::TiltedAxis);
    IBLineDistancePocaEvaluator* proc2 = new IBLineDistancePocaEvaluator();

    IBAnalyzerPoca* ap = new IBAnalyzerPoca();
    ap->SetPocaAlgorithm(processor);
    //ap->SetVariablesAlgorithm(minimizator);
    ap->SetVoxCollection(&voxels);

    IBAnalyzerPoca* ap2 = new IBAnalyzerPoca();
    ap2->SetPocaAlgorithm(proc2);
    //ap2->SetVariablesAlgorithm(minimizator);
    ap2->SetVoxCollection(&voxels2);

    std::cout << "There are " << reader->getNumberOfEvents() << " events!\n" << std::flush;

    int tot  = 0;
    int tot2 = 0;
    do {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            //std::cout << "Event : " << reader->getCurrentPosition() << "\n";
            ap->AddMuon(mu);
            ap2->AddMuon(mu);
            tot++;
            if(tot%10000 == 0 ) std::cout<<tot<<"\n"<<std::flush;
        }
        tot2++;
    } while (tot2<2500000);

    ap->Run();
    ap2->Run();
    delete ap;
    delete ap2;
    std::cout << "There are " << tot << " event actually read\n";

    voxels.ExportToVtk("poca_algorithm_noW_TAPE.vtk",0);
    voxels2.ExportToVtk("poca_algorithm_noW_LDPE.vtk",0);

    abfilt.Run();
    voxels.ExportToVtk("Trimmed_TAPE.vtk", 0);

    abfilt.SetImage(&voxels2);
    abfilt.Run();
    voxels2.ExportToVtk("Trimmed_LDPE.vtk",0);
    return 0;

}
