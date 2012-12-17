
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

using namespace uLib;

int main() {

    // errors //
    IBMuonError sigma(11.93,2.03,18.53,2.05);

    // reader //
    TFile* f = new TFile("/var/local/data/root/run_PDfit_201210/muSteel_PDfit_20121016_v10.root");
    TTree* t = (TTree*)f->Get("n");
    IBMuonEventTTreeReader* reader = IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::R3D_MC);
    reader->setTTree(t);
    reader->setError(sigma);
    reader->setMomentum(1.);
    reader->selectionCode(IBMuonEventTTreeR3DmcReader::All);

    // voxels //
    IBVoxel zero = {0,0,0};
    IBVoxCollectionCap voxels(Vector3i(160,90,60));
    voxels.SetSpacing(Vector3f(5,5,5));
    voxels.SetPosition(Vector3f(-400,-270,-150));
    voxels.InitLambda(zero);

    // filter
    IBVoxFilter_Abtrim abfilt(Vector3i(5,5,5));
    IBFilterGaussShape shape(0.2);
    abfilt.SetKernelSpherical(shape);
    abfilt.SetImage(&voxels);
    abfilt.SetABTrim(0,2);

    // tracer //
    IBVoxRaytracer* tracer = new IBVoxRaytracer(voxels);

    // variables //
    IBMinimizationVariablesEvaluator* minimizator =
            IBMinimizationVariablesEvaluator::New(IBMinimizationVariablesEvaluator::NormalPlane);
    minimizator->setRaytracer(tracer);

    //poca
    IBPocaEvaluator* processor = IBPocaEvaluator::New(IBPocaEvaluator::TiltedAxis);

    IBAnalyzerWPoca* ap = new IBAnalyzerWPoca();
    ap->SetPocaAlgorithm(processor);
    ap->SetVariablesAlgorithm(minimizator);
    ap->SetVoxCollection(&voxels);

    std::cout << "There are " << reader->getNumberOfEvents() << " events!\n" << std::flush;

    int tot  = 0;
    int tot2 = 0;
    do {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            //std::cout << "Event : " << reader->getCurrentPosition() << "\n";
            ap->AddMuon(mu);
            tot++;
            if(tot%10000 == 0 ) std::cout<<tot<<"\n"<<std::flush;
        }
        tot2++;
    } while (tot2<2500000);

    ap->Run();
    delete ap;
    std::cout << "There are " << tot << " event actually read\n";

    voxels.ExportToVtk("20121016_poca_4.vtk",0);

    abfilt.Run();
    voxels.ExportToVtk("20121016_poca_abtrim_3_4.vtk", 0);

    return 0;

}
