#include <iostream>

#include <TFile.h>
#include <TTree.h>


#include "IBPocaEvaluator.h"
#include "IBAnalyzerEM.h"
#include "IBAnalyzerEMAlgorithm.h"
#include "IBVoxRaytracer.h"
#include "IBMinimizationVariablesEvaluator.h"
#include "IBVoxCollection.h"
#include "IBMuonError.h"
#include "Detectors/MuonScatter.h"
#include "IBMuonEventTTreeReader.h"
#include "IBMuonEventTTreeR3DmcReader.h"
#include "IBVoxFilters.h"
#include "IBAnalyzerWPoca.h"
#include "IBMAPUpdateDensityAlgorithms.h"

#include <Vtk/uLibVtkViewer.h>
#include <Vtk/vtkMuonEvent.h>
#include <Vtk/vtkStructuredGrid.h>

#include "IBAnalyzerEMAlgorithmSGA.h"
#include "IBAnalyzerEMAlgorithmMGA.h"

#include "IBMuonCollection.h"

using namespace uLib;

int main(int argc, char** argv) {

    // errors //
    Voxel zero = {0};
    IBLightCollection scraps(Vector3i(140,72,60));
    scraps.SetSpacing (Vector3f(5,5,5));
    scraps.SetPosition(Vector3f(-350,-180,-150));
    scraps.InitVoxels(zero);

    for(int x=10; x < 130; ++x) {
        for (int y=10; y < 62; ++y) {
            for (int z=4; z<56; ++z) {
                Vector3i id(x,y,z);
                scraps[id].Value = 0.005127;
            }
        }
    }
    IBMuonError sigma(12.24,
                      18.85,
                      1.4); // not used if scraps image is set //
    sigma.setScrapsImage(scraps);
    sigma.azimuthalMomentumCorrection(true);


    // reader //
    TFile* f = new TFile (argv[1]);
//    TFile* f = new TFile ("/var/local/data/root/ROC_sets/201301/20130111/muSteel_PDfit_20130111_10_v12.root");
//    TFile* f = new TFile ("/var/local/data/root/ROC_sets/201212/20121223/muSteel_PDfit_20121223_10_v11.root");

    if (f->IsZombie()) {
        std::cerr << "file not found!\n";
        exit(1);
    }

    TTree* t = (TTree*)f->Get("n");
    IBMuonEventTTreeR3DmcReader *reader = new IBMuonEventTTreeR3DmcReader();
    reader->setTTree(t);
    reader->setError(sigma);
    reader->setMomentum(0.7);
    reader->selectionCode(IBMuonEventTTreeR3DmcReader::All);

    // voxels //
    IBVoxel air = {0.1E-6,0,0};
    IBVoxCollection voxels(Vector3i(140,72,60));
    voxels.SetSpacing (Vector3f(5,5,5));
    voxels.SetPosition(Vector3f(-350,-180,-150));
    voxels.InitLambda(air);

    // MAP Algorithm //
    //    IBMAPPriorTotalWeigth weight_MAP(0.19, 300E-6 * 300E-6);
    //    voxels.SetMAPAlgorithm(&weight_MAP);

    // poca //
    IBPocaEvaluator* processor = IBPocaEvaluator::New(IBPocaEvaluator::LineDistance);

    // tracer //
    IBVoxRaytracer* tracer = new IBVoxRaytracer(voxels);

    // variables //
    IBMinimizationVariablesEvaluator* minimizator =
            IBMinimizationVariablesEvaluator::New(IBMinimizationVariablesEvaluator::NormalPlane);
    minimizator->setRaytracer(tracer);

    // ML Algorithm //
    IBAnalyzerEMAlgorithmSGA *ml_algorithm = new IBAnalyzerEMAlgorithmSGA_PXTZ;

    // analyzer //
    IBAnalyzerEM* aem = new IBAnalyzerEM(voxels);
    aem->SetMLAlgorithm(ml_algorithm);
    aem->SetPocaAlgorithm(processor);
    aem->SetRayAlgorithm(tracer);
    aem->SetVarAlgorithm(minimizator);

    // filter //
    IBVoxFilter_Abtrim trim(Vector3i(3,3,3));
    IBFilterGaussShape shape(0.13);
    trim.SetKernelSpherical(shape);
    trim.SetABTrim(0,1);
    trim.SetImage(&voxels);


    reader->setAcquisitionTime(5);
    std::cout << "There are " << reader->getNumberOfEvents() << " events!\n";
    int tot=0;

    TFile v("p_stat_scaleImg_tiltC.root","RECREATE");
    TH1F h("p_in","p_in",500,0,5);
    TH1F g("p_out", "p_out", 500,0,5);

    IBMuonCollection muons;
    int ev = reader->getNumberOfEvents();
    for (int i=0; i<ev; i++) {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            muons.AddMuon(mu);
            tot++;
            h.Fill(mu.GetMomentum());
            g.Fill(mu.GetMomentumPrime());
        }
    }
    v.cd();
    h.Write();
    g.Write();
    v.Close();
    exit(1);
    std::cout << "Reader events: " << tot << "\n";

//    muons.SetHiPassAngle(0);
//    for(int i=0; i<muons.size(); ++i)
//    {
//        MuonScatter mu = muons.At(i);
//        mu.SetMomentum(0.4);
//        muons[i] = mu;
//    }
//    muons.SetLowPassAngle(0.005);
//    for(int i=0; i<muons.size(); ++i)
//    {
//        MuonScatter mu = muons.At(i);
//        mu.SetMomentum(0.9);
//        muons[i] = mu;
//    }
//    muons.SetHiPassAngle(0);

    muons.PrintSelf(std::cout);
    aem->SetMuonCollection(&muons);


    //    aem->AddVoxcollectionShift(Vector3f(2,0,0));
    //    aem->AddVoxcollectionShift(Vector3f(0,0,2));

    char file[100];

    int it   = 20;
    int drop = 20;

    aem->SijCut(60);
    std::cout << "Spared: [" << aem->Size() << "]\n";
    voxels.InitLambda(air);

    // SGA //
    std::cout << "SGA PXTZ\n";
    for (int i=1; i<=it; ++i) {
        aem->Run(drop,1);
        sprintf(file, "20130203_PXTZ_%i.vtk", i*drop);
        voxels.ExportToVtk(file,0);

        //        IBVoxCollection filtered = voxels;
        //        trim.SetImage(&filtered);
        //        trim.Run();
        //        sprintf(file, "20130203_PXTZ_p14_trim_%i.vtk", i*drop);
        //        filtered.ExportToVtk(file,0);
    }







    delete aem;
    delete minimizator;
    return 0;




}

