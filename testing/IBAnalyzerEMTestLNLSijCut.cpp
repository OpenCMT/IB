#include <iostream>

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

#include "IBPocaEvaluator.h"
#include "IBAnalyzerEM.h"
#include "IBAnalyzerEMAlgorithm.h"
#include "IBVoxRaytracer.h"
#include "IBMinimizationVariablesEvaluator.h"
#include "IBVoxCollection.h"
#include "IBMuonError.h"
#include "Detectors/MuonScatter.h"
#include "IBMuonEventTTreeReader.h"
#include "IBMuonEventTTreeLNLdataReader.h"
#include "IBMuonEventTTreeLNLmcReader.h"
#include "IBVoxFilters.h"
#include "IBAnalyzerWPoca.h"
#include "IBMAPUpdateDensityAlgorithms.h"

#include "IBAnalyzerEMAlgorithmSGA.h"
#include "IBAnalyzerEMAlgorithmMGA.h"

#include "IBMuonCollection.h"
#include "IB.h"
using namespace uLib;




int do_iterations(const char *file_in,
                  const char *file_out,
                  float min,
                  float start_min=0,
                  float vox_size=10)
{

    IB::Version::PrintSelf(std::cout);
    
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // errors //
    
    IBMuonError sigma(6.02, 7.07);
//    Voxel zero = {0};
//    IBLightCollection scraps(Vector3i(140,72,60));
//    scraps.SetSpacing (Vector3f(5,5,5));
//    scraps.SetPosition(Vector3f(-350,-180,-150));
//    scraps.InitVoxels(zero);
    
//    for(int x=10; x < 130; ++x) {
//        for (int y=10; y < 62; ++y) {
//            for (int z=4; z<56; ++z) {
//                Vector3i id(x,y,z);
//                scraps[id].Value = 1;
//            }
//        }
//    }
//    sigma.setScrapsImage(scraps);
//    sigma.averageMomentumCorrection(true);
//    sigma.azimuthalMomentumCorrection(true);

    
    
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // reader //
    
    TFile* f = new TFile (file_in);
    
//    if (f->IsZombie()) {
//        std::cerr << "file not found!\n";
//        exit(1);
//    }
    
//    TTree* t = (TTree*)f->Get("n");
    IBMuonEventTTreeLNLmcReader *reader = new IBMuonEventTTreeLNLmcReader();
    reader->setTFile(f);
    reader->setError(sigma);
    reader->setMomentum(0.7);    
    
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // voxels //
    IBVoxel air = {0.1E-6,0,0};
//    IBVoxCollection voxels(Vector3i(60,34,48));
//    voxels.SetSpacing (Vector3f(5,5,5));
//    voxels.SetPosition(Vector3f(-150,-182,-120));

    Vector3f vox_bounding(300,161,240); // centered bounding size //

    IBVoxCollection voxels(Vector3i(vox_bounding(0)/vox_size,
                                    vox_bounding(1)/vox_size,
                                    vox_bounding(2)/vox_size));
    voxels.SetSpacing (Vector3f(vox_size,
                                vox_size,
                                vox_size));
    voxels.SetPosition(Vector3f( - 150,
                                 - 172,
                                 - 120 ));

    voxels.InitLambda(air);
    
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // ALGORITHMS //
    
    // poca //
    IBPocaEvaluator* processor =
            IBPocaEvaluator::New(IBPocaEvaluator::LineDistance);
    
    // tracer //
    IBVoxRaytracer* tracer = new IBVoxRaytracer(voxels);
    
    // variables //
    IBMinimizationVariablesEvaluator* minimizator =
            IBMinimizationVariablesEvaluator::
            New(IBMinimizationVariablesEvaluator::NormalPlane);
    minimizator->setRaytracer(tracer);

    // ML Algorithm //
    //IBAnalyzerEMAlgorithmSGA *ml_algorithm = new IBAnalyzerEMAlgorithmSGA_PXTZ;
    
    // analyzer //
    IBAnalyzerEM* aem = new IBAnalyzerEM(voxels);
    //aem->SetMLAlgorithm(ml_algorithm);
    aem->SetPocaAlgorithm(processor);
    aem->SetRayAlgorithm(tracer);
    aem->SetVarAlgorithm(minimizator);
    
    
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // filter //
    IBVoxFilter_Abtrim trim(Vector3i(3,3,3));
    IBFilterGaussShape shape(0.1);
    trim.SetKernelSpherical(shape);
//    IBFilterGaussShape shape(1);
//    trim.SetKernelWeightFunction(shape);
    trim.SetABTrim(0,1);
    trim.SetImage(&voxels);

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // acquisition //
    
    reader->setAcquisitionTime(min);
    reader->setStartTime(start_min);
    std::cout << "There are " << reader->getNumberOfEvents() << " events!\n";
    int tot=0;
    
    IBMuonCollection muons;
    int ev = reader->getNumberOfEvents();
    IBAnalyzerWPoca ap;
    ap.SetPocaAlgorithm(processor);
    ap.SetVoxCollection(&voxels);
    ap.SetVariablesAlgorithm(minimizator);

    for (int i=0; i<ev; i++) {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            muons.AddMuon(mu);
            tot++;
        }
    }        
    muons.PrintSelf(std::cout);
    aem->SetMuonCollection(&muons);

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // ITERATIONS //


    int it   = 10;
    int drop = 50;

    char file[100];

    IBAnalyzerEMAlgorithmSGA_PXT ml_algorithm;
    aem->SetMLAlgorithm(&ml_algorithm);


    f->Close();
    f = new TFile("P_dump.root","RECREATE");
    TH1F *h_all = new TH1F("SijAllP","Sija",1000,0,5000);
    TH1F *h_cut = new TH1F("SijCuttedP","Sijc",1000,0,5000);

    for (int ev=0; ev < aem->Events().size(); ev++ )
        h_all->Fill(aem->Events().at(ev).header.G4P);

    aem->SijCut(60);
    voxels.InitLambda(air);

    // histogram //
    for (int ev=0; ev < aem->Events().size(); ev++ )
        h_cut->Fill(aem->Events().at(ev).header.G4P);

    f->cd();
    h_all->Write();
    h_cut->Write();
    f->Close();

//    std::cout << "SGA PXTZ ------------------------ \n";
//    for (int i=1; i<=it; ++i) {
//        aem->Run(drop,1);
//        //trim.Run();
//        sprintf(file, "%s_%i.vtk",file_out, i*drop);
//        voxels.ExportToVtk(file,0);
//    }

    
    delete reader;
    delete aem;
    delete minimizator;
}







////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// MAIN //


int main(int argc, char **argv) {
    struct {
        char  *file_in;
        char  *file_out;
        float minutes;
        float start_min;
        float vox_size;
    } parameters = {
        "/var/local/data/root/run_1363/muRadio_1363.root",
        "image_1363",
        5,
        0,
        5
    };

    if(argc == 6) {
        parameters.file_in  = argv[1];
        parameters.file_out = argv[2];
        parameters.minutes  = atof(argv[3]);
        parameters.start_min  = atof(argv[4]);
        parameters.vox_size  = atof(argv[5]);
    }
    
    do_iterations(parameters.file_in,
                  parameters.file_out,
                  parameters.minutes,
                  parameters.start_min,
                  parameters.vox_size);



    return 0;
}

