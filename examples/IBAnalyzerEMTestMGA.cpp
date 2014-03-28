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

using namespace uLib;


std::string GetFileExtension(const std::string& FileName)
{
    if(FileName.find_last_of(".") != std::string::npos)
        return FileName.substr(FileName.find_last_of(".")+1);
    return "";
}

std::string FileNameRemoveExtension(const std::string& FileName)
{
    std::string file;

    // remove path //
    if(FileName.find_last_of("/") != std::string::npos)
        file = FileName.substr(FileName.find_last_of("/")+1);
    // remove extension //
    if(file.find_last_of(".") != std::string::npos)
        return file.substr(0,file.find_last_of("."));
    return "";
}



int run_set(float ps, IBAnalyzerEMAlgorithm *ml_algorithm, const char *filename) {

    // errors //
//    IBMuonError sigma(11.93,2.03, 18.53,2.05);
    IBMuonError sigma(12.24, 0.0, 18.85, 0.0, ps );
    // reader //

    TFile *f = new TFile (filename);
    //TFile* f = new TFile ("/var/local/data/root/ROC_sets/201212/20121223/muSteel_PDfit_20121223_1_v11.root");
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


    // poca //
    IBPocaEvaluator* processor = IBPocaEvaluator::New(IBPocaEvaluator::LineDistance);

    // tracer //
    IBVoxRaytracer* tracer = new IBVoxRaytracer(voxels);

    // variables //
    IBMinimizationVariablesEvaluator* minimizator =
            IBMinimizationVariablesEvaluator::New(IBMinimizationVariablesEvaluator::NormalPlane);
    minimizator->setRaytracer(tracer);

//    // ML Algorithm //
//    IBAnalyzerEMAlgorithmSGA *ml_algorithm = new IBAnalyzerEMAlgorithmSGA_PXTZ;

    // analyzer //
    IBAnalyzerEM* aem = new IBAnalyzerEM(voxels);
    aem->SetMLAlgorithm(ml_algorithm);
    aem->SetPocaAlgorithm(processor);
    aem->SetRayAlgorithm(tracer);
    aem->SetVarAlgorithm(minimizator);


    std::cout << "There are " << reader->getNumberOfEvents() << " events!\n";
    int tot=0;


//    int ev = 1375250; // 5 min
    reader->setAcquisitionTime(5);
    for (int i=0; i<reader->getNumberOfEvents(); i++) {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            aem->AddMuon(mu);
            tot++;
        }
    }


    char file[100];

    int it   = 1;
    int drop = 500;

    aem->SijCut(60);
    std::cout << "Survived: [" << aem->Size() << "]\n";


    // SGA //

    {
        //IBAnalyzerEMAlgorithmSGA_PXTZ ml_algorithm;
        aem->SetMLAlgorithm(ml_algorithm);
        voxels.InitLambda(zero);


        std::cout << "*** PROCESSING: " << filename << " *****************\n";
        for (int i=1; i<=it; ++i) {
            aem->Run(drop,1);
            sprintf(file, "%s_%i.vtk",FileNameRemoveExtension(std::string(filename)).c_str(), i*drop);
            voxels.ExportToVtk(file,0);
        }
    }

    //MGA //

    // 1m scraps
//    float w[5] = { 0.104 , 0.286,  0.340,  0.202,  0.063  };
//    float s[5] = { 0.0006, 0.0027, 0.0082, 0.0259, 0.0718 };

    // 3m scraps
    float w[5] = { 0.093 , 0.257,  0.337,  0.234,  0.076  };
    float s[5] = { 0.0006, 0.0024, 0.0070, 0.0194, 0.0527 };

//    {
//        IBAnalyzerEMAlgorithmMGA_PXTZ<5> ml_algorithm_mga;
//        ml_algorithm_mga.SetGaussians(w,s);
//        aem->SetMLAlgorithm(&ml_algorithm_mga);
//        voxels.InitLambda(zero);

//        std::cout << "MGA PXTZ_5_\n";
//        for (int i=1; i<=it; ++i) {
//            aem->Run(drop,1);
//            sprintf(file, "20121223_MGA_5_EnoCorr_B0AZ24_PXTZ_%i.vtk", i*drop);
//            voxels.ExportToVtk(file,0);
//        }
//    }




    delete aem;
    delete minimizator;
    return 0;

}

int main()
{

    const char *path = "/var/local/data/root/";
    char filename[200];

    sprintf(filename,"%s%s",path,"muSteel_PDfit_20130204_v13.root");
    run_set(1.4,  new IBAnalyzerEMAlgorithmSGA_PXTZ, filename);

    sprintf(filename,"%s%s",path,"muSteel_PDfit_20130205_v13.root");
    run_set(1.4,  new IBAnalyzerEMAlgorithmSGA_PXTZ, filename);

    sprintf(filename,"%s%s",path,"muSteel_PDfit_20130206_v13.root");
    run_set(1.4,  new IBAnalyzerEMAlgorithmSGA_PXTZ, filename);

    sprintf(filename,"%s%s",path,"muSteel_PDfit_20130207_v13.root");
    run_set(1.4,  new IBAnalyzerEMAlgorithmSGA_PXTZ, filename);

    sprintf(filename,"%s%s",path,"muSteel_PDfit_20130208_v13.root");
    run_set(1.4,  new IBAnalyzerEMAlgorithmSGA_PXTZ, filename);


    return 0;

}
