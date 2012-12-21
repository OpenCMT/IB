#include <iostream>

#include <TFile.h>
#include <TTree.h>

#include "IBPocaEvaluator.h"
#include "IBAnalyzerEM.h"
#include "IBVoxRaytracer.h"
#include "IBMinimizationVariablesEvaluator.h"
#include "IBVoxCollectionCap.h"
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

using namespace uLib;

int main() {

    // errors //
    IBMuonError sigma(11.93,2.03, 18.53,2.05);
    // reader //
    TFile* f = new TFile
            ("/var/local/data/root/run_PDfit_201210/muSteel_PDfit_20121215_v11.root");
    TTree* t = (TTree*)f->Get("n");
    IBMuonEventTTreeReader* reader = IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::R3D_MC);
    reader->setTTree(t);
    reader->setError(sigma);
    reader->setMomentum(0.7);
    reader->selectionCode(IBMuonEventTTreeR3DmcReader::All);

    // voxels //
    IBVoxel zero = {0.1E-6,0,0};
//    IBVoxCollectionCap voxels(Vector3i(10,6,4));
//    voxels.SetSpacing(Vector3f(80,80,80));
    IBVoxCollectionCap voxels(Vector3i(60,40,60)); // BIG = 30 45 30 Portal= 80, 45, 30; Fat = 190 58 90
    voxels.SetSpacing(Vector3f(3,3,3));
    voxels.SetPosition(Vector3f(-90,-60,-90));  //BIG= -150 -270 -150 Size= -400, -270, -150
    voxels.InitLambda(zero);
//    for (int i=10; i<70; ++i) {
//        for (int k=7; k<23; ++k) {
//            for (int j=14; j<40; ++j) {
//                Vector3i idx(i,j,k);
//                voxels.SetValue(idx, (4*1.4)*1E-6); // rLen for 0.8g/cm^3 Fe
//            }
//        }
//    }



//    std::cout<<voxels.GetSpacing().transpose()<<std::endl;

//    voxels.ExportToVtk("start.vtk",0);
//    exit(0);





    // filter

    IBVoxFilter_Abtrim abfilt(Vector3i(5,5,5));
    IBFilterGaussShape shape(0.2);
    abfilt.SetKernelSpherical(shape);
    abfilt.SetImage(&voxels);
    abfilt.SetABTrim(0,2);

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


    std::cout << "There are " << reader->getNumberOfEvents() << " events!\n";
    int tot=0;

    for (int i=0; i<2500000/*reader->getNumberOfEvents()*/; i++) {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            ap->AddMuon(mu);
            aem->AddMuon(mu);
            tot++;
        }

    }


    char file[100];

    int before = aem->Size();
    std::cout << "Read " << tot << "Events \n";
    std::cout << "EM Total Events " << before << "\n";

//    aem->Run(1,1);
//    delete aem;
//    exit(0);

    aem->SijCut(60);
//    voxels.InitLambda(zero);

    float after = aem->Size();
    std::cout << "EM Events after cut " << after << "\n";
    std::cout << "Events cut: " << 100*(before-after)/before << "\%" <<std::endl;

    // RUN POCA //
//    ap->Run();
//    //abfilt.Run();
//    voxels.ExportToVtk("20121213_Poca_WphiNtht.vtk",0);
//    exit(0);
    //    voxels *= 0.0406;

    // Lambda MAP //
//    IBMAPPriorGaussianUpdateAlgorithm * lMAP = new IBMAPPriorGaussianUpdateAlgorithm(1E2);
//    IBVoxCollection prior = voxels;
//    lMAP->SetDensityPrior(&prior);
//    voxels.SetMAPAlgorithm(lMAP);



    int it   = 50;
    int drop = 10;

    std::cout << "PXTZ\n";
    for (int i=1; i<=it; ++i) {
        aem->Run(drop,1);
        sprintf(file, "20121215_v488_2PXTZ_Cut60_v11_l100_vx3_p07_Inv_%i.vtk", i*drop);
        voxels.ExportToVtk(file,0);
    }

    abfilt.Run();
    voxels.ExportToVtk("20121215_v488_2PXTZ_Cut60_v11_l100_vx3_p07_Inv_trim.vtk",0);

//    std::cout << "PX\n";
//    aem->parameters().algorithm = IBAnalyzerEM::PX;
//    voxels = prior;

//    for (int i=1; i<=it; ++i) {
//        aem->Run(drop,1);
//        sprintf(file, "20121030_v471_1PX_Cut60_poca_%i.vtk", i*drop);
//        voxels.ExportToVtk(file,0);
//    }

//    std::cout << "TZ\n";
//    aem->parameters().algorithm = IBAnalyzerEM::TZ;
//    voxels = prior;

//    for (int i=1; i<=it; ++i) {
//        aem->Run(drop,1);
//        sprintf(file, "20121030_v471_1TZ_Cut60_poca_%i.vtk", i*drop);
//        voxels.ExportToVtk(file,0);
//    }



//    std::cout << "PT\n";
//    aem->parameters().algorithm = IBAnalyzerEM::PT;
//    //voxels = prior;

//    for (int i=1; i<=it; ++i) {
//        aem->Run(drop,1);
//        sprintf(file, "20121213_v488_2PT_Cut60_v11_l100_vx3_p3G_%i.vtk", i*drop);
//        voxels.ExportToVtk(file,0);
//    }

//    std::cout << "XZ\n";
//    aem->parameters().algorithm = IBAnalyzerEM::XZ;
//    //voxels = prior;

//    for (int i=1; i<=it; ++i) {
//        aem->Run(drop,1);
//        sprintf(file, "20121213_v488_2XZ_Cut60_poca_%i.vtk", i*drop);
//        voxels.ExportToVtk(file,0);
//    }

//    std::cout << "P\n";
//    aem->parameters().algorithm = IBAnalyzerEM::P;
//    //voxels = prior;

//    for (int i=1; i<=it; ++i) {
//        aem->Run(drop,1);
//        sprintf(file, "20121213_v488_1P_Cut60_%i.vtk", i*drop);
//        voxels.ExportToVtk(file,0);
//    }

//    std::cout << "T\n";
//    aem->parameters().algorithm = IBAnalyzerEM::T;
//    //voxels = prior;

//    for (int i=1; i<=it; ++i) {
//        aem->Run(drop,1);
//        sprintf(file, "20121213_v488_1T_Cut60_%i.vtk", i*drop);
//        voxels.ExportToVtk(file,0);
//    }

//    std::cout << "X\n";
//    aem->parameters().algorithm = IBAnalyzerEM::X;
//    voxels = prior;

//    for (int i=1; i<=it; ++i) {
//        aem->Run(drop,1);
//        sprintf(file, "20121030_v471_1X_Cut60_poca_%i.vtk", i*drop);
//        voxels.ExportToVtk(file,0);
//    }

//    std::cout << "Z\n";
//    aem->parameters().algorithm = IBAnalyzerEM::Z;
//    voxels = prior;

//    for (int i=1; i<=it; ++i) {
//        aem->Run(drop,1);
//        sprintf(file, "20121030_v471_1Z_Cut60_poca_%i.vtk", i*drop);
//        voxels.ExportToVtk(file,0);
//    }



////    abfilt.Run();
////    sprintf(file, "20121016_4V_abtrim.vtk");
////    voxels.ExportToVtk(file, 0);
    delete aem;
    delete minimizator;
    return 0;

}

