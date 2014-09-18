/*//////////////////////////////////////////////////////////////////////////////
// CMT Cosmic Muon Tomography project //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

  Copyright (c) 2014, Universita' degli Studi di Padova, INFN sez. di Padova

  Coordinators: Prof. Gianni Zumerle < gianni.zumerle@pd.infn.it >
                Paolo Checchia       < paolo.checchia@pd.infn.it >

  Authors: Andrea Rigoni Garola < andrea.rigoni@pd.infn.it >
           Matteo Furlan        < nuright@gmail.com >
           Sara Vanini          < sara.vanini@pd.infn.it >

  All rights reserved
  ------------------------------------------------------------------

  This file can not be copied and/or distributed without the express
  permission of  Prof. Gianni Zumerle  < gianni.zumerle@pd.infn.it >

//////////////////////////////////////////////////////////////////////////////*/



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





void dump_chi2(IBAnalyzerEM *aem, const char *name, int index, float max_bound = 1000){

    std::cout << "chi2 dump " << name << " \n ";

    // 2 PX //
    {
        char str[100];
        sprintf(str,"%s_PX_%i",name,index);
        TH1F  *histo = new TH1F(str,str,1000,0,max_bound);
        IBAnalyzerEMAlgorithmSGA_PX ml_algorithm;

        float chi2 = 0;

        for (unsigned int i = 0; i< aem->Events().size(); ++i) {
            Matrix4f Sigma = Matrix4f::Zero();
            IBAnalyzerEM::Event &evc = aem->Events().at(i);
            ml_algorithm.ComputeSigma(Sigma,&evc);

            Matrix2f iS;
            {
                Matrix2f S;
                S << Sigma(0,0), Sigma(0,1), Sigma(1,0), Sigma(1,1);
                iS = S.inverse();
            }
            Vector2f Di(evc.header.Di(0),evc.header.Di(1));
//            Matrix2f Dn = iS * (Di * Di.transpose());
            chi2 = Di.transpose() * iS * Di;
            histo->Fill(chi2);
        }
        histo->Write();
        histo->Delete();
    }

    // 2 //
    {
        char str[100];
        sprintf(str,"%s_TZ_%i",name,index);
        TH1F  *histo = new TH1F(str,str,1000,0,max_bound);
        IBAnalyzerEMAlgorithmSGA_TZ ml_algorithm;

        float chi2 = 0;

        for (unsigned int i = 0; i< aem->Events().size(); ++i) {
            Matrix4f Sigma = Matrix4f::Zero();
            IBAnalyzerEM::Event &evc = aem->Events().at(i);
            ml_algorithm.ComputeSigma(Sigma,&evc);

            Matrix2f iS;
            {
                Matrix2f S;
                S << Sigma(2,2), Sigma(2,3), Sigma(3,2), Sigma(3,3);
                iS = S.inverse();
            }
            Vector2f Di(evc.header.Di(2),evc.header.Di(3));
//            Matrix2f Dn = iS * (Di * Di.transpose());
            chi2 = Di.transpose() * iS * Di;
            histo->Fill(chi2);
        }

        histo->Write();
        histo->Delete();
    }

    // 4 //
    {
        char str[100];
        sprintf(str,"%s_PXTZ_%i",name,index);
        TH1F  *histo = new TH1F(str,str,1000,0,max_bound);
        IBAnalyzerEMAlgorithmSGA_PXTZ ml_algorithm;

        float chi2 = 0;

        for (unsigned int i = 0; i< aem->Events().size(); ++i) {
            Matrix4f Sigma = Matrix4f::Zero();
            IBAnalyzerEM::Event &evc = aem->Events().at(i);
            ml_algorithm.ComputeSigma(Sigma,&evc);

            Vector4f Di = evc.header.Di;
            Matrix4f iS = Sigma.inverse();
            chi2 = Di.transpose() * iS * Di;
//            chi2 = Dn.trace();
            histo->Fill(chi2);
        }
        histo->Write();
        histo->Delete();
    }


}






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
    IBMuonEventTTreeLNLdataReader *reader = new IBMuonEventTTreeLNLdataReader();
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
    ap.SetVarAlgorithm(minimizator);

    std::cout << "reading mu from file: ";
    for (int i=0; i<ev; i++) {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            muons.AddMuon(mu);
            if(tot++%(ev/100) == 0) std::cout << "." << std::flush;
        }
    }        
    std::cout << "\n";

    IBAnalyzerEMAlgorithmSGA_PXT ml_algorithm;
    aem->SetMLAlgorithm(&ml_algorithm);

    muons.PrintSelf(std::cout);
    aem->SetMuonCollection(&muons);

//    std::cout << "-- TEST: removing E cross correlation\n";
//    for(int i=0; i< aem->Size(); ++i)
//    {
//        IBAnalyzerEM::Event &evc = aem->Events()[i];
//        //        evc.header.E.block<2,2>(2,0) = Matrix2f::Zero();
//        //        evc.header.E.block<2,2>(0,2) = Matrix2f::Zero();
//        //        evc.header.E(3,0) = 0;
//        //        evc.header.E(3,1) = 0;
//        //        evc.header.E(0,3) = 0;
//        //        evc.header.E(1,3) = 0;
//        //        std::cout << aem->Events()[i].header.E << "\n\n";

//        Matrix4f Sigma;

//        ml_algorithm.ComputeSigma(Sigma,evc);

//        std::cout << " -------------------------------------------------------- \n"
//                  << " -- EVENT: << " << i << "\n";
//        std::cout << " Di: " << evc.header.Di.transpose() << "\n\n"
//                  << " E: \n" << evc.header.E << "\n\n"
//                  << " Sigma: \n" << Sigma << "\n\n";

//    }




    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // CHI2 STUDY FOR CUT //

    TFile *dump_cut_file = new TFile("chi2dump.root","RECREATE");



    int it   = 10;
    int drop = 30;

    char file[100];



    // Sij CUT //
//    {
//        int size = aem->Events().size();
//        std::cout << " ////////////////////////////////////////////////////// \n"
//                  << " // Sij Cut /////////////////////////////////////////// \n"
//                  << " // total events : " << aem->Events().size(); << "\n";
//        aem->SijCut(60);
//        std::cout << " // cutted events: " << size - aem->Events().size() << "\n\n";
//    }


    voxels.InitLambda(air);
    aem->Run(9,1);

    dump_cut_file->cd();
    dump_chi2(aem,"chi2_precut",9,5000);

    // Chi2Cut //
    {
        int size = aem->Events().size();
        std::cout << " ////////////////////////////////////////////////////// \n"
                  << " // Chi2 Cut ////////////////////////////////////////// \n"
                  << " // total events : " << aem->Events().size() << "\n";
        aem->Chi2Cut(300);
        std::cout << " // cutted events: " << size - aem->Events().size() << "\n\n";
    }

    dump_cut_file->cd();
    dump_chi2(aem,"chi2_postcut",10,5000);

    voxels.InitLambda(air);
    std::cout << "SGA PXTZ ------------------------ \n";
    for (int i=1; i<=it; ++i) {
        aem->Run(drop,1);
        //trim.Run();
        sprintf(file, "%s_%i.vtk",file_out, i*drop);
        voxels.ExportToVtk(file,0);
        dump_cut_file->cd();
        dump_chi2(aem,"chi2",i*drop,5000);
    }


    
    dump_cut_file->Close();
    delete dump_cut_file;

    delete aem;
    delete minimizator;
    delete reader;
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
        30,
        0,
        10
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

