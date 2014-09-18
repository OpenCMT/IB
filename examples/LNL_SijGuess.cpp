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
#include <fstream>

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

#include "IBSubImageGrabber.h"

#include "IBAnalyzerEMAlgorithmSGA.h"
#include "IBAnalyzerEMAlgorithmMGA.h"

#include "IBMuonCollection.h"
#include "IB.h"



using namespace uLib;

////////////////////////////////////////////////////////////////////////////////
// EXPERIMENT PARAMETERS //

Vector<Vector2f> g_Sij_guess_tp_pairs;



////////////////////////////////////////////////////////////////////////////////







Vector2f do_iterations(const char *file_in,
                       const char *file_out,
                       float min,
                       float startmin=0,
                       float vox_size=10)
{


    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // errors //

    IBMuonError sigma(6.02, 7.07);


    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // reader //

    TFile* f = new TFile (file_in);
    IBMuonEventTTreeReader *reader = IBMuonEventTTreeReader::New(f);
    reader->setTFile(f);
    reader->setError(sigma);
    reader->setMomentum(1);

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // voxels //

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

    IBVoxel air = {0.1E-6,0,0};
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


    // analyzer //
    //    test::Analyzer* test_em = new test::Analyzer(voxels);
    IBAnalyzerEM *aem = new IBAnalyzerEM(voxels);
    aem->init_properties();
    aem->SetPocaAlgorithm(processor);
    aem->SetRayAlgorithm(tracer);
    aem->SetVarAlgorithm(minimizator);



    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // acquisition //

    if(!(min == 0.0))
        reader->setAcquisitionTime(min);
    reader->setStartTime(startmin);

    std::cout << "start reading file: " << file_in << "\n";
    std::cout << "There are " << reader->getNumberOfEvents() << " events!\n";
    int tot=0;

    IBMuonCollection muons;
    int ev = reader->getNumberOfEvents();

    std::cout << "Reading: ";
    for (int i=0; i<ev; i++) {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            muons.AddMuon(mu);
        }
        if(tot++%(ev/80) == 0) std::cout << "o" << std::flush;
    }
    std::cout << "\n";


    muons.PrintSelf(std::cout);
    aem->SetMuonCollection(&muons);

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // ITERATIONS //


    int it   = 300;
    int drop = 5;

    char file[100];

    IBAnalyzerEMAlgorithmSGA_PX ml_algorithm;
    aem->SetMLAlgorithm(&ml_algorithm);



    //    aem->SijCut(60);
    aem->SijGuess(g_Sij_guess_tp_pairs);

    voxels.InitLambda(air);
    std::cout << "SGA PX ------------------------ \n";
    for (int i=1; i<=it; ++i) {
        aem->Run(drop,1);
        //trim.Run();
	if(it%5 == 0) 
        sprintf(file, "%s_%i.vtk",file_out, i*drop);
        voxels.ExportToVtk(file,0);
        //sprintf(file, "%s_%i.vti",file_out, i*drop);
        //voxels.ExportToVtkXml(file,0);
    }


    delete aem;
    delete minimizator;
}







////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
// MAIN //

//RUN 20130530
// FIX ALPHA PT
//PT: RESULTS LATEST AND GREATEST muon momentum for N1*<Sij> < Sij < N2*<Sij>:
//N1=0, N2=3:    <1/p2>=1.33619, <p>=0.8651
//N1=3, N2=30:   <1/p2>=4.44864, <p>=0.474118
//N1=30, N2=100: <1/p2>=15.5233, <p>=0.25381

int main(int argc, char **argv) {
    struct {
        char  *file_in;
        char  *file_out;
        float minutes;
        float startmin;
        float vox_size;
    } parameters = {
        "NULL",
        "NULL",
        1,
        0,
        3
    };

    if(argc == 6) {
        parameters.file_in  = argv[1];
        parameters.file_out = argv[2];
        parameters.minutes  = atof(argv[3]);
        parameters.startmin = atof(argv[4]);
        parameters.vox_size  = atof(argv[5]);
    }

    IB::Version::PrintSelf(std::cout);

    g_Sij_guess_tp_pairs.resize(3);
    g_Sij_guess_tp_pairs[0] = Vector2f(0,  0.8651);
    g_Sij_guess_tp_pairs[1] = Vector2f(3,  0.474118);
    g_Sij_guess_tp_pairs[2] = Vector2f(30, 0.25381);

    do_iterations(parameters.file_in,
                  parameters.file_out,
                  parameters.minutes,
                  parameters.startmin,
                  parameters.vox_size);


    return 0;
}
