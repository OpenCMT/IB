#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TTree.h>


#include "IBPocaEvaluator.h"
#include "IBAnalyzerEM.h"
#include "IBAnalyzerEMTrim.h"
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

static struct Parameters {
    float a;
    float b;
} global_parameters = { 0,0 };


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
    // filter //
    IBVoxFilter_Abtrim trim(Vector3i(3,3,3));
    IBFilterGaussShape shape(0.3);
    trim.SetKernelWeightFunction(shape);
    trim.SetABTrim(0,1);
    trim.SetImage(&voxels);

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
    IBAnalyzerEMTrim *aem = new IBAnalyzerEMTrim(voxels);
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


    int it   = 50;
    int drop = 25;

    char file[100];

    IBAnalyzerEMAlgorithmSGA_PX ml_algorithm;
    aem->SetMLAlgorithm(&ml_algorithm);



    //    aem->SijCut(60);
    aem->SijGuess(g_Sij_guess_tp_pairs);


    aem->Run(1, 1, global_parameters.a, global_parameters.b); //
    voxels.InitLambda(air);


    std::cout << "SGA PX ----------------------- A="
              << global_parameters.a << " B=" << global_parameters.b << "\n";
    for (int i=1; i<=it; ++i) {
        if(i<3)
            aem->Run(drop,1,0,global_parameters.b);
        else
            aem->Run(drop, 1, global_parameters.a, global_parameters.b);
        //        trim.Run();
        sprintf(file, "%s_%i.vtk",file_out, i*drop);
        //        if(i%5 == 0)
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
        float start_min;
        float vox_size;
        float a;
        float b;
    } parameters = {
        "/var/local/data/root/run_1388/r1388_x.root",
        "r1388_2020_t1_v10",
        1,
        0,
        10,
        0.2, // a
        0.2  // b
    };


    if(argc >= 6) {
        parameters.file_in  = argv[1];
        parameters.file_out = argv[2];
        parameters.minutes  = atof(argv[3]);
        parameters.start_min  = atof(argv[4]);
        parameters.vox_size  = atof(argv[5]);
    }

    if(argc >= 8) {
        parameters.a = atof(argv[6]);
        parameters.b = atof(argv[7]);
    }

    IB::Version::PrintSelf(std::cout);


    global_parameters.a = parameters.a;
    global_parameters.b = parameters.b;

    g_Sij_guess_tp_pairs.resize(3);
    g_Sij_guess_tp_pairs[0] = Vector2f(0,  0.8651);
    g_Sij_guess_tp_pairs[1] = Vector2f(3,  0.474118);
    g_Sij_guess_tp_pairs[2] = Vector2f(30, 0.25381);


    do_iterations(parameters.file_in,
                  parameters.file_out,
                  parameters.minutes,
                  parameters.start_min,
                  parameters.vox_size);


    return 0;
}
