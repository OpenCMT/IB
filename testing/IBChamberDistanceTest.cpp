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
#include "IBMuonEventTTreeLNLdataReader.h"
#include "IBMuonEventTTreeLNLmcReader.h"
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

    TFile * output = new TFile("chamberdistance.root", "RECREATE");
    TTree * teeeer = new TTree("cd","cd");
    float x_in, y_in, z_in, x_out, y_out, z_out;
    float sx_in, sz_in, sx_out, sz_out;
    float pxi, pzi;

    teeeer->Branch("x_in",   &x_in,   "x_in/F");
    teeeer->Branch("y_in",   &y_in,   "y_in/F");
    teeeer->Branch("z_in",   &z_in,   "z_in/F");
    teeeer->Branch("x_out",  &x_out,  "x_out/F");
    teeeer->Branch("y_out",  &y_out,  "y_out/F");
    teeeer->Branch("z_out",  &z_out,  "z_out/F");
    teeeer->Branch("sx_in",  &sx_in,  "sx_in/F");
    teeeer->Branch("sz_in",  &sz_in,  "sz_in/F");
    teeeer->Branch("sx_out", &sx_out, "sx_out/F");
    teeeer->Branch("sz_out", &sz_out, "sz_out/F");
    teeeer->Branch("pxi",    &pxi,    "pxi/F");
    teeeer->Branch("pzi",    &pzi,    "pzi/F");

    // errors //
    Voxel zero = {0};

    IBMuonError sigma(6.02,
                      7.07); // not used if scraps image is set //

    // reader //
    TFile* f = new TFile (argv[1]);
//    TFile* f = new TFile ("/var/local/data/root/ROC_sets/201301/20130111/muSteel_PDfit_20130111_10_v12.root");
//    TFile* f = new TFile ("/var/local/data/root/ROC_sets/201212/20121223/muSteel_PDfit_20121223_10_v11.root");

    if (f->IsZombie()) {
        std::cerr << "file not found!\n";
        exit(1);
    }

    IBMuonEventTTreeLNLdataReader *reader = new IBMuonEventTTreeLNLdataReader();
    reader->setTFile(f);
    reader->setError(sigma);
    reader->setMomentum(0.7);

    // voxels //
    IBVoxel air = {0.1E-6,0,0};
    IBVoxCollection voxels(Vector3i(140,72,60));
    voxels.SetSpacing (Vector3f(5,5,5));
    voxels.SetPosition(Vector3f(-350,-180,-150));
    voxels.InitLambda(air);

    // poca //
    IBPocaEvaluator* processor = IBPocaEvaluator::New(IBPocaEvaluator::LineDistance);

    // tracer //
    IBVoxRaytracer* tracer = new IBVoxRaytracer(voxels);

    // variables //
    IBMinimizationVariablesEvaluator* minimizator =
            IBMinimizationVariablesEvaluator::New(IBMinimizationVariablesEvaluator::SimpleTwoViews);
    minimizator->setRaytracer(tracer);

    reader->setAcquisitionTime(5);
    std::cout << "There are " << reader->getNumberOfEvents() << " events!\n";
    int tot=0;

    int ev = reader->getNumberOfEvents();
    for (int i=0; i<ev; i++) {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            x_in = mu.LineIn().origin(0);
            y_in = mu.LineIn().origin(1);
            z_in = mu.LineIn().origin(2);

            x_out = mu.LineOut().origin(0);
            y_out = mu.LineOut().origin(1);
            z_out = mu.LineOut().origin(2);

            sx_in = mu.LineIn().direction(0);
            sz_in = mu.LineIn().direction(2);

            sx_out = mu.LineOut().direction(0);
            sz_out = mu.LineOut().direction(2);

            pxi = x_in - sx_in*(y_out-y_in);
            pzi = z_in - sz_in*(y_out-y_in);

            teeeer->Fill();

            tot++;
        }
    }

    output->cd();
    teeeer->Write();
    output->Close();
    std::cout << "Reader events: " << tot << "\n";

    delete minimizator;
    return 0;

}

