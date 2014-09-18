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


#include "Vtk/vtkMuonScatter.h"
#include "Vtk/vtkMuonEvent.h"

#include "Vtk/vtkVoxRaytracerRepresentation.h"


using namespace uLib;




int do_iterations(const char *file_in, const char *file_out, float min, float map_a) {
    

    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // errors //
    
    IBMuonError sigma(6.02, 7.07, 1.4);
    /* A_ph = 12.24/6.02 mrad;
     * A_th = 18.85/7.07 mrad
     * for DTChambers/CMSChambers
     */
//    Voxel zero = {0};
//    IBLightCollection scraps(Vector3i(120,50,50));
//    scraps.SetSpacing (Vector3f(5,5,5));
//    scraps.SetPosition(Vector3f(-300,-125,-125));
//    scraps.InitVoxels(zero);
//    for(int x=10; x < 130; ++x) {
//        for (int y=10; y < 62; ++y) {
//            for (int z=4; z<56; ++z) {
//                Vector3i id(x,y,z);
//                scraps[id].Value = 0.005127;
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
    
    if (f->IsZombie()) {
        std::cerr << "file not found!\n";
        exit(1);
    }
    
    TTree* t = (TTree*)f->Get("n");
    IBMuonEventTTreeLNLmcReader *reader = new IBMuonEventTTreeLNLmcReader();
    reader->setTTree(t);
    reader->setError(sigma);
    reader->setMomentum(0.7);
//    reader->readPguess(true);
//    reader->selectionCode(IBMuonEventTTreeR3DmcReader::All);
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // voxels //
    IBVoxel air = {0.1E-6,0,0};
    IBVoxCollection voxels(Vector3i(61,32,48));
    voxels.SetSpacing (Vector3f(5,5,5));
    voxels.SetPosition(Vector3f(-152.5,-171.8,-120));
    voxels.InitLambda(air);
    voxels.ExportToVtk("base.vtk",0);
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // ALGORITHMS //
    
    // poca //
    IBPocaEvaluator* processor =
            IBPocaEvaluator::New(IBPocaEvaluator::LineDistance);
    
    // tracer //
    IBVoxRaytracer tracer(voxels);
    
    
    
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // acquisition //
    
    reader->setAcquisitionTime(min);
    std::cout << "There are " << reader->getNumberOfEvents() << " events!\n";
    int tot=0;
    
    IBMuonCollection muons;
    int ev = reader->getNumberOfEvents();
    for (int i=0; i<2000; i++) {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            muons.AddMuon(mu);
            tot++;
        }
    }
    
    muons.PrintSelf(std::cout);
            
    //vtkVoxRaytracerRepresentation v_rt(tracer);

    char file[100];
    for (int i=0; i<100; ++i)
    {
        MuonScatterData mu = muons.At(i);
        vtkMuonScatter v_mu(mu);
//        if (processor->evaluate(mu))
//            v_mu.AddPocaPoint(processor->getPoca());
//        else
//            std::cout <<"error poca\n";
//        v_rt.SetMuon(v_mu);

        sprintf(file,"%s_ray_%i.vtp",file_out,i);
//        v_rt.SetRepresentationElements(vtkVoxRaytracerRepresentation::RayElements);
//        if(v_rt.GetPolyData())
            v_mu.SaveToXMLFile(file);
//        else
//            std::cout << " PROBLEMA RAY \n";

//        sprintf(file,"%s_vox_%i.vtp",file_out,i);
//        v_rt.SetRepresentationElements(vtkVoxRaytracerRepresentation::VoxelsElements);
//        v_rt.SaveToXMLFile(file);
    }
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
        float a;
    } parameters = {
        "/mnt/musteel/var/local/data/root/run_2004x/muRadio_20040.root",
        "tracer_mc",
        5,
        1
    };
    
    if(argc == 5) {
        parameters.file_in  = argv[1];
        parameters.file_out = argv[2];
        parameters.minutes  = atof(argv[3]);
        parameters.a  = atof(argv[4]);
    }
    
    do_iterations(parameters.file_in,
                  parameters.file_out,
                  parameters.minutes,
                  parameters.a);
    return 0;
}


