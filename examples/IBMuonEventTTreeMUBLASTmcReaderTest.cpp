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

#include "Core/Options.h"

#include "IBMuonEventTTreeReader.h"
#include "IBMuonEventTTreeMUBLASTmcReader.h"

#include "IBVoxCollection.h"
#include "IBVoxRaytracer.h"

#include "IBAnalyzerTrackLengths.h"
#include "IBAnalyzerWTrackLengths.h"

using namespace uLib;

static struct Parameters : Options
{
    char  *file_in;
    char  *file_out;
    float minutes;
    float start_min;
    float vox_size;
    float momentum;

    struct Image
    {
        Vector3f position;
        Vector3f size;
        Vector3f spacing;
        Vector3f origin;
        Vector3f rotation;
    } image;


    Parameters(const char *hello = "Program options") : Options(hello) {
        add_options()
                ("help", "printout help")

                // GENERAL //
                ("momentum",      &momentum,  (float)0.,    "momentum [GeV]")

                // IMAGE //
                ("image.size",    &image.size,     Vector3f(300,500,300),"image bounding size [cm]")
                ("image.spacing", &image.spacing,  Vector3f(5,5,5),"image spacing size [cm]")
                ("image.position",&image.position, Vector3f(-150,-250,-150),"image position [cm]")
                ("image.origin",  &image.origin,   Vector3f(0,0,0),"image origin [cm]")
                ("image.rotation",&image.rotation, Vector3f(0,0,0),"image YZY rotation [rad]")
                ;
    }
} p;   // <-- INSTANCE //



int main(int argc, char **argv) {

    if(argc >= 6) {
        p.file_in  = argv[1];
        p.file_out = argv[2];
        p.minutes  = atof(argv[3]);
        p.start_min  = atof(argv[4]);
        p.vox_size  = atof(argv[5]);
    }
    else {
        std::cerr << "Error in parameters try --help\n";
        exit(1);
    }

    p.parse_command_line(argc,argv);
    p.parse_config_file("cmt.config");
    p.parse_command_line(argc,argv);


    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // voxels //
    IBVoxel air = {0.1E-6,0,0};
    IBVoxel zero = {0,0,0};
    IBVoxCollection voxels;
    {
        if(p.vox_size > 0) {
            voxels.SetDims(Vector3i(p.image.size(0)/p.vox_size,
                                    p.image.size(1)/p.vox_size,
                                    p.image.size(2)/p.vox_size));
            voxels.SetSpacing (Vector3f(p.vox_size,
                                        p.vox_size,
                                        p.vox_size));
        }
        else {
            voxels.SetDims(Vector3i(p.image.size(0)/p.image.spacing(0),
                                    p.image.size(1)/p.image.spacing(1),
                                    p.image.size(2)/p.image.spacing(2)));
            voxels.SetSpacing(p.image.spacing);
        }
        voxels.SetPosition(p.image.position);
        voxels.SetOrigin(p.image.origin);
        voxels.EulerYZYRotate(p.image.rotation);
        voxels.InitLambda(zero);
    }
    std::cout << "Voxels count " << voxels.GetDims().prod()
              << " sized: " << voxels.GetSpacing().transpose()
              << "\n";



    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // reader //


    // errors //
    IBMuonError sigma(1,1); // FIX !
    // reader //

    TFile* f = new TFile(argv[1]);

    IBMuonEventTTreeMUBLASTmcReader reader;
    reader.setTFile(f);
    reader.setError(sigma);
    reader.setMomentum(1);
    //    reader.setAcquisitionTime(p.minutes);

    std::cout << "Number of Events: " << reader.getNumberOfEvents() << std::endl;


    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // ALGORITHMS //

    // tracer //
    IBVoxRaytracer* tracer = new IBVoxRaytracer(voxels);

    // analyzer for occupancy //
    IBAnalyzerTrackLengths analyzer;
    analyzer.SetVoxCollection(&voxels);
    analyzer.SetRayAlgorithm(tracer);


    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // acquisition //


    int tot=0;
    IBMuonCollection muons;

    int ev;
    ev = reader.getNumberOfEvents();


    std::cout << "Reading: ";
    for (int i=0; i<ev; i++) {
        MuonScatter mu;
        if( reader.readNext(&mu) )
            analyzer.AddMuon(mu);
        // progres bar //
        if(tot++%(ev/80) == 0) std::cout << "o" << std::flush;
    }

    analyzer.Run(1,1); // actually it does nothing

    voxels.ExportToVtk(p.file_out);


    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////////
    // purge algorithms //

    delete tracer;
}
