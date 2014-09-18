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
#include "IBVoxFilters.h"

#include "IBPocaEvaluator.h"
#include "IBAnalyzerPoca.h"
#include "IBVoxCollectionCap.h"
#include "IBMuonError.h"
#include "Detectors/MuonScatter.h"
#include "IBMuonEventTTreeReader.h"
#include "IBMuonEventTTreeR3DmcReader.h"
//#include "Vtk/vtkMuonScatter.h"
//#include "Vtk/uLibVtkViewer.h"
//#include "Vtk/vtkStructuredGrid.h"

using namespace uLib;

int main() {

    // errors //
    IBMuonError sigma(12.24,18.85);

    // reader //
//    TFile* f = new TFile("/var/local/data/root/run_PDfit_201210/muSteel_PDfit_20121112_v10.root");
    TFile* f = new TFile ("/var/local/data/root/muSteel_PDfit_20130203_v13.root");
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


    IBPocaEvaluator* processor = IBPocaEvaluator::New(IBPocaEvaluator::TiltedAxis);

    IBAnalyzerPoca ap;
    ap.SetPocaAlgorithm(processor);
    ap.SetVoxCollection(&voxels);


    //uLibVtkViewer v_iewer;
    //vtkStructuredGrid v_grid(voxels);
    //v_iewer.AddAbstractProp(v_grid);

    reader->setAcquisitionTime(5);

    char file_name[10];
    std::cout << "There are " << reader->getNumberOfEvents() << " events!\n" << std::flush;

    int tot;
    //    int ev = 1375250;
    int ev = reader->getNumberOfEvents();
    for (int i=0; i<ev; i++) {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            ap.AddMuon(mu);
            tot++;
        }
    }

    ap.Run();
    std::cout << "There are " << tot << " event actually read\n";

    voxels.ExportToVtk("20130203_poca_5min",0);

    return 0;

}
