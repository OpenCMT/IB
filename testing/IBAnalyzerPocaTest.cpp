
#include <iostream>

#include <TFile.h>
#include <TTree.h>

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
    IBMuonError sigma(11.93,2.03, 18.53,2.05);

    // reader //
    TFile* f = new TFile("/var/local/data/root/run_PDfit_201210/muSteel_PDfit_20121112_v10.root");
    TTree* t = (TTree*)f->Get("n");
    IBMuonEventTTreeReader* reader = IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::R3D_MC);
    reader->setTTree(t);
    reader->setError(sigma);
    reader->setMomentum(0.7);
    reader->selectionCode(IBMuonEventTTreeR3DmcReader::Side2Side);

    // voxels //
    IBVoxel zero = {0.1E-6,0,0};
    IBVoxCollectionCap voxels(Vector3i(160,96,60));
    voxels.SetSpacing(Vector3f(5,5,5));
    voxels.SetPosition(Vector3f(-400,-300,-150));
    voxels.InitLambda(zero);


    IBPocaEvaluator* processor = IBPocaEvaluator::New(IBPocaEvaluator::TiltedAxis);

    IBAnalyzerPoca ap;
    ap.SetPocaAlgorithm(processor);
    ap.SetVoxCollection(&voxels);


    //uLibVtkViewer v_iewer;
    //vtkStructuredGrid v_grid(voxels);
    //v_iewer.AddAbstractProp(v_grid);

    char file_name[10];
    std::cout << "There are " << reader->getNumberOfEvents() << " events!\n" << std::flush;


    int tot = 0;
    int tot2 = 0;
    do {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            //std::cout << "Event : " << reader->getCurrentPosition() << "\n";
            ap.AddMuon(mu);
            //vtkMuonScatter v_event(mu);
            //v_event.AddPocaPoint(processor->getPoca());
            //v_iewer.AddAbstractProp(v_event);
            //v_event.PrintSelf(std::cout);
            //sprintf(file_name,"muon_event_dump10_%d.vtp",i);
            //v_event.SaveToXMLFile(file_name);

            tot++;
            if(tot%1000 == 0 ) std::cout<<tot<<"\n"<<std::flush;
        }
        tot2++;
    } while (tot<96000);

    ap.Run();
    std::cout << "There are " << tot << " event actually read\n";

    voxels.ExportToVtk("latrone_orizzontale_ss.vtk",0);

    return 0;

}
