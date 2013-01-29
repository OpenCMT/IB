
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
    IBMuonError sigma(11.93,2.03, 18.53,2.05);

    // reader //
    TFile* f = new TFile("/home/eth/mustee/data/muSteel_PDfit_2012122300_v11.root");
    TTree* t = (TTree*)f->Get("n");
    IBMuonEventTTreeReader* reader = IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::R3D_MC);
    reader->setTTree(t);
    reader->setError(sigma);
    reader->setMomentum(0.7);
    reader->selectionCode(IBMuonEventTTreeR3DmcReader::Side2Side);

    // voxels //
    IBVoxel zero = {0.1E-6,0,0};
    IBVoxCollection voxels(Vector3i(140,72,60));
    voxels.SetSpacing (Vector3f(5,5,5));
    voxels.SetPosition(Vector3f(-350,-180,-150));

    voxels.InitLambda(zero);
    IBVoxCollectionCap voxelS(Vector3i(140,72,60));
    voxelS.SetSpacing(Vector3f(5,5,5));
    voxelS.SetPosition(Vector3f(-350,-180,-150));
    voxelS.InitLambda(zero);

    IBPocaEvaluator* processor = IBPocaEvaluator::New(IBPocaEvaluator::TiltedAxis);
    IBPocaEvaluator* processoR = IBPocaEvaluator::New(IBPocaEvaluator::LineDistance);

    IBAnalyzerPoca ap;
    ap.SetPocaAlgorithm(processor);
    ap.SetVoxCollection(&voxels);
    IBAnalyzerPoca aP;
    ap.SetPocaAlgorithm(processoR);
    ap.SetVoxCollection(&voxelS);


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
            ap.AddMuon(mu);
            aP.AddMuon(mu);
            tot++;
            if(tot%1000 == 0 ) std::cout<<tot<<"\n"<<std::flush;
        }
        tot2++;
    } while (tot<900000);

    ap.Run();
    aP.Run();
    std::cout << "There are " << tot << " event actually read\n";

    voxels.ExportToVtk("poca_TA.vtk",0);
    voxelS.ExportToVtk("poca_LD.vtk");


    IBVoxFilter_Abtrim abfilt(Vector3i(5,5,5));
    IBFilterGaussShape shape(0.2);
    abfilt.SetKernelSpherical(shape);
    abfilt.SetImage(&voxels);
    abfilt.SetABTrim(0,2);
    abfilt.Run();
    abfilt.SetImage(*voxelS);
    abfilt.Run();
    voxels.ExportToVtk("poca_TA_trimm.vtk",0);
    voxelS.ExportToVtk("poca_LD_trimm.vtk",0);
    return 0;

}
