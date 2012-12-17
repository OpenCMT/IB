#include <iostream>

#include <TFile.h>
#include <TTree.h>

#include "IBPocaEvaluator.h"
#include "IBMuonError.h"
#include "IBVoxCollectionCap.h"
#include "IBVoxRaytracer.h"
#include "Detectors/MuonError.h"
#include "Detectors/MuonScatter.h"
#include "IBMuonEventTTreeReader.h"
#include "IBMuonEventTTreeR3DmcReader.h"

#include <Vtk/uLibVtkViewer.h>
#include <Vtk/vtkMuonEvent.h>
#include <Vtk/vtkStructuredGrid.h>
#include <Vtk/vtkVoxRaytracerRepresentation.h>

#include "testing-prototype.h"

using namespace uLib;

int main() {
    BEGIN_TESTING(IB Raytracer);

    // errors //
    IBMuonError sigma(11.93,2.03, 18.53,2.05);

    // reader //
    TFile* f = new TFile("/var/local/data/root/run_PDfit_201210/muSteel_PDfit_20121016_v10.root");
    TTree* t = (TTree*)f->Get("n");
    IBMuonEventTTreeReader* reader = IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::R3D_MC);
    reader->setTTree(t);
    reader->setError(sigma);
    reader->setMomentum(0.7);
    reader->selectionCode(0xff);

    // voxels //
    IBVoxel zero = {0.1E-6,0,0};
    IBVoxCollectionCap voxels(Vector3i(2,2,2));
    voxels.SetSpacing(Vector3f(400,300,150));
    voxels.SetPosition(Vector3f(-400,-300,-150));
    voxels.InitLambda(zero);


    // pocal //
    IBPocaEvaluator* processor = IBPocaEvaluator::New(IBPocaEvaluator::TiltedAxis);


    // tracer //
    IBVoxRaytracer tracer(voxels);

    // visulizer //
    uLibVtkViewer v_viewer;
    vtkStructuredGrid v_voxels(voxels);
    vtkVoxRaytracerRepresentation v_tracer(tracer);
    v_viewer.AddAbstractProp(v_voxels);

    if(reader->getNumberOfEvents() >= 100)
    for (int i=18; i<100; i++) {
        MuonScatter mu;
        if(reader->readNext(&mu)) {
            std::cout << "Event : " << reader->getCurrentPosition() << "\n";
            MuonEvent mue;
            mue.LineIn() = mu.LineIn();
            mue.LineOut() = mu.LineOut();
            vtkMuonEvent v_event(mue);

            if(processor->evaluate(mu))
                v_event.AddPocaPoint(processor->getPoca());

            v_event.PrintSelf(std::cout);
            //        sprintf(file_name,"muon_event_dump8_%d.vtp",i);
            //        v_event.SaveToXMLFile(file_name);

            v_tracer.SetMuon(v_event);
            v_viewer.AddAbstractProp(v_event);
            //v_viewer.AddAbstractProp(v_tracer);

            v_viewer.Start();
            v_viewer.RemoveAbstractProp(v_event);
        }
    }

    END_TESTING;
}
