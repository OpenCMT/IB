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

int main(int argv, char** argc) {
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
    IBVoxCollectionCap voxels(Vector3i(140,72,60));
    voxels.SetSpacing(Vector3f(5,5,5));
    voxels.SetPosition(Vector3f(-350,-180,-150));
    voxels.InitLambda(zero);


    // pocal //
    IBPocaEvaluator* processor = IBPocaEvaluator::New(IBPocaEvaluator::TiltedAxis);
    IBPocaEvaluator* processoR = IBPocaEvaluator::New(IBPocaEvaluator::LineDistance);


    // tracer //
    IBVoxRaytracer tracer(voxels);
    IBVoxRaytracer traceR(voxels);
    // visulizer //
    uLibVtkViewer v_viewer;
    vtkStructuredGrid v_voxels(voxels);
    vtkVoxRaytracerRepresentation v_tracer(tracer);
    vtkVoxRaytracerRepresentation v_traceR(traceR);

    v_viewer.AddAbstractProp(v_voxels);

    if(reader->getNumberOfEvents() >= 100)
        for (int i=atoi(argc[1]); i<1000; i++) {
        MuonScatter mu;
        for (int j=0; j<=i; ++j) reader->readNext(&mu);
        if(reader->readNext(&mu)) {
            std::cout << "Event : " << reader->getCurrentPosition() << "\n";
            MuonEvent mue;
            mue.LineIn() = mu.LineIn();
            mue.LineOut() = mu.LineOut();
            vtkMuonEvent v_event(mue);
            vtkMuonEvent v_evenT(mue);
            HPoint3f poca, pocA;
            if(processor->evaluate(mu) && processoR->evaluate(mu)) {
                poca = processor->getPoca();
                pocA = processoR->getPoca();
                v_event.AddPocaPoint(poca);
                v_evenT.AddPocaPoint(pocA);
            } else continue;

            if (!voxels.IsInsideBounds(poca)&&!voxels.IsInsideBounds(pocA)) continue;

//            HPoint3f entry_pt, exit_pt;
//            tracer.GetEntryPoint(mu.LineIn(),entry_pt);
//            tracer.GetExitPoint(mu.LineOut(),exit_pt);

//            v_event.PrintSelf(std::cout);
            //        sprintf(file_name,"muon_event_dump8_%d.vtp",i);
            //        v_event.SaveToXMLFile(file_name);

            v_tracer.SetMuon(v_event);
            v_traceR.SetMuon(v_evenT);
            v_tracer.SetVoxelsColor(Vector4f(0,0.5,0,0.25));
            v_traceR.SetVoxelsColor(Vector4f(0.5,0,0,0.25));
            v_tracer.SetRayColor(Vector4f(0,1,0,1));
            v_traceR.SetRayColor(Vector4f(1,0,0,1));
            v_viewer.AddAbstractProp(v_tracer);
            v_viewer.AddAbstractProp(v_traceR);
            //v_viewer.AddAbstractProp(v_event);

            v_viewer.Start();
            //v_viewer.RemoveAbstractProp(v_event);
            v_viewer.RemoveAbstractProp(v_tracer);
            v_viewer.RemoveAbstractProp(v_traceR);

        }
    }

    END_TESTING;
}
