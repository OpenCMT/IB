#include <omp.h>
#include <TFile.h>
#include <TTree.h>
#include <Detectors/MuonError.h>
#include <Detectors/MuonEvent.h>
#include "IBPocaEvaluator.h"
#include "testing-prototype.h"
#include "Math/VoxImage.h"
#include "IBMuonEventTTreeReader.h"
#include "IBMuonEventTTreeR3DmcReader.h"
#include "Vtk/vtkMuonScatter.h"
#include "IBLineDistancePocaEvaluator.h"
#include "Vtk/uLibVtkViewer.h"
#include "Vtk/vtkVoxImage.h"

#include <IBVoxCollection.h>

using namespace uLib;

int main() {
    BEGIN_TESTING(IBPoCA);

    IBMuonError error(11.93,2.03,18.53,2.05);
    TFile* f = new TFile("/var/local/data/root/ROC_sets/201212/lead/muSteel_PDfit_2012122300_v11.root");
    TTree* t = (TTree*)f->Get("n");
    IBMuonEventTTreeReader* reader = IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::R3D_MC);
    reader->setTTree(t);
    reader->setError(error);
    reader->setMomentum(0.7);
    reader->selectionCode(IBMuonEventTTreeR3DmcReader::Whole);

    IBPocaEvaluator* pocaprocessor1 = IBPocaEvaluator::New(IBPocaEvaluator::TiltedAxis);
    IBLineDistancePocaEvaluator* pocaprocessor2 = new IBLineDistancePocaEvaluator();

    IBVoxCollection image(Vector3i(140,72,60));
    image.SetSpacing(Vector3f(5,5,5));
    image.SetPosition(Vector3f(-350,-180,-150));

    image.ExportToVtk("BBox.vtk",0);

    char fname1[20],fname2[20],fname3[20],fname4[20];
    for (int i = 0; i<40; ++i) {
        MuonScatter event;
        if(reader->readNext(&event)) {
            sprintf(fname1, "trackI_%i.vtk",i);
            sprintf(fname2, "trackO_%i.vtk",i);
            sprintf(fname3, "track1_%i.vtk",i);
            sprintf(fname4, "track2_%i.vtk",i);
            bool ok1 = pocaprocessor1->evaluate(event);
            bool ok2 = pocaprocessor2->evaluate(event);

            vtkMuonScatter v_ev1(event);
            vtkMuonScatter v_ev2(event);
            vtkMuonScatter v_ev3(event);
            vtkMuonScatter v_ev4(event);
            v_ev1.AddPocaPoint(pocaprocessor2->getInTrackPoca());
            v_ev2.AddPocaPoint(pocaprocessor2->getOutTrackPoca());
            v_ev3.AddPocaPoint(pocaprocessor1->getPoca());
            v_ev4.AddPocaPoint(pocaprocessor2->getPoca());

            std::cout << "ev: " << i << "  [ p1 = " <<
                         pocaprocessor1->getPoca().transpose() <<
                         " , p2 = " <<
                         pocaprocessor2->getPoca().transpose() <<
                         " ]\n";

            v_ev1.SaveToFile(fname1);
            v_ev2.SaveToFile(fname2);
            v_ev3.SaveToFile(fname3);
            v_ev4.SaveToFile(fname4);

        } else {
            std::cout << "Event " << reader->getCurrentPosition() << " skipped.\n";
        }
    }

    END_TESTING
}
