#include <omp.h>
#include <TFile.h>
#include <TH1F.h>
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

int main(int argc, char* argv[]) {
    BEGIN_TESTING(IBPoCA);

    IBMuonError error(12.24,18.85);
    TFile* f = new TFile(argv[1]);
    TTree* t = (TTree*)f->Get("n");
    IBMuonEventTTreeReader* reader = IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::R3D_MC);
    reader->setTTree(t);
    reader->setError(error);
    reader->setMomentum(0.7);
    reader->selectionCode(IBMuonEventTTreeR3DmcReader::Whole);

    IBPocaEvaluator* pocaprocessor1 = IBPocaEvaluator::New(IBPocaEvaluator::TiltedAxis);
    IBLineDistancePocaEvaluator* pocaprocessor2 = new IBLineDistancePocaEvaluator();

    TFile o("PocaStatistics.root","RECREATE");
    TTree * ot = new TTree("pocastat","ps");

    HPoint3f  pi, po, p1, p2, p2i, p2o;
    HVector3f si, so;
    float p1p2, p2ip2o, rdiff;
    ot->Branch("x_in",     &pi(0),  "pi0/F");
    ot->Branch("y_in",     &pi(1),  "pi1/F");
    ot->Branch("z_in",     &pi(2),  "pi2/F");
    ot->Branch("x_out",    &po(0),  "po0/F");
    ot->Branch("y_out",    &po(1),  "po1/F");
    ot->Branch("z_out",    &po(2),  "po2/F");
    ot->Branch("sx_in",    &si(0),  "si0/F");
    ot->Branch("sy_in",    &si(1),  "si1/F");
    ot->Branch("sz_in",    &si(2),  "si2/F");
    ot->Branch("sx_out",   &so(0),  "so0/F");
    ot->Branch("sy_out",   &so(1),  "so1/F");
    ot->Branch("sz_out",   &so(2),  "so2/F");
    ot->Branch("TAPx",     &p1(0),  "p10/F");
    ot->Branch("TAPy",     &p1(1),  "p11/F");
    ot->Branch("TAPz",     &p1(2),  "p12/F");
    ot->Branch("LDPx",     &p2(0),  "p20/F");
    ot->Branch("LDPy",     &p2(1),  "p21/F");
    ot->Branch("LDPz",     &p2(2),  "p22/F");
    ot->Branch("LDPx_in",  &p2i(0), "p2i0/F");
    ot->Branch("LDPy_in",  &p2i(1), "p2i1/F");
    ot->Branch("LDPz_in",  &p2i(2), "p2i2/F");
    ot->Branch("LDPx_out", &p2o(0), "p2o0/F");
    ot->Branch("LDPy_out", &p2o(1), "p2o1/F");
    ot->Branch("LDPz_out", &p2o(2), "p2o2/F");
    ot->Branch("p1_p2",    &p1p2,   "p1p2/F");
    ot->Branch("p2i_p2o",  &p2ip2o, "p2ip2o");
    ot->Branch("ldiff",    &rdiff,  "rdiff/F");

    for (int i = 0; i<200000; ++i) {
        MuonScatter event;
        if(reader->readNext(&event)) {
            pi = event.LineIn().origin;
            si = event.LineIn().direction;
            po = event.LineOut().origin;
            so = event.LineOut().direction;

            bool ok1 = pocaprocessor1->evaluate(event);
            bool ok2 = pocaprocessor2->evaluate(event);

            p1 = pocaprocessor1->getPoca();
            p2 = pocaprocessor2->getPoca();
            p2i = pocaprocessor2->getInTrackPoca();
            p2o = pocaprocessor2->getOutTrackPoca();
            p1p2 = (p2-p1).norm();
            p2ip2o = (p2o-p2i).norm();

            float np  = (pi-po).norm();
            float num = (pi-po).transpose() * (p1-p2);

            rdiff = sqrt((p1p2)*(p1p2)-pow(num/np,2));


            ot->Fill();
        } else {
            std::cout << "Event " << reader->getCurrentPosition() << " skipped.\n";
        }
    }
    o.cd();
    ot->Write();
    o.Close();

    END_TESTING
}
