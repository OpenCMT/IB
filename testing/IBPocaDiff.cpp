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

int main() {
    BEGIN_TESTING(IBPoCA);

    IBMuonError error(11.93,2.03,18.53,2.05);
    TFile* f = new TFile("/home/eth/musteel/data/muSteel_PDfit_2012122300_v11.root");
    TTree* t = (TTree*)f->Get("n");
    IBMuonEventTTreeReader* reader = IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::R3D_MC);
    reader->setTTree(t);
    reader->setError(error);
    reader->setMomentum(0.7);
    reader->selectionCode(IBMuonEventTTreeR3DmcReader::Whole);

    IBPocaEvaluator* pocaprocessor1 = IBPocaEvaluator::New(IBPocaEvaluator::TiltedAxis);
    IBLineDistancePocaEvaluator* pocaprocessor2 = new IBLineDistancePocaEvaluator();

    TFile o("PocaStatistics.root","RECREATE");
    TH1F dist("distance","d",1000,0,2000);
    for (int i = 0; i<2500000; ++i) {
        MuonScatter event;
        if(reader->readNext(&event)) {

            bool ok1 = pocaprocessor1->evaluate(event);
            bool ok2 = pocaprocessor2->evaluate(event);
            dist.Fill((pocaprocessor2->getPoca()-pocaprocessor1->getPoca()).norm());
        } else {
            std::cout << "Event " << reader->getCurrentPosition() << " skipped.\n";
        }
    }
    o.cd();
    dist.Write();
    o.Close();

    END_TESTING
}
