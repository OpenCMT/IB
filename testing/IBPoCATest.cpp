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

using namespace uLib;

int main() {
    BEGIN_TESTING(IBPoCA);

    TFile* f = new TFile("/var/local/data/root/run_PDfit_201210/muSteel_20121005_v6.root");
    TTree* t = (TTree*)f->Get("n");
    IBMuonEventTTreeReader* reader = IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::R3D_MC);
    reader->setTTree(t);
    reader->selectionCode(IBMuonEventTTreeR3DmcReader::Whole);

    IBMuonError error(Vector2f(0,0),Vector2f(0,0));
    error.Phi() = 0.0003;
    error.Theta() = 0.0007;


    IBPocaEvaluator* pocaprocessor = IBPocaEvaluator::New(IBPocaEvaluator::TiltedAxis);

    //pocaprocessor->setErrors(error);

    for (int i = 0; i<100; ++i) {
        MuonScatter event;
        if(reader->readNext(&event)) {
            bool ok = pocaprocessor->evaluate(event);
            std::cout << "PoCA is " << pocaprocessor->getPoca().transpose() <<
                         " with integrity " << ok << std::endl;
        } else {
            std::cout << "Event " << reader->getCurrentPosition() << " skipped.\n";
        }
    }

    END_TESTING
}
