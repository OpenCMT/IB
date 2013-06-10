#include <omp.h>
#include "TFile.h"
#include "TTree.h"
#include "IBMinimizationVariablesEvaluator.h"
#include "IBMuonEventTTreeReader.h"
#include "IBMuonEventTTreeR3DmcReader.h"
#include "testing-prototype.h"
#include "IBVoxRaytracer.h"
#include "IBVoxCollectionCap.h"

using namespace uLib;

//int main() {
//    BEGIN_TESTING(IBMinVar);

//    // errors //
//    IBMuonError sigma(12.24,18.85); // parameters relative to scattering angles NOT measured angles!!

//    // reader //
//    TFile* f = new TFile ("/var/local/data/root/muSteel_PDfit_20130203_v14.root");
//    TTree* t = (TTree*)f->Get("n");
//    IBMuonEventTTreeReader* reader = IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::R3D_MC);
//    reader->setTTree(t);
//    reader->setError(sigma);
//    reader->setMomentum(0.7);
//    reader->selectionCode(IBMuonEventTTreeR3DmcReader::Top2Bottom);

//    // voxels //
//    IBVoxel zero = {0.1E-6,0,0};
//    IBVoxCollection voxels(Vector3i(140,72,60));
//    voxels.SetSpacing (Vector3f(5,5,5));
//    voxels.SetPosition(Vector3f(-350,-180,-150));
//    voxels.InitLambda(zero);

//    // tracer //
//    IBVoxRaytracer* tracer = new IBVoxRaytracer(voxels);

//    // variables //
//    IBMinimizationVariablesEvaluator* minimizator =
//            IBMinimizationVariablesEvaluator::New(IBMinimizationVariablesEvaluator::NormalPlane);
//    minimizator->setRaytracer(tracer);
//    int tot  = 0;
//    int tot2 = 0;
//    reader->setAcquisitionTime(5);

//    for (int i = 0; reader->getNumberOfEvents(); ++i) {
//        MuonScatter event;
//        if (reader->readNext(&event)) {
//            if(minimizator->evaluate(event)) tot2++;
//            tot++;
//        }
//    }
//    delete minimizator;
//    std::cout << "Reader has processed " << tot  << " events" << std::endl;
//    std::cout << "MinVar has processed " << tot2 << " events" << std::endl;

//    END_TESTING
//}













////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// main
namespace ThisTest {

struct Event {
    MuonScatter mu;
    Vector4f var;
    Vector4f NPx;
    Matrix4f NPx_err;
    Vector4f s2v;
    Matrix4f s2f_err;
};

class EventCollection : public Vector<Event> {
    typedef Vector<Event> BaseClass;
public:
    EventCollection() :
        m_NPx(IBMinimizationVariablesEvaluator::New(IBMinimizationVariablesEvaluator::NormalPlane)),
        m_s2v(IBMinimizationVariablesEvaluator::New(IBMinimizationVariablesEvaluator::SimpleTwoViews))
    {
    }

    ~EventCollection() {
    }

    void Write(const char *filename) {
        TFile * tf = new TFile(filename,"RECREATE");
        TTree * tt = new TTree("NPx_vs_s2v","Scattering angles NPx vs s2v");
        //        tf->cd();

        Event ev;
        float momentum;
        tt->Branch("mu_in_dir0", &ev.mu.LineIn().direction(0));
        tt->Branch("mu_in_dir1", &ev.mu.LineIn().direction(1));
        tt->Branch("mu_in_dir2", &ev.mu.LineIn().direction(2));
        tt->Branch("mu_out_dir0", &ev.mu.LineOut().direction(0));
        tt->Branch("mu_out_dir1", &ev.mu.LineOut().direction(1));
        tt->Branch("mu_out_dir2", &ev.mu.LineOut().direction(2));
        tt->Branch("mu_in_pos0", &ev.mu.LineIn().origin(0));
        tt->Branch("mu_in_pos1", &ev.mu.LineIn().origin(1));
        tt->Branch("mu_in_pos2", &ev.mu.LineIn().origin(2));
        tt->Branch("mu_out_pos0", &ev.mu.LineOut().origin(0));
        tt->Branch("mu_out_pos1", &ev.mu.LineOut().origin(1));
        tt->Branch("mu_out_pos2", &ev.mu.LineOut().origin(2));
        tt->Branch("mu_p", &momentum);

        tt->Branch("var_0", &ev.var(0));
        tt->Branch("var_2", &ev.var(2));

        tt->Branch("NPx_0", &ev.NPx(0));
        tt->Branch("NPx_1", &ev.NPx(1));
        tt->Branch("NPx_2", &ev.NPx(2));
        tt->Branch("NPx_3", &ev.NPx(3));

        tt->Branch("s2v_0", &ev.s2v(0));
        tt->Branch("s2v_1", &ev.s2v(1));
        tt->Branch("s2v_2", &ev.s2v(2));
        tt->Branch("s2v_3", &ev.s2v(3));

        // fill cycle //
        for ( int i =0 ; i<this->size(); ++i)
        {
            ev = this->at(i);
            momentum = this->at(i).mu.GetMomentum();
            tt->Fill();
        }

        tt->Write();
        tf->Close();
    }

    void SetTracer(IBVoxRaytracer *tracer) {
        m_NPx->setRaytracer(tracer);
        m_s2v->setRaytracer(tracer);
    }

    int operator << (MuonScatter muon) {
        Event ev;
        Vector4f var;

        ev.mu = muon;
        var(0) = atan(muon.LineOut().direction(0)) - atan(muon.LineIn().direction(0));
        var(1) = 0;
        var(2) = atan(muon.LineOut().direction(2)) - atan(muon.LineIn().direction(2));
        var(3) = 0;
        ev.var = var;

        if(m_NPx->evaluate(muon) && m_s2v->evaluate(muon))
        {
            ev.NPx = m_NPx->getDataVector();
            ev.NPx_err = m_NPx->getCovarianceMatrix();
            ev.s2v = m_s2v->getDataVector();
            ev.s2f_err = m_s2v->getCovarianceMatrix();
            this->push_back(ev);
            return true;
        }
        else return false;
    }


private:
    IBMinimizationVariablesEvaluator *m_NPx;
    IBMinimizationVariablesEvaluator *m_s2v;
};

}

int main() {
    BEGIN_TESTING(IBMinVar);
    Vector3f vox_bounding(300,161,240); // centered bounding size //
    float vox_size = 10;

    // errors //
    IBMuonError sigma(6.07, 7.02); // parameters relative to scattering angles NOT measured angles!!

    // reader //
    TFile* f = new TFile ("/var/local/data/root/run_1388/r1388_x.root");
    IBMuonEventTTreeReader* reader = IBMuonEventTTreeReader::New(f);
    reader->setError(sigma);
    reader->setMomentum(1);

    // voxels //
    IBVoxel air = {0.1E-6,0,0};
    IBVoxCollection voxels(Vector3i(vox_bounding(0)/vox_size,
                                    vox_bounding(1)/vox_size,
                                    vox_bounding(2)/vox_size));
    voxels.SetSpacing (Vector3f(vox_size,
                                vox_size,
                                vox_size));
    voxels.SetPosition(Vector3f( - 150,
                                 - 172,
                                 - 120 ));
    voxels.InitLambda(air);

    // tracer //
    IBVoxRaytracer* tracer = new IBVoxRaytracer(voxels);

    int tot  = 0;
    int tot2 = 0;

    reader->setAcquisitionTime(5);

    ThisTest::EventCollection events;
    events.SetTracer(tracer);

    for (int i = 0; i < reader->getNumberOfEvents(); ++i) {
        MuonScatter mu;
        if (reader->readNext(&mu)) {
            if(i%1000==0) std::cout << "." << std::flush;
            if(events << mu) {
                tot2++;
            }
            tot++;
        }
    }

    events.Write("NPx_vs_s2v.root");

    std::cout << "Reader has processed " << tot  << " events" << std::endl;
    std::cout << "MinVar has processed " << tot2 << " events" << std::endl;

    END_TESTING
}
