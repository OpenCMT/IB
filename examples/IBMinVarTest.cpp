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

int main(int argc, char ** argv) {
    BEGIN_TESTING(IBMinVar);
    float vox_size = 2.5;

    // errors //
    IBMuonError sigma(6.07, 7.02); // parameters relative to scattering angles NOT measured angles!!

    // reader //
    TFile* f = new TFile (argv[1]);
    IBMuonEventTTreeReader* reader = IBMuonEventTTreeReader::New(f);
    reader->setError(sigma);
    reader->setMomentum(1);

    // voxels //
    Vector3f vox_bounding( 300, 63,  240   );
    Vector3f vox_pos(     -150, -183,-120  );


    IBVoxCollection voxels(Vector3i(vox_bounding(0)/vox_size,
                                    vox_bounding(1)/vox_size,
                                    vox_bounding(2)/vox_size));

    /// init density
    IBVoxel zero = {0,0,0};
    IBVoxel air = {0.1E-6,0,0};

    voxels.InitLambda(air);

    voxels.SetSpacing (Vector3f(vox_size,
                                vox_size,
                                vox_size));
    voxels.SetPosition(vox_pos);


    // tracer //
    IBVoxRaytracer* tracer = new IBVoxRaytracer(voxels);

    int tot  = 0;
    int tot2 = 0;

    reader->setAcquisitionTime(20);

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
