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



#include <stdio.h>

#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>

#include <Core/Vector.h>
#include <Core/Debug.h>

#include "IBPocaEvaluator.h"
#include "IBMinimizationVariablesEvaluator.h"
#include "IBVoxRaytracer.h"

#include "IBVoxCollectionCap.h"
#include "IBAnalyzerEM.h"

#include "IBAnalyzerEMAlgorithm.h"
#include "IBAnalyzerEMAlgorithmSGA.h"

class IBAnalyzerEMPimpl;

namespace {
typedef IBAnalyzerEM::Event Event;
//static DebugTTree trd(__FILE__);
} // namespace



////////////////////////////////////////////////////////////////////////////////
/////  PIMPL  //////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class IBAnalyzerEMPimpl {

    typedef IBAnalyzerEM::Event Event;

public:
    IBAnalyzerEMPimpl(IBAnalyzerEM *parent) :
        m_parent(parent),
        m_SijAlgorithm(NULL)
    {}


    void Project(Event *evc);

    void BackProject(Event *evc);

    void Evaluate(float muons_ratio);

    void SijCut(float threshold);

    void SijGuess(float threshold, float p);

    void Chi2Cut(float threshold);

    // members //
    IBAnalyzerEM                 *m_parent;
    IBAnalyzerEMAlgorithm        *m_SijAlgorithm;
    Vector<Event> m_Events;

};


void IBAnalyzerEMPimpl::Project(Event *evc)
{
    // compute sigma //
    Matrix4f Sigma = Matrix4f::Zero();
    m_SijAlgorithm->ComputeSigma(Sigma, evc);
    // compute sij //
    m_SijAlgorithm->evaluate(Sigma,evc);
}

void IBAnalyzerEMPimpl::BackProject(Event *evc)
{
    IBVoxel *vox;
    // sommatoria della formula 38 //
    for (unsigned int j = 0; j < evc->elements.size(); ++j) {
        vox = evc->elements[j].voxel;
#       pragma omp atomic
        vox->SijCap += evc->elements[j].Sij;
//        {
//            //            IBVoxel *v0 = &(*m_parent->GetVoxCollection()->Data().begin());
//            //            int id = (vox-v0)/sizeof(IBVoxel);
//            if( isnan(evc->elements[j].Sij) ) {
//                std::cout << "nan Sij in vox:" << vox << " mu:" << evc << "\n" << std::flush;
//            }
//            if( isnan(vox->SijCap) ) {
//                std::cout << "nan SijCap in vox:" << vox << " mu:" << evc << "\n" << std::flush;
//            }
//        }
#       pragma omp atomic
        vox->Count++;
    }
}

void IBAnalyzerEMPimpl::Evaluate(float muons_ratio)
{    
    unsigned int start = 0;
    unsigned int end = (unsigned int) (m_Events.size() * muons_ratio);

    if(m_SijAlgorithm) {
        // Projection
#       pragma omp parallel for
        for (unsigned int i = start; i < end; ++i)
            this->Project(&m_Events[i]);
#       pragma omp barrier

        // Backprojection
#       pragma omp parallel for
        for (unsigned int i = start; i < end; ++i)
            this->BackProject(&m_Events[i]);
#       pragma omp barrier
    }
    else {
        std::cerr << "Error: Lamda ML Algorithm not setted\n";
    }
}



////////////////////////////////////////////////////////////////////////////////
////// CUTS ////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// SijCut RECIPE1:  (true if Sij cut proposed) //
static bool em_test_SijCut(const Event &evc, float cut_level)
{
    int n_cuts = 0;
    for (unsigned int i = 0; i < evc.elements.size(); i++) {
        const Event::Element &el = evc.elements[i];
        if (fabs( (el.Sij * el.voxel->Count - el.voxel->SijCap)
                  / el.voxel->SijCap ) > cut_level) n_cuts++;
    }
    if (n_cuts > (int)(evc.elements.size()/3) ) return true;
    else return false;
}

void IBAnalyzerEMPimpl::SijCut(float threshold)
{
    Vector< Event >::iterator itr = this->m_Events.begin();
    while (itr != this->m_Events.end()) {
        if(em_test_SijCut(*itr, threshold))
        {
            this->m_Events.remove_element(*itr);
        }
        else ++itr;
    }
}

void IBAnalyzerEMPimpl::SijGuess(float threshold, float p)
{
    Vector< Event >::iterator itr = this->m_Events.begin();
    do {
        if(em_test_SijCut(*itr, threshold))
        {
            itr->header.InitialSqrP = m_parent->$$.nominal_momentum / p;
            itr->header.InitialSqrP *= itr->header.InitialSqrP;
        }
        itr++;
    } while (itr != this->m_Events.end());
}


void IBAnalyzerEMPimpl::Chi2Cut(float threshold)
{
    std::vector< Event >::iterator itr = this->m_Events.begin();
    do {
        Matrix4f Sigma = Matrix4f::Zero();
        Event &evc = *itr;
        this->m_SijAlgorithm->ComputeSigma(Sigma,&evc);
        Matrix4f iS = Sigma.inverse();
        Matrix4f Dn = iS * (evc.header.Di * evc.header.Di.transpose());
        if ( Dn.trace() > threshold )
            this->m_Events.remove_element(*itr);
        else ++itr;
    } while (itr != this->m_Events.end());
}





////////////////////////////////////////////////////////////////////////////////
////// UPDATE DENSITY ALGORITHM ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


class UpdateDensitySijCapAlgorithm :
        public IBInterface::IBVoxCollectionStaticUpdateAlgorithm
{
public:
    static void UpdateDensity(IBVoxCollection *voxels, unsigned int threshold)
    {
        for(unsigned int i=0; i< voxels->Data().size(); ++i) {
            IBVoxel& voxel = voxels->Data()[i];
            unsigned int tcount = voxel.Count;
            if ( voxel.Value > 0 && tcount > 0 && (threshold == 0 || tcount >= threshold) ) {
                voxel.Value += voxel.SijCap / static_cast<float>(tcount);
                if(unlikely(!isFinite(voxel.Value) || voxel.Value > 100.E-6)) {  // HARDCODED!!!
                    voxel.Value = 100.E-6;
                }
                //                 else if (unlikely(voxel.Value < 0.)) voxel.Value = 0.1E-6;
            }
            // else
            //             voxel.Value = 0;
            voxel.SijCap = 0;
        }
    }
};














////////////////////////////////////////////////////////////////////////////////
// IB ANALYZER EM  /////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

IBAnalyzerEM::IBAnalyzerEM(IBVoxCollection &voxels) :
    m_PocaAlgorithm(NULL),
    m_VarAlgorithm(NULL),
    m_RayAlgorithm(NULL),
    m_UpdateAlgorithm(NULL)
{
    BaseClass::SetVoxCollection(&voxels);
    init_properties(); // < DANGER !!! should be moved away !!
    d = new IBAnalyzerEMPimpl(this);
}

IBAnalyzerEM::~IBAnalyzerEM()
{
    delete d;
}

Vector<IBAnalyzerEM::Event> &IBAnalyzerEM::Events()
{
    return d->m_Events;
}


bool IBAnalyzerEM::AddMuon(const MuonScatterData &muon)
{
    if(unlikely(!m_RayAlgorithm || !m_VarAlgorithm)) return false;
    Event evc;

    evc.header.InitialSqrP = pow($$.nominal_momentum/muon.GetMomentum() ,2);
    if(isnan(evc.header.InitialSqrP)) std::cout << "sono in AddMuon: nominalp:" << $$.nominal_momentum << " muon.GetMomentum():" << muon.GetMomentum() <<"\n" << std::flush;
//    DBG(trd,evc.header.InitialSqrP,"invP2/F");
    if(likely(m_VarAlgorithm->evaluate(muon))) {
        evc.header.Di = m_VarAlgorithm->getDataVector();
        evc.header.E  = m_VarAlgorithm->getCovarianceMatrix();

        // HARDCODED ... ZERO CROSS CORRELATION BETWEEN VARS //
        //        evc.header.E(0,1) = 0.;
        //        evc.header.E(1,0) = 0.;
        //        evc.header.E(2,3) = 0.;
        //        evc.header.E(3,2) = 0.;
        // .................................................. //

        // HARDCODED ... ZERO CROSS CORRELATION BETWEEN VIEWS //
        //        evc.header.E.block<2,2>(2,0) = Matrix2f::Zero();
        //        evc.header.E.block<2,2>(0,2) = Matrix2f::Zero();
        // .................................................. //

        // HARDCODED ... LESS ERROR ! //
        //evc.header.E = Matrix4f::Zero();
        //evc.header.E /= 2;
        //        std::cout
        //                << " evc.header.Di " << evc.header.Di.transpose() << "\n"
        //                << " evc.header.E " << evc.header.E << "\n";
    }
    else return false;

    IBVoxRaytracer::RayData ray;
    { // Get RayTrace RayData //
        HPoint3f entry_pt,poca,exit_pt;
        if( !m_RayAlgorithm->GetEntryPoint(muon.LineIn(),entry_pt) ||
                !m_RayAlgorithm->GetExitPoint(muon.LineOut(),exit_pt) )
            return false;

        bool use_poca = false;
        if(m_PocaAlgorithm) { //TODO:  move this to poca algorithm
            use_poca = m_PocaAlgorithm->evaluate(muon);
            poca = m_PocaAlgorithm->getPoca();
//            DBG(trd,poca,"x/F:y/F:z/F:h/F");

            HVector3f in, out;
            in  = poca - muon.LineIn().origin;
            out = muon.LineOut().origin - poca;
            float poca_prj = in.transpose() * out;
//            DBG(trd,poca_prj);
            use_poca &= ( poca_prj > 0 );
        }
        if(use_poca && this->GetVoxCollection()->IsInsideBounds(poca)) {
            poca = m_PocaAlgorithm->getPoca();
            ray = m_RayAlgorithm->TraceBetweenPoints(entry_pt,poca);
            ray.AppendRay( m_RayAlgorithm->TraceBetweenPoints(poca,exit_pt) );

        }
        else {
            ray = m_RayAlgorithm->TraceBetweenPoints(entry_pt,exit_pt);
        }
    }

    Event::Element elc;
    Scalarf T = ray.TotalLength();
    for(int i=0; i<ray.Data().size(); ++i)
    {
        // voxel //
        const IBVoxRaytracer::RayData::Element &el = ray.Data().at(i);
        elc.voxel = &this->GetVoxCollection()->operator [](el.vox_id);
        // Wij   //
        Scalarf L = el.L;  T = fabs(T-L);
        elc.Wij << L ,          L*L/2 + L*T,
                   L*L/2 + L*T, L*L*L/3 + L*L*T + L*T*T;
        // pw    //
        elc.pw = evc.header.InitialSqrP;
        evc.elements.push_back(elc);
    }

    d->m_Events.push_back(evc);
//    trd.Fill();
    return true;
}

void IBAnalyzerEM::SetMuonCollection(IBMuonCollection *muons)
{
    uLibAssert(muons);
    d->m_Events.clear();
    for(int i=0; i<muons->size(); ++i)
    {
        this->AddMuon(muons->At(i));
    }
    BaseClass::SetMuonCollection(muons);
}

unsigned int IBAnalyzerEM::Size()
{
    return d->m_Events.size();
}

void IBAnalyzerEM::Run(unsigned int iterations, float muons_ratio)
{
    // performs iterations //
    for (unsigned int it = 0; it < iterations; it++) {
        fprintf(stderr,"\r[%d muons] EM -> performing iteration %i",
                (int) d->m_Events.size(), it);
        d->Evaluate(muons_ratio);          // run single iteration of proback //
        if(!m_UpdateAlgorithm)
            this->GetVoxCollection()->
                UpdateDensity<UpdateDensitySijCapAlgorithm>(10);                // HARDCODE THRESHOLD
        else
            this->m_UpdateAlgorithm->operator()(this->GetVoxCollection(),10);   // HARDCODE THRESHOLD
    }
    printf("\nEM -> done\n");
}

void IBAnalyzerEM::SetMLAlgorithm(IBAnalyzerEMAlgorithm *MLAlgorithm)
{
    d->m_SijAlgorithm = MLAlgorithm;
}

void IBAnalyzerEM::SijCut(float threshold) {
    d->Evaluate(1);
    d->SijCut(threshold);
    this->GetVoxCollection()->UpdateDensity<UpdateDensitySijCapAlgorithm>(0);   // HARDCODE THRESHOLD
}

void IBAnalyzerEM::SijGuess(Vector<Vector2f> tpv)
{
    d->Evaluate(1);
    // ATTENZIONE!! il vettore deve essere ordinato per threshold crescenti   //
    for (int i=0; i<tpv.size(); ++i)
        d->SijGuess( tpv[i](0), tpv[i](1) );
    this->GetVoxCollection()->UpdateDensity<UpdateDensitySijCapAlgorithm>(0);   // HARDCODE THRESHOLD
}

void IBAnalyzerEM::Chi2Cut(float threshold)
{
    d->Evaluate(1);
    this->GetVoxCollection()->UpdateDensity<UpdateDensitySijCapAlgorithm>(0);   // HARDCODE THRESHOLD
    d->Chi2Cut(threshold);
}


void IBAnalyzerEM::SetVoxCollection(IBVoxCollection *voxels)
{
    if(this->GetMuonCollection()) {
        BaseClass::SetVoxCollection(voxels);
        this->SetMuonCollection(BaseClass::GetMuonCollection());
    }
    else
        std::cerr << "*** Analyzer EM is unable to reset Voxels ***\n" <<
                     "*** without a defined muon collection ... ***\n";
}

void IBAnalyzerEM::SetVoxcollectionShift(Vector3f shift)
{
  if(this->GetMuonCollection()) {
    IBVoxCollection *voxels = this->GetVoxCollection();
    Vector3f pos = voxels->GetPosition();
    voxels->SetPosition(pos + shift);
    IBMuonCollection *muons = this->GetMuonCollection();
    for(int i=0; i<muons->size(); ++i)
        this->AddMuon(muons->At(i));
  }
}



void IBAnalyzerEM::DumpP(const char *filename, float x0, float x1)
{
    static int counter = 0;
    static TFile *file = new TFile(filename,"RECREATE");

    if(!filename) {
        file->Write();
        file->Close();
        delete file;
        return;
    }

    if(file) {
        char name[100];
        sprintf(name,"p_%i",counter++);
        gDirectory->cd(file->GetPath());
        TH1F *h = new TH1F(name,"1/p^2 distribution [1/GeV^2]",1000,x0,x1);
        float p0 = $$.nominal_momentum * $$.nominal_momentum;
        for(Id_t i=0; i<d->m_Events.size(); ++i)
            h->Fill(d->m_Events[i].header.InitialSqrP / p0 );
        h->Write();
        delete h;
    }
}






