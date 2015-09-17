/*//////////////////////////////////////////////////////////////////////////////
// CMT Cosmic Muon Tomography project //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

  Copyright (c) 2014, Universita' degli Studi di Padova, INFN sez. di Padova

  Coordinators: Prof. Gianni Zumerle < gianni.zumerle@pd.infn.it >
                Paolo Checchia       < paolo.checchia@pd.infn.it >

  Authors: Andrea Rigoni Garola < andrea.rigoni@pd.infn.it >
           Matteo Furlan        < nuright@gmail.com >
           Sara Vanini          < sara.vanini@pd.infn.it >
	   Joel Klinger         < klinger@pd.infn.it >

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
        m_SijAlgorithm(NULL){;}



    void Project(Event *evc);

    void BackProject(Event *evc);

    void Evaluate(float muons_ratio);

    void filterEventsVoxelMask();

    void filterEventsLineDistance(float min, float max);

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

/// filter events after voxel mask has been applied
void IBAnalyzerEMPimpl::filterEventsVoxelMask()
{
    std::cout << "\nIBAnalyzerEM: Removing frozen voxels from " << this->m_Events.size() << " muon collection.";
    Vector< Event >::iterator itr = this->m_Events.begin();
    const Vector< Event >::iterator begin = this->m_Events.begin();

    while (itr != this->m_Events.end()) {
        Event & evc = *itr;
        Vector< Event::Element >::iterator itre = evc.elements.begin();
        // create new vector and fill it with positive voxels
        Vector< Event::Element > newvelc;
        while (itre != evc.elements.end()) {
            Event::Element & elc = *itre;
            if(elc.voxel->Value <= 0){
                // add contribution to E matrix in both views
                evc.header.E.block<2,2>(2,0) += elc.Wij * fabs(elc.voxel->Value) * evc.header.InitialSqrP;
                evc.header.E.block<2,2>(0,2) += elc.Wij * fabs(elc.voxel->Value) * evc.header.InitialSqrP;
            }
            else
                newvelc.push_back(elc);
            ++itre;
        }
//        if(newvelc.size() != evc.elements.size()){
//            std::cout << newvelc.size() << " voxels left from muon of " << evc.elements.size() << " voxels " << std::endl;
//            std::cout << this->m_parent->m_MuonCollection->Data().at(itr-begin);
//        }

        evc.elements = newvelc;

        /// erase event and muon with empty voxel collection
        if(evc.elements.empty()){
            unsigned int pos = itr - begin;
            this->m_Events.remove_element(evc);
            if(this->m_parent->m_MuonCollection)
                this->m_parent->m_MuonCollection->Data().remove_element(pos);
        }
        else
            ++itr;
    }
    std::cout << " " << this->m_Events.size() << " muons left!" << std::endl;

    return;
}

////////////////////////////////////////////////////////////////////////////////
/// filter events if in-out line distance out of range
void IBAnalyzerEMPimpl::filterEventsLineDistance(float min, float max)
{
    std::cout << "\n*** Removing events with line distance out of range from " << this->m_Events.size() << " muon collection.";

    Vector< Event >::iterator itr = this->m_Events.begin();
    const Vector< Event >::iterator begin = this->m_Events.begin();

    while (itr != this->m_Events.end()) {
        Event & evc = *itr;
        unsigned int pos = itr - begin;

        MuonScatterData muon = this->m_parent->m_MuonCollection->At(pos);
        bool use_poca = this->m_parent->m_PocaAlgorithm->evaluate(muon);
        float dist = this->m_parent->m_PocaAlgorithm->getDistance();

        /// erase event and muon with distance out of range
        if(!isFinite(dist) || dist >= max || dist < min){
            this->m_Events.remove_element(evc);
            this->m_parent->m_MuonCollection->Data().remove_element(pos);
        }
        else
	  ++itr;
    }

    std::cout << " " << this->m_Events.size() << " muons left!" << std::endl;

    return;
}

////////////////////////////////////////////////////////////////////////////////
/// SijCut RECIPE1:  (true if Sij cut proposed) //
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
    const Vector< Event >::iterator begin = this->m_Events.begin();

    int count = 0;
    while (itr != this->m_Events.end()) {
        if(em_test_SijCut(*itr, threshold))
        {
            unsigned int pos = itr - begin;
            this->m_Events.remove_element(*itr);
            if(this->m_parent->m_MuonCollection)
                this->m_parent->m_MuonCollection->Data().remove_element(pos);
            count ++;
        }
        else ++itr;
    }
    std::cout << "SijCut removed muons: " << count << "\n";
}

void IBAnalyzerEMPimpl::SijGuess(float threshold, float p)
{
    Vector< Event >::iterator itr = this->m_Events.begin();
    int count = 0;
    while (itr != this->m_Events.end()) {
        if(em_test_SijCut(*itr, threshold))
        {
            itr->header.InitialSqrP = m_parent->$$.nominal_momentum / p;
            itr->header.InitialSqrP *= itr->header.InitialSqrP;
            for (unsigned int j = 0; j < itr->elements.size(); ++j) {
                itr->elements[j].pw = itr->header.InitialSqrP;
            }
            count ++;
        }
        itr++;
    }
    std::cout << "Guess class " << threshold << ", p=" << p << "   counted muons: " << count << "\n";
}


////////////////////////////////////////////////////////////////////////////////
void IBAnalyzerEMPimpl::Chi2Cut(float threshold)
{
    std::vector< Event >::iterator itr = this->m_Events.begin();
    const Vector< Event >::iterator begin = this->m_Events.begin();

    do {
        Matrix4f Sigma = Matrix4f::Zero();
        Event &evc = *itr;
        this->m_SijAlgorithm->ComputeSigma(Sigma,&evc);
        Matrix4f iS = Sigma.inverse();
        Matrix4f Dn = iS * (evc.header.Di * evc.header.Di.transpose());
        if ( Dn.trace() > threshold ){
            unsigned int pos = itr - begin;
            this->m_Events.remove_element(*itr);
            if(this->m_parent->m_MuonCollection)
                this->m_parent->m_MuonCollection->Data().remove_element(pos);
        }
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

IBAnalyzerEM::IBAnalyzerEM(IBVoxCollection &voxels, int nPath, double alpha) :
    m_PocaAlgorithm(NULL),
    m_VarAlgorithm(NULL),
    m_RayAlgorithm(NULL),
    m_UpdateAlgorithm(NULL),
    m_nPath(nPath),
    m_alpha(alpha)
{
  std::cout << "Using alpha = " << alpha << ", #path = " << nPath << std::endl;
  BaseClass::SetVoxCollection(&voxels);
  init_properties(); // < DANGER !!! should be moved away !!
  m_d = new IBAnalyzerEMPimpl(this);
}

IBAnalyzerEM::~IBAnalyzerEM()
{
    delete m_d;
}

Vector<IBAnalyzerEM::Event> &IBAnalyzerEM::Events()
{
    return m_d->m_Events;
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

    //---- Perform raytracing 
    IBVoxRaytracer::RayData ray;
    HPoint3f entry_pt,poca,exit_pt;
    bool use_poca = false;
    { // Get RayTrace RayData //

      //---- Require entry and exit points
      if( !m_RayAlgorithm->GetEntryPoint(muon.LineIn(),entry_pt) ) return false;
      if( !m_RayAlgorithm->GetExitPoint(muon.LineOut(),exit_pt)  ) return false;

      //---- Check if this track has a valid POCA
      //---- TODO:  move this to poca algorithm
      if(m_PocaAlgorithm) { 
	use_poca = m_PocaAlgorithm->evaluate(muon);
	poca = m_PocaAlgorithm->getPoca();
	
	HVector3f in, out;
	in  = poca - muon.LineIn().origin;
	out = muon.LineOut().origin - poca;
	float poca_prj = in.transpose() * out;
	use_poca &= ( poca_prj > 0 );
      }

      //---- Check if the POCA is valid
      bool valid_poca = GetVoxCollection()->IsInsideBounds(poca);
      if(!(use_poca && valid_poca)) return false;

      //---- If doing 1-path method (Entry --> Exit)
      if(m_nPath == 1){
	ray = m_RayAlgorithm->TraceBetweenPoints(entry_pt,exit_pt);	
      }
      //---- If doing 2-path method (Entry --> POCA -> Exit)
      else if(m_nPath == 2){
	ray = m_RayAlgorithm->TraceBetweenPoints(entry_pt,poca);
	ray.AppendRay( m_RayAlgorithm->TraceBetweenPoints(poca,exit_pt) );
      }
      //---- If doing 3/4-path method
      else{

	//---- Get the poca on the entry/exit tracks
	HPoint3f entry_poca = m_PocaAlgorithm->getInTrackPoca();
	HPoint3f exit_poca  = m_PocaAlgorithm->getOutTrackPoca();

	//---- Get the distance down the tracks to the inflection points
	double entry_length = m_alpha*(entry_pt - entry_poca).norm();
	double exit_length  = m_alpha*(exit_pt  - exit_poca).norm();

	//---- Get the inflection points
	HVector3f point1 = entry_pt + entry_length*muon.LineIn().direction;
	HVector3f point2 = exit_pt  - exit_length*muon.LineOut().direction;;
	
	//---- (Entry --> Point 1)
	ray = m_RayAlgorithm->TraceBetweenPoints(entry_pt,point1);
	  
	//---- 3-Path (Point 1 --> Point 2)
	if(m_nPath == 3) ray.AppendRay( m_RayAlgorithm->TraceBetweenPoints(point1,point2) );
	//---- 4-Path (Point 1 --> POCA --> Point 2)
	else{
	  ray.AppendRay( m_RayAlgorithm->TraceBetweenPoints(point1,poca) );
	  ray.AppendRay( m_RayAlgorithm->TraceBetweenPoints(poca,point2) );
	}
	//---- (Point 2 --> Exit)
	ray.AppendRay( m_RayAlgorithm->TraceBetweenPoints(point2,exit_pt) );		    
      }
    }

    //
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


        if(elc.voxel->Value <= 0){
            // add both views
            evc.header.E.block<2,2>(2,0) += elc.Wij * fabs(elc.voxel->Value) * evc.header.InitialSqrP;
            evc.header.E.block<2,2>(0,2) += elc.Wij * fabs(elc.voxel->Value) * evc.header.InitialSqrP;
        }
        else
            evc.elements.push_back(elc);
    }

    m_d->m_Events.push_back(evc);

//    trd.Fill();
    return true;
}

///
/// \note MuonCollection is syncronized with Event vector ONLY if SetMuonCollection is called
/// if AddMuons is called by the analyzer it is not!
/// \param muons
///
void IBAnalyzerEM::SetMuonCollection(IBMuonCollection *muons)
{
  std::cout << "Setting muon collection with collection " << muons << std::endl;
    uLibAssert(muons);
    std::cout << "Clearing " << std::endl;
    m_d->m_Events.clear();
    Vector<MuonScatterData>::iterator itr = muons->Data().begin();
    std::cout << "Adding muon " << std::endl;
    int iMu=0;
    while(itr != muons->Data().end()){
      iMu++;
      std::cout << "\r" <<  iMu << std::flush;
        if(!this->AddMuon(*itr))
            muons->Data().remove_element(*itr);
        else
            itr++;
    }
    std::cout << "Done, now calling base class..." << std::endl;
    BaseClass::SetMuonCollection(muons);
}

unsigned int IBAnalyzerEM::Size()
{
    return m_d->m_Events.size();
}

void IBAnalyzerEM::Run(unsigned int iterations, float muons_ratio)
{
    // performs iterations //
    for (unsigned int it = 0; it < iterations; it++) {
        fprintf(stderr,"\r[%d muons] EM -> performing iteration %i",
                (int) m_d->m_Events.size(), it);
        m_d->Evaluate(muons_ratio);          // run single iteration of proback //
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
    m_d->m_SijAlgorithm = MLAlgorithm;
}

void IBAnalyzerEM::filterEventsVoxelMask() {
    m_d->filterEventsVoxelMask();
}

void IBAnalyzerEM::filterEventsLineDistance(float min, float max) {
    m_d->filterEventsLineDistance(min, max);
}

void IBAnalyzerEM::SijCut(float threshold) {
    m_d->Evaluate(1);
    m_d->SijCut(threshold);
    this->GetVoxCollection()->UpdateDensity<UpdateDensitySijCapAlgorithm>(0);   // HARDCODE THRESHOLD
}

void IBAnalyzerEM::SijGuess(Vector<Vector2f> tpv)
{
    m_d->Evaluate(1);
    // ATTENZIONE!! il vettore deve essere ordinato per threshold crescenti   //
    for (int i=0; i<tpv.size(); ++i)
        m_d->SijGuess( tpv[i](0), tpv[i](1) );
    this->GetVoxCollection()->UpdateDensity<UpdateDensitySijCapAlgorithm>(0);   // HARDCODE THRESHOLD
}

void IBAnalyzerEM::Chi2Cut(float threshold)
{
    m_d->Evaluate(1);
    this->GetVoxCollection()->UpdateDensity<UpdateDensitySijCapAlgorithm>(0);   // HARDCODE THRESHOLD
    m_d->Chi2Cut(threshold);
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


////////////////////////////////////////////////////////////////////////////////
/////////////////// DUMP EVENTS ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

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
        gDirectory->cd(file->GetPath());
        {
            char name[100];
            sprintf(name,"inv_p_sq_%i",counter++);
            TH1F *h = new TH1F(name,"1/p^2 distribution [1/GeV^2]",1000,x0,x1);
            float p0sq = $$.nominal_momentum * $$.nominal_momentum;
            for(Id_t i=0; i<m_d->m_Events.size(); ++i)
                h->Fill(m_d->m_Events[i].header.InitialSqrP / p0sq );
            h->Write();
            delete h;
        }

        {
            char name[100];
            sprintf(name,"p_%i",counter++);
            TH1F *h = new TH1F(name,"p distribution [GeV]",1000,x0,x1);
            float p0sq = $$.nominal_momentum * $$.nominal_momentum;
            for(Id_t i=0; i<m_d->m_Events.size(); ++i)
                h->Fill( sqrt(p0sq / m_d->m_Events[i].header.InitialSqrP) );
            h->Write();
            delete h;
        }

    }
}



////////////////////////////////////////////////////////////////////////////////
/// dump events on rootuple
////////////////////////////////////////////////////////////////////////////////

void IBAnalyzerEM::dumpEventsTTree(const char *filename)
{
    /// open file, tree
    std::cout << "\n*** Dump event collection from IBAnalyzer on file " << filename << std::endl;
    static TFile *file = new TFile(filename,"update");
    gDirectory->cd(file->GetPath());

    char name[100];
    sprintf(name,"muons");
    TTree *tree = (TTree*)file->Get("muons");
    if(!tree)
        tree = new TTree(name,name);

    int ev = 0;
    float mom, sumLij, dist, DP, DX, DT, DZ;
    TBranch *bev = tree->Branch("ev",&ev,"ev/I");
    TBranch *bp = tree->Branch("p",&mom,"p/F");
    TBranch *bsumLij = tree->Branch("sumLij",&sumLij,"sumLij/F");
    TBranch *bdist = tree->Branch("dist",&dist,"dist/F");
    TBranch *bDP = tree->Branch("DP",&DP,"DP/F");
    TBranch *bDX = tree->Branch("DX",&DX,"DX/F");
    TBranch *bDT = tree->Branch("DT",&DT,"DT/F");
    TBranch *bDZ = tree->Branch("DZ",&DZ,"DZ/F");

    IBMuonCollection *muons = this->GetMuonCollection();
    /// event loop
    int pos = 0;
    Vector< Event >::iterator itr = m_d->m_Events.begin();
    while (itr != m_d->m_Events.end()) {

        Event & evc = *itr;

        /// crossed voxel loop
        sumLij = 0;
        Vector< Event::Element >::iterator itre = evc.elements.begin();
        while (itre != evc.elements.end()) {
            Event::Element & elc = *itre;
            sumLij += elc.Wij(0,0);
            ++itre;
        }
        mom = $$.nominal_momentum/sqrt(evc.header.InitialSqrP);

        DP = evc.header.Di[0];
        DX = evc.header.Di[1];
        DT = evc.header.Di[2];
        DZ = evc.header.Di[3];

        bev->Fill();
        bp->Fill();
        bsumLij->Fill();

        bDP->Fill();
        bDX->Fill();
        bDT->Fill();
        bDZ->Fill();

        pos++;
        ev++;
        itr++;
    }

    /// muon loop to add poca information
    dist = 0;

    for(int i=0; i<muons->size(); ++i){
        MuonScatterData muon = muons->At(i);
        if(m_PocaAlgorithm){
            bool use_poca = m_PocaAlgorithm->evaluate(muon);
            dist = m_PocaAlgorithm->getDistance();
            //uncomment to exclude distance when PoCA is outside voxel bounds
            //if(use_poca)
            bdist->Fill();
        }
    }

    // testing
    int sizeev = m_d->m_Events.size();
    int sizemu = muons->size();

    tree->Write("", TObject::kOverwrite);
    delete tree;

    file->Write();
    file->Close();
    delete file;

    return;
}




