#include <stdio.h>

#include <TTree.h>
#include <TFile.h>

#include <Core/Vector.h>

#include "IBPocaEvaluator.h"
#include "IBMinimizationVariablesEvaluator.h"
#include "IBVoxRaytracer.h"

#include "IBVoxCollectionCap.h"
#include "IBAnalyzerEM.h"


struct Event {
    struct Element {
        Matrix4f Wij;
        union {
            Scalarf lambda;
            Scalarf Sij;
        };
        IBVoxel *voxel;
        Scalarf pw;
    };

    struct {
        Vector4f Di;
        Matrix4f E;
        Scalarf  InitialSqrP;
    } header;
    Vector<Element> elements;
};

static int _inv_bugs = 0;


// EM ALGORITHMS DECLARATIONS //
static void four_hidden_pxtz(Matrix4f &Sigma, Event *evc);
static void two_hidden_px   (Matrix4f &Sigma, Event *evc);
static void two_hidden_tz   (Matrix4f &Sigma, Event *evc);
static void two_hidden_pt   (Matrix4f &Sigma, Event *evc);
static void two_hidden_xz   (Matrix4f &Sigma, Event *evc);
static void one_hidden_p    (Matrix4f &Sigma, Event *evc);
static void one_hidden_t    (Matrix4f &Sigma, Event *evc);
static void one_hidden_x    (Matrix4f &Sigma, Event *evc);
static void one_hidden_z    (Matrix4f &Sigma, Event *evc);

// MAP PRIOR DECLARATION //
static void gaussian_lambda_prior(float beta, Event *evc);

// PW ALGORITHM DECLARATIONS //
static void pweigth_pw(Event *evc, Scalarf nominalp);
static void pweigth_sw(Event *evc, Scalarf nominalp);
static void pweigth_cw(Event *evc, Scalarf nominalp);


////////////////////////////////////////////////////////////////////////////////
/////  PIMPL  //////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class IBAnalyzerEMPimpl {
public:
    IBAnalyzerEMPimpl(IBAnalyzerEM::Parameters *parameters) :
        m_parameters(parameters),
        m_PocaAlgorithm(NULL),
        m_VarAlgorithm(NULL),
        m_RayAlgorithm(NULL),
        m_SijSelectedFunction(four_hidden_pxtz),
        m_MAPLambdaFunction(gaussian_lambda_prior),
        m_PWeigthFunction(pweigth_pw)
    {
#ifndef NDEBUG
        ok = true;
        file = new TFile("20121030_v471_aem.root", "RECREATE");
        tree = new TTree("aem", "AnalyzerEM_Dump");
        tree->Branch("sx",        &t_dsx(0),     "sx/F");
        tree->Branch("dx",        &t_dsx(1),     "dx/F");
        tree->Branch("sz",        &t_dsx(2),     "sz/F");
        tree->Branch("dz",        &t_dsx(3),     "dz/F");
        tree->Branch("s2sx",      &t_sigma(0,0), "ssx/F");
        tree->Branch("s2dx",      &t_sigma(1,1), "sdx/F");
        tree->Branch("s2sz",      &t_sigma(2,2), "ssz/F");
        tree->Branch("s2dz",      &t_sigma(3,3), "sdz/F");
        tree->Branch("csxdx",     &t_sigma(0,1), "csxdx/F");
        tree->Branch("csxsz",     &t_sigma(0,2), "csxsz/F");
        tree->Branch("csxdz",     &t_sigma(3,3), "csxdz/F");
        tree->Branch("cdxsz",     &t_sigma(1,2), "cdxsz/F");
        tree->Branch("cdxdz",     &t_sigma(1,3), "cdxdz/F");
        tree->Branch("cszdz",     &t_sigma(2,3), "cszdz/F");
        tree->Branch("chi2_p",    &chi2.p,       "chi2p/F");
        tree->Branch("chi2_t",    &chi2.t,       "chi2t/F");
        tree->Branch("chi2_x",    &chi2.x,       "chi2x/F");
        tree->Branch("chi2_z",    &chi2.z,       "chi2z/F");
        tree->Branch("chi2_px",   &chi2.px,      "chi2px/F");
        tree->Branch("chi2_tz",   &chi2.tz,      "chi2tz/F");
        tree->Branch("chi2_pt",   &chi2.pt,      "chi2pt/F");
        tree->Branch("chi2_xz",   &chi2.xz,      "chi2xz/F");
        tree->Branch("chi2_pxtz", &chi2.pxtz,    "chi2pxtz/F");
        tree->Branch("ISx_00",    &chi2.isx_00,  "isx00/F");
        tree->Branch("ISx_01",    &chi2.isx_01,  "isx01/F");
        tree->Branch("ISx_11",    &chi2.isx_11,  "isx11/F");
        tree->Branch("Rho_X",     &chi2.rho_x,   "rhox/F");
        tree->Branch("ISz_00",    &chi2.isz_00,  "isz00/F");
        tree->Branch("ISz_01",    &chi2.isz_01,  "isz01/F");
        tree->Branch("ISz_11",    &chi2.isz_11,  "isz11/F");
        tree->Branch("Rho_Z",     &chi2.rho_z,   "rhoz/F");
#endif
    }

    void SetAlgorithm(IBAnalyzerEM::Algorithm alg);

    bool ComputeSigma(Matrix4f &Sigma, Event *evc);

    void Project(Event *evc);

    void BackProject(Event *evc);

    void Evaluate(float muons_ratio);

    void SijCut(float threshold);

    void UpdatePW(IBAnalyzerEM::PWeigthAlgorithm algorithm);

    // members //
    IBAnalyzerEM::Parameters         *m_parameters;
    IBPocaEvaluator                  *m_PocaAlgorithm;
    IBMinimizationVariablesEvaluator *m_VarAlgorithm;
    IBVoxRaytracer                   *m_RayAlgorithm;
    void (* m_SijSelectedFunction)(Matrix4f &ISigma, Event *evc);
    void (* m_MAPLambdaFunction)(float beta, Event *evc);
    void (* m_PWeigthFunction)(Event *evc, Scalarf nominalp);
    Vector<Event> m_Events;

#ifndef NDEBUG
    TFile* file;
    TTree* tree;
    Matrix4f t_sigma;
    Vector4f t_dsx;
    bool ok;
    struct CHI {
        float p, t, x, z, px, tz, pxtz, pt, xz;
        float rho_x, rho_z, isx_00, isx_11, isx_01,
                            isz_00, isz_11, isz_01;
    } chi2;
#endif

};

void IBAnalyzerEMPimpl::SetAlgorithm(IBAnalyzerEM::Algorithm alg)
{
    // when pimpl is initialized default is PXTZ .. see ctr of pimpl //
    static IBAnalyzerEM::Algorithm cache = IBAnalyzerEM::PXTZ;
    if(!(cache == alg)) {
        cache = alg;
        switch(alg){
        case IBAnalyzerEM::PXTZ:
            m_SijSelectedFunction = four_hidden_pxtz;
            break;
        case IBAnalyzerEM::PX:
            m_SijSelectedFunction = two_hidden_px;
            break;
        case IBAnalyzerEM::TZ:
            m_SijSelectedFunction = two_hidden_tz;
            break;
        case IBAnalyzerEM::PT:
            m_SijSelectedFunction = two_hidden_pt;
            break;
        case IBAnalyzerEM::XZ:
            m_SijSelectedFunction = two_hidden_xz;
            break;
        case IBAnalyzerEM::P:
            m_SijSelectedFunction = one_hidden_p;
            break;
        case IBAnalyzerEM::T:
            m_SijSelectedFunction = one_hidden_t;
            break;
        case IBAnalyzerEM::X:
            m_SijSelectedFunction = one_hidden_x;
            break;
        case IBAnalyzerEM::Z:
            m_SijSelectedFunction = one_hidden_z;
            break;
        }
    }
}

bool IBAnalyzerEMPimpl::ComputeSigma(Matrix4f &Sigma, Event *evc)
{
    for (unsigned int j = 0; j < evc->elements.size(); ++j) {
        evc->elements[j].lambda = evc->elements[j].voxel->Value;
        Sigma += evc->elements[j].Wij * evc->elements[j].lambda;
    }
    Sigma *= evc->header.InitialSqrP;
    Sigma += evc->header.E;


#ifndef NDEBUG
    t_sigma = Sigma;
    t_dsx   = evc->header.Di;
    //ok      = inv;

    Matrix4f e4_inv = Matrix4f::Zero();
    int c = 0;
    bool check = false;
    t_sigma.computeInverseWithCheck(e4_inv, check);

    c += check;

    chi2.pxtz = t_dsx.transpose() * e4_inv * t_dsx;

    // p
    chi2.p = t_dsx(0)*t_dsx(0)/t_sigma(0,0);

    // t
    chi2.t = t_dsx(2)*t_dsx(2)/t_sigma(2,2);

    // x
    chi2.x = t_dsx(1)*t_dsx(1)/t_sigma(1,1);

    // z
    chi2.z = t_dsx(3)*t_dsx(3)/t_sigma(3,3);

    // px
    Matrix2f e2_inv = Matrix2f::Zero();
    Matrix2f e2     = Matrix2f::Zero();
    Vector2f data   = Vector2f::Zero();
    e2   << t_sigma(0,0), t_sigma(0,1),
            t_sigma(1,0), t_sigma(1,1);
    data << t_dsx(0), t_dsx(1);
    e2.computeInverseWithCheck(e2_inv, check);

    c += check;
    chi2.px = data.transpose() * e2_inv * data;

    //tz
    e2_inv = Matrix2f::Zero();
    e2     = Matrix2f::Zero();
    data   = Vector2f::Zero();
    e2   << t_sigma(2,2), t_sigma(2,3),
            t_sigma(3,2), t_sigma(3,3);
    data << t_dsx(2), t_dsx(3);
    e2.computeInverseWithCheck(e2_inv, check);
    c += check;
    chi2.tz = data.transpose() * e2_inv * data;

    //pt
    e2_inv = Matrix2f::Zero();
    e2     = Matrix2f::Zero();
    data   = Vector2f::Zero();
    e2    << t_sigma(0,0), t_sigma(0,2),
             t_sigma(2,0), t_sigma(2,2);
    data <<  t_dsx(0), t_dsx(2);
    e2.computeInverseWithCheck(e2_inv, check);
    c += check;
    chi2.pt = data.transpose() * e2_inv * data;

    //xz
    e2_inv = Matrix2f::Zero();
    e2     = Matrix2f::Zero();
    data   = Vector2f::Zero();
    e2    << t_sigma(1,1), t_sigma(1,3),
             t_sigma(3,1), t_sigma(3,3);
    data <<  t_dsx(1), t_dsx(3);
    e2.computeInverseWithCheck(e2_inv, check);
    chi2.xz = data.transpose() * e2_inv * data;
    c += check;

    float den_x = (t_sigma(0,0)*t_sigma(1,1) - t_sigma(0,1)*t_sigma(0,1));
    float den_z = (t_sigma(2,2)*t_sigma(3,3) - t_sigma(2,3)*t_sigma(2,3));
    chi2.isx_00 = t_sigma(1,1) / den_x;
    chi2.isx_11 = t_sigma(0,0) / den_x;
    chi2.isx_01 = t_sigma(0,1) / den_x;
    chi2.rho_x  = t_sigma(0,1) / sqrt(t_sigma(0,0)*t_sigma(1,1));
    chi2.isz_00 = t_sigma(3,3) / den_z;
    chi2.isz_11 = t_sigma(2,2) / den_z;
    chi2.isz_01 = t_sigma(2,3) / den_z;
    chi2.rho_z  = t_sigma(2,3) / sqrt(t_sigma(2,2)*t_sigma(3,3));

    if (c) tree->Fill();
#endif
    return true;
}




void IBAnalyzerEMPimpl::Project(Event *evc)
{
    // compute sigma //
    Matrix4f Sigma = Matrix4f::Zero();
    ComputeSigma(Sigma, evc);
    // compute sij //    
    m_SijSelectedFunction(Sigma,evc);
}

void IBAnalyzerEMPimpl::BackProject(Event *evc)
{
    IBVoxel *vox;
    // sommatoria della formula 38 //
    for (unsigned int j = 0; j < evc->elements.size(); ++j) {
        vox = evc->elements[j].voxel;
        vox->SijCap += evc->elements[j].Sij;
        vox->Count++;
    }
}

void IBAnalyzerEMPimpl::Evaluate(float muons_ratio)
{
    unsigned int start = 0;
    unsigned int end = (int) (m_Events.size() * muons_ratio);

    this->SetAlgorithm(m_parameters->algorithm);
    _inv_bugs = 0;

    // Projection
    #pragma omp parallel for
    for (unsigned int i = start; i < end; ++i)
        this->Project(&m_Events[i]);
    #pragma omp barrier
    std::cout << "inv impossible at " << _inv_bugs << " tracks \n" << std::flush;    

    // Backprojection
    #pragma omp parallel for
    for (unsigned int i = start; i < end; ++i)
        this->BackProject(&m_Events[i]);
    #pragma omp barrier
}

// PWEIGTH //
void IBAnalyzerEMPimpl::UpdatePW(enum IBAnalyzerEM::PWeigthAlgorithm algorithm)
{
    static enum IBAnalyzerEM::PWeigthAlgorithm cache = IBAnalyzerEM::PWeigth_pw;
    if(cache != algorithm) {
        cache = algorithm;
        switch (algorithm) {
        case IBAnalyzerEM::PWeigth_disabled:
            this->m_PWeigthFunction = NULL;
            break;
        case IBAnalyzerEM::PWeigth_pw:
            this->m_PWeigthFunction = pweigth_pw;
            break;
        case IBAnalyzerEM::PWeigth_sw:
            this->m_PWeigthFunction = pweigth_sw;
            break;
        case IBAnalyzerEM::PWeigth_cw:
            this->m_PWeigthFunction = pweigth_cw;
            break;
        }
    }
    if (m_PWeigthFunction) {
        #pragma omp parallel for
        for (uint i = 0; i < this->m_Events.size(); ++i) {
            m_PWeigthFunction(&this->m_Events[i], m_parameters->nominal_momentum);
        }
        #pragma omp barrier
    }
}




////// CUTS //////

// true if cut proposed //
static int em_test_SijCut(Event *evc, float cut_level)
{
    int n_cuts = 0;
    for (unsigned int i = 0; i < evc->elements.size(); i++)
        if (fabs( (evc->elements[i].Sij * evc->elements[i].voxel->Count -
                   evc->elements[i].voxel->SijCap) /
                  evc->elements[i].voxel->SijCap ) > cut_level) n_cuts++;
    if (n_cuts > (int)(evc->elements.size()/3) ) return true;
    else return false;
}

void IBAnalyzerEMPimpl::SijCut(float threshold)
{
    Vector< Event >::iterator itr = this->m_Events.begin();
    do {
        if(em_test_SijCut(itr.base(), threshold))
            this->m_Events.remove_element(*itr);
        else ++itr;
    } while (itr != this->m_Events.end());
}











////////////////////////////////////////////////////////////////////////////////
// IB ANALYZER EM  /////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


IBAnalyzerEM::IBAnalyzerEM() :
    d(new IBAnalyzerEMPimpl(&parameters()))
{
    init_parameters();
}

IBAnalyzerEM::~IBAnalyzerEM()
{
#ifndef NDEBUG
    d->file->cd();
    d->tree->Write();
    d->file->Close();
#endif
    delete d;
}

void IBAnalyzerEM::AddMuon(MuonScatterData &muon)
{
    if(unlikely(!d->m_PocaAlgorithm || !d->m_RayAlgorithm || !d->m_VarAlgorithm)) return;
    Event evc;

    evc.header.InitialSqrP = parameters().nominal_momentum/muon.GetMomentum();
    evc.header.InitialSqrP *= evc.header.InitialSqrP;

    if(likely(d->m_VarAlgorithm->evaluate(muon))) {
        evc.header.Di = d->m_VarAlgorithm->getDataVector();
        evc.header.E  = d->m_VarAlgorithm->getCovarianceMatrix();
    }
    else return;

    IBVoxRaytracer::RayData ray;
    { // Get RayTrace RayData //
        HPoint3f entry_pt,poca,exit_pt;
        if( !d->m_RayAlgorithm->GetEntryPoint(muon.LineIn(),entry_pt) ||
                !d->m_RayAlgorithm->GetExitPoint(muon.LineOut(),exit_pt) )
            return;
        bool test = d->m_PocaAlgorithm->evaluate(muon);
        poca = d->m_PocaAlgorithm->getPoca();
        if(test && this->GetVoxCollection()->IsInsideBounds(poca)) {
            poca = d->m_PocaAlgorithm->getPoca();
            ray = d->m_RayAlgorithm->TraceBetweenPoints(entry_pt,poca);
            ray.AppendRay( d->m_RayAlgorithm->TraceBetweenPoints(poca,exit_pt) );
        }
        else {
            ray = d->m_RayAlgorithm->TraceBetweenPoints(entry_pt,exit_pt);
        }
    }

    Event::Element elc;
    Scalarf T = ray.TotalLength();
    for(int i=0; i<ray.Data().size(); ++i)
    {
        const IBVoxRaytracer::RayData::Element *el = &ray.Data().at(i);
        elc.voxel = &this->GetVoxCollection()->operator [](el->vox_id);
        Scalarf L = el->L;
        T -= L;

        Matrix2f wij_block;
        wij_block << L ,          L*L/2 + L*T,
                     L*L/2 + L*T, L*L*L/3 + L*L*T + L*T*T;
        elc.Wij = Matrix4f::Zero();
        elc.Wij.block<2,2>(0,0) = wij_block;
        elc.Wij.block<2,2>(2,2) = wij_block;

        evc.elements.push_back(elc);
    }
    d->m_Events.push_back(evc);

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
        this->GetVoxCollection()->UpdateDensity(10);  // update lambda         // //HARDCODE
    }
    printf("\nEM -> done\n");
}

void IBAnalyzerEM::SetPocaAlgorithm(IBPocaEvaluator *evaluator)
{
    d->m_PocaAlgorithm = evaluator;
}

void IBAnalyzerEM::SetVariablesAlgorithm(IBMinimizationVariablesEvaluator *evaluator)
{
    d->m_VarAlgorithm = evaluator;
}

void IBAnalyzerEM::SetRaytracer(IBVoxRaytracer *raytracer)
{
    d->m_RayAlgorithm = raytracer;
}

void IBAnalyzerEM::SijCut(float threshold) {
    d->Evaluate(1);
    d->SijCut(threshold);
    this->GetVoxCollection()->UpdateDensity(0); // HARDCODE THRESHOLD
}

void IBAnalyzerEM::UpdatePW(IBAnalyzerEM::PWeigthAlgorithm algorithm)
{
    d->UpdatePW(algorithm);
}



////////////////////////////////////////////////////////////////////////////////
//// Sij Evaluators ////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

// DEBUG //
static void print_matrix(Matrix4f &mat, std::ostream &o)
{
    for(int i=0; i<16; ++i) {
        o << mat.data()[i];
        if((i+1)%4==0) o << ";";
        else o << ",";
    }
}



/************************************************
 * Sij 4 HIDDEN DATA IMPLEMENTATION USING EIGEN *
 ************************************************/
static void four_hidden_pxtz(Matrix4f &Sigma, Event *evc)
 {
    bool inv = false;
    Matrix4f iS;

    Sigma.computeInverseWithCheck(iS,inv);
    iS = Sigma.inverse();


    // DUMP //
//    std::cout << "Matrice Sigma: \n";
//    print_matrix(Sigma, std::cout);
//    std::cout << "\n";
//    std::cout << "Matrice ISigma: \n";
//    print_matrix(iS, std::cout);
//    std::cout << "\n";
//    std::cout << "parametro di invertibilita: " << (iS*Sigma).sum() - 4 << " Eigen inv opinion " << inv << "\n\n";

    // if Sigma has no inverse, ISigma is computed as zero and SijCap goes 0 //
    // so we predecrement voxel Count to remove this event from mean         //
    if(unlikely(!inv)) {
       _inv_bugs++;
    }
        //        for (unsigned int j = 0; j < evc->elements.size(); ++j) {
//            evc->elements[j].voxel->Count--;
//            evc->elements[j].Sij = 0;
//        }
//        _inv_bugs++;
//        return;
//    }
    for (unsigned int j = 0; j < evc->elements.size(); ++j) {
        Matrix4f iSWij       = iS * evc->elements[j].Wij;
        float DISWISD        = evc->header.Di.transpose() * iSWij * iS * evc->header.Di;
        evc->elements[j].Sij = (DISWISD - iSWij.trace()) * evc->elements[j].lambda *
                                evc->elements[j].lambda * evc->header.InitialSqrP / 2;
    }
}


/************************************************
 * Sij 2 HIDDEN DATA (S,X) EIGEN IMPLEMENTATION *
 ************************************************/
static void two_hidden_px(Matrix4f &Sigma, Event *evc)
{
    Matrix2f iS;
    {
        Matrix2f S;
        S << Sigma(0,0), Sigma(0,1), Sigma(1,0), Sigma(1,1);
        iS = S.inverse();
    }
    Vector2f Di(evc->header.Di(0),evc->header.Di(1));

    for (unsigned int j = 0; j < evc->elements.size(); ++j) {
        Matrix2f Wij;
        Wij << evc->elements[j].Wij(0,0),evc->elements[j].Wij(0,1),
               evc->elements[j].Wij(1,0),evc->elements[j].Wij(1,1);
        Matrix2f iSWij = iS * Wij;
        float DISWISD  = Di.transpose() * iSWij * iS * Di;
        evc->elements[j].Sij = (DISWISD - iSWij.trace()) * evc->header.InitialSqrP *
                               evc->elements[j].lambda * evc->elements[j].lambda / 2.;
    }
}

static void two_hidden_tz(Matrix4f &Sigma, Event *evc)
{
    Matrix2f iS;
    {
        Matrix2f S;
        S << Sigma(2,2), Sigma(2,3), Sigma(3,2), Sigma(3,3);
        iS = S.inverse();
    }
    Vector2f Di(evc->header.Di(2),evc->header.Di(3));

    for (unsigned int j = 0; j < evc->elements.size(); ++j) {
        Matrix2f Wij;
        Wij << evc->elements[j].Wij(2,2),evc->elements[j].Wij(2,3),
               evc->elements[j].Wij(3,2),evc->elements[j].Wij(3,3);
        Matrix2f iSWij = iS * Wij;
        float DISWISD  = Di.transpose() * iSWij * iS * Di;
        evc->elements[j].Sij = (DISWISD - iSWij.trace()) * evc->header.InitialSqrP *
                               evc->elements[j].lambda * evc->elements[j].lambda / 2.;
    }
}


/******************************************
 * Sij 2 HIDDEN DATA (S,S) IMPLEMENTATION *
 ******************************************/
static void two_hidden_pt(Matrix4f &Sigma, Event *evc)
{
    Matrix2f iS;
    {
        Matrix2f S;
        S << Sigma(0,0), Sigma(0,2), Sigma(2,0), Sigma(2,2);
        iS = S.inverse();
    }
    Vector2f Di(evc->header.Di(0),evc->header.Di(2));

    for (unsigned int j = 0; j < evc->elements.size(); ++j) {
        Matrix2f Wij;
        Wij << evc->elements[j].Wij(0,0),evc->elements[j].Wij(0,2),
               evc->elements[j].Wij(2,0),evc->elements[j].Wij(2,2);
        Matrix2f iSWij = iS * Wij;
        float DISWISD  = Di.transpose() * iSWij * iS * Di;
        evc->elements[j].Sij = (DISWISD - iSWij.trace()) * evc->header.InitialSqrP *
                               evc->elements[j].lambda * evc->elements[j].lambda / 2.;
    }
}


/******************************************
* Sij 2 HIDDEN DATA (D,D) IMPLEMENTATION *
******************************************/
static void two_hidden_xz(Matrix4f &Sigma, Event *evc)
{
    Matrix2f iS;
    {
        Matrix2f S;
        S << Sigma(1,1), Sigma(1,3), Sigma(3,1), Sigma(3,3);
        iS = S.inverse();
    }
    Vector2f Di(evc->header.Di(1),evc->header.Di(3));

    for (unsigned int j = 0; j < evc->elements.size(); ++j) {
        Matrix2f Wij;
        Wij << evc->elements[j].Wij(1,1),evc->elements[j].Wij(1,3),
               evc->elements[j].Wij(3,1),evc->elements[j].Wij(3,3);
        Matrix2f iSWij = iS * Wij;
        float DISWISD  = Di.transpose() * iSWij * iS * Di;
        evc->elements[j].Sij = (DISWISD - iSWij.trace()) * evc->header.InitialSqrP *
                               evc->elements[j].lambda * evc->elements[j].lambda / 2.;
    }
}


/****************************************
 * Sij 1 HIDDEN DATA (S) IMPLEMENTATION *
 ****************************************/
static void one_hidden_p(Matrix4f &Sigma, Event *evc)
{
    float Di = evc->header.Di(0);
    float iS = 1/Sigma(0,0);
    for (unsigned int j = 0; j < evc->elements.size(); ++j) {
        float Wij = evc->elements[j].Wij(0,0);
        float DISWISD = Di * iS * Wij * iS * Di;
        evc->elements[j].Sij = (DISWISD - iS * Wij) * evc->header.InitialSqrP *
                                evc->elements[j].lambda * evc->elements[j].lambda / 2.;
    }
}
static void one_hidden_t(Matrix4f &Sigma, Event *evc)
{
    float Di = evc->header.Di(2);
    float iS = 1/Sigma(2,2);
    for (unsigned int j = 0; j < evc->elements.size(); ++j) {
        float Wij = evc->elements[j].Wij(2,2);
        float DISWISD = Di * iS * Wij * iS * Di;
        evc->elements[j].Sij = (DISWISD - iS * Wij) * evc->header.InitialSqrP *
                                evc->elements[j].lambda * evc->elements[j].lambda / 2.;
    }
}

/************************************************
* Sij 1 HIDDEN DATA (D) CLASSIC IMPLEMENTATION *
************************************************/
static void one_hidden_x(Matrix4f &Sigma, Event *evc)
{
    float Di = evc->header.Di(1);
    float iS = 1/Sigma(1,1);
    for (unsigned int j = 0; j < evc->elements.size(); ++j) {
        float Wij = evc->elements[j].Wij(1,1);
        float DISWISD = Di * iS * Wij * iS * Di;
        evc->elements[j].Sij = (DISWISD - iS * Wij) * evc->header.InitialSqrP *
                                evc->elements[j].lambda * evc->elements[j].lambda / 2.;
    }
}
static void one_hidden_z(Matrix4f &Sigma, Event *evc)
{
    float Di = evc->header.Di(3);
    float iS = 1/Sigma(3,3);
    for (unsigned int j = 0; j < evc->elements.size(); ++j) {
        float Wij = evc->elements[j].Wij(3,3);
        float DISWISD = Di * iS * Wij * iS * Di;
        evc->elements[j].Sij = (DISWISD - iS * Wij) * evc->header.InitialSqrP *
                                evc->elements[j].lambda * evc->elements[j].lambda / 2.;
    }
}


////////////////////////////////////////////////////////////////////////////////
///// LAMBDA PRIOR /////////////////////////////////////////////////////////////


static void gaussian_lambda_prior(float beta, Event *evc) {
    // foreach Element in Event
    // 1) compute beta from voxel track count
    // 2) compute new Sij_tmp from lambda and Sij_cap using Gamma function

    for(int j=0; j<evc->elements.size(); ++j)
    {
        IBVoxel *voxel = evc->elements.at(j).voxel;
        float tau = static_cast<float>(voxel->Count) / beta; //  1/bj
        float lambdaj = voxel->Value + evc->elements.at(j).Sij;
        float a = tau * lambdaj / 2;
        float b = tau * sqrt( (lambdaj * lambdaj / 4) + ( tau / 27 ) );

        // essendo lambdaj definito positivo ( vero? ) uso solo il semiasse > 0
        evc->elements[j].Sij = pow(a+b, 1/3) + pow(a-b, 1/3) - voxel->Value;    // FIXX
    }
}





////////////////////////////////////////////////////////////////////////////////
////// PWEIGTH ALGORITHMS  /////////////////////////////////////////////////////


static void pweigth_pw(Event *evc, Scalarf nominalp) {
    float pw_A = 1.42857;
    float pw_epsilon = 4.0;
    float Xres = 0;
    Matrix4f Wij;
    for (int j = evc->elements.size(); j --> 0;) //BACKWARD
    {
        Wij = evc->elements[j].Wij;

        evc->elements[j].pw = pw_A * (nominalp) *
                sqrt(pw_epsilon/(Xres + pw_epsilon));

        Xres += (evc->elements[j].voxel->Value * 1.E6 < 2.5) ?
                    Wij(0,0) * evc->elements[j].voxel->Value * 40000 :
                    Wij(0,0) * 2.5 * 0.04;

        //Xres += Wij[0] * evc->elements[j].voxel->density * 40000; //<- previous version of
                                      // p_weight: bugged but enhancing!
    }
}

static void pweigth_sw(Event *evc, Scalarf nominalp) {
    float scale = 0.8;   //this is experimental and empirical!
    float alpha = 0.436; // this is the distribution peak angle, 25 degrees
    // note: albeit NOT efficient, this procedure is fundamental in testing the algorithm capability
    TODO("Optimize PW functions");
    float b1 = 13.52; // Iron Values
    float b2 = 319.9;
    float c1 = 3.73E-4;
    float c2 = 2.55E-2;
    float d = 2.33;
    float e = 1.56;

    float sw_epsilon_1 = scale * sqrt(b1 + c1 * alpha * alpha);
    float sw_epsilon_2 = scale * sqrt(b2 + c2 * alpha * alpha);
    float sw_A = d + e * cos(alpha) * cos(alpha);

    float Xres = 0;
    Matrix4f Wij;
    for (int j = evc->elements.size(); j --> 0;) //BACKWARD
    {
        Wij = evc->elements[j].Wij;

        evc->elements[j].pw = (nominalp) * sqrt( sw_A /
                (Xres/sw_epsilon_1 + pow(Xres/sw_epsilon_2,2) + 1) );
    /*
        Xres += (evc->elements[j].voxel->density * 1.E6 < 2.5) ?
                Wij[0] * evc->elements[j].voxel->density * 40000 :
                Wij[0] * 2.5 * 0.04;
        */
        Xres += Wij(0,0) * evc->elements[j].voxel->Value * 40000; //<- previous version of
                                      // p_weight: bugged but enhancing!
    }
}

static void pweigth_cw(Event *evc, Scalarf nominalp){
    float cw_A = 1.42857;
    float cw_epsilon = 50;
    float X0_tot = 0;
    Matrix4f Wij;
    for (int j = 0; j < evc->elements.size(); ++j)
    {
        Wij = evc->elements[j].Wij;
        /*X0_tot += (evc->elements[j].voxel->density * 1.E6 < 2.5 ) ?
                 Wij[0] * evc->elements[j].voxel->density * 40000 :
                 Wij[0] * 2.5 * 0.04;
         */
        X0_tot += Wij(0,0) * evc->elements[j].voxel->Value * 40000;
    }

    float _pw = cw_A * (nominalp) *
            sqrt( cw_epsilon / ( X0_tot + cw_epsilon ));
    for (int i = 0; i < evc->elements.size(); ++i )
    {
        evc->elements[i].pw = _pw;
    }
}




