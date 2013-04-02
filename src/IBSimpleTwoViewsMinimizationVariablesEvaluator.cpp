#ifndef NDEBUG
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#endif

#include <vector>
#include "IBSimpleTwoViewsMinimizationVariablesEvaluator.h"

class IBSimpleTwoViewsMinimizationVariablesEvaluatorPimpl {

public:

    IBSimpleTwoViewsMinimizationVariablesEvaluatorPimpl() {
        m_tracer    = NULL;

#ifndef NDEBUG
        m_out  = new TFile("1363_stat_old.root", "RECREATE");
        m_tree = new TTree("stat", "Variables_Statistics");
        m_tree->Branch("dx",        &m_Data(1),          "dx/F");
        m_tree->Branch("dz",        &m_Data(3),          "dz/F");
        m_tree->Branch("sx",        &m_Data(0),          "sx/F");
        m_tree->Branch("sz",        &m_Data(2),          "sz/F");
        m_tree->Branch("S2dx",      &m_ErrorMatrix(1,1), "sdx/F");
        m_tree->Branch("S2dz",      &m_ErrorMatrix(3,3), "sdz/F");
        m_tree->Branch("S2sx",      &m_ErrorMatrix(0,0), "ssx/F");
        m_tree->Branch("S2sz",      &m_ErrorMatrix(2,2), "ssz/F");
        m_tree->Branch("Csxdx",     &m_ErrorMatrix(0,1), "sxdx/F");
        m_tree->Branch("Csxsz",     &m_ErrorMatrix(0,2), "sxsz/F");
        m_tree->Branch("Csxdz",     &m_ErrorMatrix(0,3), "sxdz/F");
        m_tree->Branch("Cszdx",     &m_ErrorMatrix(1,2), "szdx/F");
        m_tree->Branch("Cdxdz",     &m_ErrorMatrix(1,3), "dxdz/F");
        m_tree->Branch("Cszdz",     &m_ErrorMatrix(2,3), "szdz/F");
        m_tree->Branch("DisplN",    &DT.displNorm,       "dn/F");
        m_tree->Branch("PoutLinN",  &DT.poutLinNorm,     "polin/F");
        m_tree->Branch("mx_in",     &EV.mx_in,           "mx_in/F");
        m_tree->Branch("my_in",     &EV.my_in,           "my_in/F");
        m_tree->Branch("mz_in",     &EV.mz_in,           "mz_in/F");
        m_tree->Branch("chi2_p",    &chi2.p,             "chi2p/F");
        m_tree->Branch("chi2_t",    &chi2.t,             "chi2t/F");
        m_tree->Branch("chi2_x",    &chi2.x,             "chi2x/F");
        m_tree->Branch("chi2_z",    &chi2.z,             "chi2z/F");
        m_tree->Branch("chi2_px",   &chi2.px,            "chi2px/F");
        m_tree->Branch("chi2_tz",   &chi2.tz,            "chi2tz/F");
        m_tree->Branch("chi2_pt",   &chi2.pt,            "chi2pt/F");
        m_tree->Branch("chi2_xz",   &chi2.xz,            "chi2xz/F");
        m_tree->Branch("chi2_pxtz", &chi2.pxtz,          "chi2pxtz/F");
#endif

    }


    bool evaluate(MuonScatterData muon) {

        m_muon = muon;

        m_integrity = true;

        m_ErrorMatrix << 0,0,0,0,
                         0,0,0,0,
                         0,0,0,0,
                         0,0,0,0;

        m_Data        << 0,0,0,0;

        m_Data        = this->evaluateVariables();

//        if (unlikely((fabs(m_Data(0)) > 1)||(fabs(m_Data(2)) > 1))) // << HARDCODED!!!
//            m_integrity = false;

        m_ErrorMatrix = this->evaluateErrorMatrix();

#ifndef NDEBUG
        if (m_integrity) {
            EV.mx_in = m_muon.LineIn().direction(0);
            EV.my_in = m_muon.LineIn().direction(1);
            EV.mz_in = m_muon.LineIn().direction(2);

            // pxtz
            Matrix4f e4_inv = Matrix4f::Zero();
            bool check = false;
            m_ErrorMatrix.computeInverseWithCheck(e4_inv, check);
            chi2.pxtz = m_Data.transpose() * e4_inv * m_Data;

            // p
            chi2.p = m_Data(0)*m_Data(0)/m_ErrorMatrix(0,0);

            // t
            chi2.t = m_Data(2)*m_Data(2)/m_ErrorMatrix(2,2);

            // x
            chi2.x = m_Data(1)*m_Data(1)/m_ErrorMatrix(1,1);

            // z
            chi2.z = m_Data(3)*m_Data(3)/m_ErrorMatrix(3,3);

            // px
            Matrix2f e2_inv = Matrix2f::Zero();
            Matrix2f e2     = Matrix2f::Zero();
            Vector2f data   = Vector2f::Zero();
            e2   << m_ErrorMatrix(0,0), m_ErrorMatrix(0,1),
                    m_ErrorMatrix(1,0), m_ErrorMatrix(1,1);
            data << m_Data(0), m_Data(1);
            e2.computeInverseWithCheck(e2_inv, check);
            chi2.px = data.transpose() * e2_inv * data;

            //tz
            e2_inv = Matrix2f::Zero();
            e2     = Matrix2f::Zero();
            data   = Vector2f::Zero();
            e2   << m_ErrorMatrix(2,2), m_ErrorMatrix(2,3),
                    m_ErrorMatrix(3,2), m_ErrorMatrix(3,3);
            data << m_Data(2), m_Data(3);
            e2.computeInverseWithCheck(e2_inv, check);
            chi2.tz = data.transpose() * e2_inv * data;

            //pt
            e2_inv = Matrix2f::Zero();
            e2     = Matrix2f::Zero();
            data   = Vector2f::Zero();
            e2    << m_ErrorMatrix(0,0), m_ErrorMatrix(0,2),
                     m_ErrorMatrix(2,0), m_ErrorMatrix(2,2);
            data <<  m_Data(0), m_Data(2);
            e2.computeInverseWithCheck(e2_inv, check);
            chi2.pt = data.transpose() * e2_inv * data;

            //xz
            e2_inv = Matrix2f::Zero();
            e2     = Matrix2f::Zero();
            data   = Vector2f::Zero();
            e2    << m_ErrorMatrix(1,1), m_ErrorMatrix(1,3),
                     m_ErrorMatrix(3,1), m_ErrorMatrix(3,3);
            data <<  m_Data(1), m_Data(3);
            e2.computeInverseWithCheck(e2_inv, check);
            chi2.xz = data.transpose() * e2_inv * data;

            DT.displNorm = sqrt(m_Data(1)*m_Data(1)+m_Data(3)*m_Data(3));
            HVector3f n = getDirectorCosines(m_muon.LineIn().direction);
            HVector3f diff = m_muon.LineIn().origin-m_muon.LineOut().origin;
            float scal = diff.transpose()*n;
            HVector3f b = diff-scal*n;
            DT.poutLinNorm = b.head(3).norm();
            m_tree->Fill();
        }
#endif

        return m_integrity;
    }

    Vector4f evaluateVariables()
    {
        HVector3f disp;
        HVector3f dc  = getDirectorCosines(m_muon.LineIn().direction);
        HVector3f dco = getDirectorCosines(m_muon.LineOut().direction);
        float mx = (atan2(dc(0),dc(1)));
        float mz = (atan2(dc(2),dc(1)));
        float mxo = (atan2(dco(0),dco(1)));
        float mzo = (atan2(dco(2),dco(1)));
        float Lxy = m_muon.LineIn().direction.head(3).norm();
        HPoint3f in, out;
        projectionOnContainer(in, out);
        HVector3f oi = in-out;
        float dx = oi(0) - tan(mx)*oi(1);
        float dz = oi(2) - tan(mz)*oi(1);
        float dx_corr = cos(mx)*Lxy*(cos(mxo)/cos(mxo-mx));
        float dz_corr = cos(mz)*Lxy*(cos(mzo)/cos(mzo-mz));
        disp << dx*dx_corr, 0, dz*dz_corr, 0;
        float scat_x = mxo-mx;
        float scat_z = mzo-mz;
        Vector4f res(scat_x, disp(0), scat_z, disp(2));
        return res;
    }


    Matrix4f evaluateErrorMatrix()
    {
        Matrix4f covariance;
        float dy = m_muon.LineIn().origin(1) - m_muon.LineOut().origin(1);
        HVector3f dc   = getDirectorCosines(m_muon.LineIn().direction);
        float cosphi   = cos(atan2(dc(0),dc(1)));
        float costheta = cos(atan2(dc(2),dc(1)));
        float csx =  dy/cosphi;
        float csz =  dy/costheta;
        float s2x =  m_muon.ErrorIn().direction_error(0);
        float s2z = (m_muon.ErrorIn().direction_error(1)==0.f) ?
                     m_muon.ErrorIn().direction_error(2)       :
                     m_muon.ErrorIn().direction_error(1);
        s2x *= s2x;
        s2z *= s2z;
        covariance << 2*s2x,  csx*s2x,    0,  0,
                      csx*s2x,csx*csx*s2x,0,  0,
                      0,  0,              2*s2z,  csz*s2z,
                      0,  0,              csz*s2z,csz*csz*s2z;
        return covariance;
    }

    HVector3f getDirectorCosines(const HVector3f &track_direction)
    {
        return track_direction / track_direction.head(3).norm();
    }

    void projectionOnContainer(HPoint3f &in, HPoint3f &out)
    {
        if (!m_tracer->GetEntryPoint(m_muon.LineIn(),in) ||
            !m_tracer->GetExitPoint(m_muon.LineOut(),out))
            m_integrity=false;
    }

public:

    Scalarf         t_phi,
                    t_theta;
    VoxRaytracer*   m_tracer;
    MuonScatterData m_muon;
    bool            m_integrity;

    Vector4f m_Data;
    Matrix4f m_ErrorMatrix;

#ifndef NDEBUG
    TFile* m_out;
    TTree* m_tree;
    struct READ{
        float mx_in, my_in, mz_in;
    } EV;
    struct CHI {
        float p, t, x, z, px, tz, pxtz, pt, xz;
    } chi2;
    struct DaTa {
        float displNorm, poutLinNorm;
    } DT;

#endif


};

IBSimpleTwoViewsMinimizationVariablesEvaluator::IBSimpleTwoViewsMinimizationVariablesEvaluator() :
    d(new IBSimpleTwoViewsMinimizationVariablesEvaluatorPimpl) {}

IBSimpleTwoViewsMinimizationVariablesEvaluator::~IBSimpleTwoViewsMinimizationVariablesEvaluator()
{

#ifndef NDEBUG
    d->m_out->cd();
    d->m_tree->Write();
    d->m_out->Close();
#endif

    delete d;
}

bool IBSimpleTwoViewsMinimizationVariablesEvaluator::evaluate(MuonScatterData muon)
{
    d->evaluate(muon);
    return d->m_integrity;
}

Vector4f IBSimpleTwoViewsMinimizationVariablesEvaluator::getDataVector()
{
    return d->m_Data;
}

Scalarf IBSimpleTwoViewsMinimizationVariablesEvaluator::getDataVector(int i)
{
    if (unlikely(i<0||i>3)) {
        //printf("Wrong index request for MinVar Data Vector - Aborting");
        exit(1);
    }
    return d->m_Data(i);
}

Matrix4f IBSimpleTwoViewsMinimizationVariablesEvaluator::getCovarianceMatrix()
{
    return d->m_ErrorMatrix;
}

Scalarf IBSimpleTwoViewsMinimizationVariablesEvaluator::getCovarianceMatrix(int i, int j)
{
    if (unlikely(i<0||i>3||j<0||j>3)) {
        //printf("Wrong index request for MinVar Covariance Matrix - Aborting");
        exit(1);
    }
    return d->m_ErrorMatrix(i,j);
}

void IBSimpleTwoViewsMinimizationVariablesEvaluator::setRaytracer(IBVoxRaytracer *tracer)
{
    d->m_tracer = tracer;
}
