#ifndef NDEBUG
#include "TTree.h"
#include "TFile.h"
#include "TH1F.h"
#endif

#include <vector>
#include "IBNormalPlaneMinimizationVariablesEvaluator.h"

class IBNormalPlaneMinimizationVariablesEvaluatorPimpl {

public:

    IBNormalPlaneMinimizationVariablesEvaluatorPimpl() {
        m_tracer    = NULL;

#ifndef NDEBUG
        m_out  = new TFile("20040_aXY.root", "RECREATE");
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
        m_tree->Branch("alpha",     &m_alpha,            "alpha/F");
        m_tree->Branch("Phi",       &t_phi,              "phi/F");
        m_tree->Branch("Theta",     &t_theta,            "theta/F");
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

        m_Data        = this->evaluateVariables(  m_muon.LineIn(), m_muon.LineOut());

        if (unlikely((fabs(m_Data(0)) > 1)||(fabs(m_Data(2)) > 1))) // << HARDCODED!!!
            m_integrity = false;

        m_ErrorMatrix = this->evaluateErrorMatrix(m_muon.LineIn(), m_muon.LineOut());


        if (unlikely((m_ErrorMatrix(0,0)>0.03 || m_ErrorMatrix(1,1)>1500 ||
                      m_ErrorMatrix(2,2)>0.03 || m_ErrorMatrix(2,2)>1500 ))) // << HARDCODED
            m_integrity = false;

#ifndef NDEBUG
        if (m_integrity) {
            EV.mx_in = m_muon.LineIn().direction(0);
            EV.my_in = m_muon.LineIn().direction(1);
            EV.mz_in = m_muon.LineIn().direction(2);
            Matrix4f m = this->getRotationMatrix(muon.LineIn().direction);

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
            HPoint3f projected  = projectOnContainer(m_muon.LineOut());
            HVector3f diff = m_muon.LineIn().origin-projected;
            float scal = diff.transpose()*n;
            HVector3f b = diff-scal*n;
            DT.poutLinNorm = b.head(3).norm();

            m_tree->Fill();
        }
#endif

        return m_integrity;
    }

    Vector4f evaluateVariables(const HLine3f &ingoing_track, const HLine3f &outgoing_track)
    {
        Matrix4f  rotation_matrix = this->getRotationMatrix(ingoing_track.direction);

        HPoint3f  prj  = this->projectOnContainer(outgoing_track);

        HVector3f disp = rotation_matrix * (prj - ingoing_track.origin);
        HVector3f scat = rotation_matrix * this->getDirectorCosines(outgoing_track.direction);

        Scalarf scat_x = atan2(scat(0),scat(1));
        Scalarf scat_z = atan2(scat(2),scat(1));

        Vector4f out(scat_x, disp(0), scat_z, disp(2));

        return out;
    }


    Matrix4f evaluateErrorMatrix(const HLine3f &ingoing_track, const HLine3f &outgoing_track)
    {
        Matrix4f covariance_p,
                 covariance_m;
        covariance_p << 0,0,0,0,
                        0,0,0,0,
                        0,0,0,0,
                        0,0,0,0;
        if (!m_integrity) return covariance_p;
        covariance_m << 0,0,0,0,
                        0,0,0,0,
                        0,0,0,0,
                        0,0,0,0;

        //////////
        // PLUS //
        //////////
        Vector4f plus[6];
        //////////
        //  in  //
        //////////
        for (int i=0; i<3; ++i) {
            if (m_muon.ErrorIn().direction_error(i)==0.f) {
                plus[i] = m_Data;
                continue;
            }
            HLine3f in_variations = ingoing_track;
            in_variations.direction(i) += m_muon.ErrorIn().direction_error(i);
            plus[i] = evaluateVariables(in_variations, outgoing_track);
        }
        /////////
        // out //
        /////////
        for (int i=0; i<3; ++i) {
            if (m_muon.ErrorOut().direction_error(i)==0.f) {
                plus[i+3] = m_Data;
                continue;
            }
            HLine3f out_variations = outgoing_track;
            out_variations.direction(i) += m_muon.ErrorOut().direction_error(i);
            plus[i+3] = evaluateVariables(ingoing_track, out_variations);
        }
        Vector4f d_plus[6]; // Vector4f indexes run over j, array elements over i
        for (int j=0; j<6; ++j){
            d_plus[j] = plus[j]-m_Data;
        }

        for (int i=0; i<4; ++i){
            for (int k=0; k<=i; ++k){
                for (int j=0; j<6; ++j){
                    covariance_p(i,k)+=d_plus[j](i)*d_plus[j](k);
                }
            }
        }

        ///////////
        // MINUS //
        ///////////
        Vector4f minus[6];
        //////////
        //  in  //
        //////////
        for (int i=0; i<3; ++i) {
            if (m_muon.ErrorIn().direction_error(i)==0.f) {
                minus[i] = m_Data;
                continue;
            }
            HLine3f in_variations = ingoing_track;
            in_variations.direction(i) -= m_muon.ErrorIn().direction_error(i);
            minus[i] = evaluateVariables(in_variations, outgoing_track);
        }
        /////////
        // out //
        /////////
        for (int i=0; i<3; ++i) {
            if (m_muon.ErrorOut().direction_error(i)==0.f) {
                minus[i+3] = m_Data;
                continue;
            }
            HLine3f out_variations = outgoing_track;
            out_variations.direction(i) -= m_muon.ErrorOut().direction_error(i);
            minus[i+3] = evaluateVariables(ingoing_track, out_variations);
        }
        Vector4f d_minus[6]; // Vector4f indexes run over j, array elements over i
        for (int j=0; j<6; ++j){
            d_minus[j] = m_Data-minus[j];
        }

        for (int i=0; i<4; ++i){
            for (int k=0; k<=i; ++k){
                for (int j=0; j<6; ++j){
                    covariance_m(i,k)+=d_minus[j](i)*d_minus[j](k);
                }
            }
        }

        Matrix4f covariance;
        covariance = 0.5*(covariance_m + covariance_p); // p/m averaging
        for (int i=0; i<4; ++i){
            for (int k=i+1; k<4; ++k){
                covariance(i, k) = covariance(k, i); // simmetrization
            }
        }
        return covariance;
    }

    HVector3f getDirectorCosines(const HVector3f &track_direction)
    {
        return track_direction / track_direction.head(3).norm();
    }

    Matrix4f getRotationMatrix(const HVector3f &track_direction)
    {
        // directors cosines
        HVector3f dc = this->getDirectorCosines(track_direction);
        // rotation angles
        t_theta = acos(dc(1));
        t_phi   = atan2(dc(2),dc(0));
        //alpha
        this->evaluateAlpha(t_phi,t_theta);
        // YZY rotations
        Matrix4f first_y_rotation = compileYRotation(-t_phi);
        Matrix4f first_z_rotation = compileZRotation(t_theta);
        Matrix4f secnd_y_rotation = compileYRotation(m_alpha); //<<< seems legit
        Matrix4f out = secnd_y_rotation * first_z_rotation * first_y_rotation;
        return out;
    }

    Matrix4f compileYRotation(Scalarf angle)
    {
        Matrix4f out;
        out << cos(angle), 0, -sin(angle), 0,
               0,          1,  0,          0,
               sin(angle), 0,  cos(angle), 0,
               0,          0,  0,          1;
        return out;
    }

    Matrix4f compileZRotation(Scalarf angle)
    {
        Matrix4f out;
        out <<  cos(angle), -sin(angle), 0, 0,
                sin(angle),  cos(angle), 0, 0,
                0,           0,          1, 0,
                0,           0,          0, 1;
        return out;
    }

    HPoint3f projectOnContainer(const HLine3f &muon_out_track)
    {
        HPoint3f pt;
        if (!m_tracer->GetExitPoint(muon_out_track, pt)) m_integrity = false;
        return pt;
    }

    void evaluateAlpha(Scalarf& phi, Scalarf& theta) {
        // if configuration X, Z or BOTH choose combination!
        //        m_alpha = (getAlphaZ(phi,theta)+getAlphaX(phi,theta))/2.;
        m_alpha = getAlphaX(phi,theta);
        //m_alpha = 0;
    }

    Scalarf getAlphaX(Scalarf& phi, Scalarf& theta) {
        Scalarf den = sqrt(cos(phi)*cos(phi)+(sin(phi)*sin(phi))*(cos(theta)*cos(theta)));
        Scalarf a1  = (-sin(phi)*cos(theta))/den;
        Scalarf a2  = cos(phi)/den;
        return atan2(a1,a2);
    }

    Scalarf getAlphaZ(Scalarf& phi, Scalarf& theta) {
        Scalarf den = sqrt(cos(phi)*cos(phi)+(sin(phi)*sin(phi))/(cos(theta)*cos(theta)));
        Scalarf a1  = (-sin(phi)/cos(theta))/den;
        Scalarf a2  = cos(phi)/den;
        return atan2(a1,a2);
    }

public:

    Scalarf         m_alpha;
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

IBNormalPlaneMinimizationVariablesEvaluator::IBNormalPlaneMinimizationVariablesEvaluator() :
    d(new IBNormalPlaneMinimizationVariablesEvaluatorPimpl) {}

IBNormalPlaneMinimizationVariablesEvaluator::~IBNormalPlaneMinimizationVariablesEvaluator()
{

#ifndef NDEBUG
    d->m_out->cd();
    d->m_tree->Write();
    d->m_out->Close();
#endif

    delete d;
}

bool IBNormalPlaneMinimizationVariablesEvaluator::evaluate(MuonScatterData muon)
{
    d->evaluate(muon);
    return d->m_integrity;
}

Vector4f IBNormalPlaneMinimizationVariablesEvaluator::getDataVector()
{
    return d->m_Data;
}

Scalarf IBNormalPlaneMinimizationVariablesEvaluator::getDataVector(int i)
{
    if (unlikely(i<0||i>3)) {
        //printf("Wrong index request for MinVar Data Vector - Aborting");
        exit(1);
    }
    return d->m_Data(i);
}

Matrix4f IBNormalPlaneMinimizationVariablesEvaluator::getCovarianceMatrix()
{
    return d->m_ErrorMatrix;
}

Scalarf IBNormalPlaneMinimizationVariablesEvaluator::getCovarianceMatrix(int i, int j)
{
    if (unlikely(i<0||i>3||j<0||j>3)) {
        //printf("Wrong index request for MinVar Covariance Matrix - Aborting");
        exit(1);
    }
    return d->m_ErrorMatrix(i,j);
}

void IBNormalPlaneMinimizationVariablesEvaluator::setRaytracer(IBVoxRaytracer *tracer)
{
    d->m_tracer = tracer;
}
