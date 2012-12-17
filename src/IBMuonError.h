#ifndef IBMUONERROR_H
#define IBMUONERROR_H

#include <Detectors/MuonError.h>

using namespace uLib;

class IBMuonError
{
public:

    IBMuonError()
    {
        m_pX_in  = Eigen::Vector2f::Zero();
        m_pZ_in  = Eigen::Vector2f::Zero();
        m_pX_out = Eigen::Vector2f::Zero();
        m_pZ_out = Eigen::Vector2f::Zero();

    }

    IBMuonError(Scalarf p0_x, Scalarf p1_x, Scalarf p0_z, Scalarf p1_z, Scalarf ioRatio = 1.)
    {
        ioRatio  *= ioRatio;
        m_pX_in  << p0_x * p0_x, p1_x * p1_x;
        m_pZ_in  << p0_z * p0_z, p1_z * p1_x;
        m_pX_out << ioRatio * m_pX_in;
        m_pZ_out << ioRatio * m_pZ_in;
    }

    void evaluateInTrackError(Scalarf p, Scalarf tanPhi, Scalarf tanTheta)
    {
        m_In.Phi()   = 1e-3*sqrt((m_pX_in(0)/(p*p)+m_pX_in(1))/2)*(1+tanPhi*tanPhi);
        m_In.Theta() = 1e-3*sqrt((m_pZ_in(0)/(p*p)+m_pZ_in(1))/2)*(1+tanTheta*tanTheta);
    }

    void evaluateOutTrackError(Scalarf p, Scalarf tanPhi, Scalarf tanTheta, Scalarf dc = 0.f, Scalarf depth = 0.f, Scalarf rL = 0.f)
    {
        Scalarf corr = 0;
        //evaluate matter correction
        if (rL!=0.f && dc!=0.f){
          Scalarf length = fabs(depth * dc);
          corr  = (15*15/(p*p))*(length/rL);
        }
        m_Out.Phi()   = 1e-3*sqrt(((m_pX_out(0)/(p*p)+m_pX_out(1))/2)+corr)*(1+tanPhi*tanPhi);
        m_Out.Theta() = 1e-3*sqrt(((m_pZ_out(0)/(p*p)+m_pZ_out(1))/2)+corr)*(1+tanTheta*tanTheta);
    }

    uLibRefMacro(In,  MuonError)
    uLibRefMacro(Out, MuonError)
private:

    MuonError m_In;
    MuonError m_Out;
    Vector2f  m_pX_in, m_pX_out;
    Vector2f  m_pZ_in, m_pZ_out;

};



#endif // IBMUONERROR_H
