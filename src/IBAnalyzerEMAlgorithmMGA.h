/*////////////////////////////////////////////////////////////////////////////
 Copyright 2018 Istituto Nazionale di Fisica Nucleare

 Licensed under the EUPL, Version 1.2 or - as soon they will be approved by
 the European Commission - subsequent versions of the EUPL (the "Licence").
 You may not use this work except in compliance with the Licence.

 You may obtain a copy of the Licence at:

 https://joinup.ec.europa.eu/software/page/eupl

 Unless required by applicable law or agreed to in writing, software
 distributed under the Licence is distributed on an "AS IS" basis, WITHOUT
 WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
 Licence for the specific language governing permissions and limitations under
 the Licence.
////////////////////////////////////////////////////////////////////////////*/


#ifndef IBANALYZEREMALGORITHMMGA_H
#define IBANALYZEREMALGORITHMMGA_H

#include <Core/Vector.h>
#include "IBAnalyzerEMAlgorithm.h"


////////////////////////////////////////////////////////////////////////////////
/////  MGA INIT ABSTRACT CLASS /////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


template < int size = 1 >
class IBAnalyzerEMAlgorithmMGA : public IBAnalyzerEMAlgorithm {
    typedef struct IBAnalyzerEM::Event Event;
public:
    void SetGaussians(Vector<Vector2f> &WS);
    void SetGaussians(Scalarf *w, Scalarf *s);
    Scalarf GetNominalMomentum() const { return m_P; }

protected:
    virtual ~IBAnalyzerEMAlgorithmMGA() {}

    Scalarf Ki(Scalarf Chi2);
    Vector<Vector2f> normalizeScaling(Vector<Vector2f> &ws);

private:    
    Vector2f m_Kparams[size * size];
    Scalarf  m_P;
};



template < int size >
Vector<Vector2f>
IBAnalyzerEMAlgorithmMGA<size>::normalizeScaling(Vector<Vector2f> &ws)
{
    Vector<Vector2f> out = ws;

    float wsum=0;
    for(int i=0; i < size; ++i) {
        wsum += ws[i](0);
    }

    for(int i=0; i < size; ++i) {
        out[i](0) = ws[i](0)/wsum;
        out[i](1) = pow(ws[i](1),2);
    }
    return out;
}


template < int size >
void IBAnalyzerEMAlgorithmMGA<size>::SetGaussians(Scalarf *w, Scalarf *s)
{
    Vector<Vector2f> ws;
    for(int i=0 ; i<size ; ++i)
        ws.push_back(Vector2f(w[i], s[i]));
    this->SetGaussians(ws);
}

template < int size >
void IBAnalyzerEMAlgorithmMGA<size>::SetGaussians(Vector<Vector2f> &WS)
{
    Vector<Vector2f> ws = normalizeScaling(WS);
    float wsnorm = 0;
    for(int i=0; i<size; i++)
        wsnorm += ws[i].prod();

    for (int k=0 ; k<size; ++k) {
        float sk_cap = ws[k](1)/wsnorm;
        for (int l=0 ; l<size ; ++l) {
            float sl_cap = ws[l](1)/wsnorm;
            m_Kparams[k*size+l] <<
                    (sk_cap) * sqrt(sk_cap/sl_cap) * (ws[l](0)/ws[k](0)),
                    -0.5 * (1/sl_cap - 1/sk_cap);
        }
    }

}


template < int size >
Scalarf IBAnalyzerEMAlgorithmMGA<size>::Ki(Scalarf Chi2)
{
    Scalarf Ki = 0;
    float Ki_tmp = 0;
    // optimize even more: linearize for cycles //
    for (int k=0 ; k<size; ++k) {
        for (int l=0 ; l<size ; ++l) {
            Ki_tmp += this->m_Kparams[k*size+l](0) *
                    exp(Chi2 * m_Kparams[k*size+l](1));
        }
        Ki    += 1 / Ki_tmp;
        Ki_tmp = 0;
    }
    return Ki;
}







////////////////////////////////////////////////////////////////////////////////
//// ALGORITHMS  ///////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////






///// PXTZ /////////////////////////////////////////////////////////////////////

template < int size >
class IBAnalyzerEMAlgorithmMGA_PXTZ : public IBAnalyzerEMAlgorithmMGA<size> {
    typedef struct IBAnalyzerEM::Event Event;
public:
    void evaluate(Matrix4f &Sigma, Event *evc);

};

template < int size >
void IBAnalyzerEMAlgorithmMGA_PXTZ<size>::evaluate(Matrix4f &Sigma, Event *evc)
{
    Matrix4f iS = Sigma.inverse();
    Matrix4f Dn = iS * (evc->header.Di * evc->header.Di.transpose());
    Scalarf  Ki = this->Ki(Dn.trace());

    for (unsigned int j = 0; j < evc->elements.size(); ++j) {

        Matrix4f Wij = Matrix4f::Zero();
        Wij.block<2,2>(0,0) = evc->elements[j].Wij;
        Wij.block<2,2>(2,2) = evc->elements[j].Wij;
        Matrix4f Bn = iS * Wij;

        float Sij = Ki * (Bn * Dn).trace() - Bn.trace();
        evc->elements[j].Sij =  Sij * evc->elements[j].lambda *
                evc->elements[j].lambda * evc->elements[j].pw;
    }
}




///// PX   /////////////////////////////////////////////////////////////////////


template < int size >
class IBAnalyzerEMAlgorithmMGA_PX : public IBAnalyzerEMAlgorithmMGA<size> {
    typedef struct IBAnalyzerEM::Event Event;
public:
    void evaluate(Matrix4f &Sigma, Event *evc);
};


template < int size >
void IBAnalyzerEMAlgorithmMGA_PX<size>::evaluate(Matrix4f &Sigma, Event *evc) {

    Matrix2f iS;
    {
        Matrix2f S;
        S << Sigma(0,0), Sigma(0,1), Sigma(1,0), Sigma(1,1);
        iS = S.inverse();
    }
    Vector2f Di(evc->header.Di(0),evc->header.Di(1));

    Matrix2f Dn = iS * (Di * Di.transpose());
    Scalarf  Ki = this->Ki(Dn.trace());

    for (unsigned int j = 0; j < evc->elements.size(); ++j) {
        Matrix2f Bn = iS * evc->elements[j].Wij;
        float Sij = Ki * (Bn * Dn).trace() - Bn.trace();
        evc->elements[j].Sij =  Sij * evc->elements[j].lambda *
                evc->elements[j].lambda * evc->elements[j].pw;
    }
}


///// TZ   /////////////////////////////////////////////////////////////////////


template < int size >
class IBAnalyzerEMAlgorithmMGA_TZ : public IBAnalyzerEMAlgorithmMGA<size> {
    typedef struct IBAnalyzerEM::Event Event;
public:
    void evaluate(Matrix4f &Sigma, Event *evc);
};


template < int size >
void IBAnalyzerEMAlgorithmMGA_TZ<size>::evaluate(Matrix4f &Sigma, Event *evc) {

    Matrix2f iS;
    {
        Matrix2f S;
        S << Sigma(2,2), Sigma(2,3), Sigma(3,2), Sigma(3,3);
        iS = S.inverse();
    }
    Vector2f Di(evc->header.Di(2),evc->header.Di(3));


    Matrix2f Dn = iS * (Di * Di.transpose());
    Scalarf  Ki = this->Ki(Dn.trace());

    for (unsigned int j = 0; j < evc->elements.size(); ++j) {
        Matrix2f Bn = iS * evc->elements[j].Wij;
        float Sij = Ki * (Bn * Dn).trace() - Bn.trace();
        evc->elements[j].Sij =  Sij * evc->elements[j].lambda *
                evc->elements[j].lambda * evc->elements[j].pw;
    }
}



///// PT   /////////////////////////////////////////////////////////////////////



template < int size >
class IBAnalyzerEMAlgorithmMGA_PT : public IBAnalyzerEMAlgorithmMGA<size> {
    typedef struct IBAnalyzerEM::Event Event;
public:
    void evaluate(Matrix4f &Sigma, Event *evc);
};



template < int size >
void IBAnalyzerEMAlgorithmMGA_PT<size>::evaluate(Matrix4f &Sigma, Event *evc) {
    Matrix2f iS;
    {
        Matrix2f S;
        S << Sigma(0,0), Sigma(0,2), Sigma(2,0), Sigma(2,2);
        iS = S.inverse();
    }
    Vector2f Di(evc->header.Di(0),evc->header.Di(2));
    Matrix2f Dn = iS * (Di * Di.transpose());
    Scalarf  Ki = this->Ki(Dn.trace());

    for (unsigned int j = 0; j < evc->elements.size(); ++j) {
        Matrix2f Wij;
        Wij << evc->elements[j].Wij(0,0),0,
               0,evc->elements[j].Wij(0,0);
        Matrix2f Bn = iS * Wij;
        float Sij = Ki * (Bn * Dn).trace() - Bn.trace();
        evc->elements[j].Sij =  Sij *
                evc->elements[j].lambda * evc->elements[j].lambda *
                evc->elements[j].pw;
    }
}





#endif // IBANALYZEREMALGORITHMMGA_H
