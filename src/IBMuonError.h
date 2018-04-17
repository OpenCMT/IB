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



#ifndef IBMUONERROR_H
#define IBMUONERROR_H

#include <Detectors/MuonError.h>
#include <Detectors/MuonScatter.h>
#include <Math/VoxImage.h>
#include <IBVoxRaytracer.h>
#include <IBVoxel.h>

class IBPocaEvaluator;
class IBLineDistancePocaEvaluator;

using namespace uLib;

class IBMuonError
{
    friend class IBMEShader;
    friend class IBMESimpler;
public:

    IBMuonError(Scalarf xA, Scalarf zA, Scalarf xB = 0, Scalarf zB = 0, Scalarf ratio = 1);
    ~IBMuonError();

    bool evaluate(MuonScatter &event, int i, int j);
    void azimuthalMomentumCorrection(bool enable=true);
    void crossChamberErrorCorrection(bool enable=true);
    void averageMomentumCorrection(bool enable=true);
    void setOutMomentum(bool enable=true);
    void squareError(bool enable=true);

    void setScrapsImage(IBLightCollection &image);

private:
    class IBMEShader
    {
        friend class IBMuonError;
    private:
        IBMEShader(IBMuonError * ref) { d = ref; }
        bool evaluate(MuonScatter &event, int i, int j);
    private:
        IBLightCollection           *m_image;
        IBVoxRaytracer              *m_tracer;
        IBLineDistancePocaEvaluator *m_pproc;
        IBMuonError                 *d;

    };

    class IBMESimpler
    {
        friend class IBMuonError;
    private:
        IBMESimpler(IBMuonError * ref) { d = ref; }
        bool evaluate(MuonScatter &event, int i, int j);
    private:
        IBMuonError *d;
    };

    Scalarf mpdEval(Scalarf a, Scalarf p, Scalarf d);
    Scalarf mpdSquareEval(Scalarf a, Scalarf b, Scalarf p, Scalarf d);
    Scalarf      m_Ax,m_Az;
    Scalarf      m_Bx,m_Bz;
    Scalarf      m_pratio;
    bool         m_azimPcorr;
    bool         m_averPcorr;
    bool         m_usePout;
    bool         m_chamberErcorr;
    bool         m_squareError;
    IBMEShader * m_shader;
    IBMESimpler* m_simpler;

};



#endif // IBMUONERROR_H
