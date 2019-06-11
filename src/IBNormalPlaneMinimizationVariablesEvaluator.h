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



#ifndef IBNORMALPLANEMINIMIZATIONVARIABLESEVALUATOR_H
#define IBNORMALPLANEMINIMIZATIONVARIABLESEVALUATOR_H

#include "Core/Object.h"
#include "IBMinimizationVariablesEvaluator.h"

using namespace uLib;

class IBNormalPlaneMinimizationVariablesEvaluator : public IBMinimizationVariablesEvaluator
{
    uLibTypeMacro(IBNormalPlaneMinimizationVariablesEvaluator,IBMinimizationVariablesEvaluator)

public:

    IBNormalPlaneMinimizationVariablesEvaluator();
    ~IBNormalPlaneMinimizationVariablesEvaluator();

    bool evaluate(MuonScatterData muon);

    Vector4f getDataVector();
    Scalarf  getDataVector(int i);
    Matrix4f getCovarianceMatrix();
    Scalarf  getCovarianceMatrix(int i, int j);
    void setRaytracer(IBVoxRaytracer *tracer);
    void setDisplacementScatterOnly(bool,bool,bool);

    inline Scalarf GetAlphaXZ() const { return alphaXZ; }
    inline void SetAlphaXZ(Scalarf alpha) { alphaXZ = alpha; }
    inline bool GetUse_free_rotation() const { return use_free_rotation; }
    inline void SetUse_free_rotation(bool free_rot) { use_free_rotation = free_rot; }

private:

#ifndef NDEBUG
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


    Vector4f evaluateVariables(const HLine3f &ingoing_track, const HLine3f &outgoing_track);
    Matrix4f evaluateErrorMatrix(const HLine3f &ingoing_track, const HLine3f &outgoing_track);
    Matrix4f getRotationMatrix(const HVector3f &track_direction);
    HPoint3f projectOnContainer(const HLine3f &muon_out_track);

    Matrix4f compileYRotation(Scalarf angle);
    Matrix4f compileYRotation(Vector2f v);
    Matrix4f compileZRotation(Scalarf angle);
    Matrix4f compileZRotation(Vector2f v);
    inline HVector3f getDirectorCosines(const HVector3f &track_direction);

    bool            use_free_rotation;
    Scalarf         alphaXZ;
    Scalarf         m_alpha;
    Scalarf         t_phi;
    Scalarf         t_theta;
    VoxRaytracer*   m_tracer;
    MuonScatterData m_muon;
    bool            m_integrity;
    bool            m_scatterOnly;
    bool            m_displacementOnly;
    bool            m_oneD;

    Vector4f        m_Data;
    Matrix4f        m_ErrorMatrix;

#ifndef NDEBUG
    TFile*          m_out;
    TTree*          m_tree;
#endif

};


#endif // IBNORMALPLANEMINIMIZATIONVARIABLESEVALUATOR_H
