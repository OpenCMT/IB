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



#ifndef IBSIMPLETWOVIEWSMINIMIZATIONVARIABLESEVALUATOR_H
#define IBSIMPLETWOVIEWMINIMIZATIONVARIABLESEVALUATOR_H

#include "IBMinimizationVariablesEvaluator.h"

using namespace uLib;

class IBSimpleTwoViewsMinimizationVariablesEvaluator : public IBMinimizationVariablesEvaluator
{
public:

    IBSimpleTwoViewsMinimizationVariablesEvaluator();
    ~IBSimpleTwoViewsMinimizationVariablesEvaluator();

    bool evaluate(MuonScatterData muon);

    Vector4f getDataVector();
    Scalarf  getDataVector(int i);
    Matrix4f getCovarianceMatrix();
    Scalarf  getCovarianceMatrix(int i, int j);
    void setRaytracer(IBVoxRaytracer *tracer);
    void setDisplacementScatterOnly(bool,bool,bool);

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


    Vector4f evaluateVariables();
    Matrix4f evaluateErrorMatrix();
    HVector3f getDirectorCosines(const HVector3f &track_direction);
    void projectionOnContainer(HPoint3f &in, HPoint3f &out);

    Scalarf         t_phi;
    Scalarf         t_theta;
    VoxRaytracer*   m_tracer;
    MuonScatterData m_muon;
    bool            m_integrity;
    bool            m_displacementOnly;
    bool            m_scatterOnly;

    Vector4f        m_Data;
    Matrix4f        m_ErrorMatrix;

#ifndef NDEBUG
    TFile* m_out;
    TTree* m_tree;
#endif

};

#endif // IBSimpleTwoViewsMINIMIZATIONVARIABLESEVALUATOR_H
