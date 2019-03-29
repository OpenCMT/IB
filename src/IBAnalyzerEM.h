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


#ifndef IBANALYZEREM_H
#define IBANALYZEREM_H

#include "IBAnalyzer.h"
#include "IBVoxCollection.h"
#include "IBVoxRaytracer.h"
#include "IBVoxel.h"
#include <string>
#include <iomanip>

class IBPocaEvaluator;
class IBMinimizationVariablesEvaluator;
class IBAnalyzerEMAlgorithm;


class IBAnalyzerEM : public IBAnalyzer {
    uLibTypeMacro(IBAnalyzerEM,IBAnalyzer)

public:
    struct Event {
        struct Element {
            Matrix2f Wij;
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
            Scalarf  pTrue; // SV 20160921 variable to store useful quantities to study....
        } header;
        Vector<Element> elements;
    };


    properties()
    {
        Scalarf nominal_momentum;
        Scalarf SijCutEM;
    };

public:
    IBAnalyzerEM(IBVoxCollection &voxels, int nPath=2, double alpha=0., bool doRecoPath=true,
         bool oldTCalculation=false, float rankLimit=-100., IBVoxCollection* initialSqrPfromVtk=NULL, int pVoxelMean=0);
    ~IBAnalyzerEM();

    bool AddMuon(const MuonScatterData &muon);//{ return false;}
    bool AddMuonFullPath(const MuonScatterData &muon, Vector<HPoint3f>& muonPath);

    void SetMuonCollection(IBMuonCollection *muons);

    unsigned int Size();

    void Run(unsigned int iterations, float muons_ratio);

    void SetMLAlgorithm(IBAnalyzerEMAlgorithm *MLAlgorithm);

    uLibGetSetMacro(PocaAlgorithm,IBPocaEvaluator *)
    uLibGetSetMacro(VarAlgorithm,IBMinimizationVariablesEvaluator *)
    uLibGetSetMacro(RayAlgorithm,IBVoxRaytracer *)
    uLibGetSetMacro(UpdateAlgorithm,IBAbstract::IBVoxCollectionUpdateAlgorithm *)

    void filterEventsVoxelMask();

    void filterEventsLineDistance(float min, float max);

    void SijCut(float threshold);

    Vector<Event > SijCutCount(float threshold_low, float threshold_high);

    void dumpEventsSijInfo(const char *filename, Vector<float> N);

    void SijGuess(Vector<Vector2f> tpv);

    void Chi2Cut(float threshold);

    void SetVoxCollection(IBVoxCollection *voxels);

    void SetVoxcollectionShift(Vector3f shift);

    void dumpEventsTTree(const char *filename);
    void DumpP(const char *filename, float x0 = 0, float x1 = 10);
    void DumpEvent(Event *evc);

    Vector<Event> & Events();

    float SijMedian(const Event &evc);

    void SetSijMedianMomentum();

private:

    void Project(Event *evc);
    void BackProject(Event *evc);
    void Evaluate(float muons_ratio);

    IBPocaEvaluator                            *m_PocaAlgorithm;
    IBMinimizationVariablesEvaluator           *m_VarAlgorithm;
    IBVoxRaytracer                             *m_RayAlgorithm;
    IBAbstract::IBVoxCollectionUpdateAlgorithm *m_UpdateAlgorithm;
    friend class IBAnalyzerEMPimpl;
    class IBAnalyzerEMPimpl *m_d;

    int m_nPath;            //---- Int to indicate whether to build a 1, 2 or 3-path
    double m_alpha;         //---- Relative distance along the trajectory to build the 3-path
    bool m_useRecoPath;     //---- Use the true muon path from MC
    bool m_oldTCalculation; //---- Use the old method of calculating length parameter T
    float m_rankLimit;

    IBVoxCollection* m_initialSqrPfromVtk;
    int m_pVoxelMean;   //---- compute p voxel by hand
    IBVoxCollection m_imgMC;

    IBAnalyzerEMAlgorithm       *m_SijAlgorithm;
    Vector<IBAnalyzerEM::Event>  m_Events;
    bool                         m_firstIteration;

};

inline void IBAnalyzerEM::init_properties() {
    $_init();
    $$.nominal_momentum = 3;
}


#endif // IBANALYZEREM_H
