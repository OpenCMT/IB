#ifndef IBANALYZEREM_H
#define IBANALYZEREM_H

#include "IBAnalyzer.h"
#include "IBVoxRaytracer.h"
#include "IBVoxel.h"

class IBPocaEvaluator;
class IBMinimizationVariablesEvaluator;
class IBAnalyzerEMAlgorithm;


class IBAnalyzerEM : public IBAnalyzer {
    typedef IBAnalyzer BaseClass;

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
        } header;
        Vector<Element> elements;
    };


    ULIB_OBJECT_PARAMETERS(BaseClass)
    {
        Scalarf nominal_momentum;        
    };

public:
    IBAnalyzerEM(IBVoxCollection &voxels);
    ~IBAnalyzerEM();

    virtual void AddMuon(const MuonScatterData &muon);

    void SetMuonCollection(IBMuonCollection *muons);

    unsigned int Size();

    void Run(unsigned int iterations, float muons_ratio);

    void SetMLAlgorithm(IBAnalyzerEMAlgorithm *MLAlgorithm);

    uLibGetSetMacro(PocaAlgorithm,IBPocaEvaluator *)
    uLibGetSetMacro(VarAlgorithm,IBMinimizationVariablesEvaluator *)
    uLibGetSetMacro(RayAlgorithm,IBVoxRaytracer *)

    void SijCut(float threshold);

    void SetVoxCollection(IBVoxCollection *voxels);

    void AddVoxcollectionShift(Vector3f shift);

    void DumpP(const char *filename, float x0 = 0, float x1 = 10);

    Vector<Event> &Events();

private:    
    IBPocaEvaluator                  *m_PocaAlgorithm;
    IBMinimizationVariablesEvaluator *m_VarAlgorithm;
    IBVoxRaytracer                   *m_RayAlgorithm;
    friend class IBAnalyzerEMPimpl;
    class IBAnalyzerEMPimpl *d;
};


inline void IBAnalyzerEM::init_parameters() {
    Parameters &p = this->parameters();
    p.nominal_momentum = 3;
}


#endif // IBANALYZEREM_H
