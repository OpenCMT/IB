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

    IBMuonError(Scalarf xA, Scalarf zA, Scalarf ratio = 1);
    ~IBMuonError();

    bool evaluate(MuonScatter &event, int i, int j);
    void setScrapsImage(IBLightCollection &image, bool evPM = false);

private:

    class IBMEShader
    {
        friend class IBMuonError;
    private:
        IBMEShader(IBMuonError * ref) { d = ref; }
//        ~IBMEShader() { delete m_pproc; delete m_tracer; }
        bool evaluate(MuonScatter &event, int i, int j);
    private:
        bool               m_evPM;
        IBLightCollection *m_image;
        IBVoxRaytracer    *m_tracer;
        IBLineDistancePocaEvaluator   *m_pproc;
        IBMuonError       *d;

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
    Scalarf      m_Axi,m_Azi,m_Axo,m_Azo;
    IBMEShader * m_shader;
    IBMESimpler* m_simpler;

};



#endif // IBMUONERROR_H
