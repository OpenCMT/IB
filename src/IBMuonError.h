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
    void azimuthalMomentumCorrection(bool enable=true);
    void averageMomentumCorrection(bool enable=true);

private:
    void setScrapsImage(IBLightCollection &image, bool evPM = false);
    class IBMEShader
    {
        friend class IBMuonError;
    private:
        IBMEShader(IBMuonError * ref) { d = ref; }
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
    Scalarf      m_Ax,m_Az;
    Scalarf      m_pratio;
    bool         m_azimPcorr;
    bool         m_averPcorr;
    IBMEShader * m_shader;
    IBMESimpler* m_simpler;

};



#endif // IBMUONERROR_H
