#ifndef IBMUONERROR_H
#define IBMUONERROR_H

#include <Detectors/MuonError.h>
#include <Detectors/MuonScatter.h>
#include <Math/VoxImage.h>
#include <IBVoxRaytracer.h>
#include <IBPocaEvaluator.h>
#include <IBVoxel.h>

using namespace uLib;

class IBMuonError
{
public:

    IBMuonError(Scalarf xA, Scalarf zA, Scalarf ratio = 1);

    bool evaluate(MuonScatter &event, int i, int j);
    void setScrapsImage(IBLightCollection * image, bool evPM = false);

private:
    friend class IBMEShader;
    friend class IBMESimpler;

    class IBMEShader
    {
    private:
        IBMEShader(IBMuonError * ref) { d = ref; }
        void setImage(IBLightCollection* image, bool evPM)
        {
            m_image = image;
            m_evPM  = evPM;
        }
        bool evaluate(MuonScatter &event, int i, int j);
    private:
        bool               m_evPM;
        IBLightCollection *m_image;
        IBVoxRaytracer    *m_tracer;
        IBPocaEvaluator   *m_pproc;
        IBMuonError       *d;

    };

    class IBMESimpler
    {
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
