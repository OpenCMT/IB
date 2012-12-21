#ifndef IBVOXIMAGEFILTERPLASMON_HPP
#define IBVOXIMAGEFILTERPLASMON_HPP

#include "TFile.h"
#include "TH2F.h"
#include "Math/VoxImageFilter.h"

///////////////////////////////////////////////////////////////////////////////
///// VOXIMAGE FILTER OMOGENIZER //////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

using namespace uLib;

template <typename VoxelT>
class VoxImageFilterAlgorithmPlasmon :
        public VoxImageFilter<VoxelT, VoxImageFilterAlgorithmPlasmon<VoxelT> >
{

public:

    typedef VoxImageFilter<VoxelT, VoxImageFilterAlgorithmPlasmon<VoxelT> > BaseClass;

    VoxImageFilterAlgorithmPlasmon() : BaseClass(Vector3i(0,0,0)),
                                       m_Integrity(false)
    {}

    void SetMappingImage(VoxImage<VoxelT> * map) {
        m_Map = map;
        if (unlikely(m_Map->GetDims() != this->m_Image->GetDims())) {
            printf("Mapping image sizes do not match working image sizes! I'm not doing anything!");
            m_Integrity = false;
        } else m_Integrity = true;
    }

    void Run() {
        if(m_Integrity) {
            this->m_Image->operator /=(*m_Map);
        }
    }


private:

    VoxImage<VoxelT> * m_Map;
    bool m_Integrity;

    // dump useless function in order to be skipped at compile time
public:
    void SetKernelNumericXZY(const Vector<float> &numeric)    {}
    void SetKernelSpherical(float (*shape)(float))            {}
    template < typename ShapeT >
    void SetKernelSpherical( ShapeT shape )                   {}
protected:
    float Convolve(const VoxImage<VoxelT> &buffer, int index) {}
    void SetKernelOffset()                                    {}
    float Distance2(const Vector3i &v)                        {}
    // END DUMPING
};



#endif // IBVOXIMAGEFILTERPLASMON_HPP
