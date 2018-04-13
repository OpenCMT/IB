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


#ifndef IBVOXIMAGEFILTERPLASMON_HPP
#define IBVOXIMAGEFILTERPLASMON_HPP

#include "TFile.h"
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
