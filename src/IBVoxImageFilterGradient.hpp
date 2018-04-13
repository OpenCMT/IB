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


#ifndef IBVOXIMAGEFILTERGRADIENT_HPP
#define IBVOXIMAGEFILTERGRADIENT_HPP


template <typename VoxelT>
class VoxImageFilterAlgorithmGradient :
        public VoxImageFilter<VoxelT, VoxImageFilterAlgorithmGradient<VoxelT> >
{

public:

    typedef VoxImageFilter<VoxelT, VoxImageFilterAlgorithmGradient<VoxelT> > BaseClass;

    VoxImageFilterAlgorithmGradient() : BaseClass(Vector3i(0,0,0)),
        m_Factor(1)
    {}

    void SetCoeff(Scalarf m) {
        m_Factor = m;
    }

    void Run() {
        Vector3i dims = this->m_Image->GetDims();
        Vector3i index;
        for(int z=0; z < dims(2); ++z) {
            for(int y=0; y < dims(1); ++y) {
                for(int x=0; x < dims(0); ++x) {
                    index << x,y,z;
                    this->m_Image->operator [](index).Value *=
                            1 + m_Factor * ((float)y/dims(1));
                }
            }
        }
    }


private:
    Scalarf m_Factor;

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




#endif // IVBOXIMAGEFILTERGRADIENT_HPP
