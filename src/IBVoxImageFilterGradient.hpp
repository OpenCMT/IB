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
