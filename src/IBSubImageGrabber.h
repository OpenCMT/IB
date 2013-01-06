#ifndef IBSUBIMAGEGRABBER_H
#define IBSUBIMAGEGRABBER_H

#include "Math/VoxImage.h"

using namespace uLib;

struct Box { // better put this elsewhere so a "BoxSelector" can use it
    Vector3i Begins;
    Vector3i Ends;
};

template<class T> // T = VoxImage<>
class IBSubImageGrabber
{
public:
    IBSubImageGrabber() {}
    IBSubImageGrabber(T& image) : m_Image(&image) {}
    IBSubImageGrabber(T* image) : m_Image(image)  {}

    inline void SetImage(T& image) { this->m_Image = &image; }
    inline void SetImage(T* image) { this->m_Image = image ; }

    template<class K> // return VoxImage<> Type
    K GrabRegion(Vector3i start, Vector3i stop);

    template<class K>
    K GrabRegion(HPoint3f center, HVector3f size);

    template<class K>
    K GrabRegion(Box& voxBox);

private:
    T * m_Image;
};

// specializations

// TODO checks on T for type and on box border inside T
// TODO maybe move these function on the VoxImage<T> prototype and get rid of this class?
template<class T>
template<class K>
K IBSubImageGrabber<T>::GrabRegion(Box& voxBox)
{
    Vector3i sizer = (voxBox.Ends - voxBox.Begins) + Vector3i(1,1,1);
    K region(sizer);
    Vector3f spacing = this->m_Image->GetSpacing();
    region.SetSpacing(spacing);
    region.SetDataOrder(this->m_Image->GetDataOrder());
    Vector3f region_position = this->m_Image->GetPosition();
    Vector3f shift(spacing(0)*voxBox.Begins(0),
            spacing(1)*voxBox.Begins(1),
            spacing(2)*voxBox.Begins(2));
    region_position += shift;
    region.SetPosition(region_position);

    for(int i=0; i<sizer(0); ++i) {
        for (int j=0; j<sizer(1); ++j) {
            for (int k=0; k<sizer(2); ++k) {
                Vector3i box_position(i,j,k);
                Vector3i img_position = box_position + voxBox.Begins;
                region.SetValue(box_position, m_Image->GetValue(img_position));
            }
        }
    }
    return region;
}

template<class T>
template<class K>
K IBSubImageGrabber<T>::GrabRegion(Vector3i start, Vector3i stop)
{
    Box box;
    box.Begins = start;
    box.Ends   = stop;
    return this->GrabRegion<K>(box);
}

template<class T>
template<class K>
K IBSubImageGrabber<T>::GrabRegion(HPoint3f center, HVector3f size)
{
    Box box;
    box.Begins = this->m_Image->Find(center-size);
    box.Ends   = this->m_Image->Find(center+size);
    return this->GrabRegion<K>(box);
}

#endif // IBSUBIMAGEGRABBER_H
