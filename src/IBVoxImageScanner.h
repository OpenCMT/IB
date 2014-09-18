/*//////////////////////////////////////////////////////////////////////////////
// CMT Cosmic Muon Tomography project //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

  Copyright (c) 2014, Universita' degli Studi di Padova, INFN sez. di Padova

  Coordinators: Prof. Gianni Zumerle < gianni.zumerle@pd.infn.it >
                Paolo Checchia       < paolo.checchia@pd.infn.it >

  Authors: Andrea Rigoni Garola < andrea.rigoni@pd.infn.it >
           Matteo Furlan        < nuright@gmail.com >
           Sara Vanini          < sara.vanini@pd.infn.it >

  All rights reserved
  ------------------------------------------------------------------

  This file can not be copied and/or distributed without the express
  permission of  Prof. Gianni Zumerle  < gianni.zumerle@pd.infn.it >

//////////////////////////////////////////////////////////////////////////////*/



#ifndef IBVOXIMAGESCANNER_H
#define IBVOXIMAGESCANNER_H

#include "Math/VoxImage.h"

using namespace uLib;

template<class T>
class IBVoxImageScanner
{
public:
    IBVoxImageScanner() : m_Image(NULL)           {}
    IBVoxImageScanner(T * image) : m_Image(image) {}
    IBVoxImageScanner(T& image) : m_Image(&image) {}
    template<class A>
    class A::ScanData ScanImage(class A::ScanOption opt)
    {
        A scanner;
        scanner.template LoadImageData<T>(m_Image);
        return scanner.Scan(opt);
    }

    uLibGetSetMacro(Image, T*)

private:
    T * m_Image;
};

// TODO: interface for Vector<float> ImageData
// TODO: Threshold Scanner Interface
// TODO: in SubImageGrabber a mapper between V3i/int of Image and SubImage;
// move in a separate file to be included!!!

///////////////////////////////////////////////////////////////////////////////
//// SIMPLE THRESHOLD SCAN ////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
class SimpleThresholdScan
{
public:
    typedef Vector<float> ImageData;
    struct ScanData {
        float Percent;
        float Intensity;
    };
    struct ScanOption {
        float Threshold;
    };
    SimpleThresholdScan() {}

    template<class T>
    void LoadImageData(T * image);

    ScanData Scan(ScanOption opt);

    inline void ResetScanData()
    {
        m_ScanData.Percent   = 0.0;
        m_ScanData.Intensity = 0.0;
    }

    inline ScanData GetScanData() { return m_ScanData; }

protected:
    ImageData m_ImageData;
    ScanData  m_ScanData;
};

template<class T>
void SimpleThresholdScan::LoadImageData(T * image)
{
    for(int i=0; i<image->Data().size(); ++i)
        this->m_ImageData.push_back(image->At(i).Value);
}

// get rid of inline when moving to separate file!
inline SimpleThresholdScan::ScanData SimpleThresholdScan::Scan(
        SimpleThresholdScan::ScanOption opt)
{
    this->ResetScanData();
    float found = 0;
    for(int i=0; i<m_ImageData.size(); ++i) {
        if(m_ImageData.at(i)>=opt.Threshold) {
            found++;
            m_ScanData.Intensity += pow(m_ImageData.at(i)/opt.Threshold,2)-1.;
        }
    }
    m_ScanData.Intensity /= found;
    m_ScanData.Percent = 100 * ((float)found / m_ImageData.size());
    if(m_ScanData.Percent==0.f)
        m_ScanData.Intensity = 0.f;
    return m_ScanData;
}

///////////////////////////////////////////////////////////////////////////////
//// RANGE THRESHOLD SCAN /////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////

class RangeThresholdScan : public SimpleThresholdScan
{
public:
    typedef SimpleThresholdScan BaseClass;
    typedef Vector<BaseClass::ScanOption> ScanOption;
    typedef Vector<BaseClass::ScanData>   ScanData;

    RangeThresholdScan() : BaseClass() {}
    ScanData Scan(ScanOption opt);
    inline void ResetScanData() { m_ScanData.clear(); }

protected:
    ScanData m_ScanData;
};

RangeThresholdScan::ScanData RangeThresholdScan::Scan(
        RangeThresholdScan::ScanOption opt)
{
    this->ResetScanData();
    for(int i=0; i<opt.size(); ++i) {
        BaseClass::ScanData partial_data;
        partial_data = BaseClass::Scan(opt.at(i));
        m_ScanData.push_back(partial_data);
    }
    return m_ScanData;
}

#endif // IBVOXIMAGESCANNER_H
