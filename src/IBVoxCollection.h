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


/* 
 * File:   IBVoxCollection.h
 * Author: andrea
 *
 * Created on June 14, 2011, 11:43 AM
 */


#ifndef IBVOXCOLLECTION_H
#define	IBVOXCOLLECTION_H

#include <Math/Dense.h>
#include <Math/ContainerBox.h>
#include <Math/VoxImage.h>

#include "IBVoxel.h"

class IBVoxCollection;

namespace IBInterface {

struct IBVoxCollectionStaticUpdateAlgorithm {
    template <class Self> void check_structural() {
        //uLibCheckStaticFunction(Self, UpdateDensity , void ,IBVoxCollection *, unsigned int);
    }
};

}

namespace IBAbstract {

struct IBVoxCollectionUpdateAlgorithm {
    virtual void operator()(IBVoxCollection*,unsigned int) = 0;
};

struct IBVoxCollectionMAPAlgorithm {
    virtual void UpdateDensity(IBVoxCollection *voxels, unsigned int threshold) = 0;
};

}


class IBVoxCollection : public uLib::VoxImage<IBVoxel> {
    typedef uLib::VoxImage<IBVoxel> BaseClass;
    typedef IBVoxCollection VoxImageType;
public:

    IBVoxCollection();

    IBVoxCollection(const uLib::Vector3i size);

    IBVoxCollection(const IBVoxCollection &copy);

    IBVoxCollection(const BaseClass &copy);

    // templated update for analyzer specific customizations //
    template < typename StaticUpdateAlgT >
    void UpdateDensity(unsigned int threshold);

    // templated update for analyzer specific customizations //
    template < typename UpdateAlgT >
    void UpdateDensity(UpdateAlgT &algorithm, unsigned int threshold);

    void SetMAPAlgorithm(IBAbstract::IBVoxCollectionMAPAlgorithm *algorithm);

    inline void InitLambda(const IBVoxel &value);

    inline void InitCount(unsigned int count);

    int CountLambdaOverThreshold(float threshold);

    int CountLambdaOverThreshold(float threshold,uLib::Vector3i boxp1, uLib::Vector3i boxp2);

    int CountLambdaOverThreshold(float threshold,uLib::HPoint3f center, uLib::HVector3f size);    

    inline IBVoxCollection LambdaToInvLrad(float p0) { IBVoxCollection out = *this;
                                                       out*=((p0/15)*(p0/15)); // 1/cm
                                                       out*=100.;              // 1/m
                                                       return out;
                                                     }

    inline IBVoxCollection InvLradToLambda(float p0) { IBVoxCollection out = *this;
                                                       out/=((p0/15)*(p0/15)); // 1/cm
                                                       //out/=100.;              // if in 1/m
                                                       return out;
                                                     }
    IBVoxCollection getMCImage(const char* file, int nsamples);
    float getVoxelMCDensity(const char* file, uLib::HPoint3f c1, uLib::HPoint3f c2, int nrandom=1000);

private:
    IBAbstract::IBVoxCollectionMAPAlgorithm *m_MAPAlgorithm;
};


// --- inlines -------------------------------------------------------------- //

inline void IBVoxCollection::InitLambda(const IBVoxel &value)
{
    BaseClass::InitVoxels(value);
}


inline void IBVoxCollection::InitCount(unsigned int count)
{
    for(unsigned int i=0; i<this->Data().size(); ++i) {
        this->Data().operator [](i).Count = count;
    }
}


// --- Update --------------------------------------------------------------- //

template < class UpdateAlgT >
void IBVoxCollection::UpdateDensity(UpdateAlgT &algorithm,
                                    unsigned int threshold) {

    // Analyzer Update //
    algorithm(this, threshold);

    // MAP update //
    if(m_MAPAlgorithm) m_MAPAlgorithm->UpdateDensity(this, threshold);

    // Reinitialize voxels //
    this->InitCount(0);
}

template < class StaticUpdateAlgT >
void IBVoxCollection::UpdateDensity(unsigned int threshold) {

    // Analyzer Update //
    StaticUpdateAlgT::UpdateDensity(this,threshold);

    // MAP update //
    if(m_MAPAlgorithm) m_MAPAlgorithm->UpdateDensity(this, threshold);

    // Reinitialize voxels //
    this->InitCount(0);
}


#endif	/* IBVOXCOLLECTION_H */

