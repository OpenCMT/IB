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
struct IBVoxCollectionUpdateAlgorithm {
    template <class Self> void check_structural() {
        uLibCheckFunction(Self,operator(),void,IBVoxCollection *, unsigned int);
    }
};

struct IBVoxCollectionStaticUpdateAlgorithm {
    template <class Self> void check_structural() {
        //uLibCheckStaticFunction(Self, UpdateDensity , void ,IBVoxCollection *, unsigned int);
    }
};
}

namespace IBAbstract {
class IBVoxCollectionMAPAlgorithm {
public:
    virtual void UpdateDensity(IBVoxCollection *voxels, unsigned int threshold) = 0;
};
}


class IBVoxCollection : public uLib::VoxImage<IBVoxel> {
    typedef uLib::VoxImage<IBVoxel> BaseClass;
    typedef IBVoxCollection VoxImageType;
public:

    IBVoxCollection(const uLib::Vector3i size);

    IBVoxCollection(const IBVoxCollection &copy);

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

    // check structural for UpdateAlgorithm //
    uLib::Interface::IsA<UpdateAlgT, IBInterface::IBVoxCollectionUpdateAlgorithm>();

    // Analyzer Update //
    algorithm(this, threshold);

    // MAP update //
    if(m_MAPAlgorithm) m_MAPAlgorithm->UpdateDensity(this, threshold);

    // Reinitialize voxels //
    this->InitCount(0);
}

template < class StaticUpdateAlgT >
void IBVoxCollection::UpdateDensity(unsigned int threshold) {

    // check structural for UpdateAlgorithm //
    // NOT WORKING !!!!!
    // uLib::Interface::IsA<StaticUpdateAlgT, IBInterface::IBVoxCollectionStaticUpdateAlgorithm>();

    // Analyzer Update //
    StaticUpdateAlgT::UpdateDensity(this,threshold);

    // MAP update //
    if(m_MAPAlgorithm) m_MAPAlgorithm->UpdateDensity(this, threshold);

    // Reinitialize voxels //
    this->InitCount(0);
}



#endif	/* IBVOXCOLLECTION_H */

