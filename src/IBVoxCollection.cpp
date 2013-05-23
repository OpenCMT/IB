
#include "Core/StaticInterface.h"
#include "IBVoxCollection.h"


using namespace uLib;









////////////////////////////////////////////////////////////////////////////////
///// IBVOXCOLLECTION  /////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////


// --- CTR ------------------------------------------------------------------ //

IBVoxCollection::IBVoxCollection() :
    BaseClass(),
    m_MAPAlgorithm(NULL)
{}

IBVoxCollection::IBVoxCollection(const Vector3i size) :
    BaseClass(size),
    m_MAPAlgorithm(NULL)
{}

IBVoxCollection::IBVoxCollection(const IBVoxCollection &copy) :
    BaseClass(copy),
    m_MAPAlgorithm(copy.m_MAPAlgorithm)
{}

void IBVoxCollection::SetMAPAlgorithm(IBAbstract::IBVoxCollectionMAPAlgorithm *algorithm)
{
    this->m_MAPAlgorithm = algorithm;
}










// --- Threshold identification --------------------------------------------- //


int IBVoxCollection::CountLambdaOverThreshold(float threshold)
{
    return this->CountLambdaOverThreshold(threshold,
                                          Vector3i(0,0,0),
                                          this->GetDims());
}

int IBVoxCollection::CountLambdaOverThreshold(float threshold,
                                              Vector3i boxp1,
                                              Vector3i boxp2)
{
    int count = 0;
    Vector3i id;
    for(int x=boxp1(0) ; x < boxp2(0) ; ++x) {
        for(int y=boxp1(1) ; y < boxp2(1) ; ++y) {
            for(int z=boxp1(2) ; z < boxp2(2) ; ++z) {
                id << x,y,z;
                if(this->At(id).Value >= threshold) count++;
            }
        }
    }
    return count;
}

int IBVoxCollection::CountLambdaOverThreshold(float threshold,
                                              HPoint3f center,
                                              HVector3f size)
{
    Vector3i start = this->Find(center-size);
    Vector3i end = this->Find(center+size);
    return CountLambdaOverThreshold(threshold,start,end);
}






