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

IBVoxCollection::IBVoxCollection(const IBVoxCollection::BaseClass &copy) :
    BaseClass(copy),
    m_MAPAlgorithm(NULL)
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







