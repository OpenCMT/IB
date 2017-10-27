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



#ifndef IBMUONCOLLECTION_H
#define IBMUONCOLLECTION_H

#include <Core/Object.h>
#include <Core/Vector.h>

#include <Math/Dense.h>

#include "Root/RootMuonScatter.h"
#include "Detectors/MuonScatter.h"

using namespace uLib;

class IBMuonCollection : public Object {

public:

    IBMuonCollection();
    ~IBMuonCollection();

    void AddMuon(MuonScatter &mu);
    void AddMuonFullPath(Vector<HPoint3f> fullPath);
      
    Vector<MuonScatter> &Data();
    Vector<Vector<HPoint3f> > &FullPath();

    const MuonScatter &At(int i) const;

    MuonScatter &operator[](int i);

    size_t size() const;

    void SetHiPassAngle(float angle);
    void SetLowPassAngle(float angle);

    void SetHiPassMomentum(float momenutm);
    void SetLowPassMomentum(float momentum);

    void SetHiPassMomentumPrime(float momenutm);
    void SetLowPassMomentumPrime(float momentum);


    void PrintSelf(std::ostream &o);
    void DumpTTree(const char *filename);
    void DumpTxt(const char *filename);
    std::pair<HVector3f, HVector3f> GetAlignment();
    void SetAlignment(std::pair<HVector3f, HVector3f> align);



private:
    class IBMuonCollectionPimpl *d;


};


#endif // IBMUONCOLLECTION_H
