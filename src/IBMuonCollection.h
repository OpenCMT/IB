#ifndef IBMUONCOLLECTION_H
#define IBMUONCOLLECTION_H

#include <Core/Object.h>
#include <Core/Vector.h>

#include <Math/Dense.h>

#include <Detectors/MuonScatter.h>


using namespace uLib;

class IBMuonCollection : public Object {

public:

    IBMuonCollection();
    ~IBMuonCollection();

    void AddMuon(MuonScatter &mu);

    Vector<MuonScatterData> &Data();

    const MuonScatterData &At(int i) const;

    MuonScatterData &operator[](int i);

    size_t size() const;

    void SetHiPassAngle(float angle);
    void SetLowPassAngle(float angle);

    void SetHiPassMomentumPrime(float momenutm);
    void SetLowPassMomentumPrime(float momentum);


    void PrintSelf(std::ostream &o);


private:
    class IBMuonCollectionPimpl *d;


};


#endif // IBMUONCOLLECTION_H
