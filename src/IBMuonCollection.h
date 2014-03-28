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

    Vector<MuonScatter> &Data();

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

private:
    class IBMuonCollectionPimpl *d;


};


#endif // IBMUONCOLLECTION_H
