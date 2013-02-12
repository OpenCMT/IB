
#include <math.h>
#include <Math/Utils.h>

#include "IBMuonCollection.h"

class IBMuonCollectionPimpl {

    friend class IBMuonCollection;
    struct _Cmp {
        bool operator()(const MuonScatterData &data, const float value)
        {
            return MuonScatterAngle(data) <= value;
        }
    };

    static float MuonScatterAngle(const MuonScatterData &mu) {
        Vector3f in = mu.LineIn().direction.head(3);
        Vector3f out = mu.LineOut().direction.head(3);
        float a = in.transpose() * out;
        a = fabs( M_PI_2 - asin(a / (in.norm() * out.norm())) );
        if(uLib::isFinite(a)) return a;
        else return 0;
    }

public:

    IBMuonCollectionPimpl() : m_HiPass(1), m_SliceIndex(0) { }


    // members //
    bool m_HiPass;
    unsigned int m_SliceIndex;
    Vector<MuonScatterData> m_Data;
};



IBMuonCollection::IBMuonCollection() :
    d(new IBMuonCollectionPimpl)
{}

IBMuonCollection::~IBMuonCollection()
{
    delete d;
}

void IBMuonCollection::AddMuon(MuonScatter &mu)
{
    d->m_Data.push_back(mu);
    // FINIRE //
}

const MuonScatterData &IBMuonCollection::At(int i) const
{
    return d->m_Data.at(d->m_SliceIndex * d->m_HiPass + i);
}

MuonScatterData &IBMuonCollection::operator [](int i)
{    
    return d->m_Data[d->m_SliceIndex * d->m_HiPass + i];
}

size_t IBMuonCollection::size() const
{
    if(d->m_HiPass) return d->m_Data.size() - d->m_SliceIndex;
    else return d->m_SliceIndex;
}


void IBMuonCollection::SetHiPassAngle(float angle)
{
    d->m_SliceIndex =
    VectorSplice(d->m_Data.begin(),d->m_Data.end(),angle,
                 IBMuonCollectionPimpl::_Cmp());
    d->m_HiPass = 1;
}

void IBMuonCollection::SetLowPassAngle(float angle)
{
    d->m_SliceIndex =
    VectorSplice(d->m_Data.begin(),d->m_Data.end(),angle,
                 IBMuonCollectionPimpl::_Cmp());
    d->m_HiPass = 0;
}



void IBMuonCollection::PrintSelf(std::ostream &o)
{
    o << "---- Muon Collection: ----- \n";
    o << " Data size: " << d->m_Data.size() << "\n";

    //    o << " Muons over cut : " << this->size() << "\n";
}


