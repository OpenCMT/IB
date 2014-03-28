
#include "root/TFile.h"
#include "root/TTree.h"

#include <math.h>
#include <Math/Utils.h>

#include "IBMuonCollection.h"

class IBMuonCollectionPimpl {

    friend class IBMuonCollection;
    struct _Cmp {
        bool operator()(MuonScatter &data, const float value)
        {
            return MuonScatterAngle(data) <= value;
        }
    };
    struct _PCmp {
        bool operator()(const MuonScatter &data, const float value)
        {
            return data.GetMomentumPrime() <= value;
        }
    };
    struct _PPCmp {
        bool operator()(const MuonScatter &data, const float value)
        {
            return data.GetMomentum() <= value;
        }
    };

    static float MuonScatterAngle(const MuonScatter &mu) {
        Vector3f in = mu.LineIn().direction.head(3);
        Vector3f out = mu.LineOut().direction.head(3);
        float a = in.transpose() * out;
        a = fabs( acos(a / (in.norm() * out.norm())) );
        if(uLib::isFinite(a)) return a;
        else return 0;
    }

public:

    IBMuonCollectionPimpl() : m_HiPass(1), m_SliceIndex(0) { }


    // members //
    bool m_HiPass;
    unsigned int m_SliceIndex;
    Vector<MuonScatter> m_Data;
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
    // FINIRE o PENSARE perche non si puo' fare add muon dopo set Hi/LowPass //
}

Vector<MuonScatter> &IBMuonCollection::Data()
{
    return d->m_Data;
}


const MuonScatter &IBMuonCollection::At(int i) const
{
    return d->m_Data.at(d->m_SliceIndex * d->m_HiPass + i);
}

MuonScatter &IBMuonCollection::operator [](int i)
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

void IBMuonCollection::SetHiPassMomentum(float momenutm)
{
    d->m_SliceIndex =
    VectorSplice(d->m_Data.begin(),d->m_Data.end(),momenutm,
                 IBMuonCollectionPimpl::_PCmp());
    d->m_HiPass = 1;
}

void IBMuonCollection::SetLowPassMomentum(float momentum)
{
    d->m_SliceIndex =
    VectorSplice(d->m_Data.begin(),d->m_Data.end(),momentum,
                 IBMuonCollectionPimpl::_PCmp());
    d->m_HiPass = 0;
}

void IBMuonCollection::SetHiPassMomentumPrime(float momenutm)
{
    d->m_SliceIndex =
    VectorSplice(d->m_Data.begin(),d->m_Data.end(),momenutm,
                 IBMuonCollectionPimpl::_PPCmp());
    d->m_HiPass = 1;
}

void IBMuonCollection::SetLowPassMomentumPrime(float momentum)
{
    d->m_SliceIndex =
    VectorSplice(d->m_Data.begin(),d->m_Data.end(),momentum,
                 IBMuonCollectionPimpl::_PPCmp());
    d->m_HiPass = 0;
}



void IBMuonCollection::PrintSelf(std::ostream &o)
{
    o << "---- Muon Collection: ----- \n";
    o << " Data size: " << d->m_Data.size() << "\n";
    if(this->size()!=0)
        o << " Muons passed: " << this->size() << "\n";
}

void IBMuonCollection::DumpTTree(const char *filename)
{
    static int counter = 0;
    static TFile *file = new TFile(filename,"RECREATE");

    if(filename) {
        char name[100];
        sprintf(name,"p_%i",counter++);
        gDirectory->cd(file->GetPath());
        TTree *tree = new TTree(name,name);
        float p = 0;
        float angle = 0;

        tree->Branch("p",&p,"float");
        tree->Branch("angle",&angle,"float");

        for(int i=0; i < this->size(); ++i )
        {
            const MuonScatter &mu = this->At(i);
            p = mu.GetMomentum();
            angle = d->MuonScatterAngle(mu);
            tree->Fill();
        }
        tree->Write();
        delete tree;
    }
    else {
        file->Write();
        file->Close();
        delete file;
    }
        return;
}

