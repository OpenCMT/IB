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




#include "root/TFile.h"
#include "root/TTree.h"

#include <math.h>
#include <Math/Dense.h>
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

/**
 * @brief IBMuonCollection::GetAlignment
 *  compute mean of rotation and shift of muon out relative to in
 * @return Vector4f ( phi x theta z )
 */
std::pair<HVector3f,HVector3f> IBMuonCollection::GetAlignment()
{
    Vector3f direction(0,0,0);
    Vector3f position(0,0,0);
    for(int i=0; i<this->size(); ++i) {
        const MuonScatter &mu = this->At(i);
        Vector3f dir_in  = mu.LineIn().direction.head(3).normalized();
        Vector3f dir_out = mu.LineOut().direction.head(3).normalized();
        direction += dir_in.cross(dir_out);

        float y = mu.LineOut().origin(1) - mu.LineIn().origin(1);
        Vector3f shift = mu.LineOut().origin.head(3) - dir_in*(y/dir_in(1));
        position += mu.LineIn().origin.head(3) - shift;
    }
    direction /= this->size();
    position  /= this->size();
    HVector3f  pos;
    pos.head(3) = position;
    HVector3f rot;
    rot.head(3) = direction.normalized();
    rot(3) = direction.norm();
    std::cout << " ---- muons self adjust (use with care!) ----- \n";
    std::cout << "Rotation : " << rot.transpose() << "\n";
    std::cout << "Position : " << position.transpose() << "\n\n";

    return std::pair<HVector3f,HVector3f>(pos,rot);
}

void IBMuonCollection::SetAlignment(std::pair<HVector3f,HVector3f> align)
{
    const HVector3f &pos = align.first;
    const HVector3f &rot = align.second;

    Eigen::Quaternion<float> q(Eigen::AngleAxis<float>(-asin(rot(3)), rot.head(3)));
//    Eigen::Affine3f tr(Eigen::AngleAxis<float>(-rot(3), rot.head(3)));
//    std::cout << "QUATERNION: \n" << q.matrix() << "\n\nROTATION: \n" << tr.matrix() << "\n\n";

    for(int i=0; i<this->size(); ++i) {
        MuonScatter &mu = this->operator [](i);
//        if(i<10){
//            Vector3f v1 = q * mu.LineOut().direction.head<3>();
//            Vector4f v2 = tr * mu.LineOut().direction;
//            std::cout << "dif: " << v1.transpose() << " vs " << v2.transpose() << "\n";
//        }
//        mu.LineOut().direction = tr * mu.LineOut().direction;
        mu.LineOut().direction.head<3>() = q * mu.LineOut().direction.head<3>();
        mu.LineOut().direction /= fabs(mu.LineOut().direction(1)); // back to slopes
        mu.LineOut().origin += pos;
    }
}

