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

#include <math.h>
#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TROOT.h>

#include <Math/Dense.h>
#include <Math/Utils.h>
#include "Root/RootMuonScatter.h"

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
  Vector<Vector<HPoint3f> > m_FullPathData;
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

void IBMuonCollection::AddMuonFullPath(Vector<HPoint3f> fullPath)
{
    d->m_FullPathData.push_back(fullPath);
}

Vector<MuonScatter> &IBMuonCollection::Data()
{
    return d->m_Data;
}

Vector<Vector<HPoint3f> > &IBMuonCollection::FullPath()
{
    return d->m_FullPathData;
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

    for(int i=0; i < 3; ++i )
    {
        const MuonScatter &u_mu = this->At(i);
        std::cout << u_mu;
    }

}

// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////7
void IBMuonCollection::DumpTTree(const char *filename)
{ 
    std::cout << "\n\n------------- Dump muon collection on root file " << filename << std::endl;
    static TFile *file = new TFile(filename,"RECREATE");

    char name[100];
    sprintf(name,"muons");
    gDirectory->cd(file->GetPath());
    TTree *tree = new TTree(name,name);

    uLib::MuonScatter u_mu;
    ROOT::Mutom::MuonScatter mu;

    tree->Branch("mu",&mu);

    for(int i=0; i < this->size(); ++i )
    {
        const MuonScatter &u_mu = this->At(i);
        mu << u_mu;
        tree->Fill();
    }

    tree->Write();

//    file->Write();
    file->Close();
//    delete tree;
//    delete file;

    return;
}

// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////7
void IBMuonCollection::DumpSimpleTree(const char *filename)
{

    std::cout << "\n\n------------- Dump muon collection on root file " << filename << std::endl;

    static TFile *file = new TFile(filename,"RECREATE");

    /// open file, tree
    gDirectory->cd(file->GetPath());
    gROOT->ProcessLine("#include<vector>");

    // Define some simple structures
    typedef struct {Float_t x,y,z;} POINT;

    char name[100];
    sprintf(name,"muons");
    TTree *tree = (TTree*)file->Get("muons");
    if(!tree)
        tree = new TTree(name,name);

    POINT inOrigin, outOrigin;
    POINT inDir, outDir;

    tree->Branch("inOrigin",&inOrigin,"x:y:z");
    tree->Branch("inDir",&inDir,"x:y:z");
    tree->Branch("outOrigin",&outOrigin,"x:y:z");
    tree->Branch("outDir",&outDir,"x:y:z");

    /// event loop
    std::cout << "Reading " << this->size() << " muons " << std::endl;

    for(int i=0; i < this->size(); ++i )
    {
        const MuonScatter &u_mu = this->At(i);

        inOrigin.x=u_mu.LineIn().origin[0];
        inOrigin.y=u_mu.LineIn().origin[1];
        inOrigin.z=u_mu.LineIn().origin[2];

        inDir.x=u_mu.LineIn().direction[0];
        inDir.y=u_mu.LineIn().direction[1];
        inDir.z=u_mu.LineIn().direction[2];

        outOrigin.x=u_mu.LineOut().origin[0];
        outOrigin.y=u_mu.LineOut().origin[1];
        outOrigin.z=u_mu.LineOut().origin[2];

        outDir.x=u_mu.LineOut().direction[0];
        outDir.y=u_mu.LineOut().direction[1];
        outDir.z=u_mu.LineOut().direction[2];

        tree->Fill();
    }

    tree->Write();
    file->Close();

    return;
}

// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////7
void IBMuonCollection::DumpTxt(const char *filename)
{
    std::cout << "\nDUMP on txt file for Calvini 3D reconstruction\n";
    std::fstream file;
    file.open(filename, std::ios::out);

    uLib::MuonScatter mu;
    // event loop
    for(int i=0; i < this->size(); ++i )
    {
        const MuonScatter &mu = this->At(i);

        HPoint3f inPos = mu.LineIn().origin;
        HVector3f inDir = mu.LineIn().direction;
        file << inPos[0] << " " << inPos[1] << " " << inPos[2] << " " << inDir[0]/inDir[1] << " " << inDir[2]/inDir[1] << " ";

        HPoint3f outPos = mu.LineOut().origin;
        if(!std::isnan(outPos[0]))
            file << outPos[0] << " " << outPos[1] << " " << outPos[2] << "\n";
        else
            file << "\n";
    }

    file.close();
    return;
}


// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
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
//    std::cout << " ---- muons self adjust (use with care!) ----- \n";
//    std::cout << "Rotation : " << rot.transpose() << "\n";
//    std::cout << "Position : " << position.transpose() << "\n\n";

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

// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void IBMuonCollection::dataRotoTranslation(Eigen::Matrix4f  t)
{
    std::cout << "\nMuon collection: apply roto-translation matrix: \n" << t <<  "\n\n";

    for(int i=0; i<this->size(); ++i) {
        MuonScatter &mu = this->operator [](i);

        // IN muon roto-traslation
        mu.LineIn().origin = t * mu.LineIn().origin;
        mu.LineIn().direction = t * mu.LineIn().direction;
        mu.LineIn().direction /= fabs(mu.LineIn().direction(1)); // back to slopes

        // OUT muon roto-traslation
        mu.LineOut().origin = t * mu.LineOut().origin;
        mu.LineOut().direction = t * mu.LineOut().direction;
        mu.LineOut().direction /= fabs(mu.LineOut().direction(1)); // back to slopes

//        // if resulting muon isn't falling from sky.... reverse!
//        if(mu.LineIn().direction[1] > 0){
//            std::swap(mu.LineIn(),mu.LineOut());
//            mu.LineIn().direction *= -1;
//            mu.LineOut().direction *= -1;
//        }
    }

    return;
}

// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void IBMuonCollection::dataRotoTranslation(Vector3f rot, Vector3f trans)
{
    //If a standard right-handed Cartesian coordinate system is used, with the x-axis to the right and the y-axis up, the rotation R(θ) is counterclockwise
    //If a left-handed Cartesian coordinate system is used (like THIS!) , with x directed to the right but y directed down, R(θ) is clockwise

    rot = rot / 180. * M_PI;

    Eigen::Matrix4f  t = this->createAffineMatrix(rot[0], rot[1], rot[2], trans);

    std::cout << "\nMuon collection: apply roto-translation matrix: \n" << t <<  "\n\n";

    for(int i=0; i<this->size(); ++i) {
        MuonScatter &mu = this->operator [](i);

        // IN muon roto-traslation
        mu.LineIn().origin = t * mu.LineIn().origin;
        mu.LineIn().direction = t * mu.LineIn().direction;
        mu.LineIn().direction /= fabs(mu.LineIn().direction(1)); // back to slopes

        // OUT muon roto-traslation
        mu.LineOut().origin = t * mu.LineOut().origin;
        mu.LineOut().direction = t * mu.LineOut().direction;
        mu.LineOut().direction /= fabs(mu.LineOut().direction(1)); // back to slopes
    }
}

// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Eigen::Matrix4f  IBMuonCollection::createAffineMatrix(float a, float b, float c, Vector3f trans)
{
    Eigen::Transform<float, 3, Eigen::Affine> t;
    t = Eigen::AngleAxis<float>(c, Vector3f::UnitZ());
    t.prerotate(Eigen::AngleAxis<float>(b, Vector3f::UnitY()));
    t.prerotate(Eigen::AngleAxis<float>(a, Vector3f::UnitX()));
    t.pretranslate(trans);
    return t.matrix();
}
