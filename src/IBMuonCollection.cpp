/*////////////////////////////////////////////////////////////////////////////
 Copyright 2018 Istituto Nazionale di Fisica Nucleare

 Licensed under the EUPL, Version 1.2 or - as soon they will be approved by
 the European Commission - subsequent versions of the EUPL (the "Licence").
 You may not use this work except in compliance with the Licence.

 You may obtain a copy of the Licence at:

 https://joinup.ec.europa.eu/software/page/eupl

 Unless required by applicable law or agreed to in writing, software
 distributed under the Licence is distributed on an "AS IS" basis, WITHOUT
 WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. See the
 Licence for the specific language governing permissions and limitations under
 the Licence.
////////////////////////////////////////////////////////////////////////////*/
#include <math.h>
#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include <TROOT.h>

#include <Math/Dense.h>
#include <Math/Utils.h>

#include "IBMuonCollection.h"


IBMuonCollection::IBMuonCollection() :
    m_HiPass(1),
    m_SliceIndex(0)
{}

IBMuonCollection::~IBMuonCollection()
{}

void IBMuonCollection::AddMuon(MuonScatter &mu)
{
    m_Data.push_back(mu);
    // FINIRE o PENSARE perche non si puo' fare add muon dopo set Hi/LowPass //
}

void IBMuonCollection::AddMuonFullPath(std::vector<HPoint3f> fullPath)
{
    m_FullPathData.push_back(fullPath);
}

std::vector<MuonScatter> &IBMuonCollection::Data()
{
    return m_Data;
}

std::vector<std::vector<HPoint3f> > &IBMuonCollection::FullPath()
{
    return m_FullPathData;
}


const MuonScatter &IBMuonCollection::At(int i) const
{
    return m_Data.at(m_SliceIndex * m_HiPass + i);
}

MuonScatter &IBMuonCollection::operator [](int i)
{
    return m_Data[m_SliceIndex * m_HiPass + i];
}

size_t IBMuonCollection::size() const
{
    if(m_HiPass) return m_Data.size() - m_SliceIndex;
    else return m_SliceIndex;
}



void IBMuonCollection::splice(const float value,
    std::function<bool(MuonScatter&, const float)> cmp)
{
    std::vector<MuonScatter>::iterator m_it = m_Data.begin();
    std::vector<MuonScatter>::iterator m_end = m_Data.end() - 1;

    while(m_it != m_end)
    {
        if (cmp(*m_it, value))       ++m_it;
        else if(cmp(*m_end, value))  std::swap(*m_it, *m_end--);
        else                         --m_end;
    }

    m_SliceIndex = m_it - m_Data.begin();
}



const float calc_angle(const MuonScatter &data)
{
    Vector3f in = data.LineIn().direction.head(3);
    Vector3f out = data.LineOut().direction.head(3);
    float a = in.transpose() * out;
    a = fabs( acos(a / (in.norm() * out.norm())) );
    if(uLib::isFinite(a)) return a;
    return 0.;
}

void IBMuonCollection::SetHiPassAngle(float angle)
{
    splice( angle,
            [](MuonScatter &data, const float value)->bool
                {return calc_angle(data) <= value;}
            );
    m_HiPass = 1;
}

void IBMuonCollection::SetLowPassAngle(float angle)
{
    splice( angle,
            [](MuonScatter &data, const float value)->bool
                {return calc_angle(data) <= value;}
            );
    m_HiPass = 0;
}

void IBMuonCollection::SetHiPassMomentum(float momentum)
{
    splice( momentum,
            [](MuonScatter &data, const float value)->bool
                {return data.GetMomentum() <= value;}
            );
    m_HiPass = 1;
}

void IBMuonCollection::SetLowPassMomentum(float momentum)
{
    splice( momentum,
            [](MuonScatter &data, const float value)->bool
                {return data.GetMomentum() <= value;}
            );
    m_HiPass = 0;
}

void IBMuonCollection::SetHiPassMomentumPrime(float momentum)
{
    splice( momentum,
            [](MuonScatter &data, const float value)->bool
                {return data.GetMomentumPrime() <= value;}
            );
    m_HiPass = 1;
}

void IBMuonCollection::SetLowPassMomentumPrime(float momentum)
{
    splice( momentum,
            [](MuonScatter &data, const float value)->bool
                {return data.GetMomentumPrime() <= value;}
            );
    m_HiPass = 0;
}



void IBMuonCollection::PrintSelf(std::ostream &o)
{
    o << "---- Muon Collection: ----- \n";
    o << " Data size: " << m_Data.size() << "\n";
    if(this->size()!=0)
        o << " Muons passed: " << this->size() << "\n";

    o << "\n First muon: \n";
    o << this->At(0);
    o << "\n Last muon: \n";
    o << this->At(this->size() - 1);

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
    float inP, outP;

    tree->Branch("inOrigin",&inOrigin,"x:y:z");
    tree->Branch("inDir",&inDir,"x:y:z");
    tree->Branch("inP",&inP,"inP/F");
    tree->Branch("outOrigin",&outOrigin,"x:y:z");
    tree->Branch("outDir",&outDir,"x:y:z");
    tree->Branch("outP",&outP,"outP/F");

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

        inP=u_mu.GetMomentum();

        outOrigin.x=u_mu.LineOut().origin[0];
        outOrigin.y=u_mu.LineOut().origin[1];
        outOrigin.z=u_mu.LineOut().origin[2];

        outDir.x=u_mu.LineOut().direction[0];
        outDir.y=u_mu.LineOut().direction[1];
        outDir.z=u_mu.LineOut().direction[2];

        outP=u_mu.GetMomentumPrime();

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
        HPoint3f outDir = mu.LineOut().direction;

        if(!std::isnan(outPos[0]))
            file << outPos[0] << " " << outPos[1] << " " << outPos[2] << " " << outDir[0]/outDir[1] << " " << outDir[2]/outDir[1] << "\n";
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
void IBMuonCollection::PerformMuonSelfAlignment(){
    std::cout << "\n--------- DATA ALIGNMENT.... (muon self adjust)" << std::endl;
    std::pair<HVector3f,HVector3f> align = GetAlignment();
    align.first << 0,0,0,0;
    SetAlignment(align);
    align = GetAlignment();
    align.first << 0,0,0,0;
    SetAlignment(align);
    align = GetAlignment();
    align.first << 0,0,0,0;
    SetAlignment(align);
    align = GetAlignment();
    align.first << 0,0,0,0;
    SetAlignment(align);
    align = GetAlignment();
    align.first << 0,0,0,0;
    SetAlignment(align);
    align = GetAlignment();
    SetAlignment(align);
    align = GetAlignment();
    SetAlignment(align);
    align = GetAlignment();
    SetAlignment(align);
    align = GetAlignment();
    SetAlignment(align);
    align = GetAlignment();
    SetAlignment(align);
    align = GetAlignment();
    SetAlignment(align);
    align = GetAlignment();
    SetAlignment(align);

    return;
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

// //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void IBMuonCollection::AddCollection(IBMuonCollection &muonsColl){

    for(int i=0; i < muonsColl.size(); ++i ){
                MuonScatter mu = muonsColl.At(i);
                this->AddMuon(mu);
    }

    return;
}

