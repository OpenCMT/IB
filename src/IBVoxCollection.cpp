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
#include "Math/Dense.h"
#include "Math/VectorSpace.h"


#include "IBVoxCollection.h"

#include <TGeoNode.h>
#include <TGeoManager.h>
#include <TRandom.h>
#include <TFile.h>

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
{
    std::cout << "\n  ------ VOXEL COLLECTION INIT ------" << std::endl;
    std::cout << "Please cross-check that collection_dimension % voxel_dimension=0, otherwise leaks could be generated... " << std::endl;
 }

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
                                              Vector4f center,
                                              Vector4f size)
{
    Vector3i start = this->Find(center-size);
    Vector3i end = this->Find(center+size);
    return CountLambdaOverThreshold(threshold,start,end);
}

// --- Import MC truth --------------------------------------------- //
/// SV 20150525 costruzione di voxel collection con densita' montecarlo da file geant
IBVoxCollection IBVoxCollection::getMCImage(const char* file, int nsamples)
{
    // open TGeom instances
    std::cout << "Requiring MC image from file " << file << std::endl;
    TGeoManager::Import(file);
    gGeoManager->SetVerboseLevel(0);
    /// to test geometry
    //gGeoManager->CheckOverlaps(0.01);
    //gGeoManager->PrintOverlaps();
    TGeoNode *node;
    TGeoMedium *medium;
    TGeoMaterial *material;

    IBVoxCollection out = *this;
    Vector4f B = this->GetPosition().homogeneous();

    //loop over voxel collection
    Vector3i dim = this->GetDims();
    Vector3i id;
    for(int ix=0 ; ix < dim(0) ; ++ix) {
        for(int iy=0 ; iy < dim(1) ; ++iy) {
            for(int iz=0 ; iz < dim(2) ; ++iz) {

//                std::cout << "\n \n  *** NEW VOXEL " << std::endl;

                // voxel position
                Vector3i iv(ix,iy,iz);
                Vector3f v = Vector3f(iv.cast<float>().cwiseProduct(this->GetSpacing()));
                Vector4f v4 = Vector4f(v(0),v(1),v(2),1);
                Vector4f Bvox = B + v4;
                Vector4f Evox = Bvox + this->GetSpacing().homogeneous();

                // Calculating average values for different variables
                Double_t RadLen  = 0.;     // Radiation length
                Double_t invRadLen  = 0.;  // Inverse Radiation length
                Double_t IntLen  = 0.; // Interaction length
                Double_t density = 0.; // Density
                for(int ir=0; ir<nsamples; ir++) {
                    Double_t x,y,z;
                    x = gRandom->Uniform(Bvox[0],Evox[0]);
                    y = gRandom->Uniform(Bvox[1],Evox[1]);
                    z = gRandom->Uniform(Bvox[2],Evox[2]);
                    node = gGeoManager->FindNode(x,z,y);
                    medium = node->GetMedium();
                    material = medium->GetMaterial();
                    RadLen      += material->GetRadLen();
                    invRadLen   += 1/material->GetRadLen();
                    IntLen      += material->GetIntLen();
                    density     += material->GetDensity();
//                    if(1/material->GetRadLen() > 0.1)
//                    {
//                        std::cout << "Material "  << material->GetName() << std::endl;
//                        std::cout << "Radiation length [cm] "    << material->GetRadLen() << std::endl;
//                        std::cout << "Inverse radiation length [1/cm] "    << 1/material->GetRadLen() << std::endl;
//                        std::cout << "Inverse radiation length SUM "    << invRadLen << std::endl;

//                    }
                }
                RadLen     = RadLen/(double)nsamples;
                invRadLen  = invRadLen/(double)nsamples;
                IntLen     = IntLen/(double)nsamples;
                density    = density/(double)nsamples;

//                if(invRadLen > 0.024)
//                {
//                    std::cout << "  -----------  Voxel B " << Bvox[0] << ", " << Bvox[1] << ", " << Bvox[2] << ", " << std::endl;
//                    std::cout << "  -----------  Voxel E " << Evox[0] << ", " << Evox[1] << ", " << Evox[2] << ", " << std::endl;
//                    std::cout << "  -----------  Averages over " << nsamples << " samples" << std::endl;
//                    std::cout << "Radiation length [cm] "    << RadLen << std::endl;
//                    std::cout << "Inverse radiation length [1/cm] "    << invRadLen << std::endl;
//                    std::cout << "Interaction length [cm] "  << IntLen << std::endl;
//                    std::cout << "Density [g/cm3] "             << density << std::endl;
//                    // #########################################  END
//                }

                //out.SetValue(iv,invRadLen); // cm
                out.SetValue(iv,invRadLen*100); // m

            }
        }
    }

    return out;
}

/// SV 20150520 rielaborazione della funzione calculate_voxel_means.C di Gemano
float IBVoxCollection::getVoxelMCDensity(const char* file, Vector4f c1, Vector4f c2, int nrandom)
{
    if(nrandom<=0) {
        return -1;
        std::cout << "Please give a positive number for nrandom" << std::endl;
    }
    TGeoNode *node;
    TGeoMedium *medium;
    TGeoMaterial *material;
    TGeoManager::Import(file);
    gGeoManager->SetVerboseLevel(0);
    Double_t x,y,z;

    /// to test geometry
    //gGeoManager->CheckOverlaps(0.01);
    //gGeoManager->PrintOverlaps();

    // Printing information for the two corners
    // ######################################### INIT
    x = c1[0]; y = c1[1]; z = c1[2];
    node = gGeoManager->FindNode(x,z,y);
    medium = node->GetMedium();
    material = medium->GetMaterial();
//    std::cout << "First voxel corner: " << x << "," << y << "," << z << "," << std::endl;
//    std::cout << "Material " << material->GetName() << std::endl;
//    std::cout << "Density " << material->GetDensity() << " Radiation length " << material->GetRadLen() << " Interaction length " << material->GetIntLen()<< std::endl;
//    std::cout << "IsMixture " << material->IsMixture() << std::endl;
//    std::cout << "Z " << material->GetZ() << std::endl;

    x = c2[0]; y = c2[1]; z = c2[2];
    node = gGeoManager->FindNode(x,z,y);
    medium = node->GetMedium();
    material = medium->GetMaterial();
//    std::cout << "Second voxel corner: " << x << "," << y << "," << z << "," << std::endl;
//    std::cout << "Material " << material->GetName() << std::endl;
//    std::cout << "Density " << material->GetDensity() << " Radiation length " << material->GetRadLen() << " Interaction length " << material->GetIntLen()<< std::endl;
//    std::cout << "IsMixture " << material->IsMixture() << std::endl;
//    std::cout << "Z " << material->GetZ() << std::endl;
    // #########################################  END

    // Calculating average values for different variables
    // ######################################### INIT
    Double_t RadLen  = 0.;  // Radiation length
    Double_t IntLen  = 0.; // Interaction length
    Double_t density = 0.; // Density
    Double_t Z=0.;         // Z
    for(int ir=0; ir<nrandom; ir++) {
        x = gRandom->Uniform(c1[0],c2[0]);
        y = gRandom->Uniform(c1[1],c2[1]);
        z = gRandom->Uniform(c1[2],c2[2]);
        node = gGeoManager->FindNode(x,z,y);
        medium = node->GetMedium();
        material = medium->GetMaterial();
        RadLen   += material->GetRadLen();
        IntLen   += material->GetIntLen();
        density  += material->GetDensity();
    }
    RadLen     = RadLen/(double)nrandom;
    IntLen     = IntLen/(double)nrandom;
    density    = density/(double)nrandom;

//    if(density>1){
//        std::cout << "  -----------  Averages " << std::endl;
//        std::cout << "Radiation length [cm] "    << RadLen << std::endl;
//        std::cout << "Interaction length [cm] "  << IntLen << std::endl;
//        std::cout << "Density [g/cm3] "             << density << std::endl;
//        // #########################################  END
//    }
    float invrad = 1/RadLen;
    return invrad;
}

