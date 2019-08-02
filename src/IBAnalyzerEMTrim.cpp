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



#include <stdio.h>
#include <iterator>

#include "Math/Utils.h"

#include "IBAnalyzerEMTrim.h"
#include "IBAnalyzerEMAlgorithm.h"

/*
This class uses uLib Accumulator Trim to push SijCap into a trimmed set
This should be similar to median Sij effect

It is actually not very efficient due to lack of parallelization ...
... a future expansion will use boost threads or a thread safe circular buffer.
*/




using namespace uLib;


namespace IBAnalyzerEMTrimDetail {

// --- Update --------------------------------------------------------------- //


template < class StaticUpdateAlgT >
inline void IBVoxCollectionATrim::UpdateDensity(unsigned int threshold) {
    StaticUpdateAlgT::UpdateDensity(this,threshold);

    for(unsigned int i=0; i<this->Data().size(); ++i) {
        this->Data().at(i).Count = 0;
    }
}


inline IBVoxCollection& operator << (IBVoxCollection &voxels, IBVoxCollectionATrim &mdn_voxels) {
    assert(voxels.Data().size() == mdn_voxels.Data().size());
    for( int i=0; i<voxels.Data().size(); ++i )
    {
        voxels[i].Value = mdn_voxels.At(i).Value;
        voxels[i].Count = mdn_voxels.At(i).Count;
        voxels[i].SijCap = 0;
    }
    return voxels;
}

inline IBVoxCollectionATrim& operator << (IBVoxCollectionATrim &mdn_voxels, IBVoxCollection &voxels) {
    assert(voxels.Data().size() == mdn_voxels.Data().size());
    for( int i=0; i<voxels.Data().size(); ++i )
    {
        mdn_voxels[i].Value = voxels.At(i).Value;
        mdn_voxels[i].Count = voxels.At(i).Count;
    }
    return mdn_voxels;
}

inline void IBVoxCollectionATrim::SetABTrim(int a, int b) {
    for(int i=0 ; i < Data().size(); ++i)
        Data().at(i).SetABTrim(a,b);
}




////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// Update density //



class UpdateDensitySijCapATrimAlgorithm :
        public IBInterface::IBVoxCollectionStaticUpdateAlgorithm
{
public:
    static void UpdateDensity(IBAnalyzerEMTrimDetail::IBVoxCollectionATrim *voxels,
                              unsigned int threshold)
    {
        for(unsigned int i=0; i< voxels->Data().size(); ++i) {
            IBVoxelABTrim& voxel = voxels->Data()[i];
            unsigned int tcount = voxel.Count;
            if (tcount > 0 && (threshold == 0 || tcount >= threshold)) {
                voxel.Value += voxel.SijCap();
                if(unlikely(!isFinite(voxel.Value) || voxel.Value > 100.E-6)) {
                    voxel.Value = 100.E-6;
                } else if (unlikely(voxel.Value < 0.)) voxel.Value = 0.1E-6;
            }
            else
                voxel.Value = 0;
            voxel.Count = 0;
        }
    }
};


} // IBAnalyzerEMTrimDetail


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// ANALYZER //


IBAnalyzerEMTrim::IBAnalyzerEMTrim(IBVoxCollection &voxels) :
    IBAnalyzerEM(voxels),
    m_Events(this->Events()),
    m_VoxCollection(NULL),
    m_MeanMuonVoxOccupancy(0)
{}

IBAnalyzerEMTrim::~IBAnalyzerEMTrim()
{}



void IBAnalyzerEMTrim::Run(unsigned int iterations, float muons_ratio, float a, float b)
{
    //IBVoxCollection -> IBVoxCollectionMedian
    assert(this->GetVoxCollection());
    IBAnalyzerEMTrimDetail::IBVoxCollectionATrim voxels_trim(this->GetVoxCollection()->GetDims());

    // IF WE HAVE A MEAN VALUE COPUTE AB TRIM //

    int sizeA = m_MeanMuonVoxOccupancy > 100 ? 100 * a : m_MeanMuonVoxOccupancy * a;
    int sizeB = m_MeanMuonVoxOccupancy > 100 ? 100 * b : m_MeanMuonVoxOccupancy * b;

    // Setting VoxCollection for Pimpl operations
    m_VoxCollection = this->GetVoxCollection();
    m_VoxCollectionMdn = &voxels_trim;

    // performs iterations //
    for (unsigned int it = 0; it < iterations; it++) {
        // copy forward VoxCollection into median image
        voxels_trim << (*this->GetVoxCollection());
        voxels_trim.SetABTrim(sizeA, sizeB);

        fprintf(stderr,"\r[%d muons] EM Trim -> performing iteration %i  occupancy=%i,a=%i,b=%i",
                (int) m_Events.size(), it, m_MeanMuonVoxOccupancy,sizeA, sizeB);
        Evaluate(muons_ratio);
        voxels_trim.UpdateDensity<IBAnalyzerEMTrimDetail::UpdateDensitySijCapATrimAlgorithm>(2);

        // copy back VoxCollection to Parent structure
        (*this->GetVoxCollection()) << voxels_trim;
    }

    printf("\nEM Trim -> done\n");
}

void IBAnalyzerEMTrim::Project(Event *evc) {

    // compute sigma //
    Matrix4f Sigma = Matrix4f::Zero();
    m_SijAlgorithm->ComputeSigma(Sigma, evc);
    // compute sij //
    m_SijAlgorithm->evaluate(Sigma,evc);
}

void IBAnalyzerEMTrim::BackProject(Event *evc) {

    // sommatoria della formula 38 //
    //#   pragma omp parallel for
    for (unsigned int i = 0; i < evc->elements.size(); ++i) {
        IBAnalyzerEMTrimDetail::IBVoxelABTrim *vox;
        std::vector<IBVoxel>::iterator itr(evc->elements[i].voxel);
        Id_t voxid = std::distance(m_VoxCollection->Data().begin(), itr);
        vox = &m_VoxCollectionMdn->operator [](voxid);
        vox->SijCap += evc->elements[i].Sij;
        vox->Count++;
    }
    //#   pragma omp barrier
}

void IBAnalyzerEMTrim::Evaluate(float muons_ratio) {

    // TODO: Move to iterators !!! //
    unsigned int start = 0;
    unsigned int end = (unsigned int) (m_Events.size() * muons_ratio);

    if(m_SijAlgorithm) {
        // Projection
#       pragma omp parallel for
        for (unsigned int i = start; i < end; ++i)
            this->Project(&m_Events[i]);
#       pragma omp barrier

        for (int i = start; i < end; ++i)
            this->BackProject(&m_Events[i]);

    }
    else {
        std::cerr << "Error: Lamda ML Algorithm not setted\n";
    }

    // Fill Mean of MuonEvents Voxel Occupancy //
    if(unlikely(m_MeanMuonVoxOccupancy == 0))
    {
        int size = m_VoxCollectionMdn->Data().size();
        for(int i=0; i<size; i++)
            m_MeanMuonVoxOccupancy += m_VoxCollectionMdn->At(i).Count;
        m_MeanMuonVoxOccupancy /= size;
    }

}




