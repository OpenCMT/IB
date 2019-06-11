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
#include <math.h>
#include <fstream>
#include <algorithm>

#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>
#include <TROOT.h>

#include <Core/Vector.h>
#include <Core/Debug.h>

#include "IBPocaEvaluator.h"
#include "IBMinimizationVariablesEvaluator.h"
#include "IBVoxRaytracer.h"

#include "IBVoxCollectionCap.h"
#include "IBAnalyzerEM.h"

#include "IBAnalyzerEMAlgorithm.h"
#include "IBAnalyzerEMAlgorithmSGA.h"

#include <string>
#include <map>

typedef IBAnalyzerEM::Event Event;
class UpdateDensitySijCapAlgorithm;

//________________________
void IBAnalyzerEM::Project(Event *evc){
    // compute sigma //
    Matrix4f Sigma = Matrix4f::Zero();
    m_SijAlgorithm->ComputeSigma(Sigma, evc);
    // compute sij //
    m_SijAlgorithm->evaluate(Sigma,evc);
}

//________________________
void IBAnalyzerEM::BackProject(Event *evc){

    IBVoxel *vox = NULL;
    // sommatoria della formula 38 //
    for (unsigned int j = 0; j < evc->elements.size(); j++) {
      vox = evc->elements[j].voxel;
      if( vox==NULL || std::isnan(evc->elements[j].Sij) || std::isnan(vox->SijCap))
        continue;
      #pragma omp atomic
      vox->SijCap += evc->elements[j].Sij;
//        {
//            //            IBVoxel *v0 = &(*m_parent->GetVoxCollection()->Data().begin());
//            //            int id = (vox-v0)/sizeof(IBVoxel);
//            if( isnan(evc->elements[j].Sij) ) {
//                std::cout << "nan Sij in vox:" << vox << " mu:" << evc << "\n" << std::flush;
//            }
//            if( isnan(vox->SijCap) ) {
//                std::cout << "nan SijCap in vox:" << vox << " mu:" << evc << "\n" << std::flush;
//            }
//        }
	#pragma omp atomic
        vox->Count++;
    }
}

//________________________
void IBAnalyzerEM::Evaluate(float muons_ratio)
{
    unsigned int start = 0;
    // the following line gives wrong end value in some situations, don't know why......
    //unsigned int end = (unsigned int) (std::floor(m_Events.size() * muons_ratio));
    unsigned int end = (unsigned int) (m_Events.size());
    unsigned int ev = start;

    std::cout << "IBAnalyzerEM::Evaluate form start " << start << " to end " << end << " collection size " << m_Events.size() << " muons ratio " << muons_ratio << std::endl;

    if(m_SijAlgorithm) {
      // Projection
      #pragma omp parallel for
      for (unsigned int i = start; i < end; ++i){
    this->Project(&m_Events[i]);
      }
      #pragma omp barrier

      // Backprojection
      #pragma omp parallel for
      for (unsigned int i = start; i < end; ++i){
          this->BackProject(&m_Events[i]);
          ev++;
      }
      #pragma omp barrier
    }
    else {
        std::cerr << "Error: Lamda ML Algorithm not set\n";
    }

    //------------- SIJ RANK NEEDED FOR ALL SIJ RANK STUFF
    // if(m_firstIteration){
    //   m_firstIteration = false;

    //   //---- Step 1: Map voxel ('elements[j]') to sij
    //   std::map<IBVoxel*,std::vector<float> > voxelMap;
    //   for (unsigned int i = start; i < end; ++i){
    // 	for (unsigned int j = 0; j < m_Events[i].elements.size(); ++j) {
    // 	  IBVoxel* vox = m_Events[i].elements[j].voxel;
    // 	  float sij = m_Events[i].elements[j].Sij;
    // 	  if(voxelMap.find(vox)==voxelMap.end()) voxelMap[vox] = std::vector<float>();
    // 	  voxelMap[vox].push_back(sij);
    // 	}
    //   }
    //   //std::cout << "\n>> Got a voxel map of size " << voxelMap.size() << std::endl;

    //   //------------- SIJ RANK KILL
    //   //---- NEED TO ALSO SET muon_event.SetMomentum(7); in *MUBLAST*p
    //   //   //----
    // //   std::vector<Event*> killList;
    // //   for (unsigned int i = start; i < end; ++i) {
    // // 	float rankSum = 0.;
    // // 	for (unsigned int j = 0; j < m_Events[i].elements.size(); ++j) {
    // //  	  IBVoxel* vox = m_Events[i].elements[j].voxel;
    // //  	  float sij = m_Events[i].elements[j].Sij;
    // // 	  std::vector<float>& v = voxelMap[vox];
    // // 	  std::sort(v.begin(), v.end());
    // // 	  int sijPos = 0;
    // // 	  for(; sijPos < v.size() && v[sijPos] != sij; sijPos++);
    // // 	  rankSum += float(sijPos - v.size()/2)/sqrt(float(v.size()));
     // // 	}
    // // 	rankSum = rankSum/sqrt(m_Events[i].elements.size());
    // // 	if(rankSum < m_rankLimit) killList.push_back(&m_Events.at(i));
    // //   }

    // //   for(std::vector<Event*>::iterator it=killList.begin(); it!=killList.end(); it++){
    // // 	m_Events.remove_element(**it);
    // //   }
    // // }
    //   //------------- END SIJ RANK KILL

    //   //------------- END SIJ RANK FIND
    // //   //---- Step 3.5 Make TFile and TTree
    // //   gROOT->ProcessLine("#include <vector>");
    // //   TFile* tfile = new TFile("joelout.root","RECREATE");
    // //   TTree* ttree = new TTree("joeltree","joeltree");
    // //   //      std::vector<float> sigmaMean;
    // //   //      std::vector<float> sigmaMedian;
    // //   std::vector<int> position;
    // //   std::vector<int> size;
    // //   float muonMomentum;
    // //   //      ttree->Branch("sigmaMean", &sigmaMean);
    // //   //      ttree->Branch("sigmaMedian", &sigmaMedian);
    // //   ttree->Branch("position", &position);
    // //   ttree->Branch("size", &size);
    // //   ttree->Branch("muonMomentum", &muonMomentum);

    // //   //---- Step 4: For each muon i calculate the sij distribution
    // //   for (unsigned int i = start; i < end; ++i) {
    // // 	//    	sigmaMean.clear();
    // // 	//    	sigmaMedian.clear();
    // // 	position.clear();
    // // 	size.clear();
    // // 	for (unsigned int j = 0; j < m_Events[i].elements.size(); ++j) {
    // // 	  IBVoxel* vox = m_Events[i].elements[j].voxel;
    // // 	  float sij = m_Events[i].elements[j].Sij;

    // // 	  //----
    // // 	  std::vector<float>& v = voxelMap[vox];
    // // 	  std::sort(v.begin(), v.end());
    // // 	  int sijPos = 0;
    // // 	  for(; sijPos < v.size() && v[sijPos] != sij; sijPos++);
    // // 	  position.push_back(sijPos);
    // // 	  size.push_back(v.size());

    // // 	  // //----
    // // 	  // if(voxelMap_stats.find(vox) == voxelMap_stats.end()) continue;
    // // 	  // std::vector<float>& v_stats = voxelMap_stats[vox];
    // // 	  // if(v_stats[1] == 0) continue;
    // // 	  // sigmaMean.push_back((v_stats[0]-sij)/v_stats[1]);
    // // 	  // sigmaMedian.push_back((v_stats[2]-sij)/v_stats[3]);
    // // 	}
    // // 	muonMomentum = m_Events[i].header.InitialSqrP;
    // // 	std::cout << m_Events[i].header.InitialSqrP << std::endl;
    // // 	ttree->Fill();
    // //   }
    // //   ttree->Write();
    // //   tfile->Close();
    // // }
      //------------- END SIJ RANK FIND END
}


//________________________

////////////////////////////////////////////////////////////////////////////////
////// CUTS ////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

/// filter events after voxel mask has been applied
void IBAnalyzerEM::filterEventsVoxelMask()
{
  std::cout << "\nIBAnalyzerEM: Removing frozen voxels from " << this->m_Events.size() << " muon collection." << std::endl;
    Vector< Event >::iterator itr = this->m_Events.begin();
    const Vector< Event >::iterator begin = this->m_Events.begin();

    while (itr != this->m_Events.end()) {
        Event & evc = *itr;
        Vector< Event::Element >::iterator itre = evc.elements.begin();
        // create new vector and fill it with positive voxels
        Vector< Event::Element > newvelc;
        while (itre != evc.elements.end()) {
            Event::Element & elc = *itre;
            if(elc.voxel->Value <= 0){
                // add contribution to E matrix in both views
                evc.header.E.block<2,2>(2,0) += elc.Wij * fabs(elc.voxel->Value) * evc.header.InitialSqrP;
                evc.header.E.block<2,2>(0,2) += elc.Wij * fabs(elc.voxel->Value) * evc.header.InitialSqrP;
            }
            else
                newvelc.push_back(elc);
            ++itre;
        }
//        if(newvelc.size() != evc.elements.size()){
//            std::cout << newvelc.size() << " voxels left from muon of " << evc.elements.size() << " voxels " << std::endl;
//            std::cout << this->m_parent->m_MuonCollection->Data().at(itr-begin);
//        }

        evc.elements = newvelc;

        /// erase event and muon with empty voxel collection
        if(evc.elements.empty()){
            unsigned int pos = itr - begin;
            this->m_Events.remove_element(evc);
            if(m_MuonCollection)
                m_MuonCollection->Data().remove_element(pos);
        }
        else
            ++itr;
    }
    std::cout << " " << this->m_Events.size() << " muons left!" << std::endl;

    return;
}

//________________________
////////////////////////////////////////////////////////////////////////////////
/// filter events if in-out line distance out of range
void IBAnalyzerEM::filterEventsLineDistance(float min, float max)
{
  std::cout << "\n*** Removing events with line distance out of range from " << this->m_Events.size() << " muon collection." << std::endl;

    Vector< Event >::iterator itr = this->m_Events.begin();
    const Vector< Event >::iterator begin = this->m_Events.begin();

    while (itr != this->m_Events.end()) {
        Event & evc = *itr;
        unsigned int pos = itr - begin;

        MuonScatterData muon = m_MuonCollection->At(pos);
        bool use_poca = m_PocaAlgorithm->evaluate(muon);
        float dist = m_PocaAlgorithm->getDistance();

        /// erase event and muon with distance out of range
        if(!isFinite(dist) || dist >= max || dist < min){

            this->m_Events.remove_element(evc);
            m_MuonCollection->Data().remove_element(pos);
        }
        else
	  ++itr;
    }

    std::cout << " " << this->m_Events.size() << " muons left!" << std::endl;

    return;
}

//________________________
////////////////////////////////////////////////////////////////////////////////
/// SijCut RECIPE1:  (true if Sij cut proposed) //
static bool em_test_SijCut(const Event &evc, float cut_level, int &nvox_cut){
    nvox_cut = 0;
//    std::cout << "Testing treshold " << cut_level << std::endl;

   for (unsigned int i = 0; i < evc.elements.size(); i++) {
        const Event::Element &el = evc.elements[i];
        if (fabs( (el.Sij * el.voxel->Count - el.voxel->SijCap)
                  / el.voxel->SijCap ) > cut_level) nvox_cut++;
//        if(isnan(el.voxel->SijCap)){
//            std::cout << "ATTENTION: nan Sij sum!....";
//            /// this voxel should increment n_cuts, which is not at the moment. FIX it in the future
//            /// this "bug" consequence is that 2% of events is nasty and not classified with SijCutCount
//        }
    }
//    std::cout << "\n n_cuts=" << nvox_cut << ", 1/3 of the voxels=" << evc.elements.size()/3  << std::endl;
    if (nvox_cut > (int)(evc.elements.size()/3) ) return true;
    else return false;

    /// MEMORANDUM : testing this function I see that voxel Count doubles at each call ... WHY?
    /// SijCap updates, so the Sij incrementum is the same, no difference in computation.
    /// 20160408 SV stop inquiring.... no time now!
}

//________________________
float IBAnalyzerEM::SijMedian(const Event &evc)
{
    Evaluate(1);

    Vector< float > Si;
    for (unsigned int i = 0; i < evc.elements.size(); i++) {
        const Event::Element &el = evc.elements[i];
        Si.push_back(fabs( (el.Sij * el.voxel->Count - el.voxel->SijCap) / el.voxel->SijCap ));
    }

    // compute median
    float Smedian;
    std::sort(Si.begin(),Si.end());
   int nS = Si.size();
   if(nS){
       if(nS%2)
           Smedian = Si[nS / 2];
       else
           Smedian = (Si[nS / 2 - 1] + Si[nS / 2]) / 2;
    }
//    /// debug
//    if(Smedian==0){
//        for(int i=0; i<nS; i++) std::cout << Si[i] << ",";
//        std::cout << "\n     MEDIAN =" << Smedian << std::endl;
//    }
    return Smedian;
}

//________________________
Vector<Event > IBAnalyzerEM::SijCutCount(float threshold_low, float threshold_high)
{
    Evaluate(1);
    //std::cout << "Cut tresholds : " << std::dec << threshold_low << ", " << threshold_high << " ... " << std::endl;
    Vector< Event > ve;
    Vector< Event >::iterator itr = this->m_Events.begin();
    const Vector< Event >::iterator begin = this->m_Events.begin();
    int evnum = 0;
    int nvox_cut_low = 0;
    int nvox_cut_high = 0;

    while (itr != this->m_Events.end()) {
//        std::cout << "\n\n *** Event " << evnum << std::endl;
        if(em_test_SijCut(*itr, threshold_low,nvox_cut_low) && !em_test_SijCut(*itr, threshold_high,nvox_cut_high)){
            ve.push_back(*itr);
//            std::cout << "Passed!\n";
        }
        evnum++;
        ++itr;
    }
//    std::cout << "SijCutCount: muons between tresholds ["
//              << threshold_low << "," << threshold_high << "] = "
//              << ve.size() << " over " << m_Events.size() << "\n";
    return ve;
}
//________________________
void IBAnalyzerEM::dumpEventsSijInfo(const char *name, Vector<float> N)
{
    Evaluate(1);
    /// dump event Sij info on file
    std::fstream fout;
    fout.open(name, std::fstream::out | std::fstream::app);

    Vector< Event > ve;
    Vector< Event >::iterator itr = this->m_Events.begin();
    const Vector< Event >::iterator begin = this->m_Events.begin();
    int nev = 0;

    /// loop over events
    while (itr != this->m_Events.end()) {
        //std::cout << "\n\n *** Event " << nev << std::endl;
        Event & evc = *itr;
        /// momentum
        float p0sq = 3. * 3.;
        //float mom = sqrt(p0sq/evc.header.InitialSqrP);
        float mom =  evc.header.pTrue;

        /// compute and dump median
        Vector< float > Si;
        float median;
        Vector< Event::Element >::iterator itre = evc.elements.begin();
        while (itre != evc.elements.end()) {
            Event::Element & elc = *itre;
            float Nij = fabs( (elc.Sij * elc.voxel->Count - elc.voxel->SijCap) / elc.voxel->SijCap );
            Si.push_back(Nij);
            ++itre;
        }
        std::sort(Si.begin(),Si.end());
        int nS = Si.size();
        if(nS){
           if(nS%2)
               median = Si[nS / 2];
           else
               median = (Si[nS / 2 - 1] + Si[nS / 2]) / 2;
        }

        fout << nev << " " << mom << " ";
        // fout << evc.elements.size() << " ";
        //float median =  SijMedian(*itr);
        fout << median << " ";

//        /// loop over Sij tresholds
//        Vector< float >::iterator itrN = N.begin();
//        const Vector< float >::iterator beginN = N.begin();
//        while (itrN != N.end()) {
//            int nvox_cut = 0;
//            // dump number of voxels above the threshold
//            em_test_SijCut(*itr, *itrN, nvox_cut);
//            fout << nvox_cut << " ";
//            // dump N cut bin low edge
//            float threshold_low = *itrN;
//            float threshold_high = *(itrN+1);
//            if(itrN+1==N.end())
//                threshold_high = 10000;
//             if(em_test_SijCut(*itr, threshold_low,nvox_cut) && !em_test_SijCut(*itr, threshold_high,nvox_cut))
//                 fout << threshold_low << " ";

//            ++itrN;
//        }

        fout << "\n";
        nev++;
        ++itr;
    }
    fout.close();
    return;
}
//________________________
void IBAnalyzerEM::SijCut(float threshold)
{
    Evaluate(1);

    Vector< Event >::iterator itr = this->m_Events.begin();
    const Vector< Event >::iterator begin = this->m_Events.begin();

    int count = 0;
    int nvox_cut=0;
    while (itr != this->m_Events.end()) {
        if(em_test_SijCut(*itr, threshold, nvox_cut))
        {
            unsigned int pos = itr - begin;
            this->m_Events.remove_element(*itr);
            if(m_MuonCollection)
                m_MuonCollection->Data().remove_element(pos);
            count ++;
        }
        else ++itr;
    }
    std::cout << "SijCut removed muons: " << count << "\n" << std::endl;

    GetVoxCollection()->UpdateDensity<UpdateDensitySijCapAlgorithm>(0);   // HARDCODE THRESHOLD
}

//________________________
///////////////////////////////////////////////////////////////////////////////////////
void IBAnalyzerEM::SetSijMedianMomentum()
{
    Evaluate(1);

    Vector< Event >::iterator itr = this->m_Events.begin();

//    std::cout << "SetSijMedianMomentum \n"
//              << "old invp2 = " << itr->header.InitialSqrP/(m_parent->nominal_momentum * m_parent->nominal_momentum);

    while (itr != this->m_Events.end()) {
        Event & evc = (*itr);
        float m = SijMedian(evc);

        //fit function for pguess from Sij median
        /// 20160530 entire furnace
        //1/p2=13.9786 +-47.4669*1/log(x) + 78.2999*1/(x)
        //float pguess = sqrt(1/ (13.9786 - 47.4669*1/log(m) + 78.2999*1/(m)));
        /// 20160926
        float invp2guess = 0.0375 + (0.0333 *m)-(0.0002 *m*m);

        if(std::isnan(invp2guess)){
            std::cout << "ATTENTION nan invp2!!" << std::endl;
        invp2guess = 0.0004;
        }


        // cut off if p>50GeV i.e. 1/p2 < 0.0004
        if(invp2guess<0.0004)
            invp2guess = 0.0004;

        /// set InitialSqrP variable
        itr->header.InitialSqrP = nominal_momentum * nominal_momentum * invp2guess;

        /// set in addition every voxel pw, since from voxel momentum code it is used in the algorithm
        Vector< Event::Element >::iterator itre = evc.elements.begin();
        while (itre != evc.elements.end()) {
            Event::Element &elc = *itre;
            elc.pw = nominal_momentum * nominal_momentum *invp2guess;
//            std::cout << "Setting voxel pw.... " << elc.pw << std::endl;
            ++itre;
        }

//        std::cout  << ", median = " << m << ", invp2guess " << invp2guess
//                  << ", pguess = " <<  sqrt(m_parent->nominal_momentum *m_parent-> nominal_momentum/itr->header.InitialSqrP) << std::endl;
        itr++;
    }

    GetVoxCollection()->UpdateDensity<UpdateDensitySijCapAlgorithm>(0);   // HARDCODE THRESHOLD

}


//________________________
////////////////////////////////////////////////////////////////////////////////
void IBAnalyzerEM::Chi2Cut(float threshold)
{
    Evaluate(1);

    std::vector< Event >::iterator itr = this->m_Events.begin();
    const Vector< Event >::iterator begin = this->m_Events.begin();

    do {
        Matrix4f Sigma = Matrix4f::Zero();
        Event &evc = *itr;
        this->m_SijAlgorithm->ComputeSigma(Sigma,&evc);
        Matrix4f iS = Sigma.inverse();
        Matrix4f Dn = iS * (evc.header.Di * evc.header.Di.transpose());
        if ( Dn.trace() > threshold ){
            unsigned int pos = itr - begin;
            this->m_Events.remove_element(*itr);
            if(m_MuonCollection)
                m_MuonCollection->Data().remove_element(pos);
        }
        else ++itr;
    } while (itr != this->m_Events.end());

    GetVoxCollection()->UpdateDensity<UpdateDensitySijCapAlgorithm>(0);   // HARDCODE THRESHOLD
}



//________________________
////////////////////////////////////////////////////////////////////////////////
////// UPDATE DENSITY ALGORITHM ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class UpdateDensitySijCapAlgorithm :
        public IBInterface::IBVoxCollectionStaticUpdateAlgorithm
{
public:
    static void UpdateDensity(IBVoxCollection *voxels, unsigned int threshold)
    {
        for(unsigned int i=0; i< voxels->Data().size(); ++i) {
            IBVoxel& voxel = voxels->Data()[i];
            unsigned int tcount = voxel.Count;
            if ( voxel.Value > 0 && tcount > 0 && (threshold == 0 || tcount >= threshold) ) {
                voxel.Value += voxel.SijCap / static_cast<float>(tcount);
                if(unlikely(!isFinite(voxel.Value) || voxel.Value > 100.E-6)) {  // HARDCODED!!!
                    voxel.Value = 100.E-6;
                }
                //                 else if (unlikely(voxel.Value < 0.)) voxel.Value = 0.1E-6;
            }
            // else
            //             voxel.Value = 0;
            voxel.SijCap = 0;
        }
    }
};



////////////////////////////////////////////////////////////////////////////////
// IB ANALYZER EM  /////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

//___________________________
IBAnalyzerEM::IBAnalyzerEM(IBVoxCollection &voxels, int nPath, double alpha, bool useRecoPath,
               bool oldTCalculation, float rankLimit, IBVoxCollection* initialSqrPfromVtk, int pVoxelMean) :
    m_PocaAlgorithm(NULL),
    m_VarAlgorithm(NULL),
    m_RayAlgorithm(NULL),
    m_UpdateAlgorithm(NULL),
    m_nPath(nPath),
    m_alpha(alpha),
    m_useRecoPath(useRecoPath),
    m_oldTCalculation(oldTCalculation),
    m_rankLimit(rankLimit),
    m_initialSqrPfromVtk(initialSqrPfromVtk),
    m_pVoxelMean(pVoxelMean),
    m_SijAlgorithm(NULL),
    m_firstIteration(false)
{
  //---- Print the settings
  std::cout << "Using alpha = " << m_alpha << ", #path = " << m_nPath << std::endl;
  std::cout << "Reco path ("    << m_useRecoPath      << "), "
	    << "Old T ("        << m_oldTCalculation  << "), " << std::endl;
  BaseClass::SetVoxCollection(&voxels);
  nominal_momentum = 3;
}

//___________________________
IBAnalyzerEM::~IBAnalyzerEM()
{}

//___________________________
Vector<IBAnalyzerEM::Event> &IBAnalyzerEM::Events(){
    return m_Events;
}

//___________________________
void IBAnalyzerEM::DumpEvent(Event *evc){

    std::cout << "\n ----------- Dump Event: " << std::endl;
    std::cout << "InitialSqrP " << evc->header.InitialSqrP <<  ", Di " << evc->header.Di << ", E " << evc->header.E << std::endl;
    std::cout << "Voxel vector size " << evc->elements.size();

    Vector< Event::Element >::iterator itre = evc->elements.begin();
    int ivox=0;
    while (itre != evc->elements.end()) {
        Event::Element &elc = *itre;
        std::cout  << "Voxel " << ivox << ": value " << elc.voxel->Value << ", SijCap " << elc.voxel->SijCap << ", Count " << elc.voxel->Count;
        std::cout << ", Wij[0,0] " << elc.Wij(0,0) << ", lambda " << elc.lambda << ", pw " << elc.pw << ", Sij " << elc.Sij << std::endl;
        ++itre;
        ivox++;
    }

    return;
}

//___________________________
//---- **Deprecated** version of IBAnalyzerEM::AddMuon, please use IBAnalyzerEM::AddMuonFullPath
//---- Only has 2-path mode, and bugs in the calculation of L and T
bool IBAnalyzerEM::AddMuon(const MuonScatterData &muon){
  if(unlikely(!m_RayAlgorithm || !m_VarAlgorithm)) return false;
  Event evc;

  evc.header.InitialSqrP = pow(nominal_momentum/muon.GetMomentum() ,2);
  if(std::isnan(evc.header.InitialSqrP)) std::cout << "sono in AddMuon: nominalp:" << nominal_momentum << " muon.GetMomentum():" << muon.GetMomentum() <<"\n"
					      << std::flush;
  //    DBG(trd,evc.header.InitialSqrP,"invP2/F");
  if(likely(m_VarAlgorithm->evaluate(muon))) {
    evc.header.Di = m_VarAlgorithm->getDataVector();
    evc.header.E  = m_VarAlgorithm->getCovarianceMatrix();

    // HARDCODED ... ZERO CROSS CORRELATION BETWEEN VARS //
    //        evc.header.E(0,1) = 0.;
    //        evc.header.E(1,0) = 0.;
    //        evc.header.E(2,3) = 0.;
    //        evc.header.E(3,2) = 0.;
    // .................................................. //

    // HARDCODED ... ZERO CROSS CORRELATION BETWEEN VIEWS //
    //        evc.header.E.block<2,2>(2,0) = Matrix2f::Zero();
    //        evc.header.E.block<2,2>(0,2) = Matrix2f::Zero();
    // .................................................. //

    // HARDCODED ... LESS ERROR ! //
    //evc.header.E = Matrix4f::Zero();
    //evc.header.E /= 2;
    //        std::cout
    //                << " evc.header.Di " << evc.header.Di.transpose() << "\n"
    //                << " evc.header.E " << evc.header.E << "\n";


  }
  else return false;

  IBVoxRaytracer::RayData ray;
  { // Get RayTrace RayData //
    HPoint3f entry_pt,poca,exit_pt;
    if( !m_RayAlgorithm->GetEntryPoint(muon.LineIn(),entry_pt) ||
	!m_RayAlgorithm->GetExitPoint(muon.LineOut(),exit_pt) )
      return false;

    bool use_poca = false;
    if(m_PocaAlgorithm) { //TODO:  move this to poca algorithm
      use_poca = m_PocaAlgorithm->evaluate(muon);
      poca = m_PocaAlgorithm->getPoca();
      //            DBG(trd,poca,"x/F:y/F:z/F:h/F");

      HVector3f in, out;
      in  = poca - muon.LineIn().origin;
      out = muon.LineOut().origin - poca;
      float poca_prj = in.transpose() * out;
      //            DBG(trd,poca_prj);
      use_poca &= ( poca_prj > 0 );
    }
    if(use_poca && this->GetVoxCollection()->IsInsideBounds(poca)) {
      poca = m_PocaAlgorithm->getPoca();
      ray = m_RayAlgorithm->TraceBetweenPoints(entry_pt,poca);
      ray.AppendRay( m_RayAlgorithm->TraceBetweenPoints(poca,exit_pt) );
    }
    else {
      ray = m_RayAlgorithm->TraceBetweenPoints(entry_pt,exit_pt);
    }
  }

  Event::Element elc;
  Scalarf T = ray.TotalLength();
  for(int i=0; i<ray.Data().size(); ++i)
    {
      // voxel //
      const IBVoxRaytracer::RayData::Element &el = ray.Data().at(i);

      if(el.vox_id >= this->GetVoxCollection()->GetDims().prod()){
        std::cout << "ATTENTION voxel ID > size collection!! " << std::endl;
        return false;
      }
      elc.voxel = &this->GetVoxCollection()->operator [](el.vox_id);

      // Wij   //
      Scalarf L = el.L;  T = fabs(T-L);
      elc.Wij << L ,          L*L/2 + L*T,
	L*L/2 + L*T, L*L*L/3 + L*L*T + L*T*T;
      // pw    //
      elc.pw = evc.header.InitialSqrP;


      if(elc.voxel->Value <= 0){
	// add both views
	evc.header.E.block<2,2>(2,0) += elc.Wij * fabs(elc.voxel->Value) * evc.header.InitialSqrP;
	evc.header.E.block<2,2>(0,2) += elc.Wij * fabs(elc.voxel->Value) * evc.header.InitialSqrP;
      }
      else
    evc.elements.push_back(elc);
    }
  m_Events.push_back(evc);

  //    trd.Fill();
  return true;
}

//___________________________
//---- The "new" version of IBAnalyzerEM::AddMuon
//---- Compared to IBAnalyzerEM::AddMuon, this function:
//---- 1) Fixes bugs in L and T
//---- 2) Allows one to use the true muon path (provided by the argument muonPath)
//---- 3) Allows one to use any arbitrary path (see 3-path implementation)
bool IBAnalyzerEM::AddMuonFullPath(const MuonScatterData &muon, Vector<HPoint3f>& muonPath){

  bool debug = false;

  //-------------------------
  //---- STEP #1: Fill the event info

  //---- Check the ray algo (calculates the ray parameters)
  //---- and the variable algo (calculates the scattering/displacement variables)
  if(unlikely(!m_RayAlgorithm || !m_VarAlgorithm)){
      std::cout << "No RayAlgorithm or VarAlgorithm set.... EXITING!" << std::endl;
      return false;
  }

  Event evc; //<---- The event info
  evc.header.Di << NAN, NAN, NAN,NAN;
  evc.header.E << NAN, NAN, NAN,NAN,NAN, NAN, NAN,NAN,NAN, NAN, NAN,NAN,NAN, NAN, NAN,NAN;
  evc.header.InitialSqrP = NAN;
  evc.header.pTrue = NAN;
  Vector<Event::Element> voxVec(0);
  evc.elements = voxVec;
  evc.elements.clear();

  if(likely(m_VarAlgorithm->evaluate(muon))) {
    //---- Get the Data (Di) and Error (E) matrices
    evc.header.Di = m_VarAlgorithm->getDataVector();
    evc.header.E  = m_VarAlgorithm->getCovarianceMatrix();
    //---- Momentum square (the "$$" notation is a bit much...)
    evc.header.InitialSqrP = pow(nominal_momentum/muon.GetMomentum() ,2);
    // SV for Sij studies
    evc.header.pTrue = muon.GetMomentumPrime();

    if(std::isnan(evc.header.InitialSqrP)){
      std::cout << "AddMuonFullPath: nominalp:" << nominal_momentum
		<< "muon.GetMomentum():" << muon.GetMomentum() <<"\n" << std::endl;
    }
  }
  else {
      if(debug)
        std::cout << "Evaluate muon variable failed... EXITING!" << std::endl;
      return false;
  }

  //-------------------------
  //---- STEP #2: Perform raytracing

  //---------
  //---- STEP #2.1: Get all of the path points
  Vector<HPoint3f> pts;       //<---- A vector of the muon path points
  HPoint3f front_pt, back_pt; //<---- The first and last point

  //---- If reconstructing the muon's path (i.e. NOT the true path)
  if(m_useRecoPath){

    //---- Require entry and exit points
    HPoint3f entry_pt, exit_pt;
    if( !m_RayAlgorithm->GetEntryPoint(muon.LineIn(),entry_pt) ) {
        if(debug)
            std::cout << "No entry point.... EXITING!" << std::endl;
        return false;
    }
    if( !m_RayAlgorithm->GetExitPoint(muon.LineOut(),exit_pt) ) {
        if(debug)
            std::cout << "No exit point.... EXITING!" << std::endl;
        return false;
    }
    front_pt = entry_pt;
    back_pt = exit_pt;

    //---- Calculate the track length
    double trackLength = (exit_pt-entry_pt).norm();

    //---- Add the entry point
    pts.push_back(entry_pt);

    //---- If we want to build tracks with more than one line
    if(m_nPath > 1){

      //---- Evaluate the POCA and check that it is valid
      if(m_PocaAlgorithm && m_PocaAlgorithm->evaluate(muon)){;
  	HPoint3f poca = m_PocaAlgorithm->getPoca();

  	//---- Check that the POCA is valid
  	HVector3f in, out;
  	in  = poca - muon.LineIn().origin;
  	out = muon.LineOut().origin - poca;
  	float poca_prj = in.transpose() * out;
  	bool validPoca = poca_prj > 0 && GetVoxCollection()->IsInsideBounds(poca);

  	//---- If using the two-line path
  	if(m_nPath==2){
	  if(validPoca) pts.push_back(poca);
	}
  	//---- If using the three-line path
  	else if(m_nPath == 3){
  	  //---- Get the poca on the entry/exit tracks
  	  HPoint3f entry_poca = m_PocaAlgorithm->getInTrackPoca();
  	  HPoint3f exit_poca  = m_PocaAlgorithm->getOutTrackPoca();

  	  //---- Get the distance along the tracks to the inflection points
  	  double entry_length = (entry_pt - entry_poca).norm();
  	  double exit_length  = (exit_pt  - exit_poca).norm();

	  //---- Remove points which are unreasonably far away
      if(entry_length > trackLength || exit_length > trackLength){
          if(debug)
              std::cout << "Track length smaller than entry-exit distance.... EXITING!" << std::endl;
          return false;
      }

  	  //---- Get the inflection points
  	  double normIn  = muon.LineIn().direction.norm(); //<----Recently added
  	  double normOut = muon.LineOut().direction.norm();
  	  HVector3f point1 = entry_pt + m_alpha*(entry_length/normIn)*muon.LineIn().direction;
  	  HVector3f point2 = exit_pt  - m_alpha*(exit_length/normOut)*muon.LineOut().direction;

  	  //---- Add the points to the collection
  	  pts.push_back(point1);
  	  pts.push_back(point2);
  	}

  	//---- Get exit point
  	pts.push_back(exit_pt);
      }
    }
  }
  //---- If reconstructing the muon's true path
  else{
    //---- If muonPath is empty, return
    if(muonPath.size()==0){
        if(debug)
            std::cout << "Empty muon path.... EXITING!" << std::endl;
        return false;
    }
    //---- Otherwise append all the points
    for(int i=0; i<muonPath.size(); ++i) pts.push_back(muonPath[i]);
    back_pt  = pts.back();
    front_pt = pts.front();
  }

  //---------
  //---- STEP #2.2: Remove mid-voxel inflections
  std::map<int,std::vector<HPoint3f> > voxelMap; //<---- Map of voxel number to points in the voxel
  std::vector<int> voxelOrder; //<---- An ordered list of the voxels
  Scalarf totalLength = 0.;    //<---- Total length of the muon path

  //---- Loop over the points
  for(int i=0; i<pts.size()-1; ++i){
    //---- Get the points
    HPoint3f& pt1 = pts[i];
    HPoint3f& pt2 = pts[i+1];

    //---- Get the ray properties between the points
    Scalarf  rayLength = (pt2-pt1).norm();
    HPoint3f rayDir    = (pt2-pt1)/rayLength;

    IBVoxRaytracer::RayData ray = m_RayAlgorithm->TraceBetweenPoints(pt1,pt2);

    //---- Loop over the voxels in the ray
    float cumulativeLength = 0.;
    foreach(const IBVoxRaytracer::RayData::Element &el, ray.Data()){

      //---- Project the muon track to the current position
      HPoint3f pti = pt1 + cumulativeLength*rayDir;

      //---- Project the muon track to the next position
      cumulativeLength += el.L; //<---- Length of the muon track in the voxel
      HPoint3f ptj = pt1 + cumulativeLength*rayDir;

      //---- If this muon HAS NOT crossed this voxel before
      if(voxelMap.find(el.vox_id) == voxelMap.end()){
	//---- Add the voxel and points to the map
	voxelMap[el.vox_id] = std::vector<HPoint3f>();
	voxelMap[el.vox_id].push_back(pti);
	voxelMap[el.vox_id].push_back(ptj);
	//---- Add the voxel to the ordered list
	voxelOrder.push_back(el.vox_id);
      }
      //---- Otherwise, update the final point in the voxel
      voxelMap[el.vox_id][1] = ptj;
    }
    //---- Increment the total length of the muon path
    totalLength += rayLength;
  }

  //---------
  //---- STEP #2.3: Remove mid-voxel inflections

  //---- Now trace between points and calculate length parameters
  Scalarf H = muon.LineIn().direction.transpose() * (back_pt - front_pt);
  Scalarf normIn = muon.LineIn().direction.norm();

  if(!m_oldTCalculation) H = H/normIn;  //<---- The old (buggy) calculation of T needs this value of H
  Scalarf T = totalLength; //<---- Needed if using the old (buggy) calculation of T


  //std::cout << "\n *** Muon momentum " << muon.GetMomentum() << ", momentum prime " << muon.GetMomentumPrime() << std::endl;
  /// 20160731 SV voxel momentum implementation
  /// pVoxelMean=0 NO voxel momentum
  /// pVoxelMean=1 descending fixed voxel momentum
  /// pVoxelMean=2 descending fixed voxel momentum with angle dependency
  /// pVoxelMean=3 voxel momentum from real momentum
  /// pVoxelMean=4 voxel momentum from real momentum p classes
  /// pVoxelMean=5,6,7 voxel momentum from real momentum p classes I,II,II
  ///
  // parameters
  float totalLengthFurnace = 0.;
  float deltaP = 0;
  float invp2_IN = 0.;
  float invp2_OUT = 0.;
  bool noAddMuon = false;

  if(m_pVoxelMean==1){
    // compute p_in and p_out from fixed parameters, average values in angle range [60,90]
    // NB p_in 8.72177, p_out 2.2526, Dp 6.46917
    invp2_IN = 0.0131459;
    invp2_OUT = 0.197075;
  } else if (m_pVoxelMean==2){
    // compute p_in and p_out with a angle-dependent function
    float a = fabs((3.14159265359 - acos(muon.LineIn().direction[1]/normIn))/3.14159265359*180);
    invp2_IN = -0.1054 + (0.003615*a) - (0.0000269*a*a);
    invp2_OUT = -0.5815 + (0.0309*a) - (0.0002711*a*a);
    if(invp2_IN<0 || invp2_OUT<0)
      std::cout << "ATTENTION: negative 1/p2 in pVoxelMean calculation...." << std::endl;
  } else if (m_pVoxelMean==3){
      // p_in and p_out real values
      invp2_IN = 1/(muon.GetMomentum()*muon.GetMomentum());
      invp2_OUT = 1/(muon.GetMomentumPrime()*muon.GetMomentumPrime());
    } else if (m_pVoxelMean>=4){
      // use average values from 3 momentum classes. run500all. Angle range [60,90]
      // p range [0,0.5] IN <1/p2> mean : 0.0658179, OUT <1/p2> mean : 5.80614
      // p range [0.5,1] IN <1/p2> mean : 0.0570832, OUT <1/p2> mean : 1.89081
      // p range [1,10000] IN <1/p2> mean : 0.0106817, OUT <1/p2> mean : 0.0772919

      float pclass = muon.GetMomentumPrime();
//      // I class
//      if(pclass <= 0.5){
//          invp2_IN = 0.0658179;
//          invp2_OUT = 5.80614;
//          if(m_pVoxelMean != 5) noAddMuon = true;
//      }
//      // II class
//      else if(pclass > 0.5 && pclass<=1.0){
//          invp2_IN = 0.0570832;
//          invp2_OUT = 1.89081;
//          if(m_pVoxelMean !=6 ) noAddMuon = true;
//      }
//        // III class
//        else if(pclass >1.0){
//            invp2_IN = 0.0106817;
//            invp2_OUT = 0.0772919;
//            if(m_pVoxelMean !=7 ) noAddMuon = true;
//        }
        // II class
        if(pclass<=1.0){
            invp2_IN = 1/(muon.GetMomentum()*muon.GetMomentum());
            invp2_OUT = 1/(muon.GetMomentumPrime()*muon.GetMomentumPrime());
            if(m_pVoxelMean !=6 ) noAddMuon = true;
        }
      // III class
      else if(pclass >1.0){
          invp2_IN = 1./25.;
          invp2_OUT = 1./25.;
          if(m_pVoxelMean !=7 ) noAddMuon = true;
      }
      else
          std::cout << "No CLASS !" << std::endl;
      // to keep all the muons
      if(m_pVoxelMean==4)
          noAddMuon=false;
  }

  //std::cout << "\n\n Mu p_in " << sqrt(1/invp2_IN) << ", p_out " << sqrt(1/invp2_OUT) << ", invp2_IN " << invp2_IN << ", invp2_OUT " << invp2_OUT << std::endl;

  if(m_pVoxelMean){
      // loop ever voxels to find total length in furnace
      for(std::vector<int>::const_iterator it=voxelOrder.begin(); it!=voxelOrder.end(); it++){
          if( (m_imgMC.operator [](*it).Value * (1.e6)) > 0.01){
            const HPoint3f& pt1 = voxelMap[*it][0];
            const HPoint3f& pt2 = voxelMap[*it][1];
            Scalarf L = (pt2-pt1).norm();
            totalLengthFurnace += L;
        }
      }
      //std::cout << "totalLength " << totalLength << ", in furnace " << totalLengthFurnace << std::endl;
      deltaP = (sqrt(1/invp2_OUT) - sqrt(1/invp2_IN))/totalLengthFurnace;
  }

  float sumLijFurnace = 0.;

  //---- Cut muons crossing one voxel only
  if(voxelOrder.size()<2)
    return false;

  //---- Loop over the ordered list of voxels
  for(std::vector<int>::const_iterator it=voxelOrder.begin(); it!=voxelOrder.end(); it++){
    const HPoint3f& pt1 = voxelMap[*it][0];
    const HPoint3f& pt2 = voxelMap[*it][1];

    //---- Now loop over each element in the ray
    Event::Element elc;
    elc.Wij << NAN, NAN;
    elc.lambda = NAN;
    elc.Sij = NAN;
    elc.pw = NAN;
    elc.voxel = NULL;

    //---- Retrieve the voxel from the voxel collection by the ID (== *it)
    if( (*it) >= this->GetVoxCollection()->GetDims().prod()){
      std::cout << "ATTENTION voxel ID > size collection!! " << std::endl;
      return false;
    }
    elc.voxel = &this->GetVoxCollection()->operator [](*it);

    if(elc.voxel->Count != 0 || elc.voxel == NULL || elc.voxel->Value == NAN)
        continue;

    //---- Get length of ray in the voxel
    //---- (NOT == to el.L, due to mid-voxel inflections)
    Scalarf L = (pt2-pt1).norm();

    if(m_pVoxelMean){
        if( (m_imgMC.operator [](*it).Value * (1.e6)) > 0.01)
            sumLijFurnace += L;
    }

    //----> Old calculation of T
    if(m_oldTCalculation) T = fabs(T-L);
    //----> New calculation of T
    else{
      Scalarf h = (muon.LineIn().direction.transpose()*(pt2-front_pt));
      T = H - h/normIn;
    }

    //---- Negative T can arise due to precision
    if(T < 0) T = 0.;

    //---- Fill Wij (algorithm variables) and pw (momentum weight)
    elc.Wij << L, L*L/2. + L*T, L*L/2. + L*T, L*L*L/3. + L*L*T + L*T*T;
    // NB evc.header.InitialSqrP = (p0/p)^2,  if p=5, pw=0.36
    elc.pw = evc.header.InitialSqrP; //DEFAULT

    if(m_pVoxelMean){
        //std::cout << "\n*** Computing p voxel from IN <1/p2> mean : 0.0131459,         OUT <1/p2> mean : 0.197075 *** " << std::endl;
        elc.pw = 1/((deltaP * sumLijFurnace + sqrt(1/invp2_IN))*(deltaP * sumLijFurnace + sqrt(1/invp2_IN)))*  nominal_momentum *  nominal_momentum;
//        std::cout << "Voxel " <<  *it  << ", pw_hand " << elc.pw << std::endl;
//                  << ":   pw_file " <<  voxel_1op2
//                  << " MC voxel value " << imgMC.operator [](*it).Value * (1.e6) << std::endl;
    }
    else if(m_initialSqrPfromVtk){
        float voxel_1op2 = m_initialSqrPfromVtk->operator [](*it).Value * (1.e6) *  nominal_momentum *  nominal_momentum;
        if(voxel_1op2!=0.){
          elc.pw = voxel_1op2;
          //std::cout << "ATTENTION : Replacing 1/p2 " << evc.header.InitialSqrP << " with " << elc.pw << std::endl;
        }
        else
            elc.pw = 0.36;
    }


/*
    //---- Add both views to E if voxel the is "frozen"
    //---- "Frozen" means the voxel has a constant value of LSD
    if(elc.voxel->Value != NAN &&  elc.voxel->Value <= 0){
      //---- Get a 2x2 block in the top left corner of the Error (E) matrix
      evc.header.E.block<2,2>(2,0) += elc.Wij * fabs(elc.voxel->Value) * elc.pw;//evc.header.InitialSqrP;
      //---- Get a 2x2 block in the bottom right corner of the Error (E) matrix
      evc.header.E.block<2,2>(0,2) += elc.Wij * fabs(elc.voxel->Value) * elc.pw;//evc.header.InitialSqrP;
    }
    //---- If the voxel ISN'T frozen, keep the Element in the Event
    else{
      if(elc.voxel != NULL)
    }
*/
    if(!std::isnan(elc.voxel->Value))
      evc.elements.push_back(elc);
  }
  //std::cout << "=== >  Muon in FURNACE  sumLij " << sumLijFurnace << ", totalLenght " << totalLengthFurnace << std::endl;

  //---- Keep the event
  if(!noAddMuon && evc.elements.size()>1 && evc.header.Di[0]!=NAN){
    m_Events.push_back(evc);

    //---- cross check
    if(debug){
        std::cout << "\n\n Add Event to collection\n";
        this->DumpEvent(&evc);
    }
    return true;
  } else{
      if(debug)
          std::cout << "Do not add muon flag.... EXITING!" << std::endl;
      return false;
  }
}


//________________________
///
/// \note MuonCollection is syncronized with Event vector ONLY if SetMuonCollection is called
/// if AddMuons is called by the analyzer it is not!
/// \param muons
///
void IBAnalyzerEM::SetMuonCollection(IBMuonCollection *muons){
  std::cout << "Setting muon collection with collection " << muons << std::endl;
  std::cout << "'Using full path ?' = '" << !m_useRecoPath << "'" << std::endl;
  uLibAssert(muons);

  //---- Clear the event collection
  std::cout << "Clearing all events " << std::endl;
  m_Events.clear();

  //---- Iterate over the muon data, and collec the muon data, and the muon path
  Vector<MuonScatterData>::iterator itr = muons->Data().begin();
  Vector<Vector<HPoint3f> >::iterator path_itr = muons->FullPath().begin();
  std::cout << "Adding " << muons->Data().size() << " muons " << std::endl;

  if(m_pVoxelMean){
      std::cout << "\n*** Computing p voxel from linear function from <1/p2> mean IN to <1/p2> mean OUT *** " << std::endl;
      /// get MC furnace to locate voxel in furnace
      //const char *mcFurnace =  "/home/sara/workspace/experiments/radmu/mublast/analysis/20150522_imageFromMC/mcFurnace_2016-05-03_vox20_250vox.vtk";
      const char *mcFurnace =  "/mnt/mutom-gluster/data/mublast/imageFromMC/mcFurnace_2016-05-03_vox20_250vox.vtk";
      if( !m_imgMC.ImportFromVtk(mcFurnace) ){
          std::cout << "ATTENTION : error opening image from file..." << mcFurnace << std::endl;
          return;
      }
  }
  else if(m_initialSqrPfromVtk)
      std::cout << "\n*** Computing p voxel from file vtk*** " << std::endl;

  int countmu = 0;
  while(itr != muons->Data().end()){
    countmu++;
    if(countmu%100000==0)
    	std::cout << "Adding muon " << countmu << "...." << std::endl;

    //---- If the muon full path has an error, remove the muon
    if(!AddMuonFullPath(*itr, *path_itr)){
      muons->Data().remove_element(*itr);
      if(muons->FullPath().size() > 0) muons->FullPath().remove_element(*path_itr);
    }
    //---- Otherwise, iterate
    else{
      itr++;
      path_itr++;
    }
  }
  //---- Pass the muons to the base class
  std::cout << "\nDone, now calling base class... adding muon collection of " << muons->size() << " muons" << std::endl;
  BaseClass::SetMuonCollection(muons);

//  //--- cross check: dump muon collection
//  std::cout << "Dumping muon collection.... " << std::endl;
//  itr = muons->Data().begin();
//  while(itr != muons->Data().end())
//    std::cout << (*itr) << std::endl;

}

//________________________
unsigned int IBAnalyzerEM::Size(){
    return m_Events.size();
}

//________________________
void IBAnalyzerEM::Run(unsigned int iterations, float muons_ratio){
    // performs iterations //
    for (unsigned int it = 0; it < iterations; it++) {
        fprintf(stderr,"\r[%d muons] EM -> performing iteration %i",
                (int) m_Events.size(), it);
        Evaluate(muons_ratio);          // run single iteration of proback //
        if(!m_UpdateAlgorithm)
            this->GetVoxCollection()->UpdateDensity<UpdateDensitySijCapAlgorithm>(10);                // DEFAULT HARDCODE THRESHOLD
//            this->GetVoxCollection()->UpdateDensity<UpdateDensitySijCapAlgorithm>(2);                // HARDCODE THRESHOLD
        else
            this->m_UpdateAlgorithm->operator()(this->GetVoxCollection(),10);   // DEFAULT HARDCODE THRESHOLD
//            this->m_UpdateAlgorithm->operator()(this->GetVoxCollection(),2);   // HARDCODE THRESHOLD
    }
    printf("\nEM -> done\n");
}

//________________________
void IBAnalyzerEM::SetMLAlgorithm(IBAnalyzerEMAlgorithm *MLAlgorithm){
    m_SijAlgorithm = MLAlgorithm;
}

//________________________
void IBAnalyzerEM::SijGuess(Vector<Vector2f> tpv){
    Evaluate(1);
    // ATTENZIONE!! il vettore deve essere ordinato per threshold crescenti   //
    for (int i=0; i<tpv.size(); ++i) {

        Vector< Event >::iterator itr = this->m_Events.begin();
        int count = 0;
        int nvox_cut=0;
        while (itr != this->m_Events.end()) {
            if(em_test_SijCut(*itr, tpv[i](0), nvox_cut))
            {
                itr->header.InitialSqrP = nominal_momentum / tpv[i](1);
                itr->header.InitialSqrP *= itr->header.InitialSqrP;
                for (unsigned int j = 0; j < itr->elements.size(); ++j)
                    itr->elements[j].pw = itr->header.InitialSqrP;
                count ++;
            }
            itr++;
        }
        std::cout << "Guess class " << tpv[i](0) << ", p=" << tpv[i](1) << "   counted muons: " << count << "\n" << std::endl;
    }

    this->GetVoxCollection()->UpdateDensity<UpdateDensitySijCapAlgorithm>(0);   // HARDCODE THRESHOLD
}

//________________________
void IBAnalyzerEM::SetVoxCollection(IBVoxCollection *voxels){
    if(this->GetMuonCollection()) {
        BaseClass::SetVoxCollection(voxels);
        this->SetMuonCollection(BaseClass::GetMuonCollection());
    }
    else
        std::cerr << "*** Analyzer EM is unable to reset Voxels ***\n" <<
                     "*** without a defined muon collection ... ***\n";
}

//________________________
void IBAnalyzerEM::SetVoxcollectionShift(Vector3f shift){
  if(this->GetMuonCollection()) {
    IBVoxCollection *voxels = this->GetVoxCollection();
    Vector3f pos = voxels->GetPosition();
    voxels->SetPosition(pos + shift);
    IBMuonCollection *muons = this->GetMuonCollection();
    for(int i=0; i<muons->size(); ++i)
        this->AddMuon(muons->At(i));
  }
}


//________________________
////////////////////////////////////////////////////////////////////////////////
/////////////////// DUMP EVENTS ////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
void IBAnalyzerEM::DumpP(const char *filename, float x0, float x1)
{
    static int counter = 0;
    static TFile *file = new TFile(filename,"RECREATE");

    if(!filename) {
        file->Write();
        file->Close();
        delete file;
        return;
    }

    if(file) {
        gDirectory->cd(file->GetPath());
        {
            char name[100];
            sprintf(name,"inv_p_sq_%i",counter++);
            TH1F *h = new TH1F(name,"1/p^2 distribution [1/GeV^2]",1000,x0,x1);
            float p0sq = nominal_momentum * nominal_momentum;
            for(Id_t i=0; i<m_Events.size(); ++i)
                h->Fill(m_Events[i].header.InitialSqrP / p0sq );
            h->Write();
            delete h;
        }

        {
            char name[100];
            sprintf(name,"p_%i",counter++);
            TH1F *h = new TH1F(name,"p distribution [GeV]",1000,x0,x1);
            float p0sq = nominal_momentum * nominal_momentum;
            for(Id_t i=0; i<m_Events.size(); ++i)
                h->Fill( sqrt(p0sq / m_Events[i].header.InitialSqrP) );
            h->Write();
            delete h;
        }

    }
}

//________________________
////////////////////////////////////////////////////////////////////////////////
/// dump events on rootuple
////////////////////////////////////////////////////////////////////////////////
void IBAnalyzerEM::dumpEventsTTree(const char *filename)
{
    Evaluate(1);

    /// open file, tree
    std::cout << "\n*** Dump event collection from IBAnalyzer on file " << filename << std::endl;
    static TFile *file = new TFile(filename,"RECREATE");
    gDirectory->cd(file->GetPath());
    gROOT->ProcessLine("#include<vector>");

    char name[100];
    sprintf(name,"muons");
    TTree *tree = (TTree*)file->Get("muons");
    if(!tree)
        tree = new TTree(name,name);

    int ev = 0;
    float mom, sumLij, dist, DP, DX, DT, DZ;
    float Smedian;
    std::vector<float> *vpw = new std::vector<float>();

    tree->Branch("ev",&ev,"ev/I");
    tree->Branch("p",&mom,"p/F");
    tree->Branch("sumLij",&sumLij,"sumLij/F");
    tree->Branch("dist",&dist,"dist/F");
    tree->Branch("DP",&DP,"DP/F");
    tree->Branch("DX",&DX,"DX/F");
    tree->Branch("DT",&DT,"DT/F");
    tree->Branch("DZ",&DZ,"DZ/F");
    tree->Branch("Smedian",&Smedian,"Smedian/F");
    tree->Branch("pw","vector<float>",&vpw);


    /// event loop
    Vector< Event >::iterator itr = m_Events.begin();
    std::cout << "Reading " << m_Events.size() << " events " << std::endl;
    while (itr != m_Events.end()) {
        Event & evc = *itr;

        /// crossed voxel loop
        sumLij = 0;
        vpw->clear();
        Vector< float > Si;
        Vector< Event::Element >::iterator itre = evc.elements.begin();

        while (itre != evc.elements.end()) {
            Event::Element & elc = *itre;
            sumLij += elc.Wij(0,0);
            float Nij = fabs( (elc.Sij * elc.voxel->Count - elc.voxel->SijCap) / elc.voxel->SijCap );
            Si.push_back(Nij);

            vpw->push_back(nominal_momentum/sqrt(elc.pw));

            ++itre;
        }
        // momentum used in the algorithm
        mom = nominal_momentum/sqrt(evc.header.InitialSqrP);

        DP = evc.header.Di[0];
        DX = evc.header.Di[1];
        DT = evc.header.Di[2];
        DZ = evc.header.Di[3];

        // compute median
        std::sort(Si.begin(),Si.end());
       int nS = Si.size();
       if(nS){
           if(nS%2)
               Smedian = Si[nS / 2];
           else
               Smedian = (Si[nS / 2 - 1] + Si[nS / 2]) / 2;
       }

        tree->Fill();

        ev++;
        itr++;
    }

//    /// muon loop to add poca information
//    IBMuonCollection *muons = this->GetMuonCollection();
//    dist = 0;
//    for(int i=0; i<muons->size(); ++i){
//        MuonScatterData muon = muons->At(i);
//        if(m_PocaAlgorithm){
//            bool use_poca = m_PocaAlgorithm->evaluate(muon);
//            dist = m_PocaAlgorithm->getDistance();
//            //uncomment to exclude distance when PoCA is outside voxel bounds
//            //if(use_poca)
//            bdist->Fill();
//        }
//    }
//    tree->Write("", TObject::kOverwrite);

    tree->Write();
    file->Close();

    return;
}
