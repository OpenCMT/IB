/*//////////////////////////////////////////////////////////////////////////////
// CMT Cosmic Muon Tomography project //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

  Copyright (c) 2014, Universita' degli Studi di Padova, INFN sez. di Padova

  Coordinators: Prof. Gianni Zumerle < gianni.zumerle@pd.infn.it >
                Paolo Checchia       < paolo.checchia@pd.infn.it >

  Authors: Andrea Rigoni Garola < andrea.rigoni@pd.infn.it >
           Matteo Furlan        < nuright@gmail.com >
           Sara Vanini          < sara.vanini@pd.infn.it >
	   Joel Klinger         < klinger@pd.infn.it >

  All rights reserved
  ------------------------------------------------------------------

  This file can not be copied and/or distributed without the express
  permission of  Prof. Gianni Zumerle  < gianni.zumerle@pd.infn.it >

//////////////////////////////////////////////////////////////////////////////*/

#include <stdio.h>
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

class IBAnalyzerEMPimpl;

namespace {
typedef IBAnalyzerEM::Event Event;
  //static DebugTTree trd(__FILE__);
} // namespace

////////////////////////////////////////////////////////////////////////////////
/////  PIMPL  //////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

class IBAnalyzerEMPimpl {

    typedef IBAnalyzerEM::Event Event;

public:
  IBAnalyzerEMPimpl(IBAnalyzerEM *parent, float rankLimit) :
        m_parent(parent),
        m_SijAlgorithm(NULL),
	m_firstIteration(false),
	m_rankLimit(rankLimit){;}
  
  
    void Project(Event *evc);

    void BackProject(Event *evc);

    void Evaluate(float muons_ratio);

    void filterEventsVoxelMask();

    void filterEventsLineDistance(float min, float max);

    void SijCut(float threshold);

    Vector<Event > SijCutCount(float threshold_low, float threshold_high);

    float SijMedian(const Event &evc);

    void dumpEventsSijInfo(const char *filename, Vector<float> N);

    void SijGuess(float threshold, float p);

    void Chi2Cut(float threshold);

    // members //
    IBAnalyzerEM          *m_parent;
    IBAnalyzerEMAlgorithm *m_SijAlgorithm;
    Vector<Event> m_Events;

  bool m_rankLimit;      
  bool m_firstIteration;
};

//________________________
void IBAnalyzerEMPimpl::Project(Event *evc){
    // compute sigma //
    Matrix4f Sigma = Matrix4f::Zero();
    m_SijAlgorithm->ComputeSigma(Sigma, evc);
    // compute sij //
    m_SijAlgorithm->evaluate(Sigma,evc);
}

//________________________
void IBAnalyzerEMPimpl::BackProject(Event *evc){
    IBVoxel *vox;
    // sommatoria della formula 38 //
    for (unsigned int j = 0; j < evc->elements.size(); ++j) {
      vox = evc->elements[j].voxel;
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
void IBAnalyzerEMPimpl::Evaluate(float muons_ratio)
{    
    unsigned int start = 0;
    unsigned int end = (unsigned int) (m_Events.size() * muons_ratio);

    if(m_SijAlgorithm) {
      // Projection
      #pragma omp parallel for
      for (unsigned int i = start; i < end; ++i)
	this->Project(&m_Events[i]);
      #pragma omp barrier
      
      // Backprojection
      #pragma omp parallel for
      for (unsigned int i = start; i < end; ++i)
	this->BackProject(&m_Events[i]);
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
void IBAnalyzerEMPimpl::filterEventsVoxelMask()
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
            if(this->m_parent->m_MuonCollection)
                this->m_parent->m_MuonCollection->Data().remove_element(pos);
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
void IBAnalyzerEMPimpl::filterEventsLineDistance(float min, float max)
{
  std::cout << "\n*** Removing events with line distance out of range from " << this->m_Events.size() << " muon collection." << std::endl;

    Vector< Event >::iterator itr = this->m_Events.begin();
    const Vector< Event >::iterator begin = this->m_Events.begin();

    while (itr != this->m_Events.end()) {
        Event & evc = *itr;
        unsigned int pos = itr - begin;

        MuonScatterData muon = this->m_parent->m_MuonCollection->At(pos);
        bool use_poca = this->m_parent->m_PocaAlgorithm->evaluate(muon);
        float dist = this->m_parent->m_PocaAlgorithm->getDistance();

        /// erase event and muon with distance out of range
        if(!isFinite(dist) || dist >= max || dist < min){
            this->m_Events.remove_element(evc);
            this->m_parent->m_MuonCollection->Data().remove_element(pos);
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
float IBAnalyzerEMPimpl::SijMedian(const Event &evc){
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
//       // debug
//           for(int i=0; i<nS; i++) std::cout << Si[i] << ",";
//           std::cout << "\n     MEDIAN =" << Smedian << std::endl;

    return Smedian;
}

//________________________
Vector<Event > IBAnalyzerEMPimpl::SijCutCount(float threshold_low, float threshold_high)
{
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
void IBAnalyzerEMPimpl::dumpEventsSijInfo(const char *name, Vector<float> N)
{
/// dump event Sij info on file
    std::fstream fout;
    fout.open(name, std::fstream::out | std::fstream::app);

    Vector< Event > ve;
    Vector< Event >::iterator itr = this->m_Events.begin();
    const Vector< Event >::iterator begin = this->m_Events.begin();
    int nev = 0;

    /// loop over events
    while (itr != this->m_Events.end()) {
        //        std::cout << "\n\n *** Event " << evnum << std::endl;
        float p0sq = 3. * 3.;
        //float mom = sqrt(p0sq/(*itr).header.InitialSqrP);
        float mom =  (*itr).header.pTrue;
        fout << nev << " " << mom << " ";
        // fout << (*itr).elements.size() << " ";

        // just dump median
        float median =  SijMedian(*itr);
        fout << median << " ";

//        /// loop over Sij tresholds
//        Vector< float >::iterator itrN = N.begin();
//        const Vector< float >::iterator beginN = N.begin();
//        while (itrN != N.end()) {
//            int nvox_cut = 0;

////            // dump number of voxels aboce the threshold
////            em_test_SijCut(*itr, *itrN, nvox_cut);
////            fout << nvox_cut << " ";

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
void IBAnalyzerEMPimpl::SijCut(float threshold){
    Vector< Event >::iterator itr = this->m_Events.begin();
    const Vector< Event >::iterator begin = this->m_Events.begin();

    int count = 0;
    int nvox_cut=0;
    while (itr != this->m_Events.end()) {
        if(em_test_SijCut(*itr, threshold, nvox_cut))
        {
            unsigned int pos = itr - begin;
            this->m_Events.remove_element(*itr);
            if(this->m_parent->m_MuonCollection)
                this->m_parent->m_MuonCollection->Data().remove_element(pos);
            count ++;
        }
        else ++itr;
    }
    std::cout << "SijCut removed muons: " << count << "\n" << std::endl;
}

//________________________
void IBAnalyzerEMPimpl::SijGuess(float threshold, float p){
    Vector< Event >::iterator itr = this->m_Events.begin();
    int count = 0;
    int nvox_cut=0;
    while (itr != this->m_Events.end()) {
        if(em_test_SijCut(*itr, threshold, nvox_cut))
        {
            itr->header.InitialSqrP = m_parent->$$.nominal_momentum / p;
            itr->header.InitialSqrP *= itr->header.InitialSqrP;
            for (unsigned int j = 0; j < itr->elements.size(); ++j)
                itr->elements[j].pw = itr->header.InitialSqrP;
            count ++;
        }
        itr++;
    }
    std::cout << "Guess class " << threshold << ", p=" << p << "   counted muons: " << count << "\n" << std::endl;
}

//________________________
////////////////////////////////////////////////////////////////////////////////
void IBAnalyzerEMPimpl::Chi2Cut(float threshold)
{
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
            if(this->m_parent->m_MuonCollection)
                this->m_parent->m_MuonCollection->Data().remove_element(pos);
        }
        else ++itr;
    } while (itr != this->m_Events.end());
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
			   bool oldTCalculation, float rankLimit, IBVoxCollection* initialSqrPfromVtk) :
    m_PocaAlgorithm(NULL),
    m_VarAlgorithm(NULL),
    m_RayAlgorithm(NULL),
    m_UpdateAlgorithm(NULL),
    m_nPath(nPath),
    m_alpha(alpha),
    m_useRecoPath(useRecoPath),
    m_oldTCalculation(oldTCalculation),
    m_rankLimit(rankLimit),
    m_initialSqrPfromVtk(initialSqrPfromVtk)
{
  //---- Print the settings
  std::cout << "Using alpha = " << m_alpha << ", #path = " << m_nPath << std::endl;
  std::cout << "Reco path ("    << m_useRecoPath      << "), "
	    << "Old T ("        << m_oldTCalculation  << "), " << std::endl;
  BaseClass::SetVoxCollection(&voxels);
  init_properties(); // < DANGER !!! should be moved away !!
  m_d = new IBAnalyzerEMPimpl(this, m_rankLimit);
}

//___________________________
IBAnalyzerEM::~IBAnalyzerEM(){
    delete m_d;
}

//___________________________
Vector<IBAnalyzerEM::Event> &IBAnalyzerEM::Events(){
    return m_d->m_Events;
}

//___________________________
//---- **Deprecated** version of IBAnalyzerEM::AddMuon, please use IBAnalyzerEM::AddMuonFullPath
//---- Only has 2-path mode, and bugs in the calculation of L and T
bool IBAnalyzerEM::AddMuon(const MuonScatterData &muon){
  if(unlikely(!m_RayAlgorithm || !m_VarAlgorithm)) return false;
  Event evc;
  
  evc.header.InitialSqrP = pow($$.nominal_momentum/muon.GetMomentum() ,2);
  if(isnan(evc.header.InitialSqrP)) std::cout << "sono in AddMuon: nominalp:" << $$.nominal_momentum << " muon.GetMomentum():" << muon.GetMomentum() <<"\n"
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
  m_d->m_Events.push_back(evc);
  
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

  //-------------------------
  //---- STEP #1: Fill the event info
  
  //---- Check the ray algo (calculates the ray parameters)
  //---- and the variable algo (calculates the scattering/displacement variables)
  if(unlikely(!m_RayAlgorithm || !m_VarAlgorithm)) return false;

  Event evc; //<---- The event info
  if(likely(m_VarAlgorithm->evaluate(muon))) {
    //---- Get the Data (Di) and Error (E) matrices
    evc.header.Di = m_VarAlgorithm->getDataVector();
    evc.header.E  = m_VarAlgorithm->getCovarianceMatrix();
    //---- Momentum square (the "$$" notation is a bit much...)
    evc.header.InitialSqrP = pow($$.nominal_momentum/muon.GetMomentum() ,2);
    // SV for Sij studies
    evc.header.pTrue = muon.GetMomentumPrime();

    if(isnan(evc.header.InitialSqrP)){
      std::cout << "AddMuonFullPath: nominalp:" << $$.nominal_momentum
		<< "muon.GetMomentum():" << muon.GetMomentum() <<"\n" << std::endl;
    }
  }
  else return false;

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
    if( !m_RayAlgorithm->GetEntryPoint(muon.LineIn(),entry_pt) ) return false;
    if( !m_RayAlgorithm->GetExitPoint(muon.LineOut(),exit_pt) )  return false;
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
  	  if(entry_length > trackLength || exit_length > trackLength) return false;
	  
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
    if(muonPath.size()==0) return false;
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
  
  //---- Loop over the ordered list of voxels
  for(std::vector<int>::const_iterator it=voxelOrder.begin(); it!=voxelOrder.end(); it++){
    const HPoint3f& pt1 = voxelMap[*it][0];
    const HPoint3f& pt2 = voxelMap[*it][1];    

    //---- Now loop over each element in the ray
    Event::Element elc;
  
    //---- Retrieve the voxel from the voxel collection by the ID (== *it)
    elc.voxel = &this->GetVoxCollection()->operator [](*it);
  
    //---- Get length of ray in the voxel
    //---- (NOT == to el.L, due to mid-voxel inflections)
    Scalarf L = (pt2-pt1).norm();
    
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
    elc.pw = evc.header.InitialSqrP; //DEFAULT

    //---- Add both views to E if voxel the is "frozen"
    //---- "Frozen" means the voxel has a constant value of LSD
    if(elc.voxel->Value <= 0){
      //---- Get a 2x2 block in the top left corner of the Error (E) matrix
      evc.header.E.block<2,2>(2,0) += elc.Wij * fabs(elc.voxel->Value) * evc.header.InitialSqrP;
      //---- Get a 2x2 block in the bottom right corner of the Error (E) matrix
      evc.header.E.block<2,2>(0,2) += elc.Wij * fabs(elc.voxel->Value) * evc.header.InitialSqrP;
    }
    //---- If the voxel ISN'T frozen, keep the Element in the Event
    else{
      evc.elements.push_back(elc);
    }
  }

  //---- Keep the event
  m_d->m_Events.push_back(evc);
  return true;
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
  m_d->m_Events.clear();

  //---- Iterate over the muon data, and collec the muon data, and the muon path
  Vector<MuonScatterData>::iterator itr = muons->Data().begin();
  Vector<Vector<HPoint3f> >::iterator path_itr = muons->FullPath().begin();    
  std::cout << "Adding " << muons->Data().size() << " muons " << std::endl;    
  while(itr != muons->Data().end()){
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
  std::cout << "\nDone, now calling base class..." << std::endl;
  BaseClass::SetMuonCollection(muons);
}

//________________________
unsigned int IBAnalyzerEM::Size(){
    return m_d->m_Events.size();
}

//________________________
void IBAnalyzerEM::Run(unsigned int iterations, float muons_ratio){
    // performs iterations //
    for (unsigned int it = 0; it < iterations; it++) {
        fprintf(stderr,"\r[%d muons] EM -> performing iteration %i",
                (int) m_d->m_Events.size(), it);
        m_d->Evaluate(muons_ratio);          // run single iteration of proback //
        if(!m_UpdateAlgorithm)
            this->GetVoxCollection()->
                UpdateDensity<UpdateDensitySijCapAlgorithm>(10);                // HARDCODE THRESHOLD
        else
            this->m_UpdateAlgorithm->operator()(this->GetVoxCollection(),10);   // HARDCODE THRESHOLD
    }
    printf("\nEM -> done\n");
}

//________________________
void IBAnalyzerEM::SetMLAlgorithm(IBAnalyzerEMAlgorithm *MLAlgorithm){
    m_d->m_SijAlgorithm = MLAlgorithm;
}

//________________________
void IBAnalyzerEM::filterEventsVoxelMask() {
    m_d->filterEventsVoxelMask();
}

//________________________
void IBAnalyzerEM::filterEventsLineDistance(float min, float max) {
    m_d->filterEventsLineDistance(min, max);
}

//________________________
Vector<Event > IBAnalyzerEM::SijCutCount(float threshold_low, float threshold_high) {
    m_d->Evaluate(1);
    return m_d->SijCutCount(threshold_low,threshold_high);
}

//________________________
void IBAnalyzerEM::dumpEventsSijInfo(const char *name, Vector<float> N) {
    m_d->Evaluate(1);
    m_d->dumpEventsSijInfo(name, N);
    return;
}


//________________________
void IBAnalyzerEM::SijCut(float threshold) {
    m_d->Evaluate(1);
    m_d->SijCut(threshold);
    this->GetVoxCollection()->UpdateDensity<UpdateDensitySijCapAlgorithm>(0);   // HARDCODE THRESHOLD
}

//________________________
float IBAnalyzerEM::SijMedian(const Event &evc) {
    //m_d->Evaluate(1);
    m_d->SijMedian(evc);
}

//________________________
void IBAnalyzerEM::SijGuess(Vector<Vector2f> tpv){
    m_d->Evaluate(1);
    // ATTENZIONE!! il vettore deve essere ordinato per threshold crescenti   //
    for (int i=0; i<tpv.size(); ++i)
        m_d->SijGuess( tpv[i](0), tpv[i](1) );
    this->GetVoxCollection()->UpdateDensity<UpdateDensitySijCapAlgorithm>(0);   // HARDCODE THRESHOLD
}

//________________________
void IBAnalyzerEM::Chi2Cut(float threshold){
    m_d->Evaluate(1);
    this->GetVoxCollection()->UpdateDensity<UpdateDensitySijCapAlgorithm>(0);   // HARDCODE THRESHOLD
    m_d->Chi2Cut(threshold);
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
            float p0sq = $$.nominal_momentum * $$.nominal_momentum;
            for(Id_t i=0; i<m_d->m_Events.size(); ++i)
                h->Fill(m_d->m_Events[i].header.InitialSqrP / p0sq );
            h->Write();
            delete h;
        }

        {
            char name[100];
            sprintf(name,"p_%i",counter++);
            TH1F *h = new TH1F(name,"p distribution [GeV]",1000,x0,x1);
            float p0sq = $$.nominal_momentum * $$.nominal_momentum;
            for(Id_t i=0; i<m_d->m_Events.size(); ++i)
                h->Fill( sqrt(p0sq / m_d->m_Events[i].header.InitialSqrP) );
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
    m_d->Evaluate(1);

    /// open file, tree
    std::cout << "\n*** Dump event collection from IBAnalyzer on file " << filename << std::endl;
    static TFile *file = new TFile(filename,"update");
    gDirectory->cd(file->GetPath());

    char name[100];
    sprintf(name,"muons");
    TTree *tree = (TTree*)file->Get("muons");
    if(!tree)
        tree = new TTree(name,name);

    int ev = 0;
    float mom, sumLij, dist, DP, DX, DT, DZ;
    float Smedian;
    TBranch *bev = tree->Branch("ev",&ev,"ev/I");
    TBranch *bp = tree->Branch("p",&mom,"p/F");
    TBranch *bsumLij = tree->Branch("sumLij",&sumLij,"sumLij/F");
    TBranch *bdist = tree->Branch("dist",&dist,"dist/F");
    TBranch *bDP = tree->Branch("DP",&DP,"DP/F");
    TBranch *bDX = tree->Branch("DX",&DX,"DX/F");
    TBranch *bDT = tree->Branch("DT",&DT,"DT/F");
    TBranch *bDZ = tree->Branch("DZ",&DZ,"DZ/F");
    TBranch *bSmedian = tree->Branch("Smedian",&Smedian,"Smedian/F");

    /// event loop
    Vector< Event >::iterator itr = m_d->m_Events.begin();
    std::cout << "Reading " << m_d->m_Events.size() << " events " << std::endl;
    while (itr != m_d->m_Events.end()) {
        Event & evc = *itr;

        /// crossed voxel loop
        sumLij = 0;
        Vector< float > Si;
        Vector< Event::Element >::iterator itre = evc.elements.begin();

        while (itre != evc.elements.end()) {
            Event::Element & elc = *itre;
            sumLij += elc.Wij(0,0);
            float Nij = fabs( (elc.Sij * elc.voxel->Count - elc.voxel->SijCap) / elc.voxel->SijCap );
            Si.push_back(Nij);
            ++itre;
        }
        // momentum used in the algorithm
        mom = $$.nominal_momentum/sqrt(evc.header.InitialSqrP);

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

//           // debug
//           std::cout << "Event pTrue " << evc.header.pTrue << std::endl;
//           for(int i=0; i<nS; i++) std::cout << Si[i] << ",";
//           std::cout << "\n     MEDIAN =" << Smedian << std::endl;
           bSmedian->Fill();
       }

        bev->Fill();
        bp->Fill();
        bsumLij->Fill();

        bDP->Fill();
        bDX->Fill();
        bDT->Fill();
        bDZ->Fill();

        ev++;
        itr++;
    }

    /// muon loop to add poca information
    IBMuonCollection *muons = this->GetMuonCollection();
    dist = 0;

    for(int i=0; i<muons->size(); ++i){
        MuonScatterData muon = muons->At(i);
        if(m_PocaAlgorithm){
            bool use_poca = m_PocaAlgorithm->evaluate(muon);
            dist = m_PocaAlgorithm->getDistance();
            //uncomment to exclude distance when PoCA is outside voxel bounds
            //if(use_poca)
            bdist->Fill();
        }
    }

    // testing
    int sizeev = m_d->m_Events.size();
    int sizemu = muons->size();

    tree->Write("", TObject::kOverwrite);
    delete tree;

    file->Write();
    file->Close();
    delete file;

    return;
}
