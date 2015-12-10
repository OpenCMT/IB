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

#include <TTree.h>
#include <TFile.h>
#include <TH1F.h>

#include <Core/Vector.h>
#include <Core/Debug.h>

#include "IBPocaEvaluator.h"
#include "IBMinimizationVariablesEvaluator.h"
#include "IBVoxRaytracer.h"

#include "IBVoxCollectionCap.h"
#include "IBAnalyzerEM.h"

#include "IBAnalyzerEMAlgorithm.h"
#include "IBAnalyzerEMAlgorithmSGA.h"

#include "MuonProjection.hh"
#include "AlphaCalculator.hh"

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
    IBAnalyzerEMPimpl(IBAnalyzerEM *parent) :
        m_parent(parent),
        m_SijAlgorithm(NULL){;}

    void Project(Event *evc);

    void BackProject(Event *evc);

    void Evaluate(float muons_ratio);

    void filterEventsVoxelMask();

    void filterEventsLineDistance(float min, float max);

    void SijCut(float threshold);

    void SijGuess(float threshold, float p);

    void Chi2Cut(float threshold);

    // members //
    IBAnalyzerEM          *m_parent;
    IBAnalyzerEMAlgorithm *m_SijAlgorithm;
    Vector<Event> m_Events;

};


void IBAnalyzerEMPimpl::Project(Event *evc)
{
    // compute sigma //
    Matrix4f Sigma = Matrix4f::Zero();
    m_SijAlgorithm->ComputeSigma(Sigma, evc);
    // compute sij //
    m_SijAlgorithm->evaluate(Sigma,evc);
}

void IBAnalyzerEMPimpl::BackProject(Event *evc)
{
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

    // //---- Step 1: Map voxel ('elements[j]') to vector<sij + lambda??>
    // std::map<int,std::vector<float> > voxelMap;
    // for (unsigned int i = start; i < end; ++i){
    //   for (unsigned int j = 0; j < evc->elements.size(); ++j) {
    // 	int   vox = m_Events[i]->GETVOXELNUMBER;
    // 	float sij = m_Events[i]->elements[j].Sij;
    // 	if(voxelMap.find(vox)==voxelMap.end()) voxelMap[vox] = std::vector<float>();
    // 	voxelMap[vox].push_back(sij);
    //   }
    // }
    
    // //---- Step 2: Map voxel to <Sij>, sd(Sij), med(Sij), q68(Sij)/2
    // std::map<int,std::vector<float> > voxelMap2;
    // for(std::map<int,std::vector<float> >::iterator it=voxelMap.begin(); it!=voxelMap.end(); it++){
    //   std::vector<float> v = it->second;
    //   std::sort(v.begin(), v.end());
    //   float mean   = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
    //   float median = v.at(int(0.5*float(v.size())));
    //   float q68    = (v.at(int(0.84*float(v.size()))) - v.at(int(0.16*float(v.size())))) / 2.;
      
    //   std::vector<float> diff(v.size());
    //   std::transform(v.begin(), v.end(), diff.begin(),
    // 		     std::bind2nd(std::minus<float>(), mean));
    //   float sq_sum = std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0);
    //   float stdev = std::sqrt(sq_sum / v.size());
      
    //   voxelMap2[it->first] = std::vector<float>();
    //   voxelMap2[it->first].push_back(mean);
    //   voxelMap2[it->first].push_back(stdev);
    //   voxelMap2[it->first].push_back(median);
    //   voxelMap2[it->first].push_back(q68);
    // }
    
    // //---- Step 3: Delete map 1
    // voxelMap2.clear();

    // //---- Step 3.5 Make TFile and TTree
    // TFile* tfile = new TFile("joelout.root","RECREATE");
    // TTree* = make new tree;
    // add branches to tree;    
    
    // //---- Step 4: For each muon "m_Events[i]" make a pair<vector< |sij - av|/sig,> >
    // for (unsigned int i = start; i < end; ++i){
    //   std::vector<float> v1, v2;
    //   for (unsigned int j = 0; j < evc->elements.size(); ++j) {
    // 	int   vox = m_Events[i]->GETVOXELNUMBER;
    // 	float sij = m_Events[i]->elements[j].Sij;

    // 	std::vector<float> v = voxelMap2[vox];
    // 	v1.push_back(fabs(v[0]-sij)/v[1]);
    // 	v2.push_back(fabs(v[2]-sij)/v[3]);
    //   }

    //   //---- Step 5: write vectors and 1/psq to file
    //   ttree->Fill();
    // }
    // ttree->Write();
    // tfile->Close();
      
}
    


    




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

////////////////////////////////////////////////////////////////////////////////
/// SijCut RECIPE1:  (true if Sij cut proposed) //
static bool em_test_SijCut(const Event &evc, float cut_level)
{
    int n_cuts = 0;
    for (unsigned int i = 0; i < evc.elements.size(); i++) {
        const Event::Element &el = evc.elements[i];
        if (fabs( (el.Sij * el.voxel->Count - el.voxel->SijCap)
                  / el.voxel->SijCap ) > cut_level) n_cuts++;
    }
    if (n_cuts > (int)(evc.elements.size()/3) ) return true;
    else return false;
}

void IBAnalyzerEMPimpl::SijCut(float threshold)
{
    Vector< Event >::iterator itr = this->m_Events.begin();
    const Vector< Event >::iterator begin = this->m_Events.begin();

    int count = 0;
    while (itr != this->m_Events.end()) {
        if(em_test_SijCut(*itr, threshold))
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

void IBAnalyzerEMPimpl::SijGuess(float threshold, float p)
{
    Vector< Event >::iterator itr = this->m_Events.begin();
    int count = 0;
    while (itr != this->m_Events.end()) {
        if(em_test_SijCut(*itr, threshold))
        {
            itr->header.InitialSqrP = m_parent->$$.nominal_momentum / p;
            itr->header.InitialSqrP *= itr->header.InitialSqrP;
            for (unsigned int j = 0; j < itr->elements.size(); ++j) {
                itr->elements[j].pw = itr->header.InitialSqrP;
            }
            count ++;
        }
        itr++;
    }
    std::cout << "Guess class " << threshold << ", p=" << p << "   counted muons: " << count << "\n" << std::endl;
}


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

IBAnalyzerEM::IBAnalyzerEM(IBVoxCollection &voxels, int nPath, double alpha, bool useRecoPath,
			   bool oldTCalculation, bool scatterOnly, bool displacementOnly,
			   std::string projectFile, std::string alphaFile) :
    m_PocaAlgorithm(NULL),
    m_VarAlgorithm(NULL),
    m_RayAlgorithm(NULL),
    m_UpdateAlgorithm(NULL),
    m_nPath(nPath),
    m_alpha(alpha),
    m_useRecoPath(useRecoPath),
    m_oldTCalculation(oldTCalculation),
    m_scatterOnly(scatterOnly),
    m_displacementOnly(displacementOnly),
    m_project(false)
{
  //----
  if(projectFile!="" && alphaFile!=""){
    m_project = true;
    m_projector = MuonProjection(projectFile);
    m_alphaCalc = AlphaCalculator(alphaFile);
  }  

  std::cout << "Using alpha = " << m_alpha << ", #path = " << m_nPath << std::endl;
  std::cout << "Reco path ("    << m_useRecoPath      << "), "
	    << "Old T ("        << m_oldTCalculation  << "), "
	    << "Scatter only (" << m_scatterOnly      << "), "
    	    << "Displacement only (" << m_displacementOnly << "), "
    	    << "project ("      << m_project          << ")" << std::endl;
  BaseClass::SetVoxCollection(&voxels);
  init_properties(); // < DANGER !!! should be moved away !!
  m_d = new IBAnalyzerEMPimpl(this);
}

IBAnalyzerEM::~IBAnalyzerEM()
{
    delete m_d;
}

Vector<IBAnalyzerEM::Event> &IBAnalyzerEM::Events()
{
    return m_d->m_Events;
}


//___________________________
bool IBAnalyzerEM::AddMuonFullPath(const MuonScatterData &muon, Vector<HPoint3f>& muonPath){
  
  if(unlikely(!m_RayAlgorithm || !m_VarAlgorithm)) return false;
  Event evc;
  
  evc.header.InitialSqrP = pow($$.nominal_momentum/muon.GetMomentum() ,2);
  if(isnan(evc.header.InitialSqrP)){
    std::cout << "sono in AddMuon: nominalp:" << $$.nominal_momentum
	      << " muon.GetMomentum():" << muon.GetMomentum() <<"\n" << std::endl;
  }

  if(likely(m_VarAlgorithm->evaluate(muon))) {
    evc.header.Di = m_VarAlgorithm->getDataVector();
    evc.header.E  = m_VarAlgorithm->getCovarianceMatrix();
  }     
  else return false;

  //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  //---- Perform raytracing 
  //    IBVoxRaytracer::RayData ray;
  Vector<HPoint3f> pts;
  HPoint3f front_pt, back_pt;
  
  //---- If reconstructing the muon's path
  if(m_useRecoPath){
    
    //---- Require entry and exit points
    HPoint3f entry_pt, exit_pt;
    if( !m_RayAlgorithm->GetEntryPoint(muon.LineIn(),entry_pt) ) return false;
    if( !m_RayAlgorithm->GetExitPoint(muon.LineOut(),exit_pt) )  return false;
    front_pt = entry_pt;
    back_pt = exit_pt;      
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
	  //else std::cout << "Invalid poca" << std::endl;
	}
  	//---- If using the three-line path
  	else if(m_nPath == 3){
  	  //---- Get the poca on the entry/exit tracks
  	  HPoint3f entry_poca = m_PocaAlgorithm->getInTrackPoca();
  	  HPoint3f exit_poca  = m_PocaAlgorithm->getOutTrackPoca();
	  
  	  //---- Get the distance down the tracks to the inflection points
  	  double entry_length = (entry_pt - entry_poca).norm();
  	  double exit_length  = (exit_pt  - exit_poca).norm();
  	  if(entry_length > trackLength || exit_length > trackLength) return false;
	  
  	  double alpha = m_alpha;
  	  if(m_project) alpha = m_alphaCalc.Alpha(trackLength);
	  
  	  //---- Get the inflection points
  	  double normIn  = muon.LineIn().direction.norm(); //<----Recently added
  	  double normOut = muon.LineOut().direction.norm();
  	  HVector3f point1 = entry_pt + alpha*(entry_length/normIn)*muon.LineIn().direction;
  	  HVector3f point2 = exit_pt  - alpha*(exit_length/normOut)*muon.LineOut().direction;
	  
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
    if(muonPath.size()==0) return false;
    for(int i=0; i<muonPath.size(); ++i){
      pts.push_back(muonPath[i]);
    }
    if(pts.size()==0) return false;
    back_pt  = pts.back();
    front_pt = pts.front();
  }
  
  //std::cout << "Got " << pts.size() << " points " << std::endl;
  
  //---- Now loop over the points and rays in order to remove mid-voxel inflections
  std::map<int,std::vector<HPoint3f> > voxelMap;
  std::vector<int> voxelOrder; //<---- An ordered list of the voxels
  Scalarf totalLength = 0.;
  for(int i=0; i<pts.size()-1; ++i){
    HPoint3f & pt1 = pts[i];
    HPoint3f & pt2 = pts[i+1];
  
    //----
    Scalarf  rayLength = (pt2-pt1).norm();
    HPoint3f rayDir    = (pt2-pt1)/rayLength;
    IBVoxRaytracer::RayData ray = m_RayAlgorithm->TraceBetweenPoints(pt1,pt2);
  
    //----
    float cumulativeLength = 0.;
    foreach(const IBVoxRaytracer::RayData::Element &el, ray.Data()){
  
  	HPoint3f pti = pt1 + cumulativeLength*rayDir;
  	cumulativeLength += el.L;
  	HPoint3f ptj = pt1 + cumulativeLength*rayDir;
  
  	//---- If this is a new voxel for this muon
  	if(voxelMap.find(el.vox_id) == voxelMap.end()){
  	  voxelMap[el.vox_id] = std::vector<HPoint3f>();
  	  voxelOrder.push_back(el.vox_id);
  	  voxelMap[el.vox_id].push_back(pti);
  	  voxelMap[el.vox_id].push_back(ptj);
  	}
  	voxelMap[el.vox_id][1] = ptj;
    }

    //---- Increment the total length of the muon path
    totalLength += rayLength;
  }
  
  //---- Now reorder the points
  // pts.clear();
  // for(std::vector<int>::const_iterator it=voxelOrder.begin(); it!=voxelOrder.end(); it++){
  //   std::vector<HPoint3f>& voxelPts = voxelMap[*it];
  //   pts.push_back(voxelPts.front());
  //   pts.push_back(voxelPts.back());
  // }
  
  //---- Now trace between points and calculate length parameters
  m_oldTCalculation = true;
  Scalarf H = muon.LineIn().direction.transpose() * (back_pt - front_pt);
  Scalarf normIn = muon.LineIn().direction.norm();
  if(!m_oldTCalculation) H = H/normIn;  

  //----
  // m_oldTCalculation = true;
  // if(m_oldTCalculation){
  //   for(int i=0; i<pts.size()-1; ++i){
  //     HPoint3f & pt1 = pts[i];
  //     HPoint3f & pt2 = pts[i+1];
  //     IBVoxRaytracer::RayData ray = m_RayAlgorithm->TraceBetweenPoints(pt1,pt2);
  //     Scalarf T = ray.TotalLength();
  //     if(i < pts.size()-2) {
  // 	Scalarf h = muon.LineIn().direction.transpose()*(pt2-pt1);
  // 	T = T * H/h;
  //     }
  //     //---- Now loop over each element in the ray
  //     foreach(const IBVoxRaytracer::RayData::Element &el, ray.Data()){
  // 	Event::Element elc;
	
  // 	//---- Retrieve the voxel
  // 	elc.voxel = &this->GetVoxCollection()->operator [](el.vox_id);
  // 	Scalarf L = el.L;
  // 	T = fabs(T-L);
  // 	elc.Wij << L, L*L/2. + L*T, L*L/2. + L*T, L*L*L/3. + L*L*T + L*T*T;
	
  // 	//---- pw
  // 	elc.pw = evc.header.InitialSqrP;
	
  // 	//---- add both views to E if voxel is frozen 
  // 	if(elc.voxel->Value <= 0){
  // 	  evc.header.E.block<2,2>(2,0) += elc.Wij * fabs(elc.voxel->Value) * evc.header.InitialSqrP;
  // 	  evc.header.E.block<2,2>(0,2) += elc.Wij * fabs(elc.voxel->Value) * evc.header.InitialSqrP;
  // 	}
  // 	else evc.elements.push_back(elc);
  //     }
  //   }
  // }
  // else{

  Scalarf T = totalLength;
  for(std::vector<int>::const_iterator it=voxelOrder.begin(); it!=voxelOrder.end(); it++){
    const HPoint3f & pt1 = voxelMap[*it][0];
    const HPoint3f & pt2 = voxelMap[*it][1];    
    
  // for(std::map<int,std::vector<HPoint3f> >::const_iterator ptIt=voxelMap.begin(); ptIt!=voxelMap.end(); ptIt++){
  
  //   const HPoint3f & pt1 = ptIt->second[0];
  //   const HPoint3f & pt2 = ptIt->second[1];
  
    //----> New calculation of T
    //      Scalarf  rayLength = (pt2-pt1).norm();
    //      HPoint3f rayDir    = (pt2-pt1)/rayLength;
    //<---- 
    
    //----> 2nd Old calculation of T
    // Scalarf T = ray.TotalLength(); //<---- This maybe should be sum(ray)
    // if(m_oldTCalculation){
    //   if(i < pts.size()-2) {
    // 	Scalarf h = muon.LineIn().direction.transpose()*(pt2-pt1);
    // 	T = T * H/h;
    //   }
    // }
    
    
    //---- Now loop over each element in the ray
    //Scalarf cumulativeLength = 0.;
    //      foreach(const IBVoxRaytracer::RayData::Element &el, ray.Data()){
    Event::Element elc;
  
    //---- Retrieve the voxel
    //	elc.voxel = &this->GetVoxCollection()->operator [](el.vox_id);
    //    elc.voxel = &this->GetVoxCollection()->operator [](ptIt->first);
    elc.voxel = &this->GetVoxCollection()->operator [](*it);
  
    // //----> New ray tracing
    // if(voxelMap.find(el.vox_id) == voxelMap.end()){
    //   voxelMap[el.vox_id] = std::make_pair(0, std::make_pair(0.,0.));
    // }
    // voxelMap[el.vox_id].first = voxelMap[el.vox_id].first + 1;
    // //<----
    
    //---- Get length of ray in the voxel
    //Scalarf L = el.L;
    Scalarf L = (pt2-pt1).norm();
  
    // //----> New ray tracing
    // L += voxelMap[el.vox_id].second.first;
    // voxelMap[el.vox_id].second.first = L;
    // //<---- 
    
    //----> Old calculation of T
    if(m_oldTCalculation) T = fabs(T-L);
    //----> New calculation of T
    else{	  
    //cumulativeLength += el.L;
    //HPoint3f pti = (cumulativeLength*rayDir + pt1);
    //Scalarf h = (muon.LineIn().direction.transpose()*(pti-front_pt));
      Scalarf h = (muon.LineIn().direction.transpose()*(pt2-front_pt));
      T = H - h/normIn;
    }
    if(T < 0) T = 0.;
    //----> New ray tracing
    // T += voxelMap[el.vox_id].second.second;
    // voxelMap[el.vox_id].second.second = T;
    // T = (voxelMap[el.vox_id].second.second) / float(voxelMap[el.vox_id].first);
    //----
    //}
    //<----

    //---- Fill Wij
    if(m_scatterOnly)           elc.Wij << L, 0, 0, 0;
    else if(m_displacementOnly) elc.Wij << 0, 0, 0, L*L*L/3. + L*L*T + L*T*T;
    else                        elc.Wij << L, L*L/2. + L*T, L*L/2. + L*T, L*L*L/3. + L*L*T + L*T*T;
    
    //---- pw
    elc.pw = evc.header.InitialSqrP;
    
    //---- add both views to E if voxel is frozen 
    if(elc.voxel->Value <= 0){
      evc.header.E.block<2,2>(2,0) += elc.Wij * fabs(elc.voxel->Value) * evc.header.InitialSqrP;
      evc.header.E.block<2,2>(0,2) += elc.Wij * fabs(elc.voxel->Value) * evc.header.InitialSqrP;
    }
    else{
      //----> New ray tracing
      //	  elementMap[el.vox_id] = elc;
      //----> Old ray tracing
      evc.elements.push_back(elc);
      //<----
    }
  }
  //}
  
  //----> New ray tracing
  // for(std::map<int,Event::Element>::iterator it=elementMap.begin(); it!=elementMap.end(); it++){
  //   evc.elements.push_back(it->second);
  // }
  //<----
     m_d->m_Events.push_back(evc);
  return true;
  //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
}

///
/// \note MuonCollection is syncronized with Event vector ONLY if SetMuonCollection is called
/// if AddMuons is called by the analyzer it is not!
/// \param muons
///
void IBAnalyzerEM::SetMuonCollection(IBMuonCollection *muons)
{
  std::cout << "Setting muon collection with collection " << muons << std::endl;
  std::cout << "'Using full path ?' = '" << !m_useRecoPath << "'" << std::endl;
  uLibAssert(muons);
    std::cout << "Clearing " << std::endl;
    m_d->m_Events.clear();
    Vector<MuonScatterData>::iterator itr = muons->Data().begin();
    Vector<Vector<HPoint3f> >::iterator path_itr = muons->FullPath().begin();    
    std::cout << "Adding " << muons->Data().size() << " muons " << std::endl;    
    int iMu=0;
    while(itr != muons->Data().end()){
      iMu++;
      //      std::cout << "=============================" << std::endl;
      //      std::cout << "---" << iMu << " (" << path_itr->size() << ")" << std::endl;
      if(!AddMuonFullPath(*itr, *path_itr)){
	muons->Data().remove_element(*itr);
	if(muons->FullPath().size() > 0) muons->FullPath().remove_element(*path_itr);
      }
      else{
	itr++;
	path_itr++;
      }
    }
    std::cout << "\nDone, now calling base class..." << std::endl;
    BaseClass::SetMuonCollection(muons);
}

unsigned int IBAnalyzerEM::Size()
{
    return m_d->m_Events.size();
}

void IBAnalyzerEM::Run(unsigned int iterations, float muons_ratio)
{
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

void IBAnalyzerEM::SetMLAlgorithm(IBAnalyzerEMAlgorithm *MLAlgorithm)
{
    m_d->m_SijAlgorithm = MLAlgorithm;
}

void IBAnalyzerEM::filterEventsVoxelMask() {
    m_d->filterEventsVoxelMask();
}

void IBAnalyzerEM::filterEventsLineDistance(float min, float max) {
    m_d->filterEventsLineDistance(min, max);
}

void IBAnalyzerEM::SijCut(float threshold) {
    m_d->Evaluate(1);
    m_d->SijCut(threshold);
    this->GetVoxCollection()->UpdateDensity<UpdateDensitySijCapAlgorithm>(0);   // HARDCODE THRESHOLD
}

void IBAnalyzerEM::SijGuess(Vector<Vector2f> tpv)
{
    m_d->Evaluate(1);
    // ATTENZIONE!! il vettore deve essere ordinato per threshold crescenti   //
    for (int i=0; i<tpv.size(); ++i)
        m_d->SijGuess( tpv[i](0), tpv[i](1) );
    this->GetVoxCollection()->UpdateDensity<UpdateDensitySijCapAlgorithm>(0);   // HARDCODE THRESHOLD
}

void IBAnalyzerEM::Chi2Cut(float threshold)
{
    m_d->Evaluate(1);
    this->GetVoxCollection()->UpdateDensity<UpdateDensitySijCapAlgorithm>(0);   // HARDCODE THRESHOLD
    m_d->Chi2Cut(threshold);
}


void IBAnalyzerEM::SetVoxCollection(IBVoxCollection *voxels)
{
    if(this->GetMuonCollection()) {
        BaseClass::SetVoxCollection(voxels);
        this->SetMuonCollection(BaseClass::GetMuonCollection());
    }
    else
        std::cerr << "*** Analyzer EM is unable to reset Voxels ***\n" <<
                     "*** without a defined muon collection ... ***\n";
}

void IBAnalyzerEM::SetVoxcollectionShift(Vector3f shift)
{
  if(this->GetMuonCollection()) {
    IBVoxCollection *voxels = this->GetVoxCollection();
    Vector3f pos = voxels->GetPosition();
    voxels->SetPosition(pos + shift);
    IBMuonCollection *muons = this->GetMuonCollection();
    for(int i=0; i<muons->size(); ++i)
        this->AddMuon(muons->At(i));
  }
}


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



////////////////////////////////////////////////////////////////////////////////
/// dump events on rootuple
////////////////////////////////////////////////////////////////////////////////

void IBAnalyzerEM::dumpEventsTTree(const char *filename)
{
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
    TBranch *bev = tree->Branch("ev",&ev,"ev/I");
    TBranch *bp = tree->Branch("p",&mom,"p/F");
    TBranch *bsumLij = tree->Branch("sumLij",&sumLij,"sumLij/F");
    TBranch *bdist = tree->Branch("dist",&dist,"dist/F");
    TBranch *bDP = tree->Branch("DP",&DP,"DP/F");
    TBranch *bDX = tree->Branch("DX",&DX,"DX/F");
    TBranch *bDT = tree->Branch("DT",&DT,"DT/F");
    TBranch *bDZ = tree->Branch("DZ",&DZ,"DZ/F");

    IBMuonCollection *muons = this->GetMuonCollection();
    /// event loop
    int pos = 0;
    Vector< Event >::iterator itr = m_d->m_Events.begin();
    while (itr != m_d->m_Events.end()) {

        Event & evc = *itr;

        /// crossed voxel loop
        sumLij = 0;
        Vector< Event::Element >::iterator itre = evc.elements.begin();
        while (itre != evc.elements.end()) {
            Event::Element & elc = *itre;
            sumLij += elc.Wij(0,0);
            ++itre;
        }
        mom = $$.nominal_momentum/sqrt(evc.header.InitialSqrP);

        DP = evc.header.Di[0];
        DX = evc.header.Di[1];
        DT = evc.header.Di[2];
        DZ = evc.header.Di[3];

        bev->Fill();
        bp->Fill();
        bsumLij->Fill();

        bDP->Fill();
        bDX->Fill();
        bDT->Fill();
        bDZ->Fill();

        pos++;
        ev++;
        itr++;
    }

    /// muon loop to add poca information
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




