#include <typeinfo>

#include "IBAnalyzerEM3p.h"
#include "IBPocaEvaluator.h"
#include "IBLineDistancePocaEvaluator.h"
#include "IBMinimizationVariablesEvaluator.h"
#include "IBVoxCollection.h"

using namespace uLib;









void IBAnalyzerEM3p::SetPocaAlgorithm(IBPocaEvaluator *PocaAlgorithm)
{
    assert(PocaAlgorithm);
    if(typeid(*PocaAlgorithm) != typeid(IBLineDistancePocaEvaluator) ) {
        std::cerr << "Error: you passed a " << typeid(*PocaAlgorithm).name()
                  << "; only LineDistancePocaEvaluator must be used for IBAnalyzerEM3p\n";
        exit(1);
    }
    else BaseClass::SetPocaAlgorithm(PocaAlgorithm);
}




bool IBAnalyzerEM3p::AddMuon(const MuonScatterData &muon)
{
    if(unlikely(!GetRayAlgorithm() || !GetVarAlgorithm())) return false;
    Event evc;

    evc.header.InitialSqrP = $$.nominal_momentum/muon.GetMomentum();
    evc.header.InitialSqrP *= evc.header.InitialSqrP;

    if(likely(GetVarAlgorithm()->evaluate(muon))) {
        evc.header.Di = GetVarAlgorithm()->getDataVector();
        evc.header.E  = GetVarAlgorithm()->getCovarianceMatrix();

        // HARDCODED ... ZERO CROSS CORRELATION BETWEEN VIEWS //
        //        evc.header.E.block<2,2>(2,0) = Matrix2f::Zero();
        //        evc.header.E.block<2,2>(0,2) = Matrix2f::Zero();
        // .................................................. //

        // HARDCODED ... LESS ERROR ! //
        //evc.header.E /= 2;

    }
    else return false;

    IBLineDistancePocaEvaluator *PocaAlgorithm =
            static_cast<IBLineDistancePocaEvaluator *>(GetPocaAlgorithm());

    for(int i=0; i<3; ++i) {
        Event evcc = evc;
        IBVoxRaytracer::RayData ray;
        { // Get RayTrace RayData //
            HPoint3f entry_pt,poca,exit_pt;
            if( !GetRayAlgorithm()->GetEntryPoint(muon.LineIn(),entry_pt) ||
                    !GetRayAlgorithm()->GetExitPoint(muon.LineOut(),exit_pt) )
                return false;

            bool use_poca = false;
            if(PocaAlgorithm) {
                use_poca = PocaAlgorithm->evaluate(muon);
                switch(i) {
                case 0: poca = PocaAlgorithm->getPoca(); break;
                case 1: poca = PocaAlgorithm->getInTrackPoca(); break;
                case 2: poca = PocaAlgorithm->getOutTrackPoca(); break;
                }
            }
            if(use_poca && this->GetVoxCollection()->IsInsideBounds(poca)) {
                ray = GetRayAlgorithm()->TraceBetweenPoints(entry_pt,poca);
                ray.AppendRay( GetRayAlgorithm()->TraceBetweenPoints(poca,exit_pt) );
            }
            else {
                ray = GetRayAlgorithm()->TraceBetweenPoints(entry_pt,exit_pt);
            }
        }

        Event::Element elc;
        Scalarf T = ray.TotalLength();
        for(int i=0; i<ray.Data().size(); ++i)
        {
            // voxel //
            const IBVoxRaytracer::RayData::Element *el = &ray.Data().at(i);
            elc.voxel = &this->GetVoxCollection()->operator [](el->vox_id);
            // Wij   //
            Scalarf L = el->L;  T -= L;
            elc.Wij << L ,          L*L/2 + L*T,
                    L*L/2 + L*T, L*L*L/3 + L*L*T + L*T*T;
            // pw    //
            elc.pw = evcc.header.InitialSqrP;
            evcc.elements.push_back(elc);
        }
//#       pragma omp critical
        this->Events().push_back(evcc);
    }
    return true;
}








