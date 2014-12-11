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


#include <iostream>
#include <iomanip>
#include <string.h>
#include <math.h>
#include <vector>
#include "TH1.h"
#include "TH2F.h"
#include "TH2.h"
#include "TH3F.h"
#include "TF1.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include "TMath.h"
#include "TChain.h"
#include "TGraph.h"
#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TApplication.h"
#include <omp.h>
#include "TFile.h"
#include "TTree.h"
#include "IBMinimizationVariablesEvaluator.h"
#include "IBMuonEventTTreeReader.h"
#include "IBMuonEventTTreeR3DmcReader.h"
#include "testing-prototype.h"
#include "IBVoxRaytracer.h"
#include "IBVoxCollectionCap.h"

using namespace uLib;


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// main
int main(int argc, char ** argv) {

    TFile *_file0 = TFile::Open(argv[1]);
    TTree * T = (TTree*)_file0->Get("RADMU");

    ////book histograms
    ///DTBX
    // occupancy, 3 for every chamber
    TH1F * hnh[12];
    TH1F * hnh_first_hit[12];

    for(int ih=0; ih<12; ih++){
        char hname[100];
        sprintf(hname,"hnh_%i",ih);
        hnh[ih] = new TH1F(hname,"Occupancy",80,0.,80.);
        char hname_fi[100];
        sprintf(hname_fi,"hnh_fi_%i",ih);
        hnh_first_hit[ih] = new TH1F(hname_fi,"Occupancy first hit",100,0.,200.);
    }

    //DTBX block
    Int_t ihit;
    Int_t hlay[127], trtime[127], htube[127], htimein[127], htime[127], htimetube[127],
            hflag[127], hnhitstube[127];

    // associate variables to branches
    T->SetBranchAddress("DTBX_nhit", &ihit );
    T->SetBranchAddress("TRTIME", trtime );
    T->SetBranchAddress("DTBX_lay", hlay );
    T->SetBranchAddress("DTBX_tube", htube );
    T->SetBranchAddress("DTBX_time_in", htimein );
    T->SetBranchAddress("DTBX_time", htime);
    T->SetBranchAddress("DTBX_time_tube", htimetube);
    T->SetBranchAddress("DTBX_Nhit_tube", hnhitstube);

    //SEG block
    Int_t iseg;
    Int_t segN[50];
    Float_t segX[50], segS[50], segK[50], segT[50];

    // associate variables to branches
    T->SetBranchAddress( "SEG_ns", &iseg );
    T->SetBranchAddress( "SEG_sx",  segX );
    T->SetBranchAddress( "SEG_ss",  segS );
    T->SetBranchAddress( "SEG_sk",  segK );
    T->SetBranchAddress( "SEG_sn",  segN );
    T->SetBranchAddress( "SEG_t0",  segT );

    //get number of events
    Int_t evnum = T->GetEntries();
    std::cout << "Entries " << evnum << std::endl;

    //loop over the events
    for (Int_t ientry=0 ; ientry < TMath::Min(10000,evnum); ientry ++)
    {
        int shift = 0;//12000000;

        if(floor(ientry/10000)*10000==ientry)
        std::cout << "Event " << ientry << std::endl;

        //read tree
        T->GetEntry(ientry+shift);

        //occupancy of all hits
        for (int ih=0; ih<ihit; ih++){
            int ch = static_cast<int>(hlay[ih]/1000);
            int lay = hlay[ih] - 100*static_cast<int>(hlay[ih]/100);
            int sl = static_cast<int>((lay-1)/4)+1;
            int index = sl + 3*(ch-1) - 1;
            //std::cout << "ch : " << ch << ", lay " << lay << ", sl " << sl << std::endl;

            hnh[index]->Fill(htube[ih]);
        }

//        int nh[100];
//        for(int iv=0; iv<100; iv++)
//            nh[iv]=0;

//        for (int ih=0; ih<ihit; ih++)
//            if(hlay[ih]==2304)
//                nh[htube[ih]] += 1;

////        for (int it=0; it<100; it++)
////            if(nh[it]>0)
////                hnh2304_first_hit->Fill(it);
    }
    TApplication a("a", 0, 0);


    TCanvas * cHits1 = new TCanvas("cHits1","cHits",1200,900);
    cHits1->Divide(1,3);
    cHits1->cd(1);
    hnh[0]->Draw();
    cHits1->cd(2);
    hnh[1]->Draw();
    cHits1->cd(3);
    hnh[2]->Draw();

    TCanvas * cHits2 = new TCanvas("cHits2","cHits",1200,600);
    cHits2->Divide(1,3);
    cHits2->cd(1);
    hnh[3]->Draw();
    cHits2->cd(2);
    hnh[4]->Draw();
    cHits2->cd(3);
    hnh[5]->Draw();

    TCanvas * cHits3 = new TCanvas("cHits3","cHits",1200,600);
    cHits3->Divide(1,3);
    cHits3->cd(1);
    hnh[6]->Draw();
    cHits3->cd(2);
    hnh[9]->Draw();



    TPaveText *tpt = new TPaveText(0.8,0.35,0.99,0.65,"brNDC");
    tpt->SetFillColor(kWhite);
    tpt->SetTextAlign(12);
    tpt->SetBorderSize(1);
    char txt[20];
    sprintf(txt,"Run %i",2151);
    tpt->AddText(txt);
    tpt->AddText("#color[2]{First hit occupancy}");
    tpt->AddText("#color[1]{All hits occupancy}");
    tpt->Draw();

    //std::cout << "Reader has processed " << tot  << " events" << std::endl;


    a.Run(kTRUE);

    return 0;
}

