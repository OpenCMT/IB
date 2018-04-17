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
#include "TGCanvas.h"
#include "TPaveText.h"
#include "TPad.h"
#include "TFrame.h"
#include "TString.h"


#include <vector>
#include <algorithm>
#include <iostream.h>
#include <stdio.h>
#include <iomanip.h>
#include <map.h>
#include <fstream.h>


static const float schultzToInvLrad = (3./15.)*(3./15); // 1/cm

Double_t graph_surface(TGraph *g, Int_t nb)
{
   Double_t  s = 0;
   Double_t *x = g->GetX();
   Double_t *y = g->GetY();
   Int_t     n = g->GetN();
   Double_t xmin = x[0], xmax = x[0];
   Double_t ymin = y[0], ymax = y[0];
   Int_t i;

   for (i=1; i<n; i++) {
      if (xmin>x[i]) xmin = x[i];
      if (xmax<x[i]) xmax = x[i];
      if (ymin>y[i]) ymin = y[i];
      if (ymax<y[i]) ymax = y[i];
   }

   Double_t dx  = (xmax-xmin)/nb;
   Double_t dy  = (ymax-ymin)/nb;
   Double_t dx2 = dx/2;
   Double_t dy2 = dy/2;
   Double_t ds  = dx*dy;

   Double_t xi, yi, xc, yc;
   for (xi=xmin; xi<xmax; xi=xi+dx) {
      for (yi=ymin; yi<ymax; yi=yi+dy) {
        xc = xi+dx2;
        yc = yi+dy2;
        if (TMath::IsInside(xc, yc, n, x, y)) {
           s = s + ds;
        }
      }
   }
   return s;
}

void profile_different_runs() {

    bool carotaFlag = 0;

    TFile * f;
    if(carotaFlag)
        f = new TFile("carotaOutput_vox5.root");
    else
        f = new TFile("fe1Output_vox5.root");

    TCanvas *cm = new TCanvas("cm","cm",1200,600);
    if(carotaFlag)
        cm->Divide(2,1);
    cm->cd(1);

    TMultiGraph * fe2moc = new TMultiGraph("mg_2183","Density mean");
    mg_mean->Add(mg_2161);
    mg_mean->Add(mg_2162);
    mg_mean->Add(mg_2164);
    mg_mean->Draw("ALP");

    if(carotaFlag)
        mg_mean->GetXaxis()->SetTitle("CDS slices [cm]");
    else
        mg_mean->GetXaxis()->SetTitle("Fe slices [cm]");
    mg_mean->GetYaxis()->SetTitle("1/L_{rad} [1/cm]");
    mg_mean->GetYaxis()->SetTitleOffset(1.5);
    mg_mean->SetMinimum(0.);
    if(carotaFlag)
        mg_mean->SetMaximum(0.2);
    else
        mg_mean->SetMaximum(1.);

    // 2161
    mg_2161->SetMarkerStyle(20);
    mg_2161->SetMarkerSize(0.8);
    mg_2161->SetMarkerColor(1);
    mg_2161->SetLineColor(1);

    //2162
    mg_2162->SetMarkerStyle(20);
    mg_2162->SetMarkerSize(0.8);
    mg_2162->SetMarkerColor(2);
    mg_2162->SetLineColor(2);

    //2164
    mg_2164->SetMarkerStyle(20);
    mg_2164->SetMarkerSize(0.8);
    mg_2164->SetMarkerColor(3);
    mg_2164->SetLineColor(3);


    TPaveText *pt = new TPaveText(.6,.6,.95,.9,"brNDC");
    pt->SetFillColor(kWhite);
    pt->SetTextAlign(12);
    pt->SetBorderSize(1);
    if(carotaFlag)
        pt->AddText("#bf{CDS slice analysis}");
    else
        pt->AddText("#bf{Fe slice analysis}");
    pt->AddText("#bf{Image reconstruction}");
    //pt->AddText("#bf{Voxel reweight}");
    pt->AddText("#bf{NO Sij guess}");
    pt->AddText("#bf{Sij cut N=30}");
    pt->AddText("#bf{4000 iterations}");
    pt->AddText("#bf{vox 5 cm, 480 min}");
    pt->AddText("#bf{Init image: Fe=10, CDS=2, air=0.07}");
    pt->AddText("#color[1]{Run 2161 : CDS along Z}");
    pt->AddText("#color[2]{Run 2162 : CDS along X}");
    pt->AddText("#color[3]{Run 2164 : CDS along Z EMPTY}");
    pt->Draw();

    if(carotaFlag){
        TGraph * mg_diff_2161 = new TGraph();
        mg_diff_2161->SetMarkerStyle(20);
        mg_diff_2161->SetMarkerSize(0.8);
        mg_diff_2161->SetMarkerColor(6);
        mg_diff_2161->SetLineColor(6);

        TGraph * mg_diff_2162 = new TGraph();
        mg_diff_2162->SetMarkerStyle(20);
        mg_diff_2162->SetMarkerSize(0.8);
        mg_diff_2162->SetMarkerColor(7);
        mg_diff_2162->SetLineColor(7);

        for(int ip=0; ip<mg_2161->GetN(); ip++){
            Double_t x1=0;
            Double_t y1=0;
            mg_2161->GetPoint(ip,x1,y1);
            Double_t x2=0;
            Double_t y2=0;
            mg_2162->GetPoint(ip,x2,y2);
            Double_t x4=0;
            Double_t y4=0;
            mg_2164->GetPoint(ip,x4,y4);

            float y1diff = y1-y4;
            float y2diff = y2-y4;
            mg_diff_2161->SetPoint(ip,x1,y1diff);
            mg_diff_2162->SetPoint(ip,x1,y2diff);
        }

        TMultiGraph * mg_diff = new TMultiGraph("mg_diff","Density difference");
        mg_diff->Add(mg_diff_2161);
        mg_diff->Add(mg_diff_2162);

        cm->cd(2);
        mg_diff->Draw("ALP");
        mg_diff->SetMinimum(0.);
        mg_diff->SetMaximum(0.2);
        if(carotaFlag)
            mg_diff->GetXaxis()->SetTitle("CDS slices [cm]");
        else
            mg_diff->GetXaxis()->SetTitle("Fe slices [cm]");
        mg_diff->GetYaxis()->SetTitle("1/L_{rad} difference [1/cm]");
        mg_diff->GetYaxis()->SetTitleOffset(1.5);

        TPaveText *pt2 = new TPaveText(.7,.8,.95,.9,"brNDC");
        pt2->SetFillColor(kWhite);
        pt2->SetTextAlign(12);
        pt2->SetBorderSize(1);
        pt2->AddText("#color[6]{Run 2161 - 2164}");
        pt2->AddText("#color[7]{Run 2162 - 2164}");
        pt2->Draw();
    }

    return;
}

void profile_empty(int run, TGraphErrors *& mgfinal, float & norm) {

    TMultiGraph * mgcarote = new TMultiGraph();

    /// empty carota run 2208
    // default tol 10 cm radious 7.5
    TFile * fempty = new TFile("results0.95/2208/test_tol10_radius7.5/2208_carotaOutput_vox2.5.root");
    TGraphErrors * mgempty1;
    fempty->GetObject("mg_2208",mgempty1);
    mgcarote->Add(mgempty1);
    mgempty1->SetMarkerStyle(20);
    mgempty1->SetMarkerSize(0.8);
    mgempty1->SetMarkerColor(2);
    mgempty1->SetLineColor(2);
    // default tol 10 cm radious 6.5
    TFile * fempty = new TFile("results0.95/2208/test_tol10_radius6.5/2208_carotaOutput_vox2.5.root");
    TGraphErrors * mgempty2;
    fempty->GetObject("mg_2208",mgempty2);
    mgcarote->Add(mgempty2);
    mgempty2->SetMarkerStyle(20);
    mgempty2->SetMarkerSize(0.8);
    mgempty2->SetMarkerColor(3);
    mgempty2->SetLineColor(3);
    // default tol 5 cm radious 7.5
    TFile * fempty = new TFile("results0.95/2208/test_tol5/2208_carotaOutput_vox2.5.root");
    TGraphErrors * mgempty3;
    fempty->GetObject("mg_2208",mgempty3);
    mgcarote->Add(mgempty3);
    mgempty3->SetMarkerStyle(20);
    mgempty3->SetMarkerSize(0.8);
    mgempty3->SetMarkerColor(4);
    mgempty3->SetLineColor(4);
    // default tol 5 cm radious 6.5
    TFile * fempty = new TFile("results0.95/2208/test_tol5_radius6.5/2208_carotaOutput_vox2.5.root");
    TGraphErrors * mgempty4;
    fempty->GetObject("mg_2208",mgempty4);
    mgcarote->Add(mgempty4);
    mgempty4->SetMarkerStyle(20);
    mgempty4->SetMarkerSize(0.8);
    mgempty4->SetMarkerColor(6);
    mgempty4->SetLineColor(6);

    TPaveText *pt = new TPaveText(.6,.6,.95,.9,"brNDC");
    pt->SetFillColor(kWhite);
    pt->SetTextAlign(12);
    pt->SetBorderSize(1);
    pt->AddText("#bf{Run 2208}");
    pt->AddText("#bf{Empty CDS analysis}");
    pt->AddText("#color[2]{Radius 7.5 cm, XY toleance 10 cm}");
    pt->AddText("#color[3]{Radius 6.5 cm, XY toleance 10 cm}");
    pt->AddText("#color[4]{Radius 7.5 cm, XY toleance 5 cm}");
    pt->AddText("#color[6]{Radius 6.5 cm, XY toleance 5 cm}");

    TCanvas *cs = new TCanvas("cs","cs",800,600);
    cs->cd(1);
    mgcarote->Draw("ALP");
    mgcarote->SetTitle("Empty CDS Density Profile");
    mgcarote->SetMinimum(0.);
    mgcarote->SetMaximum(15.);
    mgcarote->GetXaxis()->SetTitle("CDS slices [cm]");
    mgcarote->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    mgcarote->GetYaxis()->SetTitleOffset(1.1);
    pt->Draw();

    /// compute lsd mean central 16 values
    float mean = 0;
    for(int ip=1; ip<16; ip++){
        Double_t x1=0;
        Double_t y1=0;
        mgempty2->GetPoint(ip,x1,y1);
        mean += y1;
    }
    mean /= 15;

    /// CDS mean
    std::cout << "lsd CDS mean value is " << mean <<
                 " [rad2/m], " << mean/100./schultzToInvLrad  << " schultz " << std::endl;
    return;
}
void profileFe(int run, TGraphErrors *& mgfinal, float & norm) {

    char graphname[20];
    sprintf(graphname,"mg_%i",run);
    TGraphErrors * mg;

    char file[100];
    sprintf(file, "results0.95/%i/test_normX/%i_fe1Output_vox2.5.root",run,run);
    TFile * ffeX10 = new TFile(file);
    sprintf(file, "results0.95/%i/test_normY/%i_fe1Output_vox2.5.root",run,run);
    TFile * ffeY10 = new TFile(file);
    sprintf(file, "results0.95/%i/test_normZ/%i_fe1Output_vox2.5.root",run,run);
    TFile * ffeZ10 = new TFile(file);

    sprintf(file, "results0.95/%i/test_AIR_X_tol0/%i_fe1Output_vox2.5.root",run,run);
    TFile * ffeX0 = new TFile(file);
    sprintf(file, "results0.95/%i/test_AIR_X_tol5/%i_fe1Output_vox2.5.root",run,run);
    TFile * ffeX5 = new TFile(file);
    max = 50.;

    TPaveText *pt = new TPaveText(.7,.6,.95,.9,"brNDC");
    pt->SetFillColor(kWhite);
    pt->SetTextAlign(12);
    pt->SetBorderSize(1);
    pt->AddText("#color[2]{Run 2183}");
    //pt->AddText("#bf{Slice analysis}");
    pt->AddText("#bf{Image reconstruction}");
    //pt->AddText("#bf{Voxel reweight}");
    //pt->AddText("#bf{NO Sij guess}");
    pt->AddText("#bf{Sij cut N=30}");
    pt->AddText("#bf{5000 iterations}");
    pt->AddText("#bf{vox 2.5 cm, 460 min}");
    pt->AddText("#bf{Init with real image}");
//    pt->AddText("#bf{Fe=10, Femockup=2}");
//    pt->AddText("#bf{CDS=2, air=0.07}");

    TCanvas *cfe = new TCanvas("cfe","cfe",800,600);
    cfe->Divide(3,2);
    cfe->cd(1);
    TGraphErrors * gfeX;
    ffeX10->GetObject(graphname,gfeX);
    Double_t x1=0;
    Double_t y1=0;
    gfeX->GetPoint(gfeX->GetN()-1,x1,y1);
    gfeX->RemovePoint(gfeX->GetN()-1);
    TLine * lX = new TLine(x1,y1,68,y1);
    lX->SetLineColor(3);
    lX->SetLineWidth(3.);
    lX->SetLineStyle(7);
    gfeX->SetTitle("Air block LSD profile along X");
    gfeX->GetXaxis()->SetTitle("X [cm]");
    gfeX->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfeX->GetYaxis()->SetTitleOffset(1.1);
    gfeX->SetMinimum(0.);
    gfeX->SetMaximum(max);
    gfeX->SetMarkerStyle(20);
    gfeX->SetMarkerSize(0.8);
    gfeX->SetMarkerColor(3);
    gfeX->SetLineColor(3);
    gfeX->Draw("ALP");
    lX->Draw("same");
    cfe->cd(2);
    TGraphErrors * gfeY;
    ffeY10->GetObject(graphname,gfeY);
    Double_t x2=0;
    Double_t y2=0;
    gfeY->GetPoint(gfeY->GetN()-1,x2,y2);
    TLine * lY = new TLine(x2,y2,40,y2);
    lY->SetLineColor(4);
    lY->SetLineWidth(3.);
    lY->SetLineStyle(7);
    gfeY->RemovePoint(gfeY->GetN()-1);
    gfeY->SetTitle("Air block LSD profile along Y");
    gfeY->GetXaxis()->SetTitle("Y [cm]");
    gfeY->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfeY->GetYaxis()->SetTitleOffset(1.1);
    gfeY->SetMinimum(0.);
    gfeY->SetMaximum(max);
    gfeY->SetMarkerStyle(20);
    gfeY->SetMarkerSize(0.8);
    gfeY->SetMarkerColor(4);
    gfeY->SetLineColor(4);
    gfeY->Draw("ALP");
    lY->Draw("same");
    cfe->cd(3);
    TGraphErrors * gfeZ;
    ffeZ10->GetObject(graphname,gfeZ);
    Double_t x3=0;
    Double_t y3=0;
    gfeZ->GetPoint(gfeZ->GetN()-1,x3,y3);
    TLine * lZ = new TLine(x3,y3,68,y3);
    lZ->SetLineColor(6);
    lZ->SetLineWidth(3.);
    lZ->SetLineStyle(7);
    gfeZ->RemovePoint(gfeZ->GetN()-1);
    gfeZ->SetTitle("Air block LSD profile along Z");
    gfeZ->GetXaxis()->SetTitle("Z [cm]");
    gfeZ->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfeZ->GetYaxis()->SetTitleOffset(1.1);
    gfeZ->SetMinimum(0.);
    gfeZ->SetMaximum(max);
    gfeZ->SetMarkerStyle(20);
    gfeZ->SetMarkerSize(0.8);
    gfeZ->SetMarkerColor(6);
    gfeZ->SetLineColor(6);
    gfeZ->Draw("ALP");
    lZ->Draw("same");
    //
    cfe->cd(4);
    TGraphErrors * gfeX0;
    ffeX0->GetObject(graphname,gfeX0);
    Double_t x4=0;
    Double_t y4=0;
    gfeX0->GetPoint(gfeX0->GetN()-1,x4,y4);
    TLine * lX0 = new TLine(x4,y4,68,y4);
    std::cout << "All block density " << y4 << std::endl;
    lX0->SetLineColor(6);
    lX0->SetLineWidth(3.);
    lX0->SetLineStyle(7);
    gfeX0->RemovePoint(gfeX0->GetN()-1);
    gfeX0->SetTitle("Air block LSD profile along X - tolerance 0 cm");
    gfeX0->GetXaxis()->SetTitle("X [cm]");
    gfeX0->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfeX0->GetYaxis()->SetTitleOffset(1.1);
    gfeX0->SetMinimum(0.);
    gfeX0->SetMaximum(max);
    gfeX0->SetMarkerStyle(20);
    gfeX0->SetMarkerSize(0.8);
    gfeX0->SetMarkerColor(6);
    gfeX0->SetLineColor(6);
    gfeX0->Draw("ALP");
    lX0->Draw("same");
    cfe->cd(5);
    TGraphErrors * gfeX5;
    ffeX5->GetObject(graphname,gfeX5);
    Double_t x5=0;
    Double_t y5=0;
    gfeX5->GetPoint(gfeX5->GetN()-1,x5,y5);
    TLine * lX5 = new TLine(x5,y5,68,y5);
    lX5->SetLineColor(7);
    lX5->SetLineWidth(3.);
    lX5->SetLineStyle(7);
    gfeX5->RemovePoint(gfeX5->GetN()-1);
    gfeX5->SetTitle("Air block LSD profile along X - tolerance 5 cm");
    gfeX5->GetXaxis()->SetTitle("X [cm]");
    gfeX5->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfeX5->GetYaxis()->SetTitleOffset(1.1);
    gfeX5->SetMinimum(0.);
    gfeX5->SetMaximum(max);
    gfeX5->SetMarkerStyle(20);
    gfeX5->SetMarkerSize(0.8);
    gfeX5->SetMarkerColor(7);
    gfeX5->SetLineColor(7);
    gfeX5->Draw("ALP");
    lX5->Draw("same");
    cfe->cd(6);
    gfeX->SetTitle("Air block LSD profile along X - tolerance 10 cm");
    gfeX->Draw("ALP");
    lX->Draw("same");
    pt->Draw();


    // graph integral
    std::cout << "Integral profile X = " << gfeX->Integral() << std::endl;
    std::cout << "Integral profile Y = " << gfeY->Integral() << std::endl;

    return;
}

void profileFemockupTol(int run, TGraphErrors *& mgfinal, float & norm) {

    char graphname[20];
    sprintf(graphname,"mg_%i",run);
    TGraphErrors * mg;

    char file[100];
    const char * dir1 = "_5_10";
    sprintf(file, "results0.95/%i/normX%s/%i_fe1mockupOutput_vox2.5.root",run,dir1,run);
    TFile * ffeX10 = new TFile(file);
    sprintf(file, "results0.95/%i/normY%s/%i_fe1mockupOutput_vox2.5.root",run,dir1,run);
    TFile * ffeY10 = new TFile(file);
    sprintf(file, "results0.95/%i/normZ%s/%i_fe1mockupOutput_vox2.5.root",run,dir1,run);
    TFile * ffeZ10 = new TFile(file);
    const char * dir2 = "_2_7";
    sprintf(file, "results0.95/%i/normX%s/%i_fe1mockupOutput_vox2.5.root",run,dir2,run);
    TFile * ffeX0 = new TFile(file);
    sprintf(file, "results0.95/%i/normY%s/%i_fe1mockupOutput_vox2.5.root",run,dir2,run);
    TFile * ffeY0 = new TFile(file);
    sprintf(file, "results0.95/%i/normZ%s/%i_fe1mockupOutput_vox2.5.root",run,dir2,run);
    TFile * ffeZ0 = new TFile(file);

    max = 20.;

    TPaveText *pt = new TPaveText(.7,.6,.95,.9,"brNDC");
    pt->SetFillColor(kWhite);
    pt->SetTextAlign(12);
    pt->SetBorderSize(1);
    pt->AddText("#color[2]{Run 2208}");
    //pt->AddText("#bf{Slice analysis}");
    pt->AddText("#bf{Image reconstruction}");
    //pt->AddText("#bf{Voxel reweight}");
    //pt->AddText("#bf{NO Sij guess}");
    pt->AddText("#bf{Sij cut N=30}");
    pt->AddText("#bf{5000 iterations}");
    pt->AddText("#bf{vox 2.5 cm, 460 min}");
    pt->AddText("#bf{Init with real image}");
//    pt->AddText("#bf{Fe=10, Femockup=2}");
//    pt->AddText("#bf{CDS=2, air=0.07}");

    TCanvas *cfe = new TCanvas("cfe","cfe",800,600);
    cfe->Divide(3,2);
    cfe->cd(1);
    TGraphErrors * gfeX;
    ffeX10->GetObject(graphname,gfeX);
    Double_t x1=0;
    Double_t y1=0;
    gfeX->GetPoint(gfeX->GetN()-1,x1,y1);
    gfeX->RemovePoint(gfeX->GetN()-1);
    TLine * lX = new TLine(x1,y1,68,y1);
    lX->SetLineColor(3);
    lX->SetLineWidth(3.);
    lX->SetLineStyle(7);
    gfeX->SetTitle("mockup block LSD profile - tolerance 5 10 5 cm");
    gfeX->GetXaxis()->SetTitle("X [cm]");
    gfeX->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfeX->GetYaxis()->SetTitleOffset(1.1);
    gfeX->SetMinimum(0.);
    gfeX->SetMaximum(max);
    gfeX->SetMarkerStyle(20);
    gfeX->SetMarkerSize(0.8);
    gfeX->SetMarkerColor(3);
    gfeX->SetLineColor(3);
    gfeX->Draw("ALP");
    lX->Draw("same");
    cfe->cd(2);
    TGraphErrors * gfeY;
    ffeY10->GetObject(graphname,gfeY);
    Double_t x2=0;
    Double_t y2=0;
    gfeY->GetPoint(gfeY->GetN()-1,x2,y2);
    TLine * lY = new TLine(x2,y2,60,y2);
    lY->SetLineColor(4);
    lY->SetLineWidth(3.);
    lY->SetLineStyle(7);
    gfeY->RemovePoint(gfeY->GetN()-1);
    gfeY->SetTitle("mockup block LSD profile - tolerance 5 10 5 cm");
    gfeY->GetXaxis()->SetTitle("Y [cm]");
    gfeY->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfeY->GetYaxis()->SetTitleOffset(1.1);
    gfeY->SetMinimum(0.);
    gfeY->SetMaximum(max);
    gfeY->SetMarkerStyle(20);
    gfeY->SetMarkerSize(0.8);
    gfeY->SetMarkerColor(4);
    gfeY->SetLineColor(4);
    gfeY->Draw("ALP");
    lY->Draw("same");
    cfe->cd(3);
    TGraphErrors * gfeZ;
    ffeZ10->GetObject(graphname,gfeZ);
    Double_t x3=0;
    Double_t y3=0;
    gfeZ->GetPoint(gfeZ->GetN()-1,x3,y3);
    TLine * lZ = new TLine(x3,y3,68,y3);
    lZ->SetLineColor(6);
    lZ->SetLineWidth(3.);
    lZ->SetLineStyle(7);
    gfeZ->RemovePoint(gfeZ->GetN()-1);
    gfeZ->SetTitle("mockup block LSD profile - tolerance 5 10 5 cm");
    gfeZ->GetXaxis()->SetTitle("Z [cm]");
    gfeZ->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfeZ->GetYaxis()->SetTitleOffset(1.1);
    gfeZ->SetMinimum(0.);
    gfeZ->SetMaximum(max);
    gfeZ->SetMarkerStyle(20);
    gfeZ->SetMarkerSize(0.8);
    gfeZ->SetMarkerColor(6);
    gfeZ->SetLineColor(6);
    gfeZ->Draw("ALP");
    lZ->Draw("same");
    //
    cfe->cd(4);
    TGraphErrors * gfeX0;
    ffeX0->GetObject(graphname,gfeX0);
    Double_t x4=0;
    Double_t y4=0;
    gfeX0->GetPoint(gfeX0->GetN()-1,x4,y4);
    TLine * lX0 = new TLine(x4,y4,68,y4);
    std::cout << "All block density " << y4 << std::endl;
    lX0->SetLineColor(3);
    lX0->SetLineWidth(3.);
    lX0->SetLineStyle(7);
    gfeX0->RemovePoint(gfeX0->GetN()-1);
    gfeX0->SetTitle("mockup block LSD profile - tolerance 2.5 7.5 2.5 cm");
    gfeX0->GetXaxis()->SetTitle("X [cm]");
    gfeX0->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfeX0->GetYaxis()->SetTitleOffset(1.1);
    gfeX0->SetMinimum(0.);
    gfeX0->SetMaximum(max);
    gfeX0->SetMarkerStyle(20);
    gfeX0->SetMarkerSize(0.8);
    gfeX0->SetMarkerColor(3);
    gfeX0->SetLineColor(3);
    gfeX0->Draw("ALP");
    lX0->Draw("same");
    cfe->cd(5);
    TGraphErrors * gfeY0;
    ffeY0->GetObject(graphname,gfeY0);
    Double_t x5=0;
    Double_t y5=0;
    gfeY0->GetPoint(gfeY0->GetN()-1,x5,y5);
    TLine * lX5 = new TLine(x5,y5,60,y5);
    lX5->SetLineColor(4);
    lX5->SetLineWidth(3.);
    lX5->SetLineStyle(7);
    gfeY0->RemovePoint(gfeY0->GetN()-1);
    gfeY0->SetTitle("mockup block LSD profile - tolerance 2.5 7.5 2.5 cm");
    gfeY0->GetXaxis()->SetTitle("Y [cm]");
    gfeY0->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfeY0->GetYaxis()->SetTitleOffset(1.1);
    gfeY0->SetMinimum(0.);
    gfeY0->SetMaximum(max);
    gfeY0->SetMarkerStyle(20);
    gfeY0->SetMarkerSize(0.8);
    gfeY0->SetMarkerColor(4);
    gfeY0->SetLineColor(4);
    gfeY0->Draw("ALP");
    lX5->Draw("same");
    cfe->cd(6);
    TGraphErrors * gfeZ0;
    ffeZ0->GetObject(graphname,gfeZ0);
    Double_t x5=0;
    Double_t y5=0;
    gfeZ0->GetPoint(gfeZ0->GetN()-1,x5,y5);
    TLine * lX5 = new TLine(x5,y5,68,y5);
    lX5->SetLineColor(6);
    lX5->SetLineWidth(3.);
    lX5->SetLineStyle(7);
    gfeZ0->RemovePoint(gfeZ0->GetN()-1);
    gfeZ0->SetTitle("mockup block LSD profile - tolerance 2.5 7.5 2.5 cm");
    gfeZ0->GetXaxis()->SetTitle("Z [cm]");
    gfeZ0->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfeZ0->GetYaxis()->SetTitleOffset(1.1);
    gfeZ0->SetMinimum(0.);
    gfeZ0->SetMaximum(max);
    gfeZ0->SetMarkerStyle(20);
    gfeZ0->SetMarkerSize(0.8);
    gfeZ0->SetMarkerColor(6);
    gfeZ0->SetLineColor(6);
    gfeZ0->Draw("ALP");
    lX5->Draw("same");
    pt->Draw();


    // graph integrals
    std::cout << "Tolerance 2.5 7.5 2.5"  << std::endl;
    std::cout << "Integral profile X = " << gfeX0->Integral() << std::endl;
    std::cout << "Integral profile Y = " << gfeY0->Integral() << std::endl;
    std::cout << "Integral profile Z = " << gfeZ0->Integral() << std::endl;

    std::cout << "Tolerance 5 10 5"  << std::endl;
    std::cout << "Integral profile X = " << gfeX->Integral() << std::endl;
    std::cout << "Integral profile Y = " << gfeY->Integral() << std::endl;
    std::cout << "Integral profile Z = " << gfeZ->Integral() << std::endl;

    return;
}

void profileCarotaTol(int run, TGraphErrors *& mgfinal, float & norm) {

    char graphname[20];
    sprintf(graphname,"mg_%i",run);
    TGraphErrors * mg;

    char file[100];
    const char * dir1 = "_5_10";
    sprintf(file, "results0.95/%i/normX%s/%i_carotaOutput_vox2.5.root",run,dir1,run);
    TFile * ffeX10 = new TFile(file);
    sprintf(file, "results0.95/%i/normY%s/%i_carotaOutput_vox2.5.root",run,dir1,run);
    TFile * ffeY10 = new TFile(file);
    sprintf(file, "results0.95/%i/normZ%s/%i_carotaOutput_vox2.5.root",run,dir1,run);
    TFile * ffeZ10 = new TFile(file);
    const char * dir2 = "_2_7";
    sprintf(file, "results0.95/%i/normX%s/%i_carotaOutput_vox2.5.root",run,dir2,run);
    TFile * ffeX0 = new TFile(file);
    sprintf(file, "results0.95/%i/normY%s/%i_carotaOutput_vox2.5.root",run,dir2,run);
    TFile * ffeY0 = new TFile(file);
    sprintf(file, "results0.95/%i/normZ%s/%i_carotaOutput_vox2.5.root",run,dir2,run);
    TFile * ffeZ0 = new TFile(file);

    max = 20.;

    TPaveText *pt = new TPaveText(.7,.6,.95,.9,"brNDC");
    pt->SetFillColor(kWhite);
    pt->SetTextAlign(12);
    pt->SetBorderSize(1);
    const char runlabel[10];
    sprintf(runlabel, "#color[2]{Run %i}",run);
    pt->AddText(runlabel);
    //pt->AddText("#bf{Slice analysis}");
    pt->AddText("#bf{Image reconstruction}");
    //pt->AddText("#bf{Voxel reweight}");
    //pt->AddText("#bf{NO Sij guess}");
    pt->AddText("#bf{Sij cut N=30}");
    pt->AddText("#bf{5000 iterations}");
    pt->AddText("#bf{vox 2.5 cm, 460 min}");
    pt->AddText("#bf{Init with real image}");
//    pt->AddText("#bf{Fe=10, Femockup=2}");
//    pt->AddText("#bf{CDS=2, air=0.07}");

    TCanvas *cfe = new TCanvas("cfe","cfe",800,600);
    cfe->Divide(3,2);
    //cfe->Divide(2,1);
    cfe->cd(1);
    TGraphErrors * gfeX;
    ffeX10->GetObject(graphname,gfeX);
    Double_t x1=0;
    Double_t y1=0;
    gfeX->GetPoint(gfeX->GetN()-1,x1,y1);
    gfeX->RemovePoint(gfeX->GetN()-1);
    TLine * lX = new TLine(x1,y1,30,y1);
    lX->SetLineColor(1);
    lX->SetLineWidth(3.);
    lX->SetLineStyle(7);
    gfeX->SetTitle("carota block LSD profile - tolerance 5 10 5 cm");
    gfeX->GetXaxis()->SetTitle("X [cm]");
    gfeX->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfeX->GetYaxis()->SetTitleOffset(1.1);
    gfeX->SetMinimum(0.);
    gfeX->SetMaximum(max);
    gfeX->SetMarkerStyle(20);
    gfeX->SetMarkerSize(0.8);
    gfeX->SetMarkerColor(3);
    gfeX->SetLineColor(3);
    gfeX->Draw("ALP");
    lX->Draw("same");
    cfe->cd(2);
    TGraphErrors * gfeY;
    ffeY10->GetObject(graphname,gfeY);
    Double_t x2=0;
    Double_t y2=0;
    gfeY->GetPoint(gfeY->GetN()-1,x2,y2);
    TLine * lY = new TLine(x2,y2,30,y2);
    lY->SetLineColor(1);
    lY->SetLineWidth(3.);
    lY->SetLineStyle(7);
    gfeY->RemovePoint(gfeY->GetN()-1);
    gfeY->SetTitle("carota block LSD profile - tolerance 5 10 5 cm");
    gfeY->GetXaxis()->SetTitle("Y [cm]");
    gfeY->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfeY->GetYaxis()->SetTitleOffset(1.1);
    gfeY->SetMinimum(0.);
    gfeY->SetMaximum(max);
    gfeY->SetMarkerStyle(20);
    gfeY->SetMarkerSize(0.8);
    gfeY->SetMarkerColor(4);
    gfeY->SetLineColor(4);
    gfeY->Draw("ALP");
    lY->Draw("same");
    //
    cfe->cd(4);
    TGraphErrors * gfeX0;
    ffeX0->GetObject(graphname,gfeX0);
    Double_t x4=0;
    Double_t y4=0;
    gfeX0->GetPoint(gfeX0->GetN()-1,x4,y4);
    TLine * lX0 = new TLine(x4,y4,30,y4);
    std::cout << "All block density " << y4 << std::endl;
    lX0->SetLineColor(1);
    lX0->SetLineWidth(3.);
    lX0->SetLineStyle(7);
    gfeX0->RemovePoint(gfeX0->GetN()-1);
    gfeX0->SetTitle("carota block LSD profile - tolerance 2.5 7.5 2.5 cm");
    gfeX0->GetXaxis()->SetTitle("X [cm]");
    gfeX0->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfeX0->GetYaxis()->SetTitleOffset(1.1);
    gfeX0->SetMinimum(0.);
    gfeX0->SetMaximum(max);
    gfeX0->SetMarkerStyle(20);
    gfeX0->SetMarkerSize(0.8);
    gfeX0->SetMarkerColor(3);
    gfeX0->SetLineColor(3);
    gfeX0->Draw("ALP");
    lX0->Draw("same");
    cfe->cd(5);
    TGraphErrors * gfeY0;
    ffeY0->GetObject(graphname,gfeY0);
    Double_t x5=0;
    Double_t y5=0;
    gfeY0->GetPoint(gfeY0->GetN()-1,x5,y5);
    TLine * lX5 = new TLine(x5,y5,30,y5);
    lX5->SetLineColor(1);
    lX5->SetLineWidth(3.);
    lX5->SetLineStyle(7);
    gfeY0->RemovePoint(gfeY0->GetN()-1);
    gfeY0->SetTitle("carota block LSD profile - tolerance 2.5 7.5 2.5 cm");
    gfeY0->GetXaxis()->SetTitle("Y [cm]");
    gfeY0->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfeY0->GetYaxis()->SetTitleOffset(1.1);
    gfeY0->SetMinimum(0.);
    gfeY0->SetMaximum(max);
    gfeY0->SetMarkerStyle(20);
    gfeY0->SetMarkerSize(0.8);
    gfeY0->SetMarkerColor(4);
    gfeY0->SetLineColor(4);
    gfeY0->Draw("ALP");
    lX5->Draw("same");
    cfe->cd(6);
    TGraphErrors * gfeZ0;
    ffeZ0->GetObject(graphname,gfeZ0);
    Double_t x5=0;
    Double_t y5=0;
    gfeZ0->GetPoint(gfeZ0->GetN()-1,x5,y5);
    TLine * lX5 = new TLine(15,y5,55,y5);
    lX5->SetLineColor(1);
    lX5->SetLineWidth(3.);
    lX5->SetLineStyle(7);
    gfeZ0->RemovePoint(gfeZ0->GetN()-1);
    gfeZ0->SetTitle("carota block LSD profile - tolerance 2.5 7.5 2.5 cm");
    gfeZ0->GetXaxis()->SetTitle("Z [cm]");
    gfeZ0->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfeZ0->GetYaxis()->SetTitleOffset(1.1);
    gfeZ0->SetMinimum(0.);
    gfeZ0->SetMaximum(max);
    gfeZ0->SetMarkerStyle(20);
    gfeZ0->SetMarkerSize(0.8);
    gfeZ0->SetMarkerColor(6);
    gfeZ0->SetLineColor(6);
    gfeZ0->Draw("ALP");
    lX5->Draw("same");
    cfe->cd(3);
    TGraphErrors * gfeZ;
    ffeZ10->GetObject(graphname,gfeZ);
    Double_t x3=0;
    Double_t y3=0;
    gfeZ->GetPoint(gfeZ->GetN()-1,x3,y3);
    TLine * lZ = new TLine(15.,y3,55.,y3);
    lZ->SetLineColor(1);
    lZ->SetLineWidth(3.);
    lZ->SetLineStyle(7);
    gfeZ->RemovePoint(gfeZ->GetN()-1);
    gfeZ->SetTitle("carota block LSD profile - tolerance 5 10 5 cm");
    gfeZ->GetXaxis()->SetTitle("Z [cm]");
    gfeZ->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfeZ->GetYaxis()->SetTitleOffset(1.1);
    gfeZ->SetMinimum(0.);
    gfeZ->SetMaximum(max);
    gfeZ->SetMarkerStyle(20);
    gfeZ->SetMarkerSize(0.8);
    gfeZ->SetMarkerColor(6);
    gfeZ->SetLineColor(6);
    gfeZ->Draw("ALP");
    lZ->Draw("same");
    pt->Draw();


    // graph integrals
    std::cout << "Tolerance 2.5 7.5 2.5"  << std::endl;
    std::cout << "Integral profile X = " << gfeX0->Integral() << std::endl;
    std::cout << "Integral profile Y = " << gfeY0->Integral() << std::endl;
    std::cout << "Integral profile Z = " << gfeZ0->Integral() << std::endl;

    std::cout << "Tolerance 5 10 5"  << std::endl;
    std::cout << "Integral profile X = " << gfeX->Integral() << std::endl;
    std::cout << "Integral profile Y = " << gfeY->Integral() << std::endl;
    std::cout << "Integral profile Z = " << gfeZ->Integral() << std::endl;

    return;
}

void profileFeTol(int run, TGraphErrors *& mgfinal, float & norm) {

    char graphname[20];
    sprintf(graphname,"mg_%i",run);
    TGraphErrors * mg;

    char file[100];
    const char * dir1 = "_5_10";
    sprintf(file, "results0.95/%i/normX%s/%i_fe2Output_vox2.5.root",run,dir1,run);
    TFile * ffeX10 = new TFile(file);
    sprintf(file, "results0.95/%i/normY%s/%i_fe2Output_vox2.5.root",run,dir1,run);
    TFile * ffeY10 = new TFile(file);
    sprintf(file, "results0.95/%i/normZ%s/%i_fe2Output_vox2.5.root",run,dir1,run);
    TFile * ffeZ10 = new TFile(file);
    const char * dir2 = "_2_7";
    sprintf(file, "results0.95/%i/normX%s/%i_fe2Output_vox2.5.root",run,dir2,run);
    TFile * ffeX0 = new TFile(file);
    sprintf(file, "results0.95/%i/normY%s/%i_fe2Output_vox2.5.root",run,dir2,run);
    TFile * ffeY0 = new TFile(file);
    sprintf(file, "results0.95/%i/normZ%s/%i_fe2Output_vox2.5.root",run,dir2,run);
    TFile * ffeZ0 = new TFile(file);

    max = 50.;

    TPaveText *pt = new TPaveText(.7,.6,.95,.9,"brNDC");
    pt->SetFillColor(kWhite);
    pt->SetTextAlign(12);
    pt->SetBorderSize(1);
    pt->AddText("#color[2]{Run 2208}");
    //pt->AddText("#bf{Slice analysis}");
    pt->AddText("#bf{Image reconstruction}");
    //pt->AddText("#bf{Voxel reweight}");
    //pt->AddText("#bf{NO Sij guess}");
    pt->AddText("#bf{Sij cut N=30}");
    pt->AddText("#bf{5000 iterations}");
    pt->AddText("#bf{vox 2.5 cm, 460 min}");
    pt->AddText("#bf{Init with real image}");
//    pt->AddText("#bf{Fe=10, Femockup=2}");
//    pt->AddText("#bf{CDS=2, air=0.07}");

    TCanvas *cfe = new TCanvas("cfe","cfe",800,600);
    cfe->Divide(3,2);
    cfe->cd(1);
    TGraphErrors * gfeX;
    ffeX10->GetObject(graphname,gfeX);
    Double_t x1=0;
    Double_t y1=0;
    gfeX->GetPoint(gfeX->GetN()-1,x1,y1);
    gfeX->RemovePoint(gfeX->GetN()-1);
    TLine * lX = new TLine(x1,y1,68,y1);
    lX->SetLineColor(3);
    lX->SetLineWidth(3.);
    lX->SetLineStyle(7);
    gfeX->SetTitle("Iron block LSD profile - tolerance 5 10 5 cm");
    gfeX->GetXaxis()->SetTitle("X [cm]");
    gfeX->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfeX->GetYaxis()->SetTitleOffset(1.1);
    gfeX->SetMinimum(0.);
    gfeX->SetMaximum(max);
    gfeX->SetMarkerStyle(20);
    gfeX->SetMarkerSize(0.8);
    gfeX->SetMarkerColor(3);
    gfeX->SetLineColor(3);
    gfeX->Draw("ALP");
    lX->Draw("same");
    cfe->cd(2);
    TGraphErrors * gfeY;
    ffeY10->GetObject(graphname,gfeY);
    Double_t x2=0;
    Double_t y2=0;
    gfeY->GetPoint(gfeY->GetN()-1,x2,y2);
    TLine * lY = new TLine(x2,y2,60,y2);
    lY->SetLineColor(4);
    lY->SetLineWidth(3.);
    lY->SetLineStyle(7);
    gfeY->RemovePoint(gfeY->GetN()-1);
    gfeY->SetTitle("Iron block LSD profile - tolerance 5 10 5 cm");
    gfeY->GetXaxis()->SetTitle("Y [cm]");
    gfeY->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfeY->GetYaxis()->SetTitleOffset(1.1);
    gfeY->SetMinimum(0.);
    gfeY->SetMaximum(max);
    gfeY->SetMarkerStyle(20);
    gfeY->SetMarkerSize(0.8);
    gfeY->SetMarkerColor(4);
    gfeY->SetLineColor(4);
    gfeY->Draw("ALP");
    lY->Draw("same");
    cfe->cd(3);
    TGraphErrors * gfeZ;
    ffeZ10->GetObject(graphname,gfeZ);
    Double_t x3=0;
    Double_t y3=0;
    gfeZ->GetPoint(gfeZ->GetN()-1,x3,y3);
    TLine * lZ = new TLine(x3,y3,68,y3);
    lZ->SetLineColor(6);
    lZ->SetLineWidth(3.);
    lZ->SetLineStyle(7);
    gfeZ->RemovePoint(gfeZ->GetN()-1);
    gfeZ->SetTitle("Iron block LSD profile - tolerance 5 10 5 cm");
    gfeZ->GetXaxis()->SetTitle("Z [cm]");
    gfeZ->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfeZ->GetYaxis()->SetTitleOffset(1.1);
    gfeZ->SetMinimum(0.);
    gfeZ->SetMaximum(max);
    gfeZ->SetMarkerStyle(20);
    gfeZ->SetMarkerSize(0.8);
    gfeZ->SetMarkerColor(6);
    gfeZ->SetLineColor(6);
    gfeZ->Draw("ALP");
    lZ->Draw("same");
    //
    cfe->cd(4);
    TGraphErrors * gfeX0;
    ffeX0->GetObject(graphname,gfeX0);
    Double_t x4=0;
    Double_t y4=0;
    gfeX0->GetPoint(gfeX0->GetN()-1,x4,y4);
    TLine * lX0 = new TLine(x4,y4,68,y4);
    std::cout << "All block density " << y4 << std::endl;
    lX0->SetLineColor(3);
    lX0->SetLineWidth(3.);
    lX0->SetLineStyle(7);
    gfeX0->RemovePoint(gfeX0->GetN()-1);
    gfeX0->SetTitle("Iron block LSD profile - tolerance 2.5 7.5 2.5 cm");
    gfeX0->GetXaxis()->SetTitle("X [cm]");
    gfeX0->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfeX0->GetYaxis()->SetTitleOffset(1.1);
    gfeX0->SetMinimum(0.);
    gfeX0->SetMaximum(max);
    gfeX0->SetMarkerStyle(20);
    gfeX0->SetMarkerSize(0.8);
    gfeX0->SetMarkerColor(3);
    gfeX0->SetLineColor(3);
    gfeX0->Draw("ALP");
    lX0->Draw("same");
    cfe->cd(5);
    TGraphErrors * gfeY0;
    ffeY0->GetObject(graphname,gfeY0);
    Double_t x5=0;
    Double_t y5=0;
    gfeY0->GetPoint(gfeY0->GetN()-1,x5,y5);
    TLine * lX5 = new TLine(x5,y5,60,y5);
    lX5->SetLineColor(4);
    lX5->SetLineWidth(3.);
    lX5->SetLineStyle(7);
    gfeY0->RemovePoint(gfeY0->GetN()-1);
    gfeY0->SetTitle("Iron block LSD profile - tolerance 2.5 7.5 2.5 cm");
    gfeY0->GetXaxis()->SetTitle("Y [cm]");
    gfeY0->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfeY0->GetYaxis()->SetTitleOffset(1.1);
    gfeY0->SetMinimum(0.);
    gfeY0->SetMaximum(max);
    gfeY0->SetMarkerStyle(20);
    gfeY0->SetMarkerSize(0.8);
    gfeY0->SetMarkerColor(4);
    gfeY0->SetLineColor(4);
    gfeY0->Draw("ALP");
    lX5->Draw("same");
    cfe->cd(6);
    TGraphErrors * gfeZ0;
    ffeZ0->GetObject(graphname,gfeZ0);
    Double_t x5=0;
    Double_t y5=0;
    gfeZ0->GetPoint(gfeZ0->GetN()-1,x5,y5);
    TLine * lX5 = new TLine(x5,y5,68,y5);
    lX5->SetLineColor(6);
    lX5->SetLineWidth(3.);
    lX5->SetLineStyle(7);
    gfeZ0->RemovePoint(gfeZ0->GetN()-1);
    gfeZ0->SetTitle("Iron block LSD profile - tolerance 2.5 7.5 2.5 cm");
    gfeZ0->GetXaxis()->SetTitle("Z [cm]");
    gfeZ0->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfeZ0->GetYaxis()->SetTitleOffset(1.1);
    gfeZ0->SetMinimum(0.);
    gfeZ0->SetMaximum(max);
    gfeZ0->SetMarkerStyle(20);
    gfeZ0->SetMarkerSize(0.8);
    gfeZ0->SetMarkerColor(6);
    gfeZ0->SetLineColor(6);
    gfeZ0->Draw("ALP");
    lX5->Draw("same");
    pt->Draw();


    // graph integrals
    std::cout << "Tolerance 2.5 7.5 2.5"  << std::endl;
    std::cout << "Integral profile X = " << gfeX0->Integral() << std::endl;
    std::cout << "Integral profile Y = " << gfeY0->Integral() << std::endl;
    std::cout << "Integral profile Z = " << gfeZ0->Integral() << std::endl;

    std::cout << "Tolerance 5 10 5"  << std::endl;
    std::cout << "Integral profile X = " << gfeX->Integral() << std::endl;
    std::cout << "Integral profile Y = " << gfeY->Integral() << std::endl;
    std::cout << "Integral profile Z = " << gfeZ->Integral() << std::endl;

    return;
}


float lsdMean(TGraphErrors * gr, int N){
    /// compute lsd mean over N highest lsd values
    vector<float> vm;
    for(int ip=0; ip<gr->GetN(); ip++){
        Double_t x1=0;
        Double_t y1=0;
        gr->GetPoint(ip,x1,y1);
        vm.push_back(y1);
    }
    std::sort(vm.begin(),vm.end());
    float mean = 0;
    for(int ie=1; ie<=N; ie++){
        mean += vm[vm.size()-ie];
        //std::cout << "bin " << ie << ", value " << vm[vm.size()-ie] << std::endl;
    }
    mean /= N;
    return mean;
}

void profileZ(int run, const char * dir, TGraphErrors *& mgfinal, float & norm, float & bf, float & bc, float & totlsd) {

    // empty carota run 2208
    const char * dir_empty = "normZ_5_10";
    char file[100];
    sprintf(file, "results0.95/%i/%s/%i_carotaOutput_vox2.5.root",2208,dir_empty,2208);
    TFile * fempty = new TFile(file);
    TGraphErrors * gempty;
    char graphnameempty[20];
    sprintf(graphnameempty,"mg_%i",2208);
    fempty->GetObject(graphnameempty,gempty);

    char graphname[20];
    sprintf(graphname,"mg_%i",run);
    TGraphErrors * gcarota;

    char file[100];
    // carota run
    sprintf(file, "results0.95/%i/%s/%i_carotaOutput_vox2.5.root",run,dir,run);
    std::cout << "Carota dir " << file << std::endl;
    TFile * fcarota = new TFile(file);
    // fe block from run
    sprintf(file, "results0.95/%i/%s/%i_fe1Output_vox2.5.root",run,dir,run);
    TFile * ffe = new TFile(file);
    sprintf(file, "results0.95/%i/%s/%i_fe2Output_vox2.5.root",run,dir,run);
    TFile * ffe2 = new TFile(file);
    // femockup block from run
    sprintf(file, "results0.95/%i/%s/%i_fe1mockupOutput_vox2.5.root",run,dir,run);
    TFile * ffemockup = new TFile(file);
    sprintf(file, "results0.95/%i/%s/%i_fe2mockupOutput_vox2.5.root",run,dir,run);
    TFile * ffemockup2 = new TFile(file);
//    // fe block from run 2208
//    sprintf(file, "results0.95/%i/%s/%i_fe1Output_vox2.5.root",2208,dir_empty,2208);
//    TFile * ffe = new TFile(file);
//    sprintf(file, "results0.95/%i/%s/%i_fe2Output_vox2.5.root",2208,dir_empty,2208);
//    TFile * ffe2 = new TFile(file);
//    // femockup block from run 2208
//    sprintf(file, "results0.95/%i/%s/%i_fe1mockupOutput_vox2.5.root",2208,dir_empty,2208);
//    TFile * ffemockup = new TFile(file);
//    sprintf(file, "results0.95/%i/%s/%i_fe2mockupOutput_vox2.5.root",2208,dir_empty,2208);
//    TFile * ffemockup2 = new TFile(file);


    TPaveText *pt = new TPaveText(.8,.75,.95,.94,"brNDC");
    pt->SetFillColor(kWhite);
    pt->SetTextAlign(12);
    pt->SetBorderSize(1);
    const char runlabel[10];
    sprintf(runlabel, "#color[2]{Run %i}",2208);
    pt->AddText(runlabel);
    //pt->AddText("#bf{Slice analysis}");
    pt->AddText("#bf{Image reconstruction}");
    //pt->AddText("#bf{Voxel reweight}");
    //pt->AddText("#bf{NO Sij guess}");
    pt->AddText("#bf{Sij cut N=30}");
    pt->AddText("#bf{5000 iterations}");
    pt->AddText("#bf{vox 2.5 cm, 460 min}");
    pt->AddText("#bf{Init with real image}");
//    pt->AddText("#bf{Fe=10, Femockup=2}");
//    pt->AddText("#bf{CDS=2, air=0.07}");

    TPaveText *ptr = new TPaveText(.8,.75,.95,.94,"brNDC");
    ptr->SetFillColor(kWhite);
    ptr->SetTextAlign(12);
    ptr->SetBorderSize(1);
    const char runlabel[10];
    sprintf(runlabel, "#color[2]{Run %i}",run);
    ptr->AddText(runlabel);
    //ptr->AddText("#bf{Slice analysis}");
    ptr->AddText("#bf{Image reconstruction}");
    //ptr->AddText("#bf{Voxel reweight}");
    //ptr->AddText("#bf{NO Sij guess}");
    ptr->AddText("#bf{Sij cut N=30}");
    ptr->AddText("#bf{5000 iterations}");
    ptr->AddText("#bf{vox 2.5 cm, 460 min}");
    ptr->AddText("#bf{Init with real image}");
//    ptr->AddText("#bf{Fe=10, Femockup=2}");
//    ptr->AddText("#bf{CDS=2, air=0.07}");



    TCanvas *cfe = new TCanvas("cfe","cfe",800,600);
    //c->cd(2);
    TGraphErrors * gfe1;
    ffe->GetObject(graphname,gfe1);
    //ffe->GetObject(graphnameempty,gfe1);
    gfe1->SetTitle("Fe block LSD Profile");
    gfe1->GetXaxis()->SetTitle("z [cm]");
    gfe1->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfe1->GetYaxis()->SetTitleOffset(1.1);
    gfe1->SetMinimum(0.);
    gfe1->SetMaximum(50.);
    gfe1->SetMarkerStyle(20);
    gfe1->SetMarkerSize(0.8);
    gfe1->SetMarkerColor(3);
    gfe1->SetLineColor(3);
    Double_t xf1=0;
    Double_t yf1=0;
    gfe1->GetPoint(gfe1->GetN()-1,xf1,yf1);
    gfe1->RemovePoint(gfe1->GetN()-1);
    TLine * lXf1 = new TLine(0,yf1,60,yf1);
    lXf1->SetLineColor(3);
    lXf1->SetLineWidth(3.);
    lXf1->SetLineStyle(7);
    //gfe1->Draw("ALP");
    TGraphErrors * gfe2;
    ffe2->GetObject(graphname,gfe2);
    //ffe2->GetObject(graphnameempty,gfe2);
    gfe2->SetTitle("Fe block LSD Profile");
    gfe2->GetXaxis()->SetTitle("z [cm]");
    gfe2->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gfe2->GetYaxis()->SetTitleOffset(1.1);
    gfe2->SetMinimum(0.);
    gfe2->SetMaximum(50.);
    gfe2->SetMarkerStyle(20);
    gfe2->SetMarkerSize(0.8);
    gfe2->SetMarkerColor(8);
    gfe2->SetLineColor(8);
    gfe2->SetLineStyle(9);
    Double_t xf2=0;
    Double_t yf2=0;
    gfe2->GetPoint(gfe2->GetN()-1,xf2,yf2);
    gfe2->RemovePoint(gfe2->GetN()-1);
    TLine * lXf2 = new TLine(0,yf2,60,yf2);
    lXf2->SetLineColor(8);
    lXf2->SetLineWidth(3.);
    lXf2->SetLineStyle(7);
    //gfe2->Draw("ALP");
    TMultiGraph * mgfe = new TMultiGraph();
    mgfe->Add(gfe1);
    mgfe->Add(gfe2);
    mgfe->Draw("ALP");
    mgfe->SetTitle("Fe block LSD Profile");
    mgfe->GetXaxis()->SetTitle("z [cm]");
    mgfe->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    mgfe->GetYaxis()->SetTitleOffset(1.1);
    mgfe->SetMinimum(0.);
    mgfe->SetMaximum(50.);
    lXf1->Draw("same");
    lXf2->Draw("same");
    pt->Draw();

    /// Fe mean
    std::cout << "lsd Fe mean value is " << lsdMean(gfe1,6) <<
                 " [rad2/m], " << lsdMean(gfe1,6)/100./schultzToInvLrad  << " schultz " << std::endl;
    std::cout << "Entire iron 1 mean lsd value is " << yf1 << std::endl;
    /// Fe mean
    std::cout << "lsd Fe mean value is " << lsdMean(gfe2,6) <<
                 " [rad2/m], " << lsdMean(gfe2,6)/100./schultzToInvLrad  << " schultz " << std::endl;
    std::cout << "Entire iron 2 mean lsd value is " << yf2 << std::endl;

    TCanvas *cm = new TCanvas("cm","cm",800,600);
    //c->cd(3);
    TGraphErrors * gmoc;
    ffemockup->GetObject(graphname,gmoc);
    //ffemockup->GetObject(graphnameempty,gmoc);
    gmoc->SetTitle("Mockup block LSD Profile");
    gmoc->GetXaxis()->SetTitle("Fe mockup slices [cm]");
    gmoc->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gmoc->GetYaxis()->SetTitleOffset(1.1);
    gmoc->SetMinimum(0.);
    gmoc->SetMaximum(15.);
    gmoc->SetMarkerStyle(20);
    gmoc->SetMarkerSize(0.8);
    gmoc->SetMarkerColor(4);
    gmoc->SetLineColor(4);
    Double_t xm1=0;
    Double_t ym1=0;
    gmoc->GetPoint(gmoc->GetN()-1,xm1,ym1);
    gmoc->RemovePoint(gmoc->GetN()-1);
    TLine * lXm1 = new TLine(0,ym1,60,ym1);
    lXm1->SetLineColor(4);
    lXm1->SetLineWidth(3.);
    lXm1->SetLineStyle(7);
    TGraphErrors * gmoc2;
    ffemockup2->GetObject(graphname,gmoc2);
    //ffemockup2->GetObject(graphnameempty,gmoc2);
    gmoc2->SetTitle("Mockup block LSD Profile");
    gmoc2->GetXaxis()->SetTitle("z [cm]");
    gmoc2->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    gmoc2->GetYaxis()->SetTitleOffset(1.1);
    gmoc2->SetMinimum(0.);
    gmoc2->SetMaximum(15.);
    gmoc2->SetMarkerStyle(20);
    gmoc2->SetMarkerSize(0.8);
    gmoc2->SetMarkerColor(9);
    gmoc2->SetLineColor(9);
    gmoc2->SetLineStyle(9);
    Double_t xm2=0;
    Double_t ym2=0;
    gmoc2->GetPoint(gmoc2->GetN()-1,xm2,ym2);
    gmoc2->RemovePoint(gmoc2->GetN()-1);
    TLine * lXm2 = new TLine(0,ym2,60,ym2);
    lXm2->SetLineColor(9);
    lXm2->SetLineWidth(3.);
    lXm2->SetLineStyle(7);
    TMultiGraph * mgmoc = new TMultiGraph();
    mgmoc->Add(gmoc);
    mgmoc->Add(gmoc2);
    mgmoc->Draw("ALP");
    mgmoc->SetTitle("Mockup block LSD Profile");
    mgmoc->GetXaxis()->SetTitle("z [cm]");
    mgmoc->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    mgmoc->GetYaxis()->SetTitleOffset(1.1);
    mgmoc->SetMinimum(0.);
    mgmoc->SetMaximum(15.);
    //gmoc2->Draw("ALP");
    lXm1->Draw("same");
    lXm2->Draw("same");
    pt->Draw();
    /// mean
    float norm1 = lsdMean(gmoc,6);
    std::cout << "lsd mockup 1 mean value over 6 bins is " << norm1 <<
                 " [rad2/m], " << norm1/100./schultzToInvLrad  << " schultz " << std::endl;
    std::cout << "Entire mockup 1 mean lsd value is " << ym1 << std::endl;
    float norm2 = lsdMean(gmoc2,6);
    std::cout << "lsd mockup 2 mean value over 6 bins is " << norm2 <<
                 " [rad2/m], " << norm2/100./schultzToInvLrad  << " schultz " << std::endl;
    std::cout << "Entire mockup 2 mean lsd value is " << ym2 << std::endl;
    //norm = (norm1+norm2)/2; // mean from the slices
    norm = (ym1+ym2)/2;  // mean from the entire object
    std::cout << "Normalization factor from whole mockup object " << norm << std::endl;

    /// normalize sample graphs
    fcarota->GetObject(graphname,gcarota);
    float scale = 11.7/norm;
    for (int i=0; i<gcarota->GetN();i++){
        gcarota->GetY()[i] *= scale;
        float ye = gcarota->GetErrorY(i);
        ye *= scale;
        gcarota->SetPointError(i,0.,ye);
    }
    for (int i=0; i<gempty->GetN();i++){
        gempty->GetY()[i] *= scale;
        float ye = gempty->GetErrorY(i);
        ye *= scale;
        gempty->SetPointError(i,0.,ye);
    }

    TCanvas *cs = new TCanvas("cs","cs",800,600);
    cs->cd(1);
    TMultiGraph * mgcarote = new TMultiGraph();
    mgcarote->Add(gcarota);
    mgcarote->Add(gempty);
    mgcarote->Draw("ALP");
    mgcarote->SetTitle("Filled and Empty CDS normalized LSD Profile");
    mgcarote->SetMinimum(0.);
    mgcarote->SetMaximum(15.);
    gcarota->SetMarkerStyle(20);
    gcarota->SetMarkerSize(0.8);
    gcarota->SetMarkerColor(2);
    gcarota->SetLineColor(2);
    // all sample
    Double_t x1=0;
    Double_t totlsdcarota=0;
    gcarota->GetPoint(gcarota->GetN()-1,x1,totlsdcarota);
    gcarota->RemovePoint(gcarota->GetN()-1);
    TLine * lX = new TLine(x1-25.,totlsdcarota,x1+25.,totlsdcarota);
    lX->SetLineColor(2);
    lX->SetLineWidth(3.);
    lX->SetLineStyle(7);
    lX->Draw("same");
    // only for mockup sample
    Double_t xr3=0;
    Double_t lsdIIIreg=0;
    Double_t xr2=0;
    Double_t lsdIIreg=0;
    Double_t xr1=0;
    Double_t lsdIreg=0;
    if(run==2211 || run==22092211){
        // III region mockup sample
        gcarota->GetPoint(gcarota->GetN()-1,xr3,lsdIIIreg);
        gcarota->RemovePoint(gcarota->GetN()-1);
        TLine * lXIIIr = new TLine(xr3-9,lsdIIIreg,xr3+9,lsdIIIreg);
        lXIIIr->SetLineColor(6);
        lXIIIr->SetLineWidth(3.);
        lXIIIr->SetLineStyle(7);
        lXIIIr->Draw("same");

        // II region mockup sample
        gcarota->GetPoint(gcarota->GetN()-1,xr2,lsdIIreg);
        gcarota->RemovePoint(gcarota->GetN()-1);
        TLine * lXIIr = new TLine(xr2-4,lsdIIreg,xr2+4,lsdIIreg);
        lXIIr->SetLineColor(4);
        lXIIr->SetLineWidth(3.);
        lXIIr->SetLineStyle(7);
        lXIIr->Draw("same");

        // I region mockup sample
        gcarota->GetPoint(gcarota->GetN()-1,xr1,lsdIreg);
        gcarota->RemovePoint(gcarota->GetN()-1);
        TLine * lXIr = new TLine(xr1-9,lsdIreg,xr1+9,lsdIreg);
        lXIr->SetLineColor(3);
        lXIr->SetLineWidth(3.);
        lXIr->SetLineStyle(7);
        lXIr->Draw("same");
    }
    gempty->SetMarkerStyle(20);
    gempty->SetMarkerSize(0.8);
    gempty->SetMarkerColor(7);
    gempty->SetLineColor(7);
    Double_t x2=0;
    Double_t totlsdempty=0;
    gempty->GetPoint(gempty->GetN()-1,x2,totlsdempty);
    gempty->RemovePoint(gempty->GetN()-1);
    TLine * lXe = new TLine(15,totlsdempty,55,totlsdempty);
    lXe->SetLineColor(7);
    lXe->SetLineWidth(3.);
    lXe->SetLineStyle(7);
    lXe->Draw("same");
    mgcarote->GetXaxis()->SetTitle("z [cm]");
    mgcarote->GetYaxis()->SetTitle("LSD * 11.7/LSD_{mochup} [rad^{2}/m]");
    mgcarote->GetYaxis()->SetTitleOffset(1.1);
    ptr->Draw();

    /// CDS mean
//    std::cout << "lsd CDS mean value is " << lsdMean(gcarota,18) <<
//                 " [rad2/m], " << lsdMean(gcarota,18)/100./schultzToInvLrad  << " schultz " << std::endl;
    std::cout << "Entire CDS mean lsd value is " << totlsdcarota << std::endl;

    //// final material profile
    // subtract empty CDS
    for(int ip=0; ip<gcarota->GetN(); ip++){
        Double_t x1=0;
        Double_t y1=0;
        gcarota->GetPoint(ip,x1,y1);
        float yec = gcarota->GetErrorY(ip);
        Double_t x4=0;
        Double_t y4=0;
        gempty->GetPoint(ip,x4,y4);
        float yee = gempty->GetErrorY(ip);

        //compute final error sum in quadrature
        float ye = sqrt(yec*yec + yee*yee);

        float y1diff = y1;//-y4;
        mgfinal->SetPoint(ip,x1,y1diff);
        mgfinal->SetPointError(ip,0.,ye);
    }

    TCanvas *cf = new TCanvas("cf","cf",800,600);
    mgfinal->Draw("ALP");
    mgfinal->SetTitle("CDS Material normalized LSD Profile");
    mgfinal->SetMinimum(-0.5);
    mgfinal->SetMaximum(15.);
    mgfinal->SetMarkerStyle(20);
    mgfinal->SetMarkerSize(0.8);
    mgfinal->SetMarkerColor(6);
    mgfinal->SetLineColor(6);
    mgfinal->SetLineWidth(2.);
    mgfinal->GetXaxis()->SetTitle("z [cm]");
    mgfinal->GetYaxis()->SetTitle("LSD * 11.7/LSD_{mochup} [rad^{2}/m]");
    mgfinal->GetYaxis()->SetTitleOffset(1.1);
    ptr->Draw();

    if(run==2211 || run==22092211){
        // III region mockup sample
        TLine * lXIIIr = new TLine(xr3-9,lsdIIIreg-totlsdempty,xr3+9,lsdIIIreg-totlsdempty);
        lXIIIr->SetLineColor(6);
        lXIIIr->SetLineWidth(3.);
        lXIIIr->SetLineStyle(7);
        lXIIIr->Draw("same");

        // II region mockup sample
        TLine * lXIIr = new TLine(xr2-4,lsdIIreg-totlsdempty,xr2+4,lsdIIreg-totlsdempty);
        lXIIr->SetLineColor(4);
        lXIIr->SetLineWidth(3.);
        lXIIr->SetLineStyle(7);
        lXIIr->Draw("same");

        // I region mockup sample
        TLine * lXIr = new TLine(xr1-9,lsdIreg-totlsdempty,xr1+9,lsdIreg-totlsdempty);
        lXIr->SetLineColor(3);
        lXIr->SetLineWidth(3.);
        lXIr->SetLineStyle(7);
        lXIr->Draw("same");
    }
//    /// printing
//    char outname[100];
//    sprintf(outname, "results0.95/%i/%s/%i_sample_profile_%s.png",run,dir,run,dir);
//    cs->Print(outname);
////    sprintf(outname, "results0.95/%i/%s/%i_iron_profile_%s.png",run,dir,run,dir);
////    cfe->Print(outname);
////    sprintf(outname, "results0.95/%i/%s/%i_mockup_profile_%s.png",run,dir,run,dir);
////    cm->Print(outname);
//    sprintf(outname, "results0.95/%i/%s/%i_norm_profile_%s.png",run,dir,run,dir);
//    cf->Print(outname);

//    // print all in one directory too
//    sprintf(outname, "results0.95/allCDS/%i_sample_profile_%s.png",run,dir);
//    cs->Print(outname);
//    sprintf(outname, "results0.95/allCDS/%i_norm_profile_%s.png",run,dir);
//    cf->Print(outname);

    /// compute barycentrum
//    float bfLNL = 0;
//    float bcLNL = 0;
//    barycenterLNL(run,bfLNL,bcLNL);
    bf = barycenter(gcarota);
    bc = barycenter(mgfinal);
    std::cout << "Computed barycenters:" << std::endl;
    std::cout << " ---> filled sample = " << bf << std::endl;
    std::cout << " ---> empty sample = " << barycenter(gempty) << std::endl;
    std::cout << " ---> sample content = " << bc << std::endl;

    std::cout << "Total lsd:" << std::endl;
    std::cout << " ---> filled sample = " << totlsdcarota << std::endl;
    std::cout << " ---> empty sample = " << totlsdempty << std::endl;

    totlsd = totlsdcarota;

    return;
}

void barycenterLNL(int ir, float & rbf, float & rbc){

    int nruns = 14;
    int run[14] = {2183,2184,2186,2187,2188,2190,2199,2200,2201,2202,2203,2204,2206,2207};
    int flag[14] = {204, 205, 206, 207, 101, 102, 103, 104, 105, 106, 107, 201, 202, 203};

    // find run id index
    int * p;
    p = std::find(run, run+nruns, ir);
    int ind;
    if (p != run+nruns) ind = p-run;
    else return;

    // mass and barycenter of filled sample
    float Mf[14] = {14.48,15.36,11.99,12.66,9.10,9.44,13.42,14.80,10.58,13.48,11.20,11.66,10.80,12.56};
    float bf[14] = {25.0,25.5,30.0,23.5,32.5,30.0,23.5,27.0,29.0,29.0,24.5,26.0,29.5,25.0};
    rbf = bf[ind];
    // mass and barycenter of empty sample
    float Me = 4.485;
    float be = 32.7;
    // barycenter of sample content
    rbc = (bf[ind]*Mf[ind] - be*Me) / (Mf[ind] - Me);

    std::cout << "\nPrint LNL information about run " << ir << std::endl;
    std::cout << "-------------------------------------" << std::endl;
    std::cout << "Sample identification flag = m" << flag[ind] << std::endl;
    std::cout << "Mass = " << Mf[ind] << std::endl;
    std::cout << "Measured barycenters:" << std::endl;
    std::cout << " ---> filled sample = " << rbf << std::endl;
    std::cout << " ---> empty sample = " << be << std::endl;
    std::cout << " ---> sample content = " << rbc << std::endl;

    return;
}

float barycenter(TGraphErrors * gr){
    float z_calib = 11.75;

    float sum_zd = 0;
    float sum_d = 0;

    for(int ip=0; ip<gr->GetN(); ip++){
        Double_t z=0;
        Double_t d=0;
        gr->GetPoint(ip,z,d);

        sum_zd += z*d;
        sum_d  += d;
    }

    float zcm = sum_zd / sum_d;
    return zcm - z_calib;
}

void loopProfiles(){
    // all runs
//    int nruns = 14;
//    int run[14] = {2183,2184,2186,2187,2188,2190,2199,2200,2201,2202,2203,2204,2206,2207};
    int flag[14] = {204, 205, 206, 207, 101, 102, 103, 104, 105, 106, 107, 201, 202, 203};

    // mass of filled sample
    float Mf[14] = {14.48,15.36,11.99,12.66,9.10,9.44,13.42,14.80,10.58,13.48,11.20,11.66,10.80,12.56};
    // mass, barycenter. total lsd of empty sample
    float Me = 4.485;
    float be = 32.7;
    float lsdempty = 3.57608;

    float lsdsample[14] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    float lsd_mockup[14] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    // for barycenter computation
    float bf[14] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    float bc[14] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    float bfLNL[14] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    float bcLNL[14] = {0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.};
    
//    int nruns = 2;
//    int run[2] = {2211,2211};
//    int nruns = 1;
//    int run[1] = {22092211};
    int nruns = 1;
    int run[1] = {2202};

    // for barycenter computation
    float bf[1] = {0.};
    float bc[1] = {0.};
    float bfLNL[1] = {0.};
    float bcLNL[1] = {0.};


    TMultiGraph * mgruns = new TMultiGraph();

    TGraphErrors * graph = new TGraphErrors();
    for(int ir=0; ir<nruns; ir++){

        const char * dir = "normZ_5_10";
//        // only for testing voxel mask
//        if(ir==0)
//            dir = "normZ_5_10_noMask";
//        if(ir==1)
//            dir = "normZ_5_10_voxMask";

        std::cout << "Drawing run " << run[ir] << endl;
        // nominal baricentrum
        barycenterLNL(run[ir],bfLNL[ir],bcLNL[ir]);
        // profile and computed barycentrum
        profileZ(run[ir],dir,graph,lsd_mockup[ir],bf[ir],bc[ir],lsdsample[ir]);
        graph->SetMarkerColor(ir+2);
        graph->SetLineColor(ir+2);
        mgruns->Add((TGraph*)(graph->Clone()));
    }
    // mochup lsd mean and rms
    float mean = 0;
    float rms = 0;
    for(int im=0; im<nruns; im++)
        mean += lsd_mockup[im];
    mean /= nruns;
    for(int im=0; im<nruns; im++)
        rms += (lsd_mockup[im] - mean)*(lsd_mockup[im] - mean);
    rms = sqrt(rms / (nruns));
    std::cout << "Mean mockup LSD " << mean
              << " and rms " << rms << " on " << nruns << " samples" << endl;


    TCanvas *cr = new TCanvas("cr","cr",800,600);
//    graph->Draw("ALP");
    mgruns->Draw("ALP");
//    mgruns->SetTitle("Normalized CDS Material LSD Profile");
    mgruns->SetTitle("Normalized Filled CDS LSD Profile");
    mgruns->SetMinimum(-0.5);
    mgruns->SetMaximum(15.);
//    mgruns->SetMarkerStyle(20);
//    mgruns->SetMarkerSize(0.8);
    mgruns->GetXaxis()->SetTitle("z [cm]");
    mgruns->GetYaxis()->SetTitle("LSD * 11.7/LSD_{mochup} [rad^{2}/m]");
    mgruns->GetYaxis()->SetTitleOffset(1.5);

    // label
    TPaveText *pt = new TPaveText(.7,.7,.99,.95,"brNDC");
    pt->SetFillColor(kWhite);
    pt->SetTextAlign(12);
    pt->SetBorderSize(1);
    //pt->AddText("#color[1]{LKAB CDSs analysis}");
    pt->AddText("#color[1]{Mockup CDS analysis}");
    //pt->AddText("#bf{NO Sij guess}");
    pt->AddText("#bf{Sij cut N=30}");
    pt->AddText("#bf{5000 iterations}");
    pt->AddText("#bf{vox 2.5 cm, 460 min}");
    pt->AddText("#color[2]{Init with real image}");
    pt->AddText("#color[3]{Air voxel mask and blocks real image}");
    //pt->AddText("#bf{Fe=10, Femockup=2}");
    //pt->AddText("#bf{CDS=2, air=0.07}");
    pt->Draw();

//    /// printing
//    char outname[100];
//    sprintf(outname, "results0.95/allCDS/CDSs_norm_profile_%s.png",dir);
//    cr->Print(outname);

    /// barycentrum graph
    TMultiGraph * mgbar = new TMultiGraph();
    TGraphErrors * grbf = new TGraphErrors(nruns,bfLNL,bf);
    grbf->SetMarkerStyle(20);
    grbf->SetMarkerSize(1);
    grbf->SetMarkerColor(2);
    grbf->SetLineColor(2);
    grbf->SetLineStyle(9);
    TGraphErrors * grbc = new TGraphErrors(nruns,bcLNL,bc);
    grbc->SetMarkerStyle(23);
    grbc->SetMarkerSize(1);
    grbc->SetMarkerColor(9);
    grbc->SetLineColor(9);
    grbc->SetLineStyle(9);

    mgbar->Add(grbf);
    mgbar->Add(grbc);

    TCanvas *cb = new TCanvas("cb","cb",800,600);
    mgbar->Draw("AP");
    mgbar->SetTitle("Measured vs Computed barycentrum");
    mgbar->SetMinimum(18);
    mgbar->SetMaximum(34);
    mgbar->GetXaxis()->SetTitle("Measured barycenter z [cm]");
    mgbar->GetYaxis()->SetTitle("Computed barycenter z [cm]");
    mgbar->GetYaxis()->SetTitleOffset(1.1);
    mgbar->GetXaxis()->SetRangeUser(0.,40.);

    TLine * l = new TLine(18,18,33,33);
    l->SetLineColor(1);
    l->SetLineWidth(1.);
    l->SetLineStyle(9);
    l->Draw("same");

    // label
    TPaveText *ptb = new TPaveText(.8,.8,.99,.95,"brNDC");
    ptb->SetFillColor(kWhite);
    ptb->SetTextAlign(12);
    ptb->SetBorderSize(1);
    ptb->AddText("#bf{LKAB CDSs baycenters}");
    ptb->AddText("#bf{Sij cut N=30}");
    ptb->AddText("#bf{5000 iterations}");
    ptb->AddText("#bf{vox 2.5 cm, 460 min}");
    ptb->AddText("#bf{Init with real image}");
    ptb->AddText("#bf{Tolerance [cm] 5 10 5}");
    ptb->AddText("#color[2]{Filled sample}");
    ptb->AddText("#color[9]{Sample contents}");
    ptb->Draw("same");

    /// weight graph
    TGraphErrors * grw = new TGraphErrors(nruns,Mf,lsdsample);
    grw->SetMarkerStyle(20);
    grw->SetMarkerSize(1);
    grw->SetMarkerColor(2);
    grw->SetLineColor(2);
    grw->SetLineStyle(9);

    TCanvas *cw = new TCanvas("cw","cw",800,600);
    grw->Draw("AP");
    grw->SetTitle("LSD vs Measured weight");
//    grw->SetMinimum(18);
//    grw->SetMaximum(34);
    grw->GetXaxis()->SetTitle("Measured weight [kg]");
    grw->GetYaxis()->SetTitle("LSD * 11.7/LSD_{mochup} [rad^{2}/m]");
    grw->GetYaxis()->SetTitleOffset(1.1);
    //grw->GetXaxis()->SetRangeUser(0.,40.);

//    /// dump to file for excell table
//    // open output file
//    std::ofstream of;
//    of.open("barycenter.csv", std::ofstream::app);
//    of << "CDS sample; run; weight [kg]; LSD weight [kg]; contents weight [kg]; LSD contents weight [kg]; barycenter [cm]; LSD barycenter [cm]; contents barycenter [cm]; LSD contents barycenter [cm];" << std::endl;
//    for(int is=0; is<nruns; is++)
//        of << "m" << flag[is] << ";"
//           << run[is] << ";"
//           << Mf[is] << ";"
//           << lsdsample[is] << ";"
//           << (Mf[is] - Me) << ";"
//           << (lsdsample[is] - lsdempty) << ";"
//           << bfLNL[is] << ";"
//           << bf[is] << ";"
//           << bcLNL[is] << ";"
//           << bc[is] << ";\n";
//    of << "EMPTY;"
//       << 2208 << ";"
//       << Me << ";"
//       << lsdempty << ";"
//       << 0 << ";"
//       << 0 << ";"
//       << be << ";"
//       << 29.38 << ";"
//       << 0 << ";"
//       << 0 << ";\n";

//    of.close();

    return;
}

void iterationAir(){
    char file[100];
    const char * dir = "iteration_air";
    int run = 2183;
    sprintf(file, "results0.95/%i/%s/lambda0.01_iterationsOutput.root",run,dir,run);
    TFile * fl01 = new TFile(file);
    sprintf(file, "results0.95/%i/%s/lambda0.07_iterationsOutput.root",run,dir,run);
    TFile * fl07 = new TFile(file);
    sprintf(file, "results0.95/%i/%s/lambda0.13_iterationsOutput.root",run,dir,run);
    TFile * fl13 = new TFile(file);
    sprintf(file, "results0.95/%i/%s/lambda0.2_iterationsOutput.root",run,dir,run);
    TFile * fl2 = new TFile(file);

    TCanvas *c = new TCanvas("c","c",1200,800);

    TMultiGraph * mg_it = new TMultiGraph("mg_it","Air Block LSD versus iteration");

    TGraphErrors * g01;
    fl01->cd();
    fl01->GetObject("airGraph",g01);
    mg_it->Add(g01);
    TGraphErrors * g07;
    fl07->GetObject("airGraph",g07);
    mg_it->Add(g07);
    TGraphErrors * g13;
    fl13->GetObject("airGraph",g13);
    mg_it->Add(g13);
    TGraphErrors * g2;
    fl2->GetObject("airGraph",g2);
    mg_it->Add(g2);



    mg_it->Draw("ALP");
    mg_it->GetXaxis()->SetTitle("Iterations");
    mg_it->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    mg_it->GetYaxis()->SetTitleOffset(1.5);
    mg_it->SetMinimum(0.);
    mg_it->SetMaximum(0.2);

    // lambda 0.01
    g01->SetMarkerStyle(20);
    g01->SetMarkerSize(0.8);
    g01->SetMarkerColor(1);
    g01->SetLineColor(1);
    g01->GetXaxis()->SetTitle("Iterations");
    g01->GetYaxis()->SetTitle("1/L_{rad} [1/cm]");
    g01->GetYaxis()->SetTitleOffset(1.5);

    // lambda 0.07
    g07->SetMarkerStyle(20);
    g07->SetMarkerSize(0.8);
    g07->SetMarkerColor(2);
    g07->SetLineColor(2);
    g07->GetXaxis()->SetTitle("Iterations");
    g07->GetYaxis()->SetTitle("1/L_{rad} [1/cm]");
    g07->GetYaxis()->SetTitleOffset(1.5);

    // lambda 0.13
    g13->SetMarkerStyle(20);
    g13->SetMarkerSize(0.8);
    g13->SetMarkerColor(3);
    g13->SetLineColor(3);
    g13->GetXaxis()->SetTitle("Iterations");
    g13->GetYaxis()->SetTitle("LSD [rad^{2}/m]");
    g13->GetYaxis()->SetTitleOffset(1.5);

    // lambda 0.2
    g2->SetMarkerStyle(20);
    g2->SetMarkerSize(0.8);
    g2->SetMarkerColor(4);
    g2->SetLineColor(4);
    g2->GetXaxis()->SetTitle("Iterations");
    g2->GetYaxis()->SetTitle("1/L_{rad} [1/cm]");
    g2->GetYaxis()->SetTitleOffset(1.5);

    TPaveText *pt = new TPaveText(.7,.6,.95,.9,"brNDC");
    pt->SetFillColor(kWhite);
    pt->SetTextAlign(12);
    pt->SetBorderSize(1);
    pt->AddText("#color[2]{Run 2183}");
    pt->AddText("#bf{Image reconstruction}");
    pt->AddText("#bf{Sij cut N=30}");
    pt->AddText("#bf{vox 2.5 cm, 460 min}");
    pt->AddText("#bf{Init with real image}");
    pt->AddText("#color[1]{lambda_init_air = 0.01}");
    pt->AddText("#color[2]{lambda_init_air = 0.07}");
    pt->AddText("#color[3]{lambda_init_air = 0.13}");
    pt->AddText("#color[4]{lambda_init_air = 0.20}");

    pt->Draw();

    return;
}

void iteration() {

    TFile * f0 = new TFile("results/2183/2183_iterationsOutput.root");
//    TFile * f2 = new TFile("iterationsOutput.root");

    TCanvas *c = new TCanvas("c","c",1200,800);
    c->Divide(2,2);
    c->cd(1);

    f0->cd();

    TMultiGraph * mg_it = new TMultiGraph("mg_it","Density versus iteration");
    mg_it->Add(feGraph);
    mg_it->Add(feMockupGraph);
    mg_it->Add(carotaGraph);
    mg_it->Draw("ALP");
    mg_it->GetXaxis()->SetTitle("Iterations");
    mg_it->GetYaxis()->SetTitle("1/L_{rad} [1/cm]");
    mg_it->GetYaxis()->SetTitleOffset(1.5);
    mg_it->SetMinimum(0.);
    mg_it->SetMaximum(0.3);

    // Fe
    feGraph->SetTitle("Fe density");
    feGraph->SetMarkerStyle(20);
    feGraph->SetMarkerSize(0.8);
    feGraph->SetMarkerColor(1);
    feGraph->SetLineColor(1);
    feGraph->GetXaxis()->SetTitle("Iterations");
    feGraph->GetYaxis()->SetTitle("1/L_{rad} [1/cm]");
    feGraph->GetYaxis()->SetTitleOffset(1.5);
    // Femockup
    feMockupGraph->SetTitle("Fe mockup density");
    feMockupGraph->SetMarkerStyle(20);
    feMockupGraph->SetMarkerSize(0.8);
    feMockupGraph->SetMarkerColor(2);
    feMockupGraph->SetLineColor(2);
    feMockupGraph->GetXaxis()->SetTitle("Iterations");
    feMockupGraph->GetYaxis()->SetTitle("1/L_{rad} [1/cm]");
    feMockupGraph->GetYaxis()->SetTitleOffset(1.5);
    feMockupGraph->SetMinimum(0.065);
    feMockupGraph->SetMaximum(0.1);
    // carota slice 0
    carotaGraph->SetTitle("CDS density");
    carotaGraph->SetMarkerStyle(20);
    carotaGraph->SetMarkerSize(0.8);
    carotaGraph->SetMarkerColor(3);
    carotaGraph->SetLineColor(3);
    carotaGraph->GetXaxis()->SetTitle("Iterations");
    carotaGraph->GetYaxis()->SetTitle("1/L_{rad} [1/cm]");
    carotaGraph->GetYaxis()->SetTitleOffset(1.5);
    carotaGraph->SetMinimum(0.065);
    carotaGraph->SetMaximum(0.1);


    TPaveText *pt = new TPaveText(.6,.6,.95,.9,"brNDC");
    pt->SetFillColor(kWhite);
    pt->SetTextAlign(12);
    pt->SetBorderSize(1);
    pt->AddText("#color[2]{Run 2183}");
    pt->AddText("Density versus iteration analysis");
    pt->AddText("#bf{Image reconstruction}");
    //pt->AddText("#bf{Voxel reweight}");
    pt->AddText("#bf{NO Sij guess}");
    pt->AddText("#bf{Sij cut N=30}");
    pt->AddText("#bf{40000 iterations}");
    pt->AddText("#bf{vox 2.5 cm, 480 min}");
    pt->AddText("#bf{Init image:}");
    pt->AddText("#bf{Fe=10, Femockup=2}");
    pt->AddText("#bf{CDS=2, air=0.07}");

    c->cd(2);
    feGraph->Draw("ALP");
    c->cd(3);
    feMockupGraph->Draw("ALP");
    c->cd(4);
    carotaGraph->Draw("ALP");
    pt->Draw();

    /// ratio plots
    f0->cd();
    int N = feMockupGraph->GetN();
    TGraph * feOverMockupGraph = new TGraph(N);
    feOverMockupGraph->SetTitle("Convergence of LSD Fe / LSD mockup ");
    feOverMockupGraph->SetMarkerStyle(20);
    feOverMockupGraph->SetMarkerSize(0.8);
    feOverMockupGraph->SetMarkerColor(4);
    feOverMockupGraph->SetLineColor(4);
    feOverMockupGraph->SetMinimum(3);
    feOverMockupGraph->SetMaximum(4.);

    TGraph * carotaOverMockupGraph = new TGraph(N);
    carotaOverMockupGraph->SetTitle("Convergence of LSD CDS / LSD mockup ");
    carotaOverMockupGraph->SetMarkerStyle(20);
    carotaOverMockupGraph->SetMarkerSize(0.8);
    carotaOverMockupGraph->SetMarkerColor(6);
    carotaOverMockupGraph->SetLineColor(6);
    carotaOverMockupGraph->SetMinimum(0.8);
    carotaOverMockupGraph->SetMaximum(1.);

    TGraph * carotaOverMockupGraph2 = new TGraph(N);
    carotaOverMockupGraph2->SetTitle("CDS/Femockup density");
    carotaOverMockupGraph2->SetMarkerStyle(20);
    carotaOverMockupGraph2->SetMarkerSize(0.8);
    carotaOverMockupGraph2->SetMarkerColor(7);
    carotaOverMockupGraph2->SetLineColor(7);
    carotaOverMockupGraph2->GetXaxis()->SetTitle("Iterations");
    carotaOverMockupGraph2->GetYaxis()->SetTitle("1/L_{rad} [1/cm]");
    carotaOverMockupGraph2->GetYaxis()->SetTitleOffset(1.5);
    carotaOverMockupGraph2->SetMinimum(0.8);
    carotaOverMockupGraph2->SetMaximum(1.);

    TMultiGraph * mg_it_carota = new TMultiGraph("mg_it_carota","Density versus iteration");

    for(int ip=0; ip<N; ip++){
        Double_t xm=0;
        Double_t ym=0;
        feMockupGraph->GetPoint(ip,xm,ym);
        Double_t xf=0;
        Double_t yf=0;
        feGraph->GetPoint(ip,xf,yf);
        f0->cd();
        Double_t xc=0;
        Double_t yc=0;
        carotaGraph->GetPoint(ip,xc,yc);

//        float y1diff = y1-y4;
//        float y2diff = y2-y4;
        feOverMockupGraph->SetPoint(ip,xm,yf/ym);
        carotaOverMockupGraph->SetPoint(ip,xm,yc/ym);

//        f2->cd();
//        Double_t xc2=0;
//        Double_t yc2=0;
//        carotaGraph->GetPoint(ip,xc2,yc2);
//        carotaOverMockupGraph2->SetPoint(ip,xm,yc2/ym);
    }
    feOverMockupGraph->GetXaxis()->SetTitle("Iterations");
    //feOverMockupGraph->GetYaxis()->SetTitle("1/L_{rad} [1/cm]");
    feOverMockupGraph->GetYaxis()->SetTitle("LSD_{Fe}/LSD_{mockup}");
    feOverMockupGraph->GetYaxis()->SetTitleOffset(1.2);
    carotaOverMockupGraph->GetXaxis()->SetTitle("Iterations");
    //carotaOverMockupGraph->GetYaxis()->SetTitle("1/L_{rad} [1/cm]");
    carotaOverMockupGraph->GetYaxis()->SetTitle("LSD_{CDS}/LSD_{mockup}");
    carotaOverMockupGraph->GetYaxis()->SetTitleOffset(1.5);

    mg_it_carota->Add(carotaOverMockupGraph);
//    mg_it_carota->Add(carotaOverMockupGraph2);

    TCanvas *cr = new TCanvas("cr","cr",1500,800);
    cr->Divide(2,1);
    cr->cd(1);
    gPad->SetGrid();
    feOverMockupGraph->Draw("ALP");
    cr->cd(2);
    gPad->SetGrid();
    carotaOverMockupGraph->Draw("ALP");
//    mg_it_carota->Draw("ALP");
    pt->Draw();

    mg_it_carota->SetTitle("CDS/Femockup density");
    mg_it_carota->GetXaxis()->SetTitle("Iterations");
    mg_it_carota->GetYaxis()->SetTitle("1/L_{rad} [1/cm]");
    mg_it_carota->GetYaxis()->SetTitleOffset(1.5);
    cr->Update();

    return;
}

void testPattRec(int run){

    gROOT->Reset();
    gROOT->SetStyle("Plain");
    //gStyle->SetErrorX(0);
    //gStyle->SetOptStat(1);
    gStyle->SetOptFit(1);
    gStyle->SetStatW(0.3);
    gStyle->SetStatH(0.2);

    /// open file
    char file[100];
    sprintf(file, "/mnt/mutom-gluster/data/mublast/fit/pattrecLNLdata/carote/Radmufit_r%i_0_100000000ev_PR.root",run);
    TFile *_file0 = new TFile(file);
    //TFile *_file0 =
    //TFile::Open("/mnt/mu-tom1/data/tom_data/PattRec/RootFile/Radmufit_r2151_0_25000000ev_PR.root");
    //TFile::Open("/mnt/mu-tom1/data/software/mutom/PattRec/PattRec/OUTPUT/Radmufit_r2150_0_250000ev_PR.root");
    //TFile::Open("/mnt/mutom-gluster/data/mublast/fit/pattrecLNLdata/carote/Radmufit_r2209_0_100000000ev_PR.root");
    //TFile::Open("/mnt/mutom-gluster/data/mublast/fit/pattrecLNLdata/carote/Radmufit_r2218_0_100000000ev_PR.root");
    TTree * T = (TTree*)_file0->Get("RADMU");

    ////book histograms
    ///DTBX
    // occupancy
    TH1F * hocc[24];
    TH1F * hocc_firstHit[24];
    // nhits
    TH1F * hnh[24];
    // time
    TH1F * ht[24];
    TH1F * ht_firstHit[24];

    for(int ic=0; ic<2; ic++){
        for(int il=0; il<12; il++){
            const char hname[100];
            sprintf(hname,"hocc_ch%i_lay%i",ic,il);
            hocc[ic*12 + il] = new TH1F(hname,"Occupancy",100,0.,100.);
            sprintf(hname,"hocc_firstHit_ch%i_lay%i",ic,il);
            hocc_firstHit[ic*12 + il] = new TH1F(hname,"Occupancy first hit",100,0.,100.);

            sprintf(hname,"hnh_ch%i_lay%i",ic,il);
            hnh[ic*12 + il] = new TH1F(hname,"Number of hits in channel",100,0.,100.);

            sprintf(hname,"htime_ch%i_lay%i",ic,il);
            ht[ic*12 + il] = new TH1F(hname,"TDC time in channel",150,0.,1500.);
            sprintf(hname,"htime_firstHit_ch%i_lay%i",ic,il);
            ht_firstHit[ic*12 + il] = new TH1F(hname,"TDC time first hit in channel",150,0.,1500.);
        }
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
    cout << "Entries " << evnum << endl;

    //loop over the events
    for (Int_t ientry=0 ; ientry < TMath::Min(1000000,evnum); ientry ++)
    {
        //read tree
        T->GetEntry(ientry);

        //DTBX
        //reset number of hits
        for(int ic=0; ic<2; ic++)
            for(int il=0; il<12; il++)
                hnh[ic*12 + il]->Reset();

        for (int ih=0; ih<ihit; ih++){
            int ch = int(hlay[ih]/1000);
            int lay = (hlay[ih] % 1000) % 100;

            if(ch<3){
                int ind = (ch-1)*12 + (lay-1);

                hocc[ind]->Fill(htube[ih]);

                if(hnh[ind]->GetBinContent(htube[ih])==0){
                    hocc_firstHit[ind]->Fill(htube[ih]);
                    ht_firstHit[ind]->Fill(htime[ih]);
                }

                hnh[ind]->SetBinContent(htube[ih],hnh[ind]->GetBinContent(htube[ih])+1);
                ht[ind]->Fill(htime[ih]);
            }
        }

        //SEG

    }// end event loop

    // draw
    TPaveText *tpt = new TPaveText(0.8,0.35,0.99,0.65,"brNDC");
    tpt->SetFillColor(kWhite);
    tpt->SetTextAlign(12);
    tpt->SetBorderSize(1);
    char txt[20];
    sprintf(txt,"Run %i",run);
    tpt->AddText(txt);
    tpt->AddText("#color[2]{First hit}");
    tpt->AddText("#color[1]{All hits}");

    TCanvas * cHitsCh1 = new TCanvas("cHitsCh1","cHitsCh1",1300,900);
    cHitsCh1->Divide(4,3);
    for(int il=0; il<12; il++){
        cHitsCh1->cd(1+il);
        hocc[il]->Draw();
        hocc_firstHit[il]->SetLineColor(kRed);
        hocc_firstHit[il]->Draw("same");
    }
    tpt->Draw();

    TCanvas * cHitsCh2 = new TCanvas("cHitsCh2","cHitsCh2",1300,900);
    cHitsCh2->Divide(4,3);
    for(int il=0; il<12; il++){
        cHitsCh2->cd(1+il);
        hocc[il+12]->Draw();
        hocc_firstHit[il+12]->SetLineColor(kRed);
        hocc_firstHit[il+12]->Draw("same");
    }
    tpt->Draw();

    // draw daq times
    TCanvas * cTimeCh1 = new TCanvas("cTimeCh1","cTimeCh1",1300,900);
    cTimeCh1->Divide(4,3);
    for(int il=0; il<12; il++){
        cTimeCh1->cd(1+il);
        ht[il]->Draw();
        ht_firstHit[il]->SetLineColor(kRed);
        ht_firstHit[il]->Draw("same");
    }
    tpt->Draw();

    TCanvas * cTimeCh2 = new TCanvas("cTimeCh2","cTimeCh2",1300,900);
    cTimeCh2->Divide(4,3);
    for(int il=0; il<12; il++){
        cTimeCh2->cd(1+il);
        ht[il+12]->Draw();
        ht_firstHit[il+12]->SetLineColor(kRed);
        ht_firstHit[il+12]->Draw("same");
    }
    tpt->Draw();

    return;
}


////////////////////////////////////////////////////////////////////////////////////////////////////77
/// RAW MONITOR
//for debugging
static bool DEBUG_FLAG = 0;
//define map array and function to fill it
static int nros =  1;
static int nrob =  19;
static int ntdc =  4;
static int ncha = 32;
static map<int, int> chmap;

void fillMap()
{

    ofstream myfile;
    myfile.open ("out.txt");
    if(DEBUG_FLAG)
        myfile << "reading map" << endl;

    ifstream file("legnaro2ROS25.txt");

    int ros;
    int rob;
    int tdc;
    int cha;

    int se;
    int sl;
    int lay;
    int tube;

    int ddu;
    int wh;
    int st;

    if(DEBUG_FLAG)
        cout << "map loop ... " << endl;
    while( file >> ddu >> ros >> rob >> tdc >> cha >> wh >> st >> se >> sl >> lay >> tube )
    {
        if(DEBUG_FLAG) myfile << "fillMap "
                              << ddu << " " << ros << " " << rob << " " << tdc << " " << cha  << " => "
                              << wh << " " << st << " " << se << " " << " " <<  sl << " " << lay << " " << tube << endl;

        int idTDC = getTDCid(ros,rob,tdc,cha);
        int idTube = getTubeId(se,sl,lay,tube);

        chmap.insert( std::make_pair( idTDC, idTube ) );

    }
    if(DEBUG_FLAG)
        myfile << "TDC map read" << endl;
    return;
}

void checkMap()
{

    ofstream myfile;
    myfile.open ("check.txt");
    if(DEBUG_FLAG)
        myfile << "checking map" << endl;

    ifstream file("legnaro2ROS25.txt");

    int ros;
    int rob;
    int tdc;
    int cha;

    int se;
    int sl;
    int lay;
    int tube;

    int ddu;
    int wh;
    int st;

    map<int,int>::iterator iter;

    if(DEBUG_FLAG)
        cout << "map loop ... " << endl;
    while( file >> ddu >> ros >> rob >> tdc >> cha >> wh >> st >> se >> sl >> lay >> tube )
    {
        int rawdata = getTDCid(1,rob,tdc,cha);
        //cout << "TDCid " << rawdata << endl;

        iter = chmap.find(rawdata);
        if( iter != chmap.end() )
        {
            if(DEBUG_FLAG)
            {
                cout << "SE " << getSe(iter->second) <<
                        " sl " << getSL(iter->second) <<
                        " lay " << getLay(iter->second) <<
                        " tube " << getTube(iter->second) << endl;
            }

            if(DEBUG_FLAG)
                myfile << "fillMap " << ddu << " " << ros << " " << rob << " " <<
                          tdc << " " << cha  << " => " << wh << " " << st << " " <<
                          getSe(iter->second) << " " << " " <<  getSL(iter->second) << " " <<
                          getLay(iter->second) << " " << getTube(iter->second) << endl;

        }
    }
    if(DEBUG_FLAG)
        myfile << "end map check" << endl;
    return;
}


int getTDCid( int ros, int rob, int tdc, int cha ) {
    return ( ( ( ( ( ros * nrob ) + rob ) * ntdc ) + tdc ) * ncha ) + cha;
}

int getTubeId( int se, int sl, int lay, int tube ) {
    return (se * 10000 + sl * 1000 + lay * 100 + tube);
}

int getTube(int idTube) {
    return (idTube - int(idTube/100.)*100.);
}

int getLay(int idTube) {
    return (int(idTube/100.) - int(idTube/1000.)*10);
}

int getSL(int idTube) {
    return (int(idTube/1000.) - int(idTube/10000.)*10);
}

int getSe(int idTube) {
    return int(idTube/10000.);
}

float testRawData(int run, int maxwords){
    // 20150414 SV readapted from radmuMonitorLNL 2008/01/28
    //FILE is a C-type file pointer included in stdio.h
    FILE *infile;
    //char fileName[200];
    //sprintf(fileName,"/data/radmu/flat/dataROS25/%s",runName);
    //infile = fopen(fileName,"r");
    //infile = fopen("/data/LNL/raw/MuBlast_phi_J1_J2_J3_K1_K2_K3_newrif_r2232.i0","r");
    //infile = fopen("/mnt/mutom-gluster/data/mublast/raw/MuBlast_phi_J1_J2_J3_K1_K2_K3_newrif_r2245.i9","r");
    infile = fopen("/mnt/mutom-gluster/data/mublast/raw/MuBlast_phi_A1_A2_A3_A4_newrif_r2288.i19","r");
    
    //occupancy histo cuts
    int minTime = 1200;
    int maxTime = 2200;

    //hit difference histo cuts
    int minTimeHit = 1400;
    int maxTimeHit = 2000;

    //graphic tools
    TCanvas * c1;
    TCanvas * ch1;
    TCanvas * ch2;
    TCanvas * sls;
    TCanvas * hits;

    TPad * pad1;
    TPad * pad2;
    TPad * pad3;
    TPad * pad_hits_ch1;
    TPad * pad_time_ch1;
    TPad * pad_hits_ch2;
    TPad * pad_time_ch2;
    TPad * pad_hits_sls;
    TPad * pad_time_sls;
    TPad * pad_hits_hits;

    //histos
    //TH1F * hocc_lay[32];
    TH1F * htime_sl[8];
    TH1F * hoccSL1_lay[4];
    TH1F * hoccSL2_lay[4];
    TH1F * hoccCH1_lay[12];
    TH1F * hoccCH2_lay[12];

    //TDC histograms
    TH2F * hocc  = new TH2F("hocc","Channel vs TDC",200,0,99,66,-1,32);
    TH1F * htime = new TH1F("htime","TimeBox",2000,1000,3000);
    TH2F * htimech = new TH2F("htimech","Channel vs time",2000,1000,3000,8000,0,4000);

    // nhits
    TH1F * hnh[24];
    TH1F * hhit_diff[24];
    // time
    TH1F * ht[24];
    TH1F * ht_firstHit[24];
    // time tubes ch0 layer 1
    TH1F * ht_ch0_lay1_tube[72];

    for(int ic=0; ic<2; ic++){
        for(int il=0; il<12; il++){
            const char hname[50];
            sprintf(hname,"htime_ch%i_lay%i",ic,il);
            ht[ic*12 + il] = new TH1F(hname,"TDC time in channel",200,1000,3000);
            sprintf(hname,"htime_ch%i_lay%i_firstHit",ic,il);
            ht_firstHit[ic*12 + il] = new TH1F(hname,"TDC time in channel",200,1000,3000);
            sprintf(hname,"hnh_ch%i_lay%i",ic,il);
            hnh[ic*12 + il] = new TH1F(hname,"Number of hits in channel",100,0.,100.);
            //second - first hit time histogram
            sprintf(hname,"hhit_diff_ch%i_lay%i",ic,il);
            hhit_diff[ic*12 + il] = new TH1F(hname,"TDC time_{II hit} - TDC time_{I hit} ",1250,0.,1250.);
        }
    }
    for(int it=0; it<72; it++){
        const char hname[50];
        sprintf(hname,"htime_ch0_lay1_tube%i",it);
        ht_ch0_lay1_tube[it] = new TH1F(hname,"TDC time in tubes for ch0 lay2",200,1000,3000);
    }

    //fill channel map
    fillMap();

//    if(DEBUG_FLAG)
//        checkMap();

    //map test
    /*
      int rawdata = getTDCid(1,12,0,4);
      int chdata = getTubeId(10,3,2,2);

      map<int,int>::iterator iter = chmap.find(rawdata);
      if( iter != chmap.end() ) {
            cout << "Key : " << iter->first << " entriy : " << iter->second << endl;
        cout << "se :" << getSe(iter->second) <<
            "SL :" << getSL(iter->second) <<
            "Lay :" << getLay(iter->second) <<
            "Tube :" << getTube(iter->second) << endl;

      }
    */
    /*  
    //build the MINICRATE canvas with pads
    c1 = new TCanvas("c1","RADMU monitor: TDC",1,10,800,800);
    c1->SetFillColor(10);
    //gStyle->SetFrameFillColor(18);
    c1->SetBorderSize(2);
    TPaveText * title = new TPaveText(.2,0.96,.8,.995);
    title->SetFillColor(33);
    title->AddText("RADMU MONITOR: TDC hits");
    title->Draw();
    pad1 = new TPad("pad1","The pad with hits",0.05,0.50,0.95,0.92,19,1,1);
    pad1->Draw();
    pad2 = new TPad("pad2","The pad with times",0.05,0.05,0.50,0.45,19,1,1);
    pad2->Draw();
    pad3 = new TPad("pad3","The pad with times",0.50,0.05,0.95,0.45,19,1,1);
    pad3->Draw();
    // pad1->cd();

    //Chamber 2 canvas
    ch2 = new TCanvas("ch2","RADMU monitor:CH 2",160,10,800,800);
    ch2->SetFillColor(10);
    ch2->SetBorderSize(1);

    TPaveText * title_ch2 = new TPaveText(.2,0.96,.8,.995);
    title_ch2->SetFillColor(33);
    title_ch2->AddText("RADMU MONITOR: chamber 2");
    title_ch2->Draw();

    pad_hits_ch2 = new TPad("pad_hits_ch2","The pad with hits",0.01,0.01,0.5,0.95,19,1,1);
    pad_hits_ch2->Divide(1,12);
    pad_hits_ch2->Draw();

    pad_time_ch2 = new TPad("pad_time_ch2","The pad with times",0.5,0.01,0.99,0.95,19,1,1);
    pad_time_ch2->Divide(1,3);
    pad_time_ch2->Draw();

    //Chamber 1 canvas
    ch1 = new TCanvas("ch1","RADMU monitor:CH 1",320,10,800,800);
    ch1->SetFillColor(10);
    ch1->SetBorderSize(1);

    TPaveText * title_ch1 = new TPaveText(.2,0.96,.8,.995);
    title_ch1->SetFillColor(33);
    title_ch1->AddText("RADMU MONITOR: chamber 1");
    title_ch1->Draw();

    pad_hits_ch1 = new TPad("pad_hits_ch1","The pad with hits",0.01,0.01,0.5,0.95,19,1,1);
    pad_hits_ch1->Divide(1,12);
    pad_hits_ch1->Draw();

    pad_time_ch1 = new TPad("pad_time_ch1","The pad with times",0.5,0.01,0.99,0.95,19,1,1);
    pad_time_ch1->Divide(1,3);
    pad_time_ch1->Draw();


    //SLs canvas
    sls = new TCanvas("sls","RADMU monitor:SLs",480,10,800,800);
    sls->SetFillColor(10);
    sls->SetBorderSize(1);

    TPaveText * title_sls = new TPaveText(.2,0.96,.8,.995);
    title_sls->SetFillColor(33);
    title_sls->AddText("RADMU MONITOR: SLs");
    title_sls->Draw();

    pad_hits_sls = new TPad("pad_hits_sls","The pad with hits",0.01,0.01,0.5,0.95,19,1,1);
    pad_hits_sls->Divide(1,12);
    pad_hits_sls->Draw();

    pad_time_sls = new TPad("pad_time_sls","The pad with times",0.5,0.01,0.99,0.95,19,1,1);
    pad_time_sls->Divide(1,3);
    pad_time_sls->Draw();

    //hit difference canvas
    hits = new TCanvas("hits","RADMU monitor: Hit Difference",640,10,800,800);
    hits->SetFillColor(10);
    hits->SetBorderSize(1);

    TPaveText * title_hits = new TPaveText(.2,0.96,.8,.995);
    title_hits->SetFillColor(33);
    title_hits->AddText("RADMU MONITOR: TDC-time Hit difference");
    title_hits->Draw();

    pad_hits_hits = new TPad("pad_hits_hits","The pad with hits time difference",0.01,0.01,0.95,0.95,19,1,1);
    pad_hits_hits->Divide(2,2);
    pad_hits_hits->Draw();
    */

    //occupancy histograms: all layers from bottom to top 4+4+12+12
    int nCell = 72;
    int i = 0;

    for(i=0; i<4; i++)
    {
        TString nameLayer( "hocc" );
        int layNum = i+1;
        nameLayer += "_sl1_lay_";
        nameLayer += layNum;

        hoccSL1_lay[   i] = new TH1F( nameLayer,
                                      nameLayer, nCell, 0.5, nCell + 0.5 );
    }

    for(i=0; i<4; i++)
    {
        TString nameLayer( "hocc" );
        int layNum = i+1;
        nameLayer += "_sl2_lay_";
        nameLayer += layNum;

        hoccSL2_lay[   i] = new TH1F( nameLayer,
                                      nameLayer, nCell, 0.5, nCell + 0.5 );
    }

    for(i=0; i<12; i++)
    {
        TString nameLayer( "hocc" );
        int layNum = i+1;
        nameLayer += "_ch1_lay_";
        nameLayer += layNum;

        hoccCH1_lay[   i] = new TH1F( nameLayer,
                                      nameLayer, nCell, 0.5, nCell + 0.5 );
    }

    for(i=0; i<12; i++)
    {
        TString nameLayer( "hocc" );
        int layNum = i+1;
        nameLayer += "_ch2_lay_";
        nameLayer += layNum;

        hoccCH2_lay[   i] = new TH1F( nameLayer,
                                      nameLayer, nCell, 0.5, nCell + 0.5 );
    }

    //time boxes in SLs histograms from bottom to top 1+1+3+3
    int histoWidth = 2000;
    int histoEdgeL = 1000;
    int histoEdgeH = 3000;

    for(int j=0; j<8; j++)
    {
        TString tname( "htime" );
        tname += "_sl_";
        tname += j + 1;

        htime_sl[  j] = new TH1F( tname,
                                  tname,
                                  histoWidth, histoEdgeL, histoEdgeH );
    }


    //long int is 32-bit=4-byte word
    long evbu;
    //NB size_t fread(void *ptr, size_t size, size_t nmemb, FILE *stream);
    //Reads data from the given stream into the array pointed to by ptr.
    //It reads nmemb number of elements of size size (in byte).
    //The total number of bytes read is (size*nmemb).
    //On success the number of elements read is returned.
    //On error or end-of-file the total number of elements successfully
    //read (which may be zero) is returned.

    //NB ftell(infile) returns the current file position. For binary stream,
    //then the value is the number of bytes from the beginning of the file.

    // hit array definition
    float hitArray[7000];

    //int maxwords = 100000;
    int numEvent = 0;

    //  int maxwords = 100000;
    //  if(DEBUG_FLAG)
    //	maxwords = 1000;
    {
        for(int words=1; words<maxwords; words++)
        {
            if(words/10000. == int(words/10000.))    
              cout << "Reading word " << words << endl;

            // get the 32-bit word
            if(fread(&evbu,4,1,infile)==0)
                break;

            // DMA: swap 16 lsb with 16 msb
            //long word = (evbu>>16) | (evbu<<16);
            /// SV 20150414 to function with mublast raw data....
            long word = 0x00000000;
            word = ((evbu>>16) & 0x0000FFFF) | ((evbu<<16) & 0xFFFF0000);

            // type of packet and TTCcounts
            long type = word & 0xFF000000;
            type = type >> 24;			// type of data packet: bit 24-31
            long TTCcount = word & 0x00FFFFFF;	// TTC count: bit 0-23

            // if event header reset hit array
            //cout << "Reading evbu ---> " << hex<<evbu<<dec << endl;
            //cout << "Reading word ---> " <<hex<<word<<dec << endl;
            //cout << "Reading type ---> " <<hex<<type<<dec << endl;

            if(DEBUG_FLAG)
            {
                if(type == 0x1F)
                    cout << " ---> Event header ! " << endl;

                if(type == 0x3F)
                    cout << " ---> Event trailer ! " << endl;

                if(type == 0xDF)
                    cout << " ---> Error Flag ! " << endl;

                if(type == 0xFF)
                    cout << " ---> Debugging data ! " << endl;
            }

            //variable declaration: the value is filled at the first occurance
            long Event_Id, Bunch_Id, TDC_Id, channel, time, ROB_Id;
            map<int,int>::iterator iter;

            if(type == 0x1F) // Event header
            {
                if(DEBUG_FLAG)
                    cout << "Resetting hit array " << endl;
                for(int i=0; i<7000;i++)
                    hitArray[i]=0.;

                //reset number of hits
                for(int ic=0; ic<2; ic++)
                    for(int il=0; il<12; il++)
                        hnh[ic*12 + il]->Reset();

                //reset variables
                Event_Id = 999.;
                Bunch_Id = 999.;
                ROB_Id = 999.;

                numEvent +=1;
            }

            if(type == 0x3F) // Event trailer
            {
                if(DEBUG_FLAG)
                    cout << "Resetting hit array " << endl;
                for(int i=0; i<7000;i++)
                    hitArray[i]=0.;

                //reset variables
                Event_Id = 999.;
                Bunch_Id = 999.;
                ROB_Id = 999.;
            }

            if(type != 0x1F && type != 0x3F && type != 0xDF && type != 0xFF )
            {
                long group = type & 0xE0;
                group = group >> 5;

                //variable declaration
                //long Event_Id, Bunch_Id, TDC_Id, channel, time, ROB_Id;

                if(DEBUG_FLAG)
                    cout << "group " << group << endl;
                switch(group)
                {
                case 5:
                    if(DEBUG_FLAG)
                        cout << " Trailing measurment " << endl;
                    break;
                case 6:
                    if(DEBUG_FLAG)
                        cout << " Errors " << endl;
                    break;
                case 7:
                    if(DEBUG_FLAG)
                        cout << " Debugging data " << endl;
                    break;
                case 0:
                    Event_Id = TTCcount & 0xFFF000;
                    Event_Id  = Event_Id >> 12;
                    Bunch_Id = TTCcount & 0xFFF;
                    ROB_Id = type & 0x1F;
                    if(DEBUG_FLAG)
                    {
                        cout << "Group header: ";
                        cout << " Event " << setbase(10) << Event_Id <<
                                " bunch " << Bunch_Id <<
                                " rob " << ROB_Id << endl;
                    }

                    break;
                case 1:
                    Event_Id = TTCcount & 0xFFF000;
                    Event_Id  = Event_Id >> 12;
                    Bunch_Id = TTCcount & 0xFFF;
                    if(DEBUG_FLAG)
                    {
                        cout << "Group trailer: ";
                        cout << " Event " << setbase(10) << Event_Id <<
                                " bunch " << Bunch_Id <<
                                " rob " << ROB_Id << endl;
                    }

                    //reset variables
                    Event_Id = 999.;
                    Bunch_Id = 999.;
                    ROB_Id = 999.;

                    break;
                case 2:
                    Event_Id = TTCcount & 0xFFF000;
                    Event_Id  = Event_Id >> 12;
                    Bunch_Id = TTCcount & 0xFFF;
                    TDC_Id = ROB_Id & 0x3;
                    if(DEBUG_FLAG)
                    {
                        cout << "TDC Id " << setbase(10) << TDC_Id << " header: ";
                        cout << " Event " << setbase(10) << Event_Id <<
                                " bunch " << Bunch_Id <<
                                " rob " << ROB_Id << endl;
                    }
                    break;
                case 3:
                    Event_Id = TTCcount & 0xFFF000;
                    Event_Id  = Event_Id >> 12;
                    Bunch_Id = TTCcount & 0xFFF;
                    TDC_Id = ROB_Id & 0x3;
                    if(DEBUG_FLAG)
                    {
                        cout << "TDC Id " << setbase(10) << TDC_Id << " trailer: ";
                        cout << " Event " << setbase(10) << Event_Id <<
                                " bunch " << Bunch_Id <<
                                " rob " << ROB_Id << endl;
                    }
                    break;
                case 4:
                    if(ROB_Id==999.)
                        break;
                    TDC_Id = type & 0x3;
                    channel = TTCcount & 0xF80000;
                    channel = channel >> 19;
                    time = TTCcount & 0x7FFFF;


                    //htime->Fill(time/4.);

                    // ROB: 0..24;  TDC:0..3;  ch:0..31
                    int tdcflag = TDC_Id + 4*ROB_Id;
                    int chflag = channel + 32*tdcflag;
                    //if(DEBUG_FLAG)
                    if(chflag>3200)
                        cout <<	"Event_Id " << Event_Id << " tdc " << TDC_Id <<
                                " rob " << ROB_Id << " ch " << channel <<
                                " tdcflag " << tdcflag << " chflag " << chflag <<
                                " numEvent " << numEvent << endl;

//                    //fill tdc histograms
//                    if(time/4. > minTime && time/4. < maxTime)
//                        hocc->Fill(tdcflag,channel);
//                    htimech->Fill(time/4.,chflag);

                    //fill ch histograms
                    int rawdata = getTDCid(1,ROB_Id,TDC_Id,channel);


                    //cout << "TDCid " << rawdata << endl;
                    int lay;
                    int sl;
                    iter = chmap.find(rawdata);

                    if(DEBUG_FLAG && iter != chmap.end())
                    {
                        cout << "SE " << getSe(iter->second) <<
                                " sl " << getSL(iter->second) <<
                                " lay " << getLay(iter->second) <<
                                " tube " << getTube(iter->second) << endl;
                    }

                    // 20150414 SV fill htime histos for every layer
                    int ch =  abs(getSe(iter->second) -8 - 3);
                    if(ch < 2){
                        int lay = getLay(iter->second) + 4*(getSL(iter->second)-1);
                        int ind = ch*12 + (lay - 1);
                        int tube = getTube(iter->second);
                        if(DEBUG_FLAG)
                            cout << "SE " << getSe(iter->second) <<
                                " ch " << ch <<
                                " sl " << getSL(iter->second) <<
                                " lay " << getLay(iter->second) <<
                                " tube " << tube <<
                                " ind " << ind <<
                                " fill time/4. " << time/4. << endl;
                        ht[ind]->Fill(time/4.);

                        if(hnh[ind]->GetBinContent(tube)==0)
                            ht_firstHit[ind]->Fill(time/4.);

                        if(time/4. > minTimeHit && time/4. < maxTimeHit){
                            if(hnh[ind]->GetBinContent(tube)==0)
                                hitArray[rawdata] = time/4.;

                            if(hnh[ind]->GetBinContent(tube)==1)
                                hhit_diff[ind]->Fill(time/4. - hitArray[rawdata]);

                            hnh[ind]->SetBinContent(tube,hnh[ind]->GetBinContent(tube)+1);
                        }
                        if(ch == 0 && (lay-1)==2)
                                ht_ch0_lay1_tube[tube]->Fill(time/4.);
                    }

//                    if( iter != chmap.end() )
//                    {
//                        if( getSe(iter->second) == 8 || getSe(iter->second) == 9
//                                ||getSe(iter->second) == 10 || getSe(iter->second) == 11 ){

//                            if(time/4. > minTimeHit && time/4. < maxTimeHit)
//                            {
//                                //fill hit vector
//                                if(DEBUG_FLAG)
//                                    cout << "hitArray[" << rawdata << "] = " << hitArray[rawdata] << endl;
//                                if(hitArray[rawdata]==0.)
//                                    hitArray[rawdata] = time/4.;
//                                else
//                                {
//                                    hhit_diff[getSe(iter->second)-8]->Fill(time/4. - hitArray[rawdata]);
//                                    if(DEBUG_FLAG)
//                                        cout << "Filling hhit_diff with " << time/4. - hitArray[rawdata] << endl;
//                                }
//                            }

//                            if(getSe(iter->second) == 8){
//                                if(time/4. > minTime && time/4. < maxTime)
//                                    hoccSL1_lay[getLay(iter->second)-1]->Fill(getTube(iter->second));
//                                sl = 0;
//                            }
//                            if(getSe(iter->second) == 9){
//                                if(time/4. > minTime && time/4. < maxTime)
//                                    hoccSL2_lay[getLay(iter->second)-1]->Fill(getTube(iter->second));
//                                sl = 1;
//                            }
//                            if(getSe(iter->second) == 10){
//                                if(time/4. > minTime && time/4. < maxTime){
//                                    hoccCH1_lay[getLay(iter->second)-1 + (getSL(iter->second)-1)*4]
//                                            ->Fill(getTube(iter->second));
//                                }
//                                sl = getSL(iter->second) + 1;
//                            }
//                            if(getSe(iter->second) == 11){
//                                if(time/4. > minTime && time/4. < maxTime)
//                                    hoccCH2_lay[getLay(iter->second)-1 + (getSL(iter->second)-1)*4]
//                                            ->Fill(getTube(iter->second));
//                                sl = getSL(iter->second) + 4;
//                            }

//                            htime_sl[sl]->Fill(time/4.);

//                        }}
//                    if(DEBUG_FLAG)
//                        cout 	<< " Leading measurment " << " TDC " << setbase(10) << TDC_Id
//                                << " ch " << channel << " time " << time << endl;
                    //reset
                    TDC_Id = 999.;
                    channel = 999.;
                    time = 999.;

                    break;

                }//end switch


            }//end no ROS type


            //update TDC plots
//            if(words/100. == int(words/100.))
//            {
//                pad1->cd();
//                hocc->Draw();
//                pad1->Update();
//                pad2->cd();
//                htimech->Draw();
//                pad2->Update();
//                pad3->cd();
//                htime->Draw();
//                pad3->Update();

//                //chamber 2 plots
//                ch2->cd();
//                pad_hits_ch2->cd(1);
//                hoccCH2_lay[11]->Draw();
//                pad_hits_ch2->cd(2);
//                hoccCH2_lay[10]->Draw();
//                pad_hits_ch2->cd(3);
//                hoccCH2_lay[9]->Draw();
//                pad_hits_ch2->cd(4);
//                hoccCH2_lay[8]->Draw();
//                pad_hits_ch2->cd(5);
//                hoccCH2_lay[7]->Draw();
//                pad_hits_ch2->cd(6);
//                hoccCH2_lay[6]->Draw();
//                pad_hits_ch2->cd(7);
//                hoccCH2_lay[5]->Draw();
//                pad_hits_ch2->cd(8);
//                hoccCH2_lay[4]->Draw();
//                pad_hits_ch2->cd(9);
//                hoccCH2_lay[3]->Draw();
//                pad_hits_ch2->cd(10);
//                hoccCH2_lay[2]->Draw();
//                pad_hits_ch2->cd(11);
//                hoccCH2_lay[1]->Draw();
//                pad_hits_ch2->cd(12);
//                hoccCH2_lay[0]->Draw();

//                pad_time_ch2->cd(1);
//                htime_sl[7]->Draw();
//                pad_time_ch2->cd(2);
//                htime_sl[6]->Draw();
//                pad_time_ch2->cd(3);
//                htime_sl[5]->Draw();

//                ch2->Update();

//                //chamber 1 plots
//                ch1->cd();
//                pad_hits_ch1->cd(1);
//                hoccCH1_lay[11]->Draw();
//                pad_hits_ch1->cd(2);
//                hoccCH1_lay[10]->Draw();
//                pad_hits_ch1->cd(3);
//                hoccCH1_lay[9]->Draw();
//                pad_hits_ch1->cd(4);
//                hoccCH1_lay[8]->Draw();
//                pad_hits_ch1->cd(5);
//                hoccCH1_lay[7]->Draw();
//                pad_hits_ch1->cd(6);
//                hoccCH1_lay[6]->Draw();
//                pad_hits_ch1->cd(7);
//                hoccCH1_lay[5]->Draw();
//                pad_hits_ch1->cd(8);
//                hoccCH1_lay[4]->Draw();
//                pad_hits_ch1->cd(9);
//                hoccCH1_lay[3]->Draw();
//                pad_hits_ch1->cd(10);
//                hoccCH1_lay[2]->Draw();
//                pad_hits_ch1->cd(11);
//                hoccCH1_lay[1]->Draw();
//                pad_hits_ch1->cd(12);
//                hoccCH1_lay[0]->Draw();
//                pad_hits_ch1->Update();

//                pad_time_ch1->cd(1);
//                htime_sl[4]->Draw();
//                pad_time_ch1->cd(2);
//                htime_sl[3]->Draw();
//                pad_time_ch1->cd(3);
//                htime_sl[2]->Draw();

//                ch1->Update();

//                //SLs plots
//                sls->cd();
//                pad_hits_sls->cd(1);
//                hoccSL2_lay[3]->Draw();
//                pad_hits_sls->cd(2);
//                hoccSL2_lay[2]->Draw();
//                pad_hits_sls->cd(3);
//                hoccSL2_lay[1]->Draw();
//                pad_hits_sls->cd(4);
//                hoccSL2_lay[0]->Draw();
//                pad_hits_sls->cd(9);
//                hoccSL1_lay[3]->Draw();
//                pad_hits_sls->cd(10);
//                hoccSL1_lay[2]->Draw();
//                pad_hits_sls->cd(11);
//                hoccSL1_lay[1]->Draw();
//                pad_hits_sls->cd(12);
//                hoccSL1_lay[0]->Draw();
//                pad_hits_sls->Update();

//                pad_time_sls->cd(1);
//                htime_sl[1]->Draw();
//                pad_time_sls->cd(3);
//                htime_sl[0]->Draw();

//                sls->Update();

//                // hit difference plot
//                pad_hits_hits->cd(1);
//                hhit_diff[3]->Draw();
//                pad_hits_hits->cd(2);
//                hhit_diff[2]->Draw();
//                pad_hits_hits->cd(3);
//                hhit_diff[1]->Draw();
//                pad_hits_hits->cd(4);
//                hhit_diff[0]->Draw();
//                pad_hits_hits->Update();

//                hits->Update();
//            }
        }//end loop on words

    }//end while

    // draw
    TPaveText *tpt = new TPaveText(0.8,0.35,0.99,0.65,"brNDC");
    tpt->SetFillColor(kWhite);
    tpt->SetTextAlign(12);
    tpt->SetBorderSize(1);
    char txt[20];
    sprintf(txt,"Run %i",run);
    tpt->AddText(txt);
    tpt->AddText("#color[2]{First hit}");
    tpt->AddText("#color[1]{All hits}");
    tpt->AddText("Raw time box");

    // draw daq times
    TCanvas * cTimeCh1 = new TCanvas("cTimeCh1","cTimeCh1",1300,900);
    cTimeCh1->Divide(4,3);
    for(int il=0; il<12; il++){
        cTimeCh1->cd(1+il);
        ht[il]->Draw();
        ht_firstHit[il]->SetLineColor(kRed);
        ht_firstHit[il]->Draw("same");
    }
    tpt->Draw();

    // draw daq times ch0 lay1
    TCanvas * cTimeTubes1 = new TCanvas("cTimeTubes1","cTimeTubes1",1300,900);
    TCanvas * cTimeTubes2 = new TCanvas("cTimeTubes2","cTimeTubes2",1300,900);
    TCanvas * cTimeTubes3 = new TCanvas("cTimeTubes3","cTimeTubes3",1300,900);
    cTimeTubes1->Divide(6,4);
    cTimeTubes3->Divide(6,4);
    cTimeTubes2->Divide(6,4);
    for(int it=0; it<24; it++){
        cTimeTubes1->cd(1+it);
        ht_ch0_lay1_tube[it]->Draw();
        cTimeTubes2->cd(1+it);
        ht_ch0_lay1_tube[it+24]->Draw();
        cTimeTubes3->cd(1+it);
        ht_ch0_lay1_tube[it+48]->Draw();
    }
    tpt->Draw();

    TCanvas * cTimeCh2 = new TCanvas("cTimeCh2","cTimeCh2",1300,900);
    cTimeCh2->Divide(4,3);
    for(int il=0; il<12; il++){
        cTimeCh2->cd(1+il);
        ht[il+12]->Draw();
        ht_firstHit[il+12]->SetLineColor(kRed);
        ht_firstHit[il+12]->Draw("same");    }
    tpt->Draw();


    // draw hit time difference
    TCanvas * cdiffTimeCh1 = new TCanvas("cdiffTimeCh1","cdiffTimeCh1",1300,900);
    cdiffTimeCh1->Divide(4,3);
    for(int il=0; il<12; il++){
        cdiffTimeCh1->cd(1+il);
        hhit_diff[il]->Draw();
    }
    tpt->Draw();

    TCanvas * cdiffTimeCh2 = new TCanvas("cdiffTimeCh2","cdiffTimeCh2",1300,900);
    cdiffTimeCh2->Divide(4,3);
    for(int il=0; il<12; il++){
        cdiffTimeCh2->cd(1+il);
        hhit_diff[il+12]->Draw();
    }
    tpt->Draw();

    ///print desidered canvas
//    TString nameOut("hit_diff_run_");
//    nameOut += runName;
//    nameOut += ".eps";
//    hits->Print(nameOut,"eps");
    char outname[100];
    sprintf(outname, "~/workspace/experiments/radmu/mublast/analysis/20160308_variazione_LSD/%i_raw_timebox_ch1.png",run,run);
    cTimeCh1->Print(outname);
    sprintf(outname, "~/workspace/experiments/radmu/mublast/analysis/20160308_variazione_LSD/%i_raw_timebox_ch2.png",run,run);
    cTimeCh2->Print(outname);
//    sprintf(outname, "results0.95/%i/%i_raw_timebox_ch1.root",run,run);
//    cTimeCh1->Print(outname);
//    sprintf(outname, "results0.95/%i/%i_raw_timebox_ch2.root",run,run);
//    cTimeCh2->Print(outname);

    return;
}

void countAllFilesEvents(){

    ofstream myfile;
    myfile.open ("NumberEvents.txt");

    myfile << "\n *** Run " << 2326 << endl;
    for(int i=0; i<20;i++){
        char fileName[200];
        //sprintf(fileName,"/data/experiments/radmu/mublast/raw/MuBlast_vuoto_2400_2400_4600_mezzoferro_newrif_r2326.i%i",i);
        sprintf(fileName,"/mudata/LNL/muonegrafia/MuBlast_vuoto_2400_2400_4600_mezzoferro_newrif_r2326.i%i",i);
        int n = countRawEvents(fileName);
        myfile << n << endl;
    }

    myfile << "\n *** Run " << 2328 << endl;
    for(int i=0; i<20;i++){
        char fileName[200];
        //sprintf(fileName,"/data/experiments/radmu/mublast/raw/MuBlast_vuoto_2400_2400_4600_mezzoferro_newrif_2HVB_1HVC_GIF_r2328.i%i",i);
        sprintf(fileName,"/mudata/LNL/muonegrafia/MuBlast_vuoto_2400_2400_4600_mezzoferro_newrif_2HVB_1HVC_GIF_r2328.i%i",i);
        int n = countRawEvents(fileName);
        myfile << n << endl;
    }

    myfile << "\n *** Run " << 2331 << endl;
    for(int i=0; i<20;i++){
        char fileName[200];
        //sprintf(fileName,"/data/experiments/radmu/mublast/raw/MuBlast_vuoto_2400_2400_4600_mezzoferro_newrif_2HVB_1HVC_GIF_r2328_r2331.i%i",i);
        sprintf(fileName,"/mudata/LNL/muonegrafia/MuBlast_vuoto_2400_2400_4600_mezzoferro_newrif_2HVB_1HVC_GIF_r2328_r2331.i%i",i);
        int n = countRawEvents(fileName);
        myfile << n << endl;
    }

    myfile << "\n *** Run " << 2333 << endl;
    for(int i=0; i<100;i++){
        char fileName[200];
        //sprintf(fileName,"/data/experiments/radmu/mublast/raw/MuBlast_vuoto_2400_2400_4600_mezzoferro_newrif_2HVB_1HVC_GIF_r2333.i%i",i);
        sprintf(fileName,"/mudata/LNL/muonegrafia/MuBlast_vuoto_2400_2400_4600_mezzoferro_newrif_2HVB_1HVC_GIF_r2333.i%i",i);
        int n = countRawEvents(fileName);
        myfile << n << endl;
    }

    return;
}

int countRawEvents(char * fileName){
    // 20150414 SV readapted from radmuMonitorLNL 2008/01/28
    //FILE is a C-type file pointer included in stdio.h
    FILE *infile;

//    char fileName[200];
//    sprintf(fileName,"/data/experiments/radmu/mublast/raw/MuBlast_vuoto_2400_2400_4600_mezzoferro_newrif_r%i.i%i",run,num_file);
//    sprintf(fileName,"/data/experiments/radmu/mublast/raw/MuBlast_vuoto_2400_2400_4600_mezzoferro_newrif_r%i.i%i",run,num_file);
    infile = fopen(fileName,"r");

    //fill channel map
    fillMap();

    //long int is 32-bit=4-byte word
    long evbu;
    //NB size_t fread(void *ptr, size_t size, size_t nmemb, FILE *stream);
    //Reads data from the given stream into the array pointed to by ptr.
    //It reads nmemb number of elements of size size (in byte).
    //The total number of bytes read is (size*nmemb).
    //On success the number of elements read is returned.
    //On error or end-of-file the total number of elements successfully
    //read (which may be zero) is returned.

    //NB ftell(infile) returns the current file position. For binary stream,
    //then the value is the number of bytes from the beginning of the file.

    // hit array definition
    float hitArray[7000];

    //int maxwords = 100000;
    int numEvent = 0;

    while(!feof(infile))
    {
        // get the 32-bit word
        if(fread(&evbu,4,1,infile)==0)
            break;

        // DMA: swap 16 lsb with 16 msb
        //long word = (evbu>>16) | (evbu<<16);
        /// SV 20150414 to function with mublast raw data....
        long word = 0x00000000;
        word = ((evbu>>16) & 0x0000FFFF) | ((evbu<<16) & 0xFFFF0000);

        // type of packet and TTCcounts
        long type = word & 0xFF000000;
        type = type >> 24;			// type of data packet: bit 24-31
        long TTCcount = word & 0x00FFFFFF;	// TTC count: bit 0-23

        // if event header reset hit array
        //cout << "Reading evbu ---> " << hex<<evbu<<dec << endl;
        //cout << "Reading word ---> " <<hex<<word<<dec << endl;
        //cout << "Reading type ---> " <<hex<<type<<dec << endl;

        if(DEBUG_FLAG)
        {
            if(type == 0x1F)
                cout << " ---> Event header ! " << endl;

            if(type == 0x3F)
                cout << " ---> Event trailer ! " << endl;

            if(type == 0xDF)
                cout << " ---> Error Flag ! " << endl;

            if(type == 0xFF)
                cout << " ---> Debugging data ! " << endl;
        }

        //variable declaration: the value is filled at the first occurance
        long Event_Id, Bunch_Id, TDC_Id, channel, time, ROB_Id;
        map<int,int>::iterator iter;

        if(type == 0x1F) // Event header
        {
            if(DEBUG_FLAG)
                cout << "Resetting hit array " << endl;
            for(int i=0; i<7000;i++)
                hitArray[i]=0.;

            //reset variables
            Event_Id = 999.;
            Bunch_Id = 999.;
            ROB_Id = 999.;

            numEvent +=1;
        }

        if(type == 0x3F) // Event trailer
        {
            if(DEBUG_FLAG)
                cout << "Resetting hit array " << endl;
            for(int i=0; i<7000;i++)
                hitArray[i]=0.;

            //reset variables
            Event_Id = 999.;
            Bunch_Id = 999.;
            ROB_Id = 999.;
        }

        if(type != 0x1F && type != 0x3F && type != 0xDF && type != 0xFF )
        {
            long group = type & 0xE0;
            group = group >> 5;

            //variable declaration
            //long Event_Id, Bunch_Id, TDC_Id, channel, time, ROB_Id;

            if(DEBUG_FLAG)
                cout << "group " << group << endl;
            switch(group)
            {
            case 5:
                if(DEBUG_FLAG)
                    cout << " Trailing measurment " << endl;
                break;
            case 6:
                if(DEBUG_FLAG)
                    cout << " Errors " << endl;
                break;
            case 7:
                if(DEBUG_FLAG)
                    cout << " Debugging data " << endl;
                break;
            case 0:
                Event_Id = TTCcount & 0xFFF000;
                Event_Id  = Event_Id >> 12;
                Bunch_Id = TTCcount & 0xFFF;
                ROB_Id = type & 0x1F;
                if(DEBUG_FLAG)
                {
                    cout << "Group header: ";
                    cout << " Event " << setbase(10) << Event_Id <<
                            " bunch " << Bunch_Id <<
                            " rob " << ROB_Id << endl;
                }

                break;
            case 1:
                Event_Id = TTCcount & 0xFFF000;
                Event_Id  = Event_Id >> 12;
                Bunch_Id = TTCcount & 0xFFF;
                if(DEBUG_FLAG)
                {
                    cout << "Group trailer: ";
                    cout << " Event " << setbase(10) << Event_Id <<
                            " bunch " << Bunch_Id <<
                            " rob " << ROB_Id << endl;
                }

                //reset variables
                Event_Id = 999.;
                Bunch_Id = 999.;
                ROB_Id = 999.;

                break;
            case 2:
                Event_Id = TTCcount & 0xFFF000;
                Event_Id  = Event_Id >> 12;
                Bunch_Id = TTCcount & 0xFFF;
                TDC_Id = ROB_Id & 0x3;
                if(DEBUG_FLAG)
                {
                    cout << "TDC Id " << setbase(10) << TDC_Id << " header: ";
                    cout << " Event " << setbase(10) << Event_Id <<
                            " bunch " << Bunch_Id <<
                            " rob " << ROB_Id << endl;
                }
                break;
            case 3:
                Event_Id = TTCcount & 0xFFF000;
                Event_Id  = Event_Id >> 12;
                Bunch_Id = TTCcount & 0xFFF;
                TDC_Id = ROB_Id & 0x3;
                if(DEBUG_FLAG)
                {
                    cout << "TDC Id " << setbase(10) << TDC_Id << " trailer: ";
                    cout << " Event " << setbase(10) << Event_Id <<
                            " bunch " << Bunch_Id <<
                            " rob " << ROB_Id << endl;
                }
                break;
            case 4:
                if(ROB_Id==999.)
                    break;
                TDC_Id = type & 0x3;
                channel = TTCcount & 0xF80000;
                channel = channel >> 19;
                time = TTCcount & 0x7FFFF;

                // ROB: 0..24;  TDC:0..3;  ch:0..31
                int tdcflag = TDC_Id + 4*ROB_Id;
                int chflag = channel + 32*tdcflag;
                //if(DEBUG_FLAG)
                if(chflag>3200)
                    cout <<	"Event_Id " << Event_Id << " tdc " << TDC_Id <<
                            " rob " << ROB_Id << " ch " << channel <<
                            " tdcflag " << tdcflag << " chflag " << chflag <<
                            " numEvent " << numEvent << endl;

                //fill ch histograms
                int rawdata = getTDCid(1,ROB_Id,TDC_Id,channel);


                //cout << "TDCid " << rawdata << endl;
                int lay;
                int sl;
                iter = chmap.find(rawdata);

                if(DEBUG_FLAG && iter != chmap.end())
                {
                    cout << "SE " << getSe(iter->second) <<
                            " sl " << getSL(iter->second) <<
                            " lay " << getLay(iter->second) <<
                            " tube " << getTube(iter->second) << endl;
                }

                //reset
                TDC_Id = 999.;
                channel = 999.;
                time = 999.;

                break;
            }//end switch
        }//end no ROS type
    }//end while

//    cout << "Total numbers of events in file = "
    cout << numEvent << endl;

    return numEvent;
}


void testIBMuons(){

    gSystem->Load("/home/vanini/mutom/uLib/trunk/build/lib/libmutomRoot.so");

    // *** root style
    gROOT->SetStyle("Plain");
    //gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    gStyle->SetStatH(0.1);
    gStyle->SetStatW(0.35);

    // *** histograms definition
    TH1F*  hXin    = new TH1F("hXin", "In track X", 400, -200, 200);
    TH1F*  hYin    = new TH1F("hYin", "In track Y", 400, -200, 200);
    TH1F*  hZin    = new TH1F("hZin", "In track Z", 400, -200, 200);
    TH1F*  hXout    = new TH1F("hXout", "out track X", 400, -200, 200);
    TH1F*  hYout    = new TH1F("hYout", "out track Y", 400, -200, 200);
    TH1F*  hZout    = new TH1F("hZout", "out track Z", 400, -200, 200);

    // ***  open rootuples
    char filename[100];
    sprintf(filename,"2218_dump.root");
    cout << "Processing file " << filename << endl;

    TFile * file = new TFile(filename);
    TTree * tree = (TTree*)file->Get("muons");

    ROOT::Mutom::MuonScatter *mu;
    tree->GetBranch("mu")->SetAutoDelete(kFALSE);
    tree->SetBranchAddress("mu",&mu);

    int nentries = tree->GetEntries();

    // *** loop over events
    for (int i=0;i<TMath::Min(5000000,nentries);i++) {
      tree->GetEntry(i);

      hXin->Fill(mu->in.X);
      hYin->Fill(mu->in.Y);
      hZin->Fill(mu->in.Z);

      hXout->Fill(mu->out.X);
      hYout->Fill(mu->out.Y);
      hZout->Fill(mu->out.Z);

    } // end event loop

    file->Close();

    // *** draw canvas
    TCanvas* c = new TCanvas("c", "Track Position", 1600,600);
    c->Divide(2,1);
    c->cd(1);
    hXin->Draw();
    hXout->SetLineColor(kRed);
    hXout->Draw("sames");
    c->cd(2);
    hZin->Draw();
    hZout->SetLineColor(kRed);
    hZout->Draw("sames");

    return;
}

void dumpFromHistosFile(){

    gSystem->Load("/home/vanini/mutom/uLib/trunk/build/lib/libmutomRoot.so");

    // *** root style
    gROOT->SetStyle("Plain");
    //gStyle->SetOptStat(0);
    gStyle->SetOptFit(1);
    gStyle->SetStatH(0.1);
    gStyle->SetStatW(0.35);

    // *** histograms definition
    TH1F*  hXin    = new TH1F("hXin", "In track X", 400, -200, 200);
    TH1F*  hYin    = new TH1F("hYin", "In track Y", 400, -200, 200);
    TH1F*  hZin    = new TH1F("hZin", "In track Z", 400, -200, 200);
    TH1F*  hXout    = new TH1F("hXout", "out track X", 400, -200, 200);
    TH1F*  hYout    = new TH1F("hYout", "out track Y", 400, -200, 200);
    TH1F*  hZout    = new TH1F("hZout", "out track Z", 400, -200, 200);

    // ***  open rootuples
    char file_2326[100];
    sprintf(file_2326,"/data/experiments/radmu/mublast/fit/prova_HIT_r2326_0_1000000ev.root");
    cout << "Processing file " << file_2326 << endl;

    char file_2328[100];
    sprintf(file_2328,"/data/experiments/radmu/mublast/fit/prova_HIT_r2328_0_1000000ev.root");
    cout << "Processing file " << file_2328 << endl;

    char file_2331[100];
    sprintf(file_2331,"/data/experiments/radmu/mublast/fit/prova_HIT_r2331_0_1000000ev.root");
    cout << "Processing file " << file_2331 << endl;

    TFile * f26 = new TFile(file_2326);
    TFile * f28 = new TFile(file_2328);
    TFile * f31 = new TFile(file_2331);


    // *** draw canvas MP
    TCanvas* c = new TCanvas("c", "Track Position", 1600,600);
    c->Divide(2,1);
    c->cd(1);
    f26->cd();
    h_MP_CH1_glo->Draw();
    f28->cd();
    h_MP_CH1_glo->SetLineColor(kRed);
    h_MP_CH1_glo->Draw("sames");
    f31->cd();
    h_MP_CH1_glo->SetLineColor(kGreen);
    h_MP_CH1_glo->Draw("sames");

    c->cd(2);
    f26->cd();
    h_MP_CH2_glo->Draw();
    f28->cd();
    h_MP_CH2_glo->SetLineColor(kRed);
    h_MP_CH2_glo->Draw("sames");
    f31->cd();
    h_MP_CH2_glo->SetLineColor(kGreen);
    h_MP_CH2_glo->Draw("sames");


    // *** draw canvas SLOPES
    TCanvas* cs = new TCanvas("cs", "Track Position", 1600,600);
    cs->Divide(2,1);
    cs->cd(1);
    //gPad->SetLogy();
    f26->cd();
    h_dphi->Draw();
    cs->Update();
    TPaveStats *p = (TPaveStats*)h_dphi->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_dphi->SetLineColor(kRed);
    h_dphi->Draw("sames");
    cs->Update();
    TPaveStats *p = (TPaveStats*)h_dphi->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_dphi->SetLineColor(kGreen);
    h_dphi->Draw("sames");
    cs->Update();
    TPaveStats *p = (TPaveStats*)h_dphi->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    cs->cd(2);
    //gPad->SetLogy();
    f26->cd();
    h_dthe->Draw();
    cs->Update();
    TPaveStats *p = (TPaveStats*)h_dthe->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_dthe->SetLineColor(kRed);
    h_dthe->Draw("sames");
    cs->Update();
    TPaveStats *p = (TPaveStats*)h_dthe->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_dthe->SetLineColor(kGreen);
    h_dthe->Draw("sames");
    cs->Update();
    TPaveStats *p = (TPaveStats*)h_dthe->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    TPaveText *tpt = new TPaveText(0.8,0.35,0.99,0.65,"brNDC");
    tpt->SetFillColor(kWhite);
    tpt->SetTextAlign(12);
    tpt->SetBorderSize(1);
    tpt->AddText("#color[1]{Run 2326}");
    tpt->AddText("#color[2]{Run 2328}");
    tpt->AddText("#color[3]{Run 2331}");
    tpt->Draw();

//    const char outname[100];
//    sprintf(outname, "~/workspace/experiments/radmu/mublast/analysis/20160308_variazione_LSD/MP_comparison.png");
//    c->Print(outname);
//    sprintf(outname, "~/workspace/experiments/radmu/mublast/analysis/20160308_variazione_LSD/Scattering_comparison.png");
//    cs->Print(outname);

    // *** draw canvas SLOPE distribution
    TCanvas* csd = new TCanvas("csd", "Slopes", 1600,600);
    csd->Divide(3,2);
    csd->cd(1);
    //gPad->SetLogy();
    f26->cd();
    h_SLOPE_0->Draw();
    h_SLOPE_0->Rebin(8);
    csd->Update();
    TPaveStats *p = (TPaveStats*)h_SLOPE_0->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_SLOPE_0->SetLineColor(kRed);
    h_SLOPE_0->Draw("sames");
    h_SLOPE_0->Rebin(8);
    csd->Update();
    TPaveStats *p = (TPaveStats*)h_SLOPE_0->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_SLOPE_0->SetLineColor(kGreen);
    h_SLOPE_0->Draw("sames");
    h_SLOPE_0->Rebin(8);
    csd->Update();
    TPaveStats *p = (TPaveStats*)h_SLOPE_0->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    csd->cd(2);
    //gPad->SetLogy();
    f26->cd();
    h_SLOPE_1->Draw();
    h_SLOPE_1->Rebin(8);
    csd->Update();
    TPaveStats *p = (TPaveStats*)h_SLOPE_1->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_SLOPE_1->SetLineColor(kRed);
    h_SLOPE_1->Draw("sames");
    h_SLOPE_1->Rebin(8);
    csd->Update();
    TPaveStats *p = (TPaveStats*)h_SLOPE_1->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_SLOPE_1->SetLineColor(kGreen);
    h_SLOPE_1->Draw("sames");
    h_SLOPE_1->Rebin(8);
    csd->Update();
    TPaveStats *p = (TPaveStats*)h_SLOPE_1->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    csd->cd(3);
    //gPad->SetLogy();
    f26->cd();
    h_SLOPE_2->Draw();
    h_SLOPE_2->Rebin(8);
    csd->Update();
    TPaveStats *p = (TPaveStats*)h_SLOPE_2->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_SLOPE_2->SetLineColor(kRed);
    h_SLOPE_2->Draw("sames");
    h_SLOPE_2->Rebin(8);
    csd->Update();
    TPaveStats *p = (TPaveStats*)h_SLOPE_2->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_SLOPE_2->SetLineColor(kGreen);
    h_SLOPE_2->Draw("sames");
    h_SLOPE_2->Rebin(8);
    csd->Update();
    TPaveStats *p = (TPaveStats*)h_SLOPE_2->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    csd->cd(4);
    //gPad->SetLogy();
    f26->cd();
    h_SLOPE_3->Draw();
    h_SLOPE_3->Rebin(8);
    csd->Update();
    TPaveStats *p = (TPaveStats*)h_SLOPE_3->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_SLOPE_3->SetLineColor(kRed);
    h_SLOPE_3->Draw("sames");
    h_SLOPE_3->Rebin(8);
    csd->Update();
    TPaveStats *p = (TPaveStats*)h_SLOPE_3->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_SLOPE_3->SetLineColor(kGreen);
    h_SLOPE_3->Draw("sames");
    h_SLOPE_3->Rebin(8);
    csd->Update();
    TPaveStats *p = (TPaveStats*)h_SLOPE_3->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    csd->cd(5);
    //gPad->SetLogy();
    f26->cd();
    h_SLOPE_4->Draw();
    h_SLOPE_4->Rebin(8);
    csd->Update();
    TPaveStats *p = (TPaveStats*)h_SLOPE_4->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_SLOPE_4->SetLineColor(kRed);
    h_SLOPE_4->Draw("sames");
    h_SLOPE_4->Rebin(8);
    csd->Update();
    TPaveStats *p = (TPaveStats*)h_SLOPE_4->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_SLOPE_4->SetLineColor(kGreen);
    h_SLOPE_4->Draw("sames");
    h_SLOPE_4->Rebin(8);
    csd->Update();
    TPaveStats *p = (TPaveStats*)h_SLOPE_4->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    csd->cd(6);
    //gPad->SetLogy();
    f26->cd();
    h_SLOPE_5->Draw();
    h_SLOPE_5->Rebin(8);
    csd->Update();
    TPaveStats *p = (TPaveStats*)h_SLOPE_5->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_SLOPE_5->SetLineColor(kRed);
    h_SLOPE_5->Draw("sames");
    h_SLOPE_5->Rebin(8);
    csd->Update();
    TPaveStats *p = (TPaveStats*)h_SLOPE_5->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_SLOPE_5->SetLineColor(kGreen);
    h_SLOPE_5->Draw("sames");
    h_SLOPE_5->Rebin(8);
    csd->Update();
    TPaveStats *p = (TPaveStats*)h_SLOPE_5->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    // *** draw canvas NHITS distribution
    TCanvas* ch = new TCanvas("ch", "Slopes", 1600,1200);
    ch->Divide(3,3);
    ch->cd(1);
    gPad->SetLogy();
    f26->cd();
    h_Nhit_xSL0->Draw();
    h_Nhit_xSL0->GetXaxis()->SetRangeUser(0, 20);
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xSL0->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_Nhit_xSL0->SetLineColor(kRed);
    h_Nhit_xSL0->Draw("sames");
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xSL0->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_Nhit_xSL0->SetLineColor(kGreen);
    h_Nhit_xSL0->Draw("sames");
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xSL0->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    ch->cd(2);
    gPad->SetLogy();
    f26->cd();
    h_Nhit_xSL1->Draw();
    h_Nhit_xSL1->GetXaxis()->SetRangeUser(0, 20);
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xSL1->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_Nhit_xSL1->SetLineColor(kRed);
    h_Nhit_xSL1->Draw("sames");
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xSL1->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_Nhit_xSL1->SetLineColor(kGreen);
    h_Nhit_xSL1->Draw("sames");
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xSL1->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    ch->cd(3);
    gPad->SetLogy();
    f26->cd();
    h_Nhit_xSL2->Draw();
    h_Nhit_xSL2->GetXaxis()->SetRangeUser(0, 20);
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xSL2->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_Nhit_xSL2->SetLineColor(kRed);
    h_Nhit_xSL2->Draw("sames");
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xSL2->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_Nhit_xSL2->SetLineColor(kGreen);
    h_Nhit_xSL2->Draw("sames");
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xSL2->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    ch->cd(4);
    //gPad->SetLogy();
    f26->cd();
    h_Nhit_xSL3->Draw();
    h_Nhit_xSL3->GetXaxis()->SetRangeUser(0, 20);
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xSL3->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_Nhit_xSL3->SetLineColor(kRed);
    h_Nhit_xSL3->Draw("sames");
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xSL3->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_Nhit_xSL3->SetLineColor(kGreen);
    h_Nhit_xSL3->Draw("sames");
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xSL3->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    ch->cd(5);
    //gPad->SetLogy();
    f26->cd();
    h_Nhit_xSL4->Draw();
    h_Nhit_xSL4->GetXaxis()->SetRangeUser(0, 20);
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xSL4->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_Nhit_xSL4->SetLineColor(kRed);
    h_Nhit_xSL4->Draw("sames");
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xSL4->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_Nhit_xSL4->SetLineColor(kGreen);
    h_Nhit_xSL4->Draw("sames");
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xSL4->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    ch->cd(6);
    //gPad->SetLogy();
    f26->cd();
    h_Nhit_xSL5->Draw();
    h_Nhit_xSL5->GetXaxis()->SetRangeUser(0, 20);
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xSL5->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_Nhit_xSL5->SetLineColor(kRed);
    h_Nhit_xSL5->Draw("sames");
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xSL5->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_Nhit_xSL5->SetLineColor(kGreen);
    h_Nhit_xSL5->Draw("sames");
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xSL5->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    ch->cd(7);
    //gPad->SetLogy();
    f26->cd();
    h_Nhit_xSL6->Draw();
    h_Nhit_xSL6->GetXaxis()->SetRangeUser(0, 20);
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xSL6->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_Nhit_xSL6->SetLineColor(kRed);
    h_Nhit_xSL6->Draw("sames");
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xSL6->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_Nhit_xSL6->SetLineColor(kGreen);
    h_Nhit_xSL6->Draw("sames");
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xSL6->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    ch->cd(8);
    //gPad->SetLogy();
    f28->cd();
    h_Nhit_xSL7->Draw();
    h_Nhit_xSL7->GetXaxis()->SetRangeUser(0, 20);
    h_Nhit_xSL7->SetLineColor(kRed);
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xSL7->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f26->cd();
    h_Nhit_xSL7->SetLineColor(kBlack);
    h_Nhit_xSL7->Draw("sames");
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xSL7->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_Nhit_xSL7->SetLineColor(kGreen);
    h_Nhit_xSL7->Draw("sames");
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xSL7->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    ch->cd(9);
    //gPad->SetLogy();
    f26->cd();
    h_Nhit->Draw();
    h_Nhit->GetXaxis()->SetRangeUser(0, 100);
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_Nhit->SetLineColor(kRed);
    h_Nhit->Draw("sames");
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_Nhit->SetLineColor(kGreen);
    h_Nhit->Draw("sames");
    ch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    TPaveText *tpt = new TPaveText(0.8,0.35,0.99,0.65,"brNDC");
    tpt->SetFillColor(kWhite);
    tpt->SetTextAlign(12);
    tpt->SetBorderSize(1);
    tpt->AddText("#color[1]{Run 2326}");
    tpt->AddText("#color[2]{Run 2328}");
    tpt->AddText("#color[3]{Run 2331}");
    tpt->Draw();

//    const char outname[100];
//    sprintf(outname, "~/workspace/experiments/radmu/mublast/analysis/20160308_variazione_LSD/MP_comparison.png");
//    c->Print(outname);
//    sprintf(outname, "~/workspace/experiments/radmu/mublast/analysis/20160308_variazione_LSD/Scattering_comparison.png");
//    csd->Print(outname);

    TCanvas* cch = new TCanvas("cch", "Hits", 1200,600);
    cch->Divide(2,1);
    cch->cd(1);
    //gPad->SetLogy();
    f26->cd();
    h_Nhit_xCH0->Draw();
    h_Nhit_xCH0->GetXaxis()->SetRangeUser(0, 50);
    cch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xCH0->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_Nhit_xCH0->SetLineColor(kRed);
    h_Nhit_xCH0->Draw("sames");
    cch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xCH0->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_Nhit_xCH0->SetLineColor(kGreen);
    h_Nhit_xCH0->Draw("sames");
    cch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xCH0->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    cch->cd(2);
    //gPad->SetLogy();
    f26->cd();
    h_Nhit_xCH1->Draw();
    h_Nhit_xCH1->GetXaxis()->SetRangeUser(0, 50);
    cch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xCH1->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_Nhit_xCH1->SetLineColor(kRed);
    h_Nhit_xCH1->Draw("sames");
    cch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xCH1->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_Nhit_xCH1->SetLineColor(kGreen);
    h_Nhit_xCH1->Draw("sames");
    cch->Update();
    TPaveStats *p = (TPaveStats*)h_Nhit_xCH1->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);
    tpt->Draw();

    TCanvas* cr = new TCanvas("cr", "Residuals", 1200,600);
    cr->Divide(2,2);
    cr->cd(1);
    //gPad->SetLogy();
    f26->cd();
    h_resCH1Phi_glo->Draw();
    cr->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1Phi_glo->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_resCH1Phi_glo->SetLineColor(kRed);
    h_resCH1Phi_glo->Draw("sames");
    cr->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1Phi_glo->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_resCH1Phi_glo->SetLineColor(kGreen);
    h_resCH1Phi_glo->Draw("sames");
    cr->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1Phi_glo->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    cr->cd(2);
    //gPad->SetLogy();
    f26->cd();
    h_resCH2Phi_glo->Draw();
    cr->Update();
    TPaveStats *p = (TPaveStats*)h_resCH2Phi_glo->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_resCH2Phi_glo->SetLineColor(kRed);
    h_resCH2Phi_glo->Draw("sames");
    cr->Update();
    TPaveStats *p = (TPaveStats*)h_resCH2Phi_glo->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_resCH2Phi_glo->SetLineColor(kGreen);
    h_resCH2Phi_glo->Draw("sames");
    cr->Update();
    TPaveStats *p = (TPaveStats*)h_resCH2Phi_glo->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    cr->cd(3);
    //gPad->SetLogy();
    f26->cd();
    h_resCH1The_glo->Draw();
    cr->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1The_glo->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_resCH1The_glo->SetLineColor(kRed);
    h_resCH1The_glo->Draw("sames");
    cr->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1The_glo->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_resCH1The_glo->SetLineColor(kGreen);
    h_resCH1The_glo->Draw("sames");
    cr->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1The_glo->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    cr->cd(4);
    //gPad->SetLogy();
    f26->cd();
    h_resCH2The_glo->Draw();
    cr->Update();
    TPaveStats *p = (TPaveStats*)h_resCH2The_glo->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_resCH2The_glo->SetLineColor(kRed);
    h_resCH2The_glo->Draw("sames");
    cr->Update();
    TPaveStats *p = (TPaveStats*)h_resCH2The_glo->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_resCH2The_glo->SetLineColor(kGreen);
    h_resCH2The_glo->Draw("sames");
    cr->Update();
    TPaveStats *p = (TPaveStats*)h_resCH2The_glo->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    TCanvas* crc1 = new TCanvas("crc1", "Residuals CH1", 1200,600);
    crc1->Divide(4,2);
    crc1->cd(1);
    //gPad->SetLogy();
    f26->cd();
    h_resCH1phi_glo_0->Draw();
    h_resCH1phi_glo_0->Rebin(4);
    h_resCH1phi_glo_0->GetXaxis()->SetRangeUser(-1000.,1000.);
    crc1->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1phi_glo_0->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_resCH1phi_glo_0->SetLineColor(kRed);
    h_resCH1phi_glo_0->Draw("sames");
    h_resCH1phi_glo_0->Rebin(4);
    crc1->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1phi_glo_0->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_resCH1phi_glo_0->SetLineColor(kGreen);
    h_resCH1phi_glo_0->Draw("sames");
    h_resCH1phi_glo_0->Rebin(4);
    crc1->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1phi_glo_0->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    crc1->cd(2);
    //gPad->SetLogy();
    f26->cd();
    h_resCH1phi_glo_1->Draw();
    h_resCH1phi_glo_1->Rebin(4);
    h_resCH1phi_glo_1->GetXaxis()->SetRangeUser(-1000.,1000.);
    crc1->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1phi_glo_1->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_resCH1phi_glo_1->SetLineColor(kRed);
    h_resCH1phi_glo_1->Draw("sames");
    h_resCH1phi_glo_1->Rebin(4);
    crc1->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1phi_glo_1->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_resCH1phi_glo_1->SetLineColor(kGreen);
    h_resCH1phi_glo_1->Draw("sames");
    h_resCH1phi_glo_1->Rebin(4);
    crc1->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1phi_glo_1->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    crc1->cd(3);
    //gPad->SetLogy();
    f26->cd();
    h_resCH1phi_glo_2->Draw();
    h_resCH1phi_glo_2->Rebin(4);
    h_resCH1phi_glo_2->GetXaxis()->SetRangeUser(-1000.,1000.);
    crc1->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1phi_glo_2->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_resCH1phi_glo_2->SetLineColor(kRed);
    h_resCH1phi_glo_2->Draw("sames");
    h_resCH1phi_glo_2->Rebin(4);
    crc1->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1phi_glo_2->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_resCH1phi_glo_2->SetLineColor(kGreen);
    h_resCH1phi_glo_2->Draw("sames");
    h_resCH1phi_glo_2->Rebin(4);
    crc1->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1phi_glo_2->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    crc1->cd(4);
    //gPad->SetLogy();
    f26->cd();
    h_resCH1phi_glo_3->Draw();
    h_resCH1phi_glo_3->Rebin(4);
    h_resCH1phi_glo_3->GetXaxis()->SetRangeUser(-1000.,1000.);
    crc1->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1phi_glo_3->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_resCH1phi_glo_3->SetLineColor(kRed);
    h_resCH1phi_glo_3->Draw("sames");
    h_resCH1phi_glo_3->Rebin(4);
    crc1->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1phi_glo_3->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_resCH1phi_glo_3->SetLineColor(kGreen);
    h_resCH1phi_glo_3->Draw("sames");
    h_resCH1phi_glo_3->Rebin(4);
    crc1->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1phi_glo_3->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    crc1->cd(5);
    //gPad->SetLogy();
    f26->cd();
    h_resCH1phi_glo_4->Draw();
    h_resCH1phi_glo_4->Rebin(4);
    h_resCH1phi_glo_4->GetXaxis()->SetRangeUser(-1000.,1000.);
    crc1->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1phi_glo_4->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_resCH1phi_glo_4->SetLineColor(kRed);
    h_resCH1phi_glo_4->Draw("sames");
    h_resCH1phi_glo_4->Rebin(4);
    crc1->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1phi_glo_4->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_resCH1phi_glo_4->SetLineColor(kGreen);
    h_resCH1phi_glo_4->Draw("sames");
    h_resCH1phi_glo_4->Rebin(4);
    crc1->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1phi_glo_4->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    crc1->cd(6);
    //gPad->SetLogy();
    f26->cd();
    h_resCH1phi_glo_5->Draw();
    h_resCH1phi_glo_5->Rebin(4);
    h_resCH1phi_glo_5->GetXaxis()->SetRangeUser(-1000.,1000.);
    crc1->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1phi_glo_5->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_resCH1phi_glo_5->SetLineColor(kRed);
    h_resCH1phi_glo_5->Draw("sames");
    h_resCH1phi_glo_5->Rebin(4);
    crc1->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1phi_glo_5->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_resCH1phi_glo_5->SetLineColor(kGreen);
    h_resCH1phi_glo_5->Draw("sames");
    h_resCH1phi_glo_5->Rebin(4);
    crc1->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1phi_glo_5->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    crc1->cd(7);
    //gPad->SetLogy();
    f26->cd();
    h_resCH1phi_glo_6->Draw();
    h_resCH1phi_glo_6->Rebin(4);
    h_resCH1phi_glo_6->GetXaxis()->SetRangeUser(-1000.,1000.);
    crc1->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1phi_glo_6->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_resCH1phi_glo_6->SetLineColor(kRed);
    h_resCH1phi_glo_6->Draw("sames");
    h_resCH1phi_glo_6->Rebin(4);
    crc1->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1phi_glo_6->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_resCH1phi_glo_6->SetLineColor(kGreen);
    h_resCH1phi_glo_6->Draw("sames");
    h_resCH1phi_glo_6->Rebin(4);
    crc1->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1phi_glo_6->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    crc1->cd(8);
    //gPad->SetLogy();
    f26->cd();
    h_resCH1phi_glo_7->Draw();
    h_resCH1phi_glo_7->Rebin(4);
    h_resCH1phi_glo_7->GetXaxis()->SetRangeUser(-1000.,1000.);
    crc1->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1phi_glo_7->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_resCH1phi_glo_7->SetLineColor(kRed);
    h_resCH1phi_glo_7->Draw("sames");
    h_resCH1phi_glo_7->Rebin(4);
    crc1->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1phi_glo_7->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_resCH1phi_glo_7->SetLineColor(kGreen);
    h_resCH1phi_glo_7->Draw("sames");
    h_resCH1phi_glo_7->Rebin(4);
    crc1->Update();
    TPaveStats *p = (TPaveStats*)h_resCH1phi_glo_7->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);


    TCanvas* crc2 = new TCanvas("crc2", "Residuals CH2", 1200,600);
    crc2->Divide(2,2);
    crc2->cd(1);
    //gPad->SetLogy();
    f26->cd();
    h_resCH2phi_glo_0->Draw();
    h_resCH2phi_glo_0->Rebin(8);
    crc2->Update();
    TPaveStats *p = (TPaveStats*)h_resCH2phi_glo_0->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_resCH2phi_glo_0->SetLineColor(kRed);
    h_resCH2phi_glo_0->Draw("sames");
    h_resCH2phi_glo_0->Rebin(8);
    crc2->Update();
    TPaveStats *p = (TPaveStats*)h_resCH2phi_glo_0->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_resCH2phi_glo_0->SetLineColor(kGreen);
    h_resCH2phi_glo_0->Draw("sames");
    h_resCH2phi_glo_0->Rebin(8);
    crc2->Update();
    TPaveStats *p = (TPaveStats*)h_resCH2phi_glo_0->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    crc2->cd(2);
    //gPad->SetLogy();
    f26->cd();
    h_resCH2phi_glo_1->Draw();
    h_resCH2phi_glo_1->Rebin(8);
    crc2->Update();
    TPaveStats *p = (TPaveStats*)h_resCH2phi_glo_1->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_resCH2phi_glo_1->SetLineColor(kRed);
    h_resCH2phi_glo_1->Draw("sames");
    h_resCH2phi_glo_1->Rebin(8);
    crc2->Update();
    TPaveStats *p = (TPaveStats*)h_resCH2phi_glo_1->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_resCH2phi_glo_1->SetLineColor(kGreen);
    h_resCH2phi_glo_1->Draw("sames");
    h_resCH2phi_glo_1->Rebin(8);
    crc2->Update();
    TPaveStats *p = (TPaveStats*)h_resCH2phi_glo_1->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    crc2->cd(3);
    //gPad->SetLogy();
    f26->cd();
    h_resCH2phi_glo_2->Draw();
    h_resCH2phi_glo_2->Rebin(8);
    crc2->Update();
    TPaveStats *p = (TPaveStats*)h_resCH2phi_glo_2->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_resCH2phi_glo_2->SetLineColor(kRed);
    h_resCH2phi_glo_2->Draw("sames");
    h_resCH2phi_glo_2->Rebin(8);
    crc2->Update();
    TPaveStats *p = (TPaveStats*)h_resCH2phi_glo_2->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_resCH2phi_glo_2->SetLineColor(kGreen);
    h_resCH2phi_glo_2->Draw("sames");
    h_resCH2phi_glo_2->Rebin(8);
    crc2->Update();
    TPaveStats *p = (TPaveStats*)h_resCH2phi_glo_2->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    crc2->cd(4);
    //gPad->SetLogy();
    f26->cd();
    h_resCH2phi_glo_3->Draw();
    h_resCH2phi_glo_3->Rebin(8);
    crc2->Update();
    TPaveStats *p = (TPaveStats*)h_resCH2phi_glo_3->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(1);
    f28->cd();
    h_resCH2phi_glo_3->SetLineColor(kRed);
    h_resCH2phi_glo_3->Draw("sames");
    h_resCH2phi_glo_3->Rebin(8);
    crc2->Update();
    TPaveStats *p = (TPaveStats*)h_resCH2phi_glo_3->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(2);
    f31->cd();
    h_resCH2phi_glo_3->SetLineColor(kGreen);
    h_resCH2phi_glo_3->Draw("sames");
    h_resCH2phi_glo_3->Rebin(8);
    crc2->Update();
    TPaveStats *p = (TPaveStats*)h_resCH2phi_glo_3->GetListOfFunctions()->FindObject("stats");
    p->SetTextColor(3);

    TCanvas* ct = new TCanvas("ct", "Drift times", 1200,600);
    ct->Divide(3,2);
    ct->cd(1);
    //gPad->SetLogy();
    f26->cd();
    h_tempo_0->Draw();
    h_tempo_0->GetXaxis()->SetRangeUser(-100.,600.);
    f28->cd();
    h_tempo_0->SetLineColor(kRed);
    h_tempo_0->Draw("sames");
    f31->cd();
    h_tempo_0->SetLineColor(kGreen);
    h_tempo_0->Draw("sames");

    ct->cd(2);
    //gPad->SetLogy();
    f26->cd();
    h_tempo_1->Draw();
    h_tempo_1->GetXaxis()->SetRangeUser(-100.,600.);
    f28->cd();
    h_tempo_1->SetLineColor(kRed);
    h_tempo_1->Draw("sames");
    f31->cd();
    h_tempo_1->SetLineColor(kGreen);
    h_tempo_1->Draw("sames");

    ct->cd(3);
    //gPad->SetLogy();
    f26->cd();
    h_tempo_2->Draw();
    h_tempo_2->GetXaxis()->SetRangeUser(-100.,600.);
    f28->cd();
    h_tempo_2->SetLineColor(kRed);
    h_tempo_2->Draw("sames");
    f31->cd();
    h_tempo_2->SetLineColor(kGreen);
    h_tempo_2->Draw("sames");

    ct->cd(4);
    //gPad->SetLogy();
    f26->cd();
    h_tempo_3->Draw();
    h_tempo_3->GetXaxis()->SetRangeUser(-100.,600.);
    f28->cd();
    h_tempo_3->SetLineColor(kRed);
    h_tempo_3->Draw("sames");
    f31->cd();
    h_tempo_3->SetLineColor(kGreen);
    h_tempo_3->Draw("sames");

    ct->cd(5);
    //gPad->SetLogy();
    f26->cd();
    h_tempo_4->Draw();
    h_tempo_4->GetXaxis()->SetRangeUser(-100.,600.);
    f28->cd();
    h_tempo_4->SetLineColor(kRed);
    h_tempo_4->Draw("sames");
    f31->cd();
    h_tempo_4->SetLineColor(kGreen);
    h_tempo_4->Draw("sames");

    ct->cd(6);
    //gPad->SetLogy();
    f26->cd();
    h_tempo_5->Draw();
    h_tempo_5->GetXaxis()->SetRangeUser(-100.,600.);
    f28->cd();
    h_tempo_5->SetLineColor(kRed);
    h_tempo_5->Draw("sames");
    f31->cd();
    h_tempo_5->SetLineColor(kGreen);
    h_tempo_5->Draw("sames");
    tpt->Draw();

    return;
}

TGraphErrors * createTGraph(const Vector<Vector2f> &v, const Vector<float> &e, int mcolor, int mstyle=7) {
    if(mcolor==10)
        mcolor+=1;

    /// check that the error are flat over density values
    //return verifyErrors(v,e,3);

    // 20150129 assume the same error in every point, computed as the mean of the relative errors
    float err_rel_mean = 0;
    for(int i=0 ; i<v.size(); ++i){
        float err_rel = e[i]/v[i](1);
        err_rel_mean += err_rel;
    }
    err_rel_mean /= v.size();


    TGraphErrors *gr = new TGraphErrors();
    for(int i=0 ; i<v.size(); ++i) {
        gr->SetPoint(i,v[i](0),v[i](1));
        //gr->SetPointError(i,0,e[i]);
        float e = err_rel_mean*v[i](1);
        gr->SetPointError(i,0,e);
    }
    gr->SetMarkerColor(mcolor);
    gr->SetLineColor(mcolor);
    gr->SetLineWidth(2.4);
    gr->SetMarkerStyle(7);
    gr->SetMarkerSize(1.7);
    //mstyle++;
    mcolor++;
    return gr;
}
