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



#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <stdio.h>
#include "IBVoxCollectionCap.h"
#include "IBSubImageGrabber.h"
#include "IBVoxImageScanner.h"
#include <IBVoxFilters.h>

void PrintRData(RangeThresholdScan::ScanData d);
void PrintData(SimpleThresholdScan::ScanData d);

int main(int argc, char** argv)
{
    std::cout << "Initializing containers..." << std::flush;
    TFile f(argv[4],"RECREATE");

    TTree t("ROC", argv[4]);
    RangeThresholdScan::ScanOption opt;
    for(int i=0; i<200; ++i) {
        SimpleThresholdScan::ScanOption o;
        o.Threshold = (0.1e-6)*i;
        opt.push_back(o);
    }

    float perc_b_t[12][opt.size()];
    float inte_b_t[12][opt.size()];
    float iden_b_t[12][opt.size()];

    for (int i=0; i<12; ++i) {
        for (int j=0; j<opt.size(); ++j) {
            perc_b_t[i][j] = 0.0;
            inte_b_t[i][j] = 0.0;
            iden_b_t[i][j] = 0.0;
        }
    }

    float perc[12],inte[12],iden[12];
    float thres;

    t.Branch("threshold", &thres, "thres/F");
    char bnamePerc[12][40],bnameInte[12][40],bnameIden[12][40];
    for(int i=0; i<12; ++i) {
        perc[i] = 0.0;
        inte[i] = 0.0;
        iden[i] = 0.0;
        sprintf(bnamePerc[i], "AveragePercent_%i",        i);
        sprintf(bnameInte[i], "AverageIntensity_%i",      i);
        sprintf(bnameIden[i], "IdentificationPercent_%i", i);
        t.Branch(bnamePerc[i], &perc[i], "perc/F");
        t.Branch(bnameInte[i], &inte[i], "inte/F");
        t.Branch(bnameIden[i], &iden[i], "iden/F");
    }

    IBVoxCollectionCap image(Vector3i(0,0,0));
    std::cout << "done!\n" << std::flush;
    std::cout << "Scanning the 9 areas of POSITIVE situation...\n" << std::flush;
    // lead
    int fbulk = 0;
    for(int ii=1; ii<=atoi(argv[3]); ++ii) {
        fbulk++;
        char fname[100];
        sprintf(fname, "%s%i.vtk",argv[1],ii);
        image.ImportFromVtk(fname);
        // filter
        IBVoxFilter_Abtrim trim(Vector3i(5,5,5));
        IBFilterGaussShape shape(0.2);
        trim.SetKernelSpherical(shape);
        trim.SetABTrim(0,2);
        trim.SetImage(&image);
        trim.Run();
//        Vector<float> kern;
//        for(int i=0; i<8; ++i) kern.push_back(static_cast<float>(1));
//        IBVoxFilter_Linear average(Vector3i(2,2,2));
//        average.SetImage(&image);
//        average.SetKernelNumericXZY(kern);
//        average.Run();
        // grabbing
        Vector3i z(0,0,0);
        IBSubImageGrabber<IBVoxCollectionCap> grabber(image);
        IBLightCollection imgLR_B(z),imgLR_M(z), imgLR_F(z);
        imgLR_B = grabber.GrabRegion<IBLightCollection>(Vector3i(21,14,9),Vector3i(28,22,16));
        imgLR_M = grabber.GrabRegion<IBLightCollection>(Vector3i(46,14,9),Vector3i(53,22,16));
        imgLR_F = grabber.GrabRegion<IBLightCollection>(Vector3i(76,14,9),Vector3i(83,22,16));
        IBLightCollection imgMR_B(z), imgMR_M(z), imgMR_F(z);
        imgMR_B = grabber.GrabRegion<IBLightCollection>(Vector3i(36,31,26),Vector3i( 43,39,33));
        imgMR_M = grabber.GrabRegion<IBLightCollection>(Vector3i(66,31,26),Vector3i( 73,39,33));
        imgMR_F = grabber.GrabRegion<IBLightCollection>(Vector3i(96,31,26),Vector3i(103,39,33));
        IBLightCollection imgHR_B(z), imgHR_M(z), imgHR_F(z);
        imgHR_B = grabber.GrabRegion<IBLightCollection>(Vector3i( 56,49,43),Vector3i( 63,57,50));
        imgHR_M = grabber.GrabRegion<IBLightCollection>(Vector3i( 86,49,43),Vector3i( 93,57,50));
        imgHR_F = grabber.GrabRegion<IBLightCollection>(Vector3i(111,49,43),Vector3i(118,57,50));
        std::vector<IBLightCollection> boxC;
        boxC.push_back(imgLR_B);
        boxC.push_back(imgLR_M);
        boxC.push_back(imgLR_F);
        boxC.push_back(imgMR_B);
        boxC.push_back(imgMR_M);
        boxC.push_back(imgMR_F);
        boxC.push_back(imgHR_B);
        boxC.push_back(imgHR_M);
        boxC.push_back(imgHR_F);
        // scanning
        IBVoxImageScanner<IBLightCollection> scanner;
        for(int i=0; i<boxC.size(); ++i) {
            scanner.SetImage(&boxC.at(i));
            RangeThresholdScan::ScanData res = scanner.ScanImage<RangeThresholdScan>(opt);
            for(int j=0; j<opt.size(); ++j) {
                perc_b_t[i][j] += res.at(j).Percent;
                inte_b_t[i][j] += res.at(j).Intensity;
                iden_b_t[i][j] += (res.at(j).Percent>0.f) ? 1.0 : 0.0;
            }
        }
        std::cout << "\rProcessing " << (float)100*ii/atoi(argv[3]) << "\% complete." << std::flush;
    }
    std::cout << "\n...done!\n" << std::flush;

    // no lead
    std::cout << "Scanning the whole images of NEGATIVE situation...\n" << std::flush;
    fbulk = 0;
    for(int ii=1; ii<=atoi(argv[3]); ++ii) {
        fbulk++;
        char fname[100];
        sprintf(fname, "%s%i.vtk",argv[2], ii);
        image.ImportFromVtk(fname);
        //filtering
        IBVoxFilter_Abtrim trim(Vector3i(5,5,5));
        IBFilterGaussShape shape(0.2);
        trim.SetKernelSpherical(shape);
        trim.SetABTrim(0,2);
        trim.SetImage(&image);
        trim.Run();
//        Vector<float> kern;
//        for(int i=0; i<8; ++i) kern.push_back(static_cast<float>(1));
//        IBVoxFilter_Linear average(Vector3i(2,2,2));
//        average.SetImage(&image);
//        average.SetKernelNumericXZY(kern);
//        average.Run();
        // grabbing
        IBSubImageGrabber<IBVoxCollectionCap> grabber(image);
        IBLightCollection imgContBott(Vector3i(0,0,0)), imgContMidd(Vector3i(0,0,0)), imgContTop(Vector3i(0,0,0));
        imgContBott = grabber.GrabRegion<IBLightCollection>(Vector3i(10,10,4),Vector3i(129,26,55));
	imgContMidd = grabber.GrabRegion<IBLightCollection>(Vector3i(10,27,4),Vector3i(129,44,55));
	imgContTop  = grabber.GrabRegion<IBLightCollection>(Vector3i(10,45,4),Vector3i(129,61,55));
	std::vector<IBLightCollection> boxS;
	boxS.push_back(imgContBott);
	boxS.push_back(imgContMidd);
	boxS.push_back(imgContTop );
	IBVoxImageScanner<IBLightCollection> scanner;
	for(int i=0; i<boxS.size(); ++i) {
	    scanner.SetImage(&boxS.at(i));
	    RangeThresholdScan::ScanData res = scanner.ScanImage<RangeThresholdScan>(opt);
	    for(int j=0; j<opt.size(); ++j) {
	        perc_b_t[9+i][j] += res.at(j).Percent;
			inte_b_t[9+i][j] += res.at(j).Intensity;
			iden_b_t[9+i][j] += (res.at(j).Percent>0.f) ? 1.0 : 0.0;
	    }
	}
        std::cout << "\rProcessing " << (float)100*ii/atoi(argv[3]) << "\% complete." << std::flush;
    }

    std::cout << "\n...done!\n" << std::flush;

    std::cout << "Finalizing Data and Saving..." << std::flush;

    for(int j=0; j<opt.size(); ++j) {
        for(int i=0; i<12; ++i) {
            perc[i] = perc_b_t[i][j] / fbulk;
            inte[i] = inte_b_t[i][j] / fbulk;
            iden[i] = iden_b_t[i][j] / fbulk;
        }
        thres = opt.at(j).Threshold;
        t.Fill();
    }

    t.Write();
    f.Close();
    std::cout << "done!\nExiting!\n" << std::flush;
    return 0;
}
