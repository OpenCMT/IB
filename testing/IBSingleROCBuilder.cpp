#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <stdio.h>
#include "IBVoxCollectionCap.h"
#include "IBSubImageGrabber.h"
#include "IBVoxImageScanner.h"
#include <IBVoxFilters.h>
int main(int argc, char** argv)
{
    std::cout << "Initializing containers..." << std::flush;
    TFile f(argv[3],"RECREATE");

    TTree t("ROC", argv[3]);
    RangeThresholdScan::ScanOption opt;
    for(int i=0; i<200; ++i) {
        SimpleThresholdScan::ScanOption o;
        o.Threshold = (0.1e-6)*i;
        opt.push_back(o);
    }

    float perc_b_t[2][opt.size()];
    float inte_b_t[2][opt.size()];
    float iden_b_t[2][opt.size()];

    for (int i=0; i<2; ++i) {
        for (int j=0; j<opt.size(); ++j) {
            perc_b_t[i][j] = 0.0;
            inte_b_t[i][j] = 0.0;
            iden_b_t[i][j] = 0.0;
        }
    }

    float perc[2],inte[2],iden[2];
    float thres;
    t.Branch("threshold", &thres, "thres/F");
    char bnamePerc[2][40],bnameInte[2][40],bnameIden[2][40];
    for(int i=0; i<2; ++i) {
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
    std::cout << "Scanning Leaded Samples...\n" << std::flush;
    // lead
    int fbulk = 0;
    for(int y=0; y<3; ++y) {
        char fname[200];
        sprintf(fname, "%s%i_200.vtk", argv[1], y);
        image.ImportFromVtk(fname);
        IBVoxFilter_Abtrim trim(Vector3i(5,5,5));
        IBFilterGaussShape shape(0.2);
        trim.SetKernelSpherical(shape);
        trim.SetABTrim(0,2);
        trim.SetImage(&image);
        //trim.Run();
        Vector3i zero(0,0,0);
        IBSubImageGrabber<IBVoxCollectionCap> grabber(image);

        for (int x=0; x<6; ++x) {
            for (int z=0; z<3; ++z) {
                fbulk++;
                IBLightCollection img(zero);
                img = grabber.GrabRegion<IBLightCollection>(HPoint3f(-250+x*100,
                                                                     -100+y*100,
                                                                     -100+z*100),
                                                            HVector3f(40,40,40));
                IBVoxImageScanner<IBLightCollection> scanner;
                scanner.SetImage(&img);
                RangeThresholdScan::ScanData res = scanner.ScanImage<RangeThresholdScan>(opt);
                for(int j=0; j<opt.size(); ++j) {
                    perc_b_t[0][j] += res.at(j).Percent;
                    inte_b_t[0][j] += res.at(j).Intensity;
                    iden_b_t[0][j] += (res.at(j).Percent>0.f) ? 1.0 : 0.0;
                }
            }
        }
        std::cout << "\rProcessing " << (float)100*(y+1)/3 << "\% complete." << std::flush;
    }
    std::cout << "\n...done!\n" << std::flush;

    // no lead
    std::cout << "Scanning Unleaded Samples...\n" << std::flush;
    int fBulk = 0;
    for(int ii=0; ii<54; ++ii) {
        fBulk++;
        char fname[200];
        sprintf(fname, "%s%i_200.vtk",argv[2], ii);
        image.ImportFromVtk(fname);
        //filtering
        IBVoxFilter_Abtrim trim(Vector3i(5,5,5));
        IBFilterGaussShape shape(0.2);
        trim.SetKernelSpherical(shape);
        trim.SetABTrim(0,2);
        trim.SetImage(&image);
        //trim.Run();
        IBVoxImageScanner<IBVoxCollectionCap> scanner;
        scanner.SetImage(&image);
        RangeThresholdScan::ScanData res = scanner.ScanImage<RangeThresholdScan>(opt);
        for(int j=0; j<opt.size(); ++j) {
            perc_b_t[1][j] += res.at(j).Percent;
            inte_b_t[1][j] += res.at(j).Intensity;
            iden_b_t[1][j] += (res.at(j).Percent>0.f) ? 1.0 : 0.0;
        }
        std::cout << "\rProcessing " << (float)100*(ii+1)/54 << "\% complete." << std::flush;
    }

    std::cout << "\n...done!\n" << std::flush;

    std::cout << "Finalizing Data and Saving..." << std::flush;

    for(int j=0; j<opt.size(); ++j) {
        perc[0] = perc_b_t[0][j] / fbulk;
        inte[0] = inte_b_t[0][j] / fbulk;
        iden[0] = iden_b_t[0][j] / fbulk;
        perc[1] = perc_b_t[1][j] / fBulk;
        inte[1] = inte_b_t[1][j] / fBulk;
        iden[1] = iden_b_t[1][j] / fBulk;
        thres = opt.at(j).Threshold;
        t.Fill();
    }

    t.Write();
    f.Close();
    std::cout << "done!\nExiting!\n" << std::flush;
    return 0;
}

