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


#include <TFile.h>
#include <TTree.h>
#include <TGraph.h>
#include <stdio.h>
#include <fstream>
#include "IBVoxCollectionCap.h"
#include "IBSubImageGrabber.h"
#include "IBVoxImageScanner.h"
#include <IBVoxFilters.h>

#include "testing-prototype.h"

#define CSV_SEPARATOR ";"


namespace Recipes {

//struct Gauss3 {
//    static bool Run(IBVoxCollection *image) {
//        // RECIPE // -------------------------------------------------------- //
////        IBVoxFilter_Abtrim trim(Vector3i(3,3,3));
////        IBVoxFilter_Linear trim(Vector3i(3,3,3));
//        IBVoxFilter_Median trim(Vector3i(3,3,3));
//        IBFilterGaussShape shape(0.2);
//        std::vector <float> values;
//        for (int i=0; i<trim.GetKernelData().GetDims().prod(); ++i) {
//            values.push_back(1.);
//        }
//        trim.SetKernelNumericXZY(values);
////        trim.SetKernelSpherical(shape);
////        trim.SetABTrim(0,1);
//        trim.SetImage(&image);
//        trim.Run();
//        // ------------------------------------------------------------------ //
//        return true;
//    }

struct NoFilter {
    static const char *name() { return "NoFilter"; }
    static bool Run(IBVoxCollection *image) {
        // RECIPE // -------------------------------------------------------- //
        // ------------------------------------------------------------------ //
        return true;
    }
};

struct Gauss3 {
    static const char *name() { return "Gauss3"; }
    static bool Run(IBVoxCollection *image) {
        // RECIPE // -------------------------------------------------------- //
        IBVoxFilter_Linear trim(Vector3i(3,3,3));
        IBFilterGaussShape shape(0.7);
        trim.SetKernelWeightFunction(shape);
        trim.SetImage(image);
        trim.Run();
        // ------------------------------------------------------------------ //
        return true;
    }
};

struct Gauss5 {
    static const char *name() { return "Gauss5"; }
    static bool Run(IBVoxCollection *image) {
        // RECIPE // -------------------------------------------------------- //
        IBVoxFilter_Linear trim(Vector3i(5,5,5));
        IBFilterGaussShape shape(0.7);
        trim.SetKernelWeightFunction(shape);
        trim.SetImage(image);
        trim.Run();
        // ------------------------------------------------------------------ //
        return true;
    }
};


struct Avg {
    static const char *name() { return "Avg"; }
    static bool Run(IBVoxCollection *image) {
        // RECIPE // -------------------------------------------------------- //
        IBVoxFilter_Linear trim(Vector3i(3,3,3));
        std::vector <float> values;
        for (int i=0; i<trim.GetKernelData().GetDims().prod(); ++i) {
            values.push_back(1.);
        }
        trim.SetKernelNumericXZY(values);
        trim.SetImage(image);
        trim.Run();
        // ------------------------------------------------------------------ //
        return true;
    }
};

struct Median {
    static const char *name() { return "Median"; }
    static bool Run(IBVoxCollection *image) {
        // RECIPE // -------------------------------------------------------- //
        IBVoxFilter_Median trim(Vector3i(3,3,3));
        std::vector <float> values;
        for (int i=0; i<trim.GetKernelData().GetDims().prod(); ++i) {
            values.push_back(1.);
        }
        trim.SetKernelNumericXZY(values);
        trim.SetImage(image);
        trim.Run();
        // ------------------------------------------------------------------ //
        return true;
    }
};


struct Trim3 {
    static const char *name() { return "Trim3"; }
    static bool Run(IBVoxCollection *image) {
        // RECIPE // -------------------------------------------------------- //
        IBVoxFilter_Abtrim trim(Vector3i(3,3,3));
        IBFilterGaussShape shape(0.7);
        trim.SetKernelWeightFunction(shape);
        trim.SetABTrim(0,1);
        trim.SetImage(image);
        trim.Run();
        // ------------------------------------------------------------------ //
        return true;
    }
};

struct Trim3u {
    static const char *name() { return "Trim3u"; }
    static bool Run(IBVoxCollection *image) {
        // RECIPE // -------------------------------------------------------- //
        IBVoxFilter_Abtrim trim(Vector3i(3,3,3));
        std::vector <float> values;
        for (int i=0; i<trim.GetKernelData().GetDims().prod(); ++i) {
            values.push_back(1.);
        }
        trim.SetKernelNumericXZY(values);
        trim.SetABTrim(0,1);
        trim.SetImage(image);
        trim.Run();
        // ------------------------------------------------------------------ //
        return true;
    }
};

struct Trim5 {
    static const char *name() { return "Trim5"; }
    static bool Run(IBVoxCollection *image) {
        // RECIPE // -------------------------------------------------------- //
        IBVoxFilter_Abtrim trim(Vector3i(5,5,5));
        IBFilterGaussShape shape(0.7);
        trim.SetKernelWeightFunction(shape);
        trim.SetABTrim(0,2);
        trim.SetImage(image);
        trim.Run();
        // ------------------------------------------------------------------ //
        return true;
    }
};



}





template < class RecipeT >
int process_ROC(int argc, char** argv, int sequence_number=-1)
{
    std::cout << "Initializing containers..." << std::flush;


    RangeThresholdScan::ScanOption opt;
    float start = atof(argv[5]);
    float stop = atof(argv[6]);
    float step = atof(argv[7]);
    int howmany = floor((stop-start)/step);
    for(int i=0; i<howmany; ++i) {
        SimpleThresholdScan::ScanOption o;
        o.Threshold = (1E-6)*(start+i*step);
        opt.push_back(o);
    }
    std::cout << "Scanning thresholds from " << (1E-6)*(start) << " to "<< (1E-6)*(start+howmany*step) << std::endl;
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
    char bnamePerc[2][40],bnameInte[2][40],bnameIden[2][40];
    for(int i=0; i<2; ++i) {
        perc[i] = 0.0;
        inte[i] = 0.0;
        iden[i] = 0.0;
        sprintf(bnamePerc[i], "AveragePercent_%i",        i);
        sprintf(bnameInte[i], "AverageIntensity_%i",      i);
        sprintf(bnameIden[i], "IdentificationPercent_%i", i);
    }

    IBVoxCollection image(Vector3i::Zero());


    std::cout << "done!\n" << std::flush;
    std::cout << "Scanning Leaded Samples...\n" << std::flush;


    ////////////////////////////////////
    ///  L E A D     D A T A S E T   ///
    ////////////////////////////////////
    int fbulk = 0;
    float iron_average_accumulator_1 = 0.;

    char fname[200];
    sprintf(fname, "%s%i_%i.vtk", argv[1], 0, atoi(argv[4]));

    while (image.ImportFromVtk(fname)) {
        fbulk++;
        sprintf(fname, "%s%i_%i.vtk", argv[1], fbulk, atoi(argv[4]));

        //        bool read_status = image.ImportFromVtk(fname);
        //        if(!read_status) {
        //            std::cerr << "ERROR file not found !!";
        //            exit(1);
        //        }

        // FILTER RECIPE //
        RecipeT::Run(&image);

        // GRAB IMAGE //
        IBSubImageGrabber<IBVoxCollectionCap> grabber(image);
        IBLightCollection img(Vector3i::Zero());
        IBLightCollection ref_1(Vector3i::Zero());
        IBLightCollection ref_2(Vector3i::Zero());

        img = grabber.GrabRegion<IBLightCollection>(HPoint3f(-5,-95,-5),//65-95
                                                    HVector3f(55,70,50));//40-70
        ref_1 = grabber.GrabRegion<IBLightCollection>(HPoint3f(-96,-162,-69),
                                                    HVector3f(20,10,20));
        ref_2 = grabber.GrabRegion<IBLightCollection>(HPoint3f(95,-162,65),
                                                      HVector3f(20,10,20));

        // IRON AVERAGE //
        float refdens_1 = 0.;
        float refdens_2 = 0.;
        for (int i = 0; i<ref_1.Data().size(); ++i) {
            refdens_1 += 1E6*ref_1.GetValue(i);
        }
        for (int i = 0; i<ref_2.Data().size(); ++i) {
            refdens_2 += 1E6*ref_2.GetValue(i);
        }

        float voxvol = 4000 / image.GetSpacing().prod();
        //        std::cout << "SS image spacing: " << image.GetSpacing().transpose() << " RefVol: " << voxvol;

        iron_average_accumulator_1 += ((refdens_1 + refdens_2) / (2*voxvol));

        // SCAN THRESHOLD //
        IBVoxImageScanner<IBLightCollection> scanner;
        scanner.SetImage(&img); //
	    RangeThresholdScan::ScanData res = scanner.ScanImage<RangeThresholdScan>(opt);
	    for(int j=0; j<opt.size(); ++j) {
	        perc_b_t[0][j] += res.at(j).Percent;
	        inte_b_t[0][j] += res.at(j).Intensity;
	        iden_b_t[0][j] += (res.at(j).Percent>0.f) ? 1.0 : 0.0;
        }


        std::cout << "\rProcessing image: " << fbulk << " complete." << std::flush;
    }
    iron_average_accumulator_1 /= fbulk;


    float iron_average_accumulator_2 = 0.;
    ////////////////////////////////////
    ///  E M P T Y    D A T A S E T  ///
    ////////////////////////////////////
    std::cout << "Scanning Unleaded Samples...\n" << std::flush;
    int fBulk = 0;


    sprintf(fname, "%s%i_%i.vtk", argv[2], 0, atoi(argv[4]));

    while ( image.ImportFromVtk(fname) ){
        fBulk++;
        sprintf(fname, "%s%i_%i.vtk", argv[2], fBulk, atoi(argv[4]));

        // FILTER RECIPE //
        RecipeT::Run(&image);

        // GRAB IMAGE //
        IBSubImageGrabber<IBVoxCollectionCap> grabber(image);
        IBLightCollection img(Vector3i::Zero());
        IBLightCollection ref_1(Vector3i::Zero());
        IBLightCollection ref_2(Vector3i::Zero());
        img = grabber.GrabRegion<IBLightCollection>(HPoint3f(-5,-95,-5),
                                                    HVector3f(55,70,50));
        ref_1 = grabber.GrabRegion<IBLightCollection>(HPoint3f(-96,-162,-69),
                                                      HVector3f(20,10,20));
        ref_2 = grabber.GrabRegion<IBLightCollection>(HPoint3f(95,-162,65),
                                                      HVector3f(20,10,20));
        // REF DENSITY //
        float refdens_1 = 0.;
        float refdens_2 = 0.;
        for (int i = 0; i<ref_1.Data().size(); ++i) {
            refdens_1 += 1E6*ref_1.GetValue(i);
        }
        for (int i = 0; i<ref_2.Data().size(); ++i) {
            refdens_2 += 1E6*ref_2.GetValue(i);
        }

        float voxvol = 4000 / image.GetSpacing().prod();
        iron_average_accumulator_2 += ((refdens_1 + refdens_2) / (2*voxvol));

        // THRESHOLD //
        IBVoxImageScanner<IBLightCollection> scanner;
        scanner.SetImage(&img);
        RangeThresholdScan::ScanData res = scanner.ScanImage<RangeThresholdScan>(opt);
        for(int j=0; j<opt.size(); ++j) {
            perc_b_t[1][j] += res.at(j).Percent;
            inte_b_t[1][j] += res.at(j).Intensity;
            iden_b_t[1][j] += (res.at(j).Percent>0.f) ? 1.0 : 0.0;
        }

        std::cout << "\rProcessing image: " << fBulk << " complete." << std::flush;
    }
    iron_average_accumulator_2 /= fBulk;


    std::cout << "\n...done!\nExamined " << fBulk << " samples\n\n" << std::flush;
    float norm = 14.2/((iron_average_accumulator_1 + iron_average_accumulator_2) / 2);
    std::cout << "Threshold Normalization Factor: " << norm << "\n\n" << std::flush;




    // FILE SAVE ------------------------------------------------------------ //
    std::cout << "Finalizing Data and Saving..." << std::flush;
    std::ofstream fout;

//    char fname[200];
    sprintf(fname,argv[3],RecipeT::name());
    fout.open( fname );



    //    image_mean[0].ExportToVtkXml( std::string(FileNameRemoveExtension(fname)+"SS.vti").c_str() );
    //    image_mean[1].ExportToVtkXml( std::string(FileNameRemoveExtension(fname)+"NS.vti").c_str() );

    // csv header //
    fout << fname << CSV_SEPARATOR << "Awo" << CSV_SEPARATOR << "Owa\n";
    std::cout << fname << CSV_SEPARATOR << "Awo" << CSV_SEPARATOR << "Owa\n";
    for(int j=0; j<opt.size(); ++j) {
        perc[0] = perc_b_t[0][j] / (float)fbulk;
        inte[0] = inte_b_t[0][j] / (float)fbulk;
        iden[0] = 100*(iden_b_t[0][j] / (float)fbulk);
        perc[1] = perc_b_t[1][j] / (float)fBulk;
        inte[1] = inte_b_t[1][j] / (float)fBulk;
        iden[1] = 100*(iden_b_t[1][j] / (float)fBulk);
        thres = norm * opt.at(j).Threshold;

        fout << 1E6*thres << CSV_SEPARATOR <<  iden[1] << CSV_SEPARATOR << 100-iden[0]   << "\n";
//        std::cout << 1E6*thres << CSV_SEPARATOR <<  (float)iden[1] << CSV_SEPARATOR << (float)(100-iden[0])   << "\n";
        printf("%f == %f == %f \n",1E6*thres , iden[1] ,100-iden[0]);
    }
	fout.close();

    // ---------------------------------------------------------------------- //


    std::cout << "done!\nExiting!\n" << std::flush;
    return 0;
}



int main(int argc, char **argv)
{

    process_ROC<Recipes::NoFilter>(argc,argv);
    process_ROC<Recipes::Gauss3>(argc,argv);
//    process_ROC<Recipes::Gauss5>(argc,argv);
    process_ROC<Recipes::Avg>(argc,argv);
//    process_ROC<Recipes::Median>(argc,argv);
//    process_ROC<Recipes::Trim3s2>(argc,argv);
    process_ROC<Recipes::Trim3u>(argc,argv);
//    process_ROC<Recipes::Trim3>(argc,argv);
//    process_ROC<Recipes::Trim5>(argc,argv);
    return 0;
}
