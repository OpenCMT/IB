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
#include <fstream>
#include <string>
#include <algorithm>


#include "IBVoxCollection.h"
#include "IBVoxFilters.h"
#include "IBSubImageGrabber.h"
#include "IBROC.h"
#include "testing-prototype.h"

#define CSV_SEPARATOR ';'

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// GLOBAL PARAMETERS


struct Parameters {
    char *file_inTp;
    char *file_inFp;
    char *file_out;
    int  samples;
    int iteration;
    std::string algs;
} p = {
    (char*)"/var/local/data/vtk/ROC/ROC_LNL_1388/Source_12l/vtk_R7_10cm/min_0.5/image_%i_%i.vtk",
    (char*)"/var/local/data/vtk/ROC/ROC_LNL_1388/NOSource/vtk_R7_10cm/min_0.5/image_%i_%i.vtk",
    (char*)"ROC_r1388_%s_v10_t0.5.csv",
    1000,
    300,
    (char*)"NoFilter"
};



HPoint3f  img_center(-5,-95,-5);
HVector3f img_hsize(55,70,50);




////////////////////////////////////////////////////////////////////////////////
// RECIPES //



namespace Recipes {

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
        Vector <float> values;
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
        Vector <float> values;
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
        //        Vector <float> values;
        //        for (int i=0; i<trim.GetKernelData().GetDims().prod(); ++i) {
        //            values.push_back(1.);
        //        }
        //        trim.SetKernelNumericXZY(values);
        IBFilterGaussShape shape(1);
        trim.SetKernelWeightFunction(shape);
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

struct BilateralTrim {
    static const char *name() { return "BilateralTrim"; }
    static bool Run(IBVoxCollection *image) {
        // RECIPE // -------------------------------------------------------- //
        IBVoxFilter_BilateralTrim trim(Vector3i(5,5,5));
        IBFilterGaussShape shape(0.7);
        trim.SetKernelWeightFunction(shape);
        trim.SetIntensitySigma(30);
        trim.SetABTrim(0,2);
        trim.SetImage(image);
        trim.Run();
        // ------------------------------------------------------------------ //
        return true;
    }
};



}






////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// ROC builder



float max_of_image(IBVoxCollection &image) {
    Vector<IBVoxel>::Iterator itr;
    float max = 0;
    for(itr = image.Data().begin(); itr != image.Data().end(); ++itr)
        if(itr->Value > max) max = itr->Value;
    return max;
}

template < class RecipeT >
int process_ROC(int argc, char** argv, int sequence_number=-1)
{
    std::cout << "Initializing containers..." << std::flush;
    char fname[200];


    IBVoxCollection image(Vector3i::Zero());

    ////////////////////////////////////
    // Auto Range                     //
    ////////////////////////////////////
    float max = 0;
    {
        int index = 0;
        sprintf(fname, p.file_inTp, index++, p.iteration);
        while (image.ImportFromVtk(fname)) {
            sprintf(fname, p.file_inTp, index++, p.iteration);
            IBSubImageGrabber<IBVoxCollection> grabber(image);
            IBVoxCollection img = grabber.GrabRegion<IBVoxCollection>(img_center,img_hsize);
            RecipeT::Run(&img);
            float imgmax = max_of_image(img);
            if (imgmax > max) max = imgmax;
        }
    }
    max += max;

    uLib::IBROC roc(p.samples);
    for(int i=0; i < p.samples; ++i )
    {
        roc[i].X() = max/p.samples * i;
        roc[i].Awo() = 0;
        roc[i].Owa() = 0;
    }


    ////////////////////////////////////
    //   L E A D     D A T A S E T    //
    ////////////////////////////////////
    int fbulk = 0;
    sprintf(fname, p.file_inTp, fbulk, p.iteration);
    while (image.ImportFromVtk(fname)) {
        sprintf(fname, p.file_inTp, fbulk++, p.iteration);

        // FILTER RECIPE //
        RecipeT::Run(&image);

        // SCAN THRESHOLD //
        for(uLib::IBROC::Iterator itr = roc.begin(); itr != roc.end(); itr++)
            itr->Owa() += (image.CountLambdaOverThreshold(itr->X(),
                                                         img_center,
                                                         img_hsize)) > 0 ? 1 : 0;

        std::cout << "\rProcessing image: " << fbulk << std::flush;
    }

    ////////////////////////////////////
    //   E M P T Y    D A T A S E T   //
    ////////////////////////////////////

    int fBulk = 0;
    sprintf(fname, p.file_inFp, fBulk, p.iteration);
    while ( image.ImportFromVtk(fname) ){
        sprintf(fname, p.file_inFp, fBulk++, p.iteration);

        // FILTER RECIPE //
        RecipeT::Run(&image);

        // SCAN THRESHOLD //
        for(uLib::IBROC::Iterator itr = roc.begin(); itr != roc.end(); itr++)
            itr->Awo() += (image.CountLambdaOverThreshold(itr->X(),
                                                          img_center,
                                                          img_hsize) > 0) ? 1 : 0;

        std::cout << "\rProcessing image: " << fBulk << std::flush;
    }


    for(uLib::IBROC::Iterator itr = roc.begin(); itr != roc.end(); itr++) {
        itr->X()   *= 1E6;
        itr->Awo() *= 100./fBulk;
        itr->Owa()  = 100 * (1 - itr->Owa()/fbulk);
    }    

    roc.Samples() << fBulk , fbulk;

    // FILE SAVE ------------------------------------------------------------ //
    std::cout << "Finalizing Data and Saving..." << std::flush;

    std::ofstream fout;
    sprintf(fname,p.file_out,RecipeT::name());
    fout.open( fname );
    fout << roc;
    fout.close();

    // ---------------------------------------------------------------------- //

    std::cout << "done!\nExiting!\n" << std::flush;
    return 0;
}







////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// MAIN

#include "boost/tokenizer.hpp"

#include "Core/Options.h"


int main(int argc, char **argv)
{



    Options opt("usage: sourcefile_%i_%i, nosourcefile_%i_%i, fileoutname  ( %i_%i is: imagesample_iteration )");
    opt.add_options()
            ("help", "get this help and exit")
            ("samples", &p.samples, "samples to quantize ROC threshold")
            ("iteration", &p.iteration, "iteration to run")
            ("algs",&p.algs, "Algorithms to use");

    opt.parse_command_line(argc,argv);

    boost::tokenizer<> tok(p.algs);
    for(boost::tokenizer<>::iterator beg = tok.begin(); beg!=tok.end(); ++beg){
        std::cout << *beg << "\n";
    }

    if(argc >= 4) {
        p.file_inTp  = argv[1];
        p.file_inFp  = argv[2];
        p.file_out   = argv[3];
    }
    else {
        std::cerr << "error parsing command line: use --help \n";
        exit(1);
    }

    exit(0);

    // used
    process_ROC<Recipes::NoFilter>(argc,argv);
    process_ROC<Recipes::Avg>(argc,argv);
    process_ROC<Recipes::Trim3u>(argc,argv);
    process_ROC<Recipes::Trim3>(argc,argv);
    process_ROC<Recipes::BilateralTrim>(argc,argv);

    // not used
    //    process_ROC<Recipes::Gauss3>(argc,argv);
    //    process_ROC<Recipes::Gauss5>(argc,argv);
    //    process_ROC<Recipes::Median>(argc,argv);
    //    process_ROC<Recipes::Trim3s2>(argc,argv);
    //    process_ROC<Recipes::Trim5>(argc,argv);
    return 0;
}
