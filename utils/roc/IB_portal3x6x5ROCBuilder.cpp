
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
        Vector <float> values;
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
// GLOBAL PARAMETERS


struct Parameters {
    char *file_inTp;
    char *file_inFp;
    char *file_out;
    int  samples;
    int iteration;
} p = {
    "/var/local/data/pubs/IEEE_P/src/m20130220/vtk_R7_5cm/min_1.5/image_%i_tr_%i.vtk",
    "/var/local/data/pubs/IEEE_P/src/m20130214/vtk_R7_5cm/min_1.5/image_%i_tr_%i.vtk",
    "test_ROC.csv",
    500,
    120
};



HPoint3f  img_center(-5,-95,-5);
HVector3f img_hsize(55,70,50);






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
        sprintf(fname, p.file_inTp, index, p.iteration);
        std::cout << "starting autorange with file: " << fname << "\n";
        while (image.ImportFromVtk(fname)) {
            sprintf(fname, p.file_inTp, ++index, p.iteration);
            //IBSubImageGrabber<IBVoxCollection> grabber(image);
            //IBVoxCollection img = grabber.GrabRegion<IBVoxCollection>(img_center,img_hsize);
            RecipeT::Run(&image);
            float imgmax = max_of_image(image);
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
    int y=0;
    sprintf(fname, p.file_inTp, y, p.iteration);
    std::cout << "starting Lead dataset with file: " << fname << "\n";
    while (image.ImportFromVtk(fname)) {
        // FILTER RECIPE //
        RecipeT::Run(&image);
        IBVoxCollection img;

        IBSubImageGrabber<IBVoxCollection> grabber(image);

        for (int x=0; x<6; ++x) {
            for (int z=0; z<3; ++z) {
                fbulk++;
                img = grabber.GrabRegion<IBVoxCollection>(HPoint3f(-250+x*100,
                                                                     -100+y*50,
                                                                     -100+z*100),
                                                            HVector3f(30,25,30));
                for(uLib::IBROC::Iterator itr = roc.begin(); itr != roc.end(); itr++)
                    itr->Owa() += img.CountLambdaOverThreshold(itr->X()) > 0 ? 1 : 0;
                //#               pragma omp parallel for
//                for(uint i=0; i < roc.size(); ++i)
//                    roc[i].Owa() += img.CountLambdaOverThreshold(roc[i].X()) > 0 ? 1 : 0;
            }
        }
        sprintf(fname, p.file_inTp, ++y, p.iteration);
        std::cout << "\rProcessing LEAD: " << y << std::flush;
    }
    std::cout << std::endl;

    ////////////////////////////////////
    //   E M P T Y    D A T A S E T   //
    ////////////////////////////////////

    int fBulk = 0;
    sprintf(fname, p.file_inFp, fBulk, p.iteration);
    std::cout << "starting Scraps dataset with file: " << fname << "\n";
    while ( image.ImportFromVtk(fname) ){
        sprintf(fname, p.file_inFp, ++fBulk, p.iteration);

        // FILTER RECIPE //
        RecipeT::Run(&image);

        // SCAN THRESHOLD //
        for(uLib::IBROC::Iterator itr = roc.begin(); itr != roc.end(); itr++)
            itr->Awo() += (image.CountLambdaOverThreshold(itr->X()) > 0) ? 1 : 0;
        //#       pragma omp parallel for
        //        for(uint i=0; i < roc.size(); ++i)
        //            roc[i].Owa() += image.CountLambdaOverThreshold(roc[i].X()) > 0 ? 1 : 0;

        std::cout << "\rProcessing SCRAP: " << fBulk << std::flush;
    }


    for(uLib::IBROC::Iterator itr = roc.begin(); itr != roc.end(); itr++) {
        itr->X()   *= 1E6;
        itr->Awo() *= 100./fBulk;
        itr->Owa()  = 100 * (1 - itr->Owa()/fbulk);
    }

    // FILE SAVE ------------------------------------------------------------ //
    std::cout << "Finalizing Data and Saving..." << std::flush;

    ofstream fout;
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


int main(int argc, char **argv)
{

    if(argc == 1) {
        std::cout << "working in test mode with default values \n";
    }
    else if(argc == 6) {
        p.file_inTp  = argv[1];
        p.file_inFp  = argv[2];
        p.file_out   = argv[3];
        p.samples    = atof(argv[4]);
        p.iteration  = atof(argv[5]);
    }
    else {
        std::cout << "invalid command. please use this args: sourcefile_%i_%i, nosourcefile_%i_%i, fileoutname, nsamples, iteration\n";
        exit(1);
    }

    // used
//    process_ROC<Recipes::NoFilter>(argc,argv);
//    process_ROC<Recipes::Avg>(argc,argv);
    process_ROC<Recipes::Trim3u>(argc,argv);
//    process_ROC<Recipes::Trim3>(argc,argv);
    process_ROC<Recipes::BilateralTrim>(argc,argv);

    // not used
    //    process_ROC<Recipes::Gauss3>(argc,argv);
    //    process_ROC<Recipes::Gauss5>(argc,argv);
    //    process_ROC<Recipes::Median>(argc,argv);
    //    process_ROC<Recipes::Trim3s2>(argc,argv);
    //    process_ROC<Recipes::Trim5>(argc,argv);
    return 0;
}

