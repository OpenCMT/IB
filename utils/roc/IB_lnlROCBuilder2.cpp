
#include <stdio.h>
#include <fstream>
#include <string>
#include <algorithm>


#include "IBVoxCollection.h"
#include "IBVoxFilters.h"
#include "IBSubImageGrabber.h"

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
// ROC


struct ROCElement : public Vector3f {
    ROCElement (){}
    ROCElement(float X, float Awo, float Owa) : Vector3f(X,Awo,Owa) {}
    inline float &X() { return this->operator ()(0); }
    inline float &Awo() { return this->operator ()(1); }
    inline float &Owa() { return this->operator ()(2); }
};

typedef Vector<ROCElement> ROC;

inline ROC read_roc_with_header(std::ifstream &file) {
    ROC roc;
    std::string line;

    std::getline(file,line); // header

    while ( std::getline(file, line) ) {
        std::istringstream csvStream(line);

        std::string col;
        ROCElement elemt;
        int i=0;
        while( i<3 && std::getline(csvStream, col, CSV_SEPARATOR ) )
            elemt(i++) = atof(col.c_str());
        roc.push_back(elemt);
    }
    return roc;
}


inline void shift_roc(ROC &roc, float shift ) {
    for(ROC::Iterator itr = roc.begin(); itr<roc.end(); itr++)
        itr->X() += shift;
}

inline void scale_roc(ROC &roc, float scale ) {
    for(ROC::Iterator itr = roc.begin(); itr<roc.end(); itr++)
        itr->X() *= scale;
}


inline float match_midpt(ROC &roc, float th, float match_pt = 50) {
    // first point match_pt [%]
    ROC::Iterator itr = roc.begin();
    while (itr != roc.end() && itr->Awo() > match_pt ) itr++;
    float begin = itr->X();

    // second point match_pt [%]
    itr = roc.begin();
    while (itr != roc.end() && itr->Owa() < match_pt ) itr++;
    float end = itr->X();

    float scale_factor = th / (begin + end) * 2;
    scale_roc( roc, scale_factor );
    return scale_factor;
}


inline std::fstream&
operator<< (std::fstream& stream, ROC &roc) {
    stream << "X" << CSV_SEPARATOR << "Awo" << CSV_SEPARATOR << "Owa\n";
    for (ROC::Iterator itr = roc.begin(); itr < roc.end(); itr++)
        stream << itr->X() << CSV_SEPARATOR
               << itr->Awo() << CSV_SEPARATOR
               << itr->Owa() << "\n";
    return stream;
}

inline std::ostream&
operator<< (std::ostream& stream, ROC &roc) {
    stream << "X" << CSV_SEPARATOR << "Awo" << CSV_SEPARATOR << "Owa\n";
    for (ROC::Iterator itr = roc.begin(); itr < roc.end(); itr++)
        stream << itr->X() << CSV_SEPARATOR
               << itr->Awo() << CSV_SEPARATOR
               << itr->Owa() << "\n";
    return stream;
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
    "/var/local/data/vtk/ROC/ROC_LNL_1388/Source_12l/vtk_R7_10cm/min_0.5/image_%i_%i.vtk",
    "/var/local/data/vtk/ROC/ROC_LNL_1388/NOSource/vtk_R7_10cm/min_0.5/image_%i_%i.vtk",
    "ROC_r1388_%s_v10_t0.5.csv",
    1000,
    300
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

    ROC roc(p.samples);
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
        for(ROC::Iterator itr = roc.begin(); itr != roc.end(); itr++)
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
        for(ROC::Iterator itr = roc.begin(); itr != roc.end(); itr++)
            itr->Awo() += (image.CountLambdaOverThreshold(itr->X(),
                                                          img_center,
                                                          img_hsize) > 0) ? 1 : 0;

        std::cout << "\rProcessing image: " << fBulk << std::flush;
    }


    for(ROC::Iterator itr = roc.begin(); itr != roc.end(); itr++) {
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

    if(argc == 6) {
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
