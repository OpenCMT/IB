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

#include <Math/VectorSpace.h>

#include "IBVoxCollectionCap.h"
#include "IBVoxFilters.h"
#include "IBSubImageGrabber.h"
#include "IBVoxImageScanner.h"

#include "IBROCBuilder.h"

using namespace uLib;




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
        IBFilterGaussShape shape(0.5);
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
        IBFilterGaussShape shape(0.5);
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
        IBFilterGaussShape shape(0.5);
        trim.SetKernelWeightFunction(shape);
        trim.SetABTrim(0,1);
        trim.SetImage(image);
        trim.Run();
        // ------------------------------------------------------------------ //
        return true;
    }
};

struct Trim3s2 {
    static const char *name() { return "Trim3s2"; }
    static bool Run(IBVoxCollection *image) {
        // RECIPE // -------------------------------------------------------- //
        IBVoxFilter_Abtrim trim(Vector3i(3,3,3));
        IBFilterGaussShape shape(0.2);
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
        IBFilterGaussShape shape(0.5);
        trim.SetKernelWeightFunction(shape);
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
// CLASS

ROCBuilder::ROCBuilder()
{}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// IBROC BUILDER 2



float max_of_image(IBVoxCollection &image) {
    Vector<IBVoxel>::Iterator itr;
    float max = 0;
    for(itr = image.Data().begin(); itr != image.Data().end(); ++itr)
        if(itr->Value > max) max = itr->Value;
    return max;
}





template < class RecipeT >
IBROC ROCBuilder::BuildRoc(Vector<IBVoxCollection> Owa, Vector<IBVoxCollection> Awo)
{
    typedef Vector<IBVoxCollection> ImgSequence;
    static const Vector4f  img_center = HPoint3f(-5,-95,-5);
    static const Vector4f img_hsize = HVector3f(55,70,50);


    std::cout << "Initializing containers..." << std::flush;


    IBVoxCollection image(Vector3i::Zero());

    ////////////////////////////////////
    // Auto Range                     //
    ////////////////////////////////////
    float max = 0;
    {
        for (ImgSequence::Iterator Itr = Owa.begin(); Itr < Owa.end(); Itr++) {
            image = *Itr;
            IBSubImageGrabber<IBVoxCollection> grabber(image);
            IBVoxCollection img = grabber.GrabRegion<IBVoxCollection>(img_center,img_hsize);
            RecipeT::Run(&img);
            float imgmax = max_of_image(img);
            if (imgmax > max) max = imgmax;
        }
    }
    max += max/2;

    IBROC roc($$.samples);
    for(int i=0; i < $$.samples; ++i )
    {
        roc[i].X() = max/$$.samples * i;
        roc[i].Awo() = 0;
        roc[i].Owa() = 0;
    }


    ////////////////////////////////////
    //   L E A D     D A T A S E T    //
    ////////////////////////////////////
    int fbulk = 0;
    for (ImgSequence::Iterator Itr = Owa.begin(); Itr < Owa.end(); Itr++) {
        fbulk ++;
        image = *Itr;

        // FILTER RECIPE //
        RecipeT::Run(&image);

        // SCAN THRESHOLD //
        for(IBROC::Iterator itr = roc.begin(); itr != roc.end(); itr++)
            itr->Owa() += (int)(image.CountLambdaOverThreshold(itr->X(),
                                                         img_center,
                                                         img_hsize) > 0);

        std::cout << "\rProcessing image: " << fbulk << std::flush;
    }

    ////////////////////////////////////
    //   E M P T Y    D A T A S E T   //
    ////////////////////////////////////

    int fBulk = 0;
    for (ImgSequence::Iterator Itr = Awo.begin(); Itr < Awo.end(); Itr++) {
        fBulk ++;
        image = *Itr;

        // FILTER RECIPE //
        RecipeT::Run(&image);

        // SCAN THRESHOLD //
        for(IBROC::Iterator itr = roc.begin(); itr != roc.end(); itr++)
            itr->Awo() += (int)(image.CountLambdaOverThreshold(itr->X(),
                                                         img_center,
                                                         img_hsize) > 0);

        std::cout << "\rProcessing image: " << fBulk << std::flush;
    }


    for(IBROC::Iterator itr = roc.begin(); itr != roc.end(); itr++) {
        itr->X()   *= 1E6;
        itr->Awo() *= 100./fBulk;
        itr->Owa()  = 100 * (1 - itr->Owa()/fbulk);
    }

    return roc;
}





IBROC ROCBuilder::BuildRoc(Vector<IBVoxCollection> Owa,
                         Vector<IBVoxCollection> Awo,
                         ROCBuilder::ROCRecipeEnum recipe)
{
    IBROC ret;
    switch(recipe) {
    case ROCBuilder::Gauss3:
        ret = BuildRoc<Recipes::Gauss3>(Owa,Awo);
        break;
    case ROCBuilder::Gauss5:
        ret = BuildRoc<Recipes::Gauss5>(Owa,Awo);
        break;
    case ROCBuilder::Avg:
        ret = BuildRoc<Recipes::Avg>(Owa,Awo);
        break;
    case ROCBuilder::Median:
        ret = BuildRoc<Recipes::Median>(Owa,Awo);
        break;
    case ROCBuilder::Trim3u:
        ret = BuildRoc<Recipes::Trim3u>(Owa,Awo);
        break;
    case ROCBuilder::Trim3:
        ret = BuildRoc<Recipes::Trim3>(Owa,Awo);
        break;
    case ROCBuilder::Trim5:
        ret = BuildRoc<Recipes::Trim5>(Owa,Awo);
        break;
    case ROCBuilder::NoFilter:
        ret = BuildRoc<Recipes::NoFilter>(Owa,Awo);
        break;
    }

    return ret;
}



float ROCBuilder::Ratio(IBROC roc, float y)
{
    // first point y [%] Owa
    IBROC::Iterator itr = roc.begin();
    while (itr != roc.end() && itr->Awo() > y ) itr++;
    float begin = itr->X();

    // second point y [%] Awo
    itr = roc.end()-1;
    while (itr != roc.begin() && itr->Owa() > y ) itr--;
    float end = itr->X();

    return end / begin;
}


float ROCBuilder::FSI(IBROC roc, float y)
{
    // first point y [%] Owa
    IBROC::Iterator itr = roc.begin();
    while (itr != roc.end() && itr->Awo() > y ) itr++;
    float begin = itr->X();

    // second point y [%] Awo
    itr = roc.end()-1;
    while (itr != roc.begin() && itr->Owa() > y ) itr--;
    float end = itr->X();

    std::cout << "end " << end <<  ", begin " << begin << std::endl;
    return end - begin;
}

float ROCBuilder::AUC(IBROC &roc)
{
    float AUC = 0;
    for(IBROC::Iterator itr=roc.begin(); itr<roc.end()-1; itr++) {
        AUC += (itr->Awo() - (itr+1)->Awo()) * (100 - itr->Owa() + 100 - (itr+1)->Owa()) / 2E4;
    }
    return AUC;
}

