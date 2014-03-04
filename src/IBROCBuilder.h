#ifndef IBROCBUILDER_H
#define IBROCBUILDER_H

#include <string>
#include <algorithm>
#include <fstream>

#include "Core/Object.h"
#include "IBVoxCollection.h"

#include "IBROC.h"

using namespace uLib;


class ROCBuilder : public Object
{
    uLibTypeMacro(ROCBuilder,Object)
public:
    properties()
    {
        float start; // actually not used
        float stop;  // actually not used
        unsigned int samples;
    };

public:
    ROCBuilder();

    enum ROCRecipeEnum {
        NoFilter = 0,
        Gauss3,
        Gauss5,
        Avg,
        Median,
        Trim3u,
        Trim3,
        Trim5
    };

    template < class RecipeT >
    IBROC BuildRoc(Vector<IBVoxCollection> Owa, Vector<IBVoxCollection> Awo);

    IBROC BuildRoc(Vector<IBVoxCollection> Owa, Vector<IBVoxCollection> Awo, ROCRecipeEnum recipe = NoFilter);

    float Ratio(IBROC roc, float y);

    float FSI(IBROC roc, float y);

    float AUC(IBROC &roc);

};


inline void ROCBuilder::init_properties() {
    $_init()
    $$.start = 0;
    $$.stop  = 100;
    $$.samples  = 1000;
}

#endif // IBROCBUILDER_H
