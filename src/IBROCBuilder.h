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
public:
    float start; // actually not used
    float stop;  // actually not used
    unsigned int samples;

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
    IBROC BuildRoc(std::vector<IBVoxCollection> Owa, std::vector<IBVoxCollection> Awo);

    IBROC BuildRoc(std::vector<IBVoxCollection> Owa, std::vector<IBVoxCollection> Awo, ROCRecipeEnum recipe = NoFilter);

    float Ratio(IBROC roc, float y);

    float FSI(IBROC roc, float y);

    float AUC(IBROC &roc);

};


#endif // IBROCBUILDER_H
