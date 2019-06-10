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



/*
 * File:   IB.h
 * Author: andrea
 *
 * Created on June 8, 2011, 3:13 PM
 */


// list of IB Headers for fast inclusion //

#include <Core/Object.h>
#include <Detectors/MuonScatter.h>


namespace IB {

class Version {
public:
    static const char *PackageName;
    static const char *VersionNumber;
    static const char *Release;

    static void PrintSelf(std::ostream &o);
};


//class TestObject : public Object {
//    uLibTypeMacro(TestObject,Object)
//public:
//    properties() {
//        int a;
//        int b;
//    };

//};

//inline void TestObject::init_properties() {
//    $_init();
//    $$.a = 123;
//    $$.b = 5552368;
//}


}
