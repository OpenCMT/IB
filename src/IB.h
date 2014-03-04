/* 
 * File:   IB.h
 * Author: andrea
 *
 * Created on June 8, 2011, 3:13 PM
 */


// list of IB Headers for fast inclusion //

#include <Core/Macros.h>
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
