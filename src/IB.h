/* 
 * File:   IB.h
 * Author: andrea
 *
 * Created on June 8, 2011, 3:13 PM
 */


// list of IB Headers for fast inclusion //

#include <Core/Macros.h>
#include <Core/Object.h>

using namespace uLib;

namespace IB {

class Version {
public:
    static const char *PackageName;
    static const char *VersionNumber;
    static const char *Release;
};


}
