
#include <config.h>
#include "IB.h"


const char *IB::Version::PackageName = PACKAGE_NAME;
const char *IB::Version::VersionNumber = PACKAGE_VERSION;
const char *IB::Version::Release = SVN_REVISION;



void IB::Version::PrintSelf(std::ostream &o)
{

    o << " --------------------------------- \n"
      << " INFN - Muon Tomography Project    \n"
      << " IMAGE BUILDER Ver. " << VersionNumber << "\n"
      << " Release: " << Release << "\n"
      << " --------------------------------- \n";
}
