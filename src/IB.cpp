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



#include <config.h>
//#include <RConfigure.h>
#include "IB.h"


//const char *IB::Version::PackageName = PACKAGE_NAME;
//const char *IB::Version::VersionNumber = PACKAGE_VERSION;
const char *IB::Version::PackageName = "mutom";
const char *IB::Version::VersionNumber = "0.2";
const char *IB::Version::Release = "x"; // SVN_REVISION;



void IB::Version::PrintSelf(std::ostream &o)
{

    o << " --------------------------------- \n"
      << " INFN - Muon Tomography Project    \n"
      << " IMAGE BUILDER Ver. " << VersionNumber << "\n"
      << " Release: " << Release << "\n"
      << " --------------------------------- \n";
}
