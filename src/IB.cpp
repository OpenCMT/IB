/*//////////////////////////////////////////////////////////////////////////////
// CMT Cosmic Muon Tomography project //////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

  Copyright (c) 2014, Universita' degli Studi di Padova, INFN sez. di Padova

  Coordinators: Prof. Gianni Zumerle < gianni.zumerle@pd.infn.it >
                Paolo Checchia       < paolo.checchia@pd.infn.it >

  Authors: Andrea Rigoni Garola < andrea.rigoni@pd.infn.it >
           Matteo Furlan        < nuright@gmail.com >
           Sara Vanini          < sara.vanini@pd.infn.it >

  All rights reserved
  ------------------------------------------------------------------

  This file can not be copied and/or distributed without the express
  permission of  Prof. Gianni Zumerle  < gianni.zumerle@pd.infn.it >

//////////////////////////////////////////////////////////////////////////////*/




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
