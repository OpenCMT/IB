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



#include <iostream>
#include <assert.h>

#include <root/TFile.h>
#include <root/TTree.h>


#include <Detectors/MuonScatter.h>
#include "IBMuonTTreeReader.h"


#include "testing-prototype.h"


using namespace uLib;

int main()
{

  BEGIN_TESTING(IBMuon);

  TFile* hbfile = new TFile("/var/local/data/root/run_1363/muRadio_1363.root");
  if (hbfile->IsZombie()) {
      exit(-1);
  }
  TTree* tree = (TTree*) hbfile->Get("RADMU");
  assert(tree);

  IBMuonReader reader(tree,Demo_LNL);

  std::cout << "read->num of events = " << reader.GetNumberOfEvents() << "\n";

  MuonScatter *muon = reader.GetNext();





  END_TESTING;
}




