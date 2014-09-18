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




#include <TFile.h>
#include <TTree.h>

#include <iostream>


#include "testing-prototype.h"
#include "IBMuonEventTTreeReader.h"

using namespace uLib;

enum Chamber {
    ChUp    = 1 << 3,
    ChDown  = 1 << 2,
    ChLeft  = 1 << 1,
    ChRight = 1 << 0
};


class R3DmReaderFlags : public uLib::Flags4B {
    typedef uLib::Flags4B BaseClass;
public:
    R3DmReaderFlags() : BaseClass() {}
    R3DmReaderFlags(int i) : BaseClass(i) {}
    R3DmReaderFlags(const R3DmReaderFlags &f) : BaseClass(f) {}

    inline static int In(int f)   { return (f & 0xf) << 4; }
    inline static int Out(int f)  { return f & 0xf; }
    inline static int Both(int f) { f &= 0xf; return f | (f<<4); }

    void  PrintSelf(std::ostream &o)  {
        int f = this->operator ()() & 0xff;
        o << " In/Out | U D L R U D L R |\n";
        o << " Code:  [";
        for ( int i = 8 ; i-->0 ;) {
            o << " " << ( (f & (1<<i) ) > 0);
        }
        o << " ]\n";
    }
};



int main()
{
    BEGIN_TESTING(IBReader);

    // some test on Flags for Track selection //
    typedef R3DmReaderFlags Flags;
    R3DmReaderFlags fla;
    fla = Flags::In(ChUp) | Flags::Out(ChDown);
    std::cout << fla() << "\n";
    fla.PrintSelf(std::cout);


    TEST1( fla.testFlag(ChUp) );
    TEST1( fla.testFlag(fla.Out(ChDown)));
    TEST0( fla.testFlag(ChLeft));
    TEST0( fla.testFlag(ChRight));

    fla = fla.Both(ChUp | ChDown);
    std::cout << fla() << "\n";
    fla.PrintSelf(std::cout);


    std::cout << fla.testFlag(fla.In(ChLeft)) << "\n";



    TFile f("/var/local/data/root/run_PDfit_201210/muSteel_PDfit_20121005_v2.root");
    TTree * t = (TTree*)f.Get("n");
    TFile f2("/var/local/data/root/run_10000x/muSteel_20120924.root");
    TTree * t2 = (TTree*)f2.Get("n");

    IBMuonEventTTreeReader* reader = IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::R3D_MC);
    IBMuonEventTTreeReader* reader_old = IBMuonEventTTreeReader::New(IBMuonEventTTreeReader::R2D_MC);

    reader->setTTree(t);
    std::cout << "TTree R3D initialized!\n";
    reader_old->setTTree(t2);
    reader_old->setHitCuts(4,4);
    std::cout << "TTree R2D initialized!\n";
    std::cout << "There are " << reader->getNumberOfEvents() << " events inside R3D! \n";
    std::cout << "There are " << reader_old->getNumberOfEvents() << " events inside R2D! \n";


    MuonScatter* the_event = new MuonScatter();
    for (int i = 0; i < 10; ++i) {
        if (reader->readNext(the_event)) {
            std::cout << "Event Read: " << reader->getCurrentPosition() << "\n" \
                         " -> IN \n" \
                         "     [Origin]   : " << the_event->LineIn().origin.transpose() << "\n" <<
                         "     [Direction]: " << the_event->LineIn().direction.transpose() << "\n" <<
                         " -> OUT \n" \
                         "     [Origin]   : " << the_event->LineOut().origin.transpose() << "\n" <<
                         "     [Direction]: " << the_event->LineOut().direction.transpose() << "\n" ; //<<
//                         " -> Momentum : " << the_event->Momentum() << "\n";
        }
    }
    std::cout << "---------------------------------------------------\n";
    MuonScatter* the_event2 = new MuonScatter();
    for (int i = 0; i < 10; ++i) {
        if (reader_old->readNext(the_event2)) {
            std::cout << "Event Read: " << reader_old->getCurrentPosition() << "\n" \
                         " -> IN \n" \
                         "     [Origin]   : " << the_event2->LineIn().origin.transpose() << "\n" <<
                         "     [Direction]: " << the_event2->LineIn().direction.transpose() << "\n" <<
                         " -> OUT \n" \
                         "     [Origin]   : " << the_event2->LineOut().origin.transpose() << "\n" <<
                         "     [Direction]: " << the_event2->LineOut().direction.transpose() << "\n";
//                         " -> Momentum : " << the_event2->Momentum() << "\n";
        }
    }

    END_TESTING;
}



