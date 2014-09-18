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
#include <fstream>
#include <algorithm>

#include <Core/Vector.h>
#include <Math/Dense.h>

using namespace std;
using namespace uLib;


struct ROCElement : public Vector3f {
    inline float &X() { return this->operator ()(0); }
    inline float &Awo() { return this->operator ()(1); }
    inline float &Owa() { return this->operator ()(2); }
};

typedef Vector<ROCElement> ROC;

#define CSVSEP ';'


ROC read_roc_with_header(ifstream &file) {
    ROC roc;
    string line;

    std::getline(file,line); // header

    while ( std::getline(file, line) ) {
        istringstream csvStream(line);

        string col;
        ROCElement elemt;
        int i=0;
        while( i<3 && std::getline(csvStream, col, CSVSEP ) )
            elemt(i++) = atof(col.c_str());
        roc.push_back(elemt);
    }
    return roc;
}


void shift_roc(ROC &roc, float shift ) {
    for(ROC::Iterator itr = roc.begin(); itr<roc.end(); itr++)
        itr->X() += shift;
}

void scale_roc(ROC &roc, float scale ) {
    for(ROC::Iterator itr = roc.begin(); itr<roc.end(); itr++)
        itr->X() *= scale;
}


int main (int argc,char **argv)
{
    struct {
        char *file_in1;
        char *file_out;
        float bulk;
    } arg;

    if(argc==4) {
        arg.file_in1 = argv[1];
        arg.file_out = argv[2];
        arg.bulk = atof(argv[3]);
    }
    else {
        std::cerr << "requires 3 argumets\n";
        exit(1);
    }

    // open file 1
    std::ifstream fin;

    fin.open(arg.file_in1);
    if (!fin.is_open()) exit (1);
    ROC roc1 = read_roc_with_header(fin);
    fin.close();

    // first point 50% Owa
    ROC::Iterator itr = roc1.begin();
    while (itr != roc1.end() && itr->Awo() > 50 ) itr++;
    float begin = itr->X();

    // second point 50% Awo
    itr = roc1.begin();
    while (itr != roc1.end() && itr->Owa() < 50 ) itr++;
    float end = itr->X();


    float scale_factor = arg.bulk / (begin + end) * 2;
    scale_roc( roc1, scale_factor );


    std::ofstream fout;
    fout.open(arg.file_out);


    // Header
    fout << "Scaled " << scale_factor << CSVSEP
         << "Awo" << CSVSEP
         << "Owa" << "\n";


    for(itr = roc1.begin(); itr < roc1.end(); ++itr) {
        fout
                << itr->X() << CSVSEP
                << itr->Awo() << CSVSEP
                << itr->Owa() << "\n";
    }



}

