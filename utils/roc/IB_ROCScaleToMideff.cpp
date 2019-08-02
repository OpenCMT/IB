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



#include <iostream>
#include <fstream>
#include <algorithm>

#include <Math/Dense.h>

using namespace std;
using namespace uLib;


struct ROCElement : public Vector3f {
    inline float &X() { return this->operator ()(0); }
    inline float &Awo() { return this->operator ()(1); }
    inline float &Owa() { return this->operator ()(2); }
};

typedef std::vector<ROCElement> ROC;

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
    for(ROC::iterator itr = roc.begin(); itr<roc.end(); itr++)
        itr->X() += shift;
}

void scale_roc(ROC &roc, float scale ) {
    for(ROC::iterator itr = roc.begin(); itr<roc.end(); itr++)
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

    // first point
    ROC::iterator itr = roc1.begin();
    while (itr != roc1.end() && itr->Awo() > 0 ) itr++;

    float eff = 0;
    while (eff < 20 && itr >= roc1.begin()) eff += itr--->Awo();
    float begin = itr->X();

    // second point 50%
    itr = roc1.begin();
    while (itr != roc1.end() && itr->Owa() <= 0 ) itr++;

    eff = 0;
    while (eff < 20 && itr < roc1.end()) eff += itr++->Owa();
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

