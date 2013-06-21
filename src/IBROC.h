#ifndef IBROC_H
#define IBROC_H


#include <stdio.h>
#include <fstream>
#include <string>
#include <algorithm>
#include "Math/Dense.h"
#include "Core/Vector.h"

#define CSV_SEPARATOR ';'

namespace uLib {

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// ROC


struct ROCElement : public Vector3f {
    ROCElement (){}
    ROCElement(float X, float Awo, float Owa) : Vector3f(X,Awo,Owa) {}
    inline float &X() { return this->operator ()(0); }
    inline float &Awo() { return this->operator ()(1); }
    inline float &Owa() { return this->operator ()(2); }
};


class IBROC : public Vector<ROCElement>{
  public:

    virtual void Update() {}

    virtual void read_roc_with_header(const char *file_name) {
        std::ifstream file;
        file.open(file_name);
        this->read_roc_with_header(file);
    }

    virtual void read_roc_with_header(std::ifstream &file) {
        std::string line;
        std::string col;

        std::getline(file, line);
        std::istringstream csvStream(line);

        // header
        int j=0;
        while( j<3 && std::getline(csvStream, col, CSV_SEPARATOR ) )
            m_labels[j++] = col;

        while ( std::getline(file, line) ) {
            std::istringstream csvStream(line);
            ROCElement elemt;
            int i=0;
            while( i<3 && std::getline(csvStream, col, CSV_SEPARATOR ) )
                elemt(i++) = atof(col.c_str());
            this->push_back(elemt);
        }

        this->Update();
    }

    Vector2f GetRange(float y = 100 ) {
        IBROC::Iterator itr = this->begin();
        while (itr != this->end() && itr->Awo() >= y ) itr++;
        float begin = itr->X();
        itr = this->end()-1;
        while (itr != this->begin() && itr->Owa() >= y ) itr--;
        float end = itr->X();
        return Vector2f(begin,end);
    }


    inline void shift_roc(float shift ) {
        for(IBROC::Iterator itr = this->begin(); itr<this->end(); itr++)
            itr->X() += shift;
        this->Update();
    }

    inline void scale_roc(float scale ) {
        for(IBROC::Iterator itr = this->begin(); itr<this->end(); itr++)
            itr->X() *= scale;
        this->Update();
    }

    inline float match_midpt( float th, float match_pt = 50) {
        // first point match_pt [%]
        Vector2f range = this->GetRange(match_pt);
        float scale_factor = th / range.sum() * 2;
        scale_roc(scale_factor);
        return scale_factor;
    }

public:
    std::string m_labels[3];
};

inline std::fstream&
operator<< (std::fstream& stream, IBROC &roc) {
    stream << "X" << CSV_SEPARATOR << "Awo" << CSV_SEPARATOR << "Owa\n";
    for (IBROC::Iterator itr = roc.begin(); itr < roc.end(); itr++)
        stream << itr->X() << CSV_SEPARATOR
               << itr->Awo() << CSV_SEPARATOR
               << itr->Owa() << "\n";
    return stream;
}

inline std::ostream&
operator<< (std::ostream& stream, IBROC &roc) {
    stream << "X" << CSV_SEPARATOR << "Awo" << CSV_SEPARATOR << "Owa\n";
    for (IBROC::Iterator itr = roc.begin(); itr < roc.end(); itr++)
        stream << itr->X() << CSV_SEPARATOR
               << itr->Awo() << CSV_SEPARATOR
               << itr->Owa() << "\n";
    return stream;
}


} // uLib

#endif // IBROC_H
