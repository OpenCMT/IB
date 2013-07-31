#ifndef IBANALYZEREMTRIM_H
#define IBANALYZEREMTRIM_H

#include "IBAnalyzerEM.h"

using namespace uLib;

class IBAnalyzerEMTrim : public IBAnalyzerEM {

    typedef IBAnalyzerEM BaseClass;
public:

    IBAnalyzerEMTrim(IBVoxCollection &voxels);
    ~IBAnalyzerEMTrim();

    void Run(unsigned int iterations, float muons_ratio, float a=0, float b=0);

    void SetMLAlgorithm(IBAnalyzerEMAlgorithm *MLAlgorithm);

private:
    friend class IBAnalyzerEMTrimPimpl;
    class IBAnalyzerEMTrimPimpl *d;
};



#endif // IBANALYZEREMTRIM_H
