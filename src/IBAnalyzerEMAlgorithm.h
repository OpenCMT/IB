#ifndef IBANALYZEREMALGORITHM_H
#define IBANALYZEREMALGORITHM_H

#include "IBAnalyzerEM.h"

using namespace uLib;





class IBAnalyzerEMAlgorithm {
protected:
    typedef struct IBAnalyzerEM::Event Event;
public:

    virtual void evaluate(Matrix4f &Sigma, Event *evc) = 0;

    virtual bool ComputeSigma(Matrix4f &Sigma, Event *evc);

protected:
    virtual ~IBAnalyzerEMAlgorithm() {}

};



#endif // IBANALYZEREMALGORITHM_H
