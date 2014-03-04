#ifndef IBANALYZEREMALGORITHM_H
#define IBANALYZEREMALGORITHM_H

#include "Core/Object.h"
#include "IBAnalyzerEM.h"

using namespace uLib;



class IBAnalyzerEMAlgorithm : public Object {
protected:
    typedef struct IBAnalyzerEM::Event Event;
    uLibTypeMacro(IBAnalyzerEMAlgorithm, uLib::Object)
public:
    properties() {
        Scalarf inertia;
    };

    IBAnalyzerEMAlgorithm() { init_properties(); }

    virtual void evaluate(Matrix4f &Sigma, Event *evc) = 0;

    virtual bool ComputeSigma(Matrix4f &Sigma, Event *evc);

protected:
    virtual ~IBAnalyzerEMAlgorithm() {}

};


inline void IBAnalyzerEMAlgorithm::init_properties() {
    $_init();
    $$.inertia = 1;
}

#endif // IBANALYZEREMALGORITHM_H
