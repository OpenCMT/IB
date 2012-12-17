#ifndef IBTILTEDAXISPOCAEVALUATOR_H
#define IBTILTEDAXISPOCAEVALUATOR_H

#include "IBPocaEvaluator.h"

class IBTiltedAxisPocaEvaluator : public IBPocaEvaluator
{
public:
    IBTiltedAxisPocaEvaluator();
    ~IBTiltedAxisPocaEvaluator();

    bool     evaluate(MuonScatterData muon);
    HPoint3f getPoca();

private:
    friend class IBTiltedAxisPocaEvaluatorPimpl;
    class IBTiltedAxisPocaEvaluatorPimpl *d;
};

#endif // IBTILTEDAXISPOCAEVALUATOR_H
