#ifndef IBLINEDISTANCEPOCAEVALUATOR_H
#define IBLINEDISTANCEPOCAEVALUATOR_H

#include "IBPocaEvaluator.h"

class IBLineDistancePocaEvaluator : public IBPocaEvaluator
{
public:
    IBLineDistancePocaEvaluator();
    ~IBLineDistancePocaEvaluator();

    bool evaluate(MuonScatterData muon);
    HPoint3f getPoca();
    HPoint3f getOutTrackPoca();
    HPoint3f getInTrackPoca();
    void setDistanceCut(Scalarf length);
private:
    friend class IBLineDistancePocaEvaluatorPimpl;
    class IBLineDistancePocaEvaluatorPimpl *d;
};

#endif // IBLINEDISTANCEPOCAEVALUATOR_H
