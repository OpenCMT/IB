#include "IBPocaEvaluator.h"
#include "IBTiltedAxisPocaEvaluator.h"
#include "IBLineDistancePocaEvaluator.h"

using namespace uLib;

IBPocaEvaluator *IBPocaEvaluator::New(IBPocaEvaluator::IBPocaEvaluationAlgorithms S)
{
    switch (S) {
    case TiltedAxis:
        return new IBTiltedAxisPocaEvaluator;
    case LineDistance:
        return new IBLineDistancePocaEvaluator;
    }
}
