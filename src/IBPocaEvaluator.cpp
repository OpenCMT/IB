#include "IBPocaEvaluator.h"
#include "IBTiltedAxisPocaEvaluator.h"

using namespace uLib;

IBPocaEvaluator *IBPocaEvaluator::New(IBPocaEvaluator::IBPocaEvaluationAlgorithms S)
{
    switch (S) {
    case TiltedAxis:
        return new IBTiltedAxisPocaEvaluator;
    }
}
