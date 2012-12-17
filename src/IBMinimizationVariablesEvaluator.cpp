#include "IBMinimizationVariablesEvaluator.h"
#include "IBNormalPlaneMinimizationVariablesEvaluator.h"

IBMinimizationVariablesEvaluator *IBMinimizationVariablesEvaluator::New(IBMinimizationVariablesEvaluator::IBMinVarEvaluatorAlgorithm S)
{
    switch(S){
    case IBMinimizationVariablesEvaluator::NormalPlane:
        return new IBNormalPlaneMinimizationVariablesEvaluator;

    }
}
