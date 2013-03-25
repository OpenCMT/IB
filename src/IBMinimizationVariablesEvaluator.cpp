#include "IBMinimizationVariablesEvaluator.h"
#include "IBNormalPlaneMinimizationVariablesEvaluator.h"
#include "IBSimpleTwoViewsMinimizationVariablesEvaluator.h"

IBMinimizationVariablesEvaluator *IBMinimizationVariablesEvaluator::New(IBMinimizationVariablesEvaluator::IBMinVarEvaluatorAlgorithm S)
{
    switch(S){
    case IBMinimizationVariablesEvaluator::NormalPlane:
        return new IBNormalPlaneMinimizationVariablesEvaluator;
    case IBMinimizationVariablesEvaluator::SimpleTwoViews:
        return new IBSimpleTwoViewsMinimizationVariablesEvaluator;
    }
}
