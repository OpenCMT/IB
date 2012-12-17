
#include "IBAnalyzer.h"

#include "IBAnalyzerPoca.h"
#include "IBAnalyzerEM.h"
#include "IBAnalyzerWPoca.h"



IBAnalyzer *IBAnalyzer::New(IBAnalyzerFactoryType id)
{
    switch(id) {
    case IBAnalyzer::Poca:
        return new IBAnalyzerPoca;
        break;
    case IBAnalyzer::EM:
        return new IBAnalyzerEM;
        break;
    case IBAnalyzer::WPoca:
        return new IBAnalyzerWPoca;

    }
}
