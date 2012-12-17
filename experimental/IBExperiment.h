#ifndef IBEXPERIMENT_H
#define IBEXPERIMENT_H

#include <Core/Object.h>
#include <Core/Vector.h>
#include <Detectors/ChamberDetector.h>

using namespace uLib;

class IBExperiment : public Object {
public:

    uLibRefMacro(Chambers,Vector<ChamberDetector>)

private:
    Vector<ChamberDetector> m_Chambers;
};


#endif // IBEXPERIMENT_H
