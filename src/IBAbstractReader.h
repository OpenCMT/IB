#ifndef IBABSTRACTREADER_H
#define IBABSTRACTREADER_H

#include <Core/Object.h>

template <typename T>
class IBAbstractReader : public uLib::Object
{
    typedef uLib::Object BaseClass;
public:
    virtual ~IBAbstractReader() {}

    virtual unsigned long GetNumberOfEvents() =0;
    virtual T *GetNext() =0;

protected:
    IBAbstractReader() {}
    IBAbstractReader(IBAbstractReader &copy) {}
};










#endif // IBABSTRACTREADER_H
