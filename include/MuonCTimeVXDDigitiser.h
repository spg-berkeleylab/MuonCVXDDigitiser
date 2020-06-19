#ifndef MuonCTimeVXDDigitiser_h
#define MuonCTimeVXDDigitiser_h 1

#include "EVENT/LCIO.h"
#include "MuonCVXDDigitiser.h"

class MuonCTimeVXDDigitiser : MuonCVXDDigitiser
{
public:

    virtual Processor* newProcessor() { return new MuonCTimeVXDDigitiser; }

    MuonCTimeVXDDigitiser();

    virtual void processEvent(LCEvent* evt);
protected:
    float _tclick;
    float _window_size;
};

#endif //MuonCTimeVXDDigitiser_h

