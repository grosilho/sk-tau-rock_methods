#ifndef R_TAU_ROCK_H
#define R_TAU_ROCK_H

#include "StabilizedTauLeapMethod.h"

class RtauROCK: public StabilizedTauLeapMethod
{
public:
    RtauROCK(Parameters param_, ChemicalSystem* cs_);
    ~RtauROCK();
     

protected:
    void step();
    bool ensure_stability();
  
    Real alpha;
};

#endif /* R_TAU_ROCK_H */

