#ifndef TAU_ROCK_H
#define TAU_ROCK_H

#include "StabilizedTauLeapMethod.h"

class tauROCK: public StabilizedTauLeapMethod
{
public:
    tauROCK(Parameters param_, ChemicalSystem* cs_);
    ~tauROCK();
     

protected:
    void step();
    bool ensure_stability();
    
};

#endif /* TAU_ROCK_H */

