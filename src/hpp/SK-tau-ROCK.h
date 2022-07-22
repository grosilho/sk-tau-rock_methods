#ifndef SK_TAU_ROCK_H
#define SK_TAU_ROCK_H

#include "StabilizedTauLeapMethod.h"

class SKtauROCK: public StabilizedTauLeapMethod
{
public:
    SKtauROCK(Parameters param_, ChemicalSystem* cs_);
    ~SKtauROCK();
     

protected:
    void step();
    bool ensure_stability();
    void post_process();
    
    Real alpha;
    
};

#endif /* SK_TAU_ROCK_H */

