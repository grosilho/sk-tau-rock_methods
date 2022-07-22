#ifndef TAULEAP_H
#define TAULEAP_H

#include "TauLeapMethod.h"

class TauLeap: public TauLeapMethod
{
public:
    TauLeap(Parameters param_, ChemicalSystem* cs_);
    ~TauLeap();
     

protected:
    void step();
    bool ensure_stability();
};

#endif /* TAULEAP_H */

