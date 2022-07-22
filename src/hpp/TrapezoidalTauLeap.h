#ifndef TRAPEZOIDALTAULEAP_H
#define TRAPEZOIDALTAULEAP_H

#include "TauLeapMethod.h"

class TrapezoidalTauLeap: public TauLeapMethod
{
public:
    TrapezoidalTauLeap(Parameters param_, ChemicalSystem* cs_);
    ~TrapezoidalTauLeap();
        
    
protected:
    void step();
    
    MatrixXd J, I;
    VectorXd atmp[2];
};

#endif /* TRAPEZOIDALTAULEAP_H */

