#ifndef STABILIZEDTAULEAPMETHOD_H
#define STABILIZEDTAULEAPMETHOD_H

#include "TauLeapMethod.h"

class StabilizedTauLeapMethod: public TauLeapMethod
{
public:
    StabilizedTauLeapMethod(Parameters param_, ChemicalSystem* cs_);
    ~StabilizedTauLeapMethod();
     

protected:
    void print_step_info();
    
    vector<Real> mu;
    vector<Real> nu;
    vector<Real> kappa;
            
    Real damping;
    Real beta;
    unsigned int s_old;
    
    unsigned int coeff_recomputation;
};

#endif /* STABILIZEDTAULEAPMETHOD_H */

