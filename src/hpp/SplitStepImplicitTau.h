#ifndef SPLITSTEPIMPLICITTAU_H
#define SPLITSTEPIMPLICITTAU_H

#include "TauLeapMethod.h"

// This scheme comes from the reference below, but here it is not multi level monte carlo.
// It is similar to Implicit-tau but it first does an implicit Euler step (only drift, no noise) 
// obtaining a solution Y and then to Y it adds the zero-mean noise evaluated at Y.
// Hammouda, C. B., Moraes, A., and Tempone, R. (2017). Multilevel hybrid split-step implicit tau-leap. Numer. Algorithms , 74(2):527-560
// We added as well a post-processing procedure, the same as implicit tau.
class SplitStepImplicitTau: public TauLeapMethod
{
public:
    SplitStepImplicitTau(Parameters param_, ChemicalSystem* cs_);
    ~SplitStepImplicitTau();
         
    
protected:
    void step();
    
    MatrixXd J, I;
    VectorXd atmp[2];
    
    // Options for interlaced tau
    unsigned int n_ET_steps; // number of explicit tau steps in the interlaced implicit tau method
    Real size_EE_steps;      // relative size of the steps, compared to maximal step size 2/rho
};

#endif /* SPLITSTEPIMPLICITTAU_H */

