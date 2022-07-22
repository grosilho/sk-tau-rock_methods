#ifndef IMPLICITTAULEAP_H
#define IMPLICITTAULEAP_H

#include "TauLeapMethod.h"

// If the post processing option is enabled then this scheme is the interlaced implicit tau, where 
// a few small steps resolving the fast variables are performed at the end of the simulation.
// Rathinam, M., Petzold, L. R., Cao, Y., and Gillespie, D. T. (2003). Stiffness in stochastic chemically reacting systems: The implicit tau-leaping method. Journal of Chemical Physics, 119(24):12784-12794.
class ImplicitTauLeap: public TauLeapMethod
{
public:
    ImplicitTauLeap(Parameters param_, ChemicalSystem* cs_);
    ~ImplicitTauLeap();
         
    
protected:
    void step();
    
    MatrixXd J, I;
    VectorXd atmp[2];
    
    // Options for interlaced tau
    bool ET_postprocess;     // if true then postprocessing is done with explicit tau, else with implicit tau
    unsigned int n_ET_steps; // number of postprocessing steps
    Real size_ET_steps;      // relative size of the steps, compared to 2/rho
};

#endif /* IMPLICITTAULEAP_H */

