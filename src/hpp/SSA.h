#ifndef SSA_H
#define SSA_H

#include "Solver.h"

class SSA: public Solver
{
public:
    SSA(Parameters param_, ChemicalSystem* cs_);
    ~SSA();
        
    bool run();
    bool run_once_save_path();
    
protected:
    void print_step_info();
    
    int j;         // next reaction index

    Real psum_a; // temporary var for the partial sum of prop fun
    Real r;     //temporary variable holding a sampling
    uniform_real_distribution<Real> unif_dist;
};

#endif /* SSA_H */

