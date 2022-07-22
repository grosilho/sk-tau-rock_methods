#ifndef TAULEAPMETHOD_H
#define TAULEAPMETHOD_H

#include "Solver.h"

class TauLeapMethod: public Solver
{
public:
    TauLeapMethod(Parameters param_, ChemicalSystem* cs_);
    virtual ~TauLeapMethod();

    bool run();
    bool run_once_save_path();
    
protected:
    virtual void step() =0;
    VectorXd f(VectorXd& Y);
    VectorXd Qdt(VectorXd& Y, Real dt);
    
    virtual bool ensure_stability();
    void sample_K(const VectorXd& means);
    void compute_Jacobian(MatrixXd& J, VectorXd& Xeval, 
                          VectorXd& Xtmp, VectorXd atmp[2]);
    bool rho();
    virtual void print_step_info();
    
    void init_eigvec_tmp_vec(bool init_eigvec);
    bool isNan(const VectorXd& vec);
    
    unsigned int n_tmp_vec;
    vector<VectorXd> tmp_vec;
    VectorXd eigvec;
    Real eigmax;
    unsigned int rho_or_Newton_iter;
    
    unsigned int s;
    
    const Real tau_max; // upper limit for tau, defined by the user
    VectorXd K; // vector holding the Poisson samplings
    
    bool last; // indicator that tells if it is the last step
    
    poisson_distribution<int> poisson_dist;
        
};

#endif /* TAULEAPMETHOD_H */

