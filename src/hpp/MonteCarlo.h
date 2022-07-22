#ifndef MONTECARLO_H
#define MONTECARLO_H

#include "headers.h"
#include "Parameters.h"
#include "Solver.h"

class MonteCarlo
{
public:
    MonteCarlo(Parameters param_);
    ~MonteCarlo();
    
    void run();
    void import_process_data();
   
protected:
    void process_data(vector<bool>& succeed,vector<Real>& tau_mean_vec,
                      vector<Real>& s_mean_vec, vector<Real>& rho_Nit_mean_vec,
                      vector<Real>& damping_mean_vec);
    void print_data(string solver_name, bool post_proc, Real Ntol);
    void write_data(vector<bool>& succeed);
    void compute_print_errors();
    
    void import_data(string filename, vector<VectorXd>& Y);
    void read_solver_info(string filename, Real& time, unsigned int& iter,
                 unsigned int& succeed_iter, string& solver_name_str,
                 Real& tau_mean, Real& s_mean, Real& damping_mean, Real& rho_Nit_mean,
                 bool& post_proc, Real& Ntol);
    void read_statistics(string filename, ArrayXd& avg, ArrayXd& std);
    
    void print_solver_info(Real time, unsigned int iter, unsigned int succ_iter, 
                           string solver_name, Real tau_max, Real s_mean, Real damping_mean,
                           Real rho_Nit_mean, bool postproc, Real Ntol);
    void print_statistics(const ArrayXd& avg, const ArrayXd& std);
    
    void compute_density_distance_area(vector<unsigned int>& n_bins, ArrayXd& dda,
        vector<VectorXd>& X1, vector<VectorXd>& X2);
    void estimate_self_distances(vector<unsigned int>& n_bins, unsigned int N1, unsigned int N2,
                                 ArrayXd& sd1, ArrayXd& sd2);
    
    Real get_cpu_time();
    
    Parameters param;
    unsigned int MCiter;
    vector<VectorXd> X;
    Real tau_mean;
    Real s_mean;
    Real rho_Nit_mean;
    Real damping_mean;
    
    ArrayXd mean;
    ArrayXd std_dev;
    unsigned int succeed_MCiter;
    
    Real elapsed_time;
};

#endif /* MONTECARLO_H */
