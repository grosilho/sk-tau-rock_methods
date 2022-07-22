#include "Solver.h"

#ifdef _OPENMP
#include <omp.h>
#endif

Solver::Solver(Parameters param_, ChemicalSystem* cs_)
:param(param_), cs(cs_), T(cs_->get_final_time()), Nu(cs_->get_nu()), a(cs_->get_num_reactions())
{
#ifdef _OPENMP
    unsigned seed = chrono::system_clock::now().time_since_epoch().count()
                    +omp_get_thread_num();
#else
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
#endif
    generator.seed(seed);
}

Solver::~Solver()
{
    
}

VectorXd Solver::get_sol()
{
    return X;
}

void Solver::save_solution(const string& filename)
{
    unsigned int N = X.size();
    
    ofstream ofile(filename+string(".bin"), ios::out | ios::binary | ios::app);

    ofile.write((char*)&t, sizeof(double));
    ofile.write((char*)&X(0), N*sizeof(double));
    
    ofile.close();
}

void Solver::post_process()
{
}

Real Solver::get_tau_mean()
{
    return tau_mean;
}

Real Solver::get_s_mean()
{
    return s_mean;
}

Real Solver::get_rho_or_Nit_mean()
{
    return rho_or_Nit_mean;
}

Real Solver::get_damping_mean()
{
    return damping_mean;
}