/*
    Copyright (C) 2022 Giacomo Rosilho de Souza

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

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
