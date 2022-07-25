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

#ifndef SOLVER_H
#define SOLVER_H

#include "headers.h"
#include "Parameters.h"
#include "ChemicalSystem.h"

class Solver
{
public:
    Solver(Parameters param_, ChemicalSystem* cs_);
    virtual ~Solver();
    
    virtual bool run() =0;
    virtual bool run_once_save_path() =0;
    
    VectorXd get_sol();
    Real get_tau_mean();
    Real get_s_mean();
    Real get_damping_mean();
    Real get_rho_or_Nit_mean();
    
protected:
    void save_solution(const string& filename);
    virtual void print_step_info() =0;
    virtual void post_process();
    
    Parameters param;
    ChemicalSystem* cs;
    Real t;         // current time 
    VectorXd X;     // state-vector
    
    unsigned int nsteps;
    Real tau;       // waiting time or leap time
    Real tau_mean;  // the average tau used in the simulation
    Real s_mean;    // the average number of stages
    Real damping_mean; //average damping parameter
    Real rho_or_Nit_mean; //the average number of iterations to compute rho or to solve the nonlinear system
    
    // auxiliary variables
    const Real T;   // final time
    const MatrixXd& Nu;    // soichiometrix matrix 
    VectorXd a;     // propensity functions, evaluated
    Real a0;        // sum of propensity functions
    
    default_random_engine generator;
};

#endif /* SOLVER_H */

