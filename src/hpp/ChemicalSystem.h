#ifndef CHEMICALSYSTEM_H
#define CHEMICALSYSTEM_H

#include "headers.h"

class ChemicalSystem
{
public:
    ChemicalSystem(const int N_, const int M_);
    virtual ~ChemicalSystem();
    
    virtual void eval_propensity_fun(const VectorXd& X, VectorXd& a) =0;
    Real get_final_time();
    VectorXd get_X0();
    const MatrixXd& get_nu();
    int get_num_species();
    int get_num_reactions();
    string get_problem_name();
    
protected:
    const int N;             // number of species
    const int M;             // number of reactions
    Real T;         // final time
    VectorXd X0;          // initial values of state-vector
    MatrixXd nu;    // stoichiometric matrix
    VectorXd c;     // rate constants
    
    string problem_name;
};

#endif /* CHEMICALSYSTEM_H */

