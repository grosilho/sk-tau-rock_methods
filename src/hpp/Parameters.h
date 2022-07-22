#ifndef INITIALIZER_H
#define INITIALIZER_H

#include "headers.h"

class Solver;
class ChemicalSystem;

class Parameters
{
public:
    Parameters(int argc, char** argv, int prob_num_, string solver_name_, string filename_, 
               const unsigned int MCiter_, Real tau_max_, int nout_,
               bool post_proc_, unsigned int n_bins_, string refsol_, Real Newton_tol_);
    ~Parameters();
    
    void init_ChemicalSystem(ChemicalSystem*& cs);
    void init_Solver_ChemicalSystem(Solver*& sol, ChemicalSystem*& cs);
    
    string get_filename();
    unsigned int get_num_species();
    unsigned int get_MCiter();
    unsigned int get_nout();
    bool get_post_proc();
    unsigned int get_n_bins();
    string get_refsol();
    Real get_tau_max();
    string get_solver();
    Real get_final_time();
    Real get_Newton_tol();
    
protected:
    void read_command_line(int argc, char** argv);
    
    unsigned int prob_num;
    unsigned int N_species;
    Real final_time;
    unsigned int MCiter;
    Real tau_max;
    unsigned int nout;
    string solver_name;
    string problem_name;
    string filename;
    bool post_proc;
    unsigned int n_bins;
    Real Newton_tol;
    string refsol;
};

#endif /* INITIALIZER_H */

