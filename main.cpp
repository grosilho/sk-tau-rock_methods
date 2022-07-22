#include "headers.h"
#include "MonteCarlo.h"
#include "ProblemsList.h"
#include "SSA.h"

int main(int argc, char** argv) 
{
    
    string filename = "sol";
    string solver_name = "tl";
    int prob_num = 2;
    unsigned int MCiter = 1e3;
    Real tau_max = 0.005;
    unsigned int nout = 1e4;
    bool post_process = false;
    unsigned int n_bins = 0;
    string refsol = "";
    Real Newton_tol = 1e-2;
    
    Parameters param(argc, argv, prob_num, solver_name, filename, MCiter, 
                     tau_max, nout, post_process, n_bins, refsol, Newton_tol);
    
    if(param.get_solver()==string("none")) // we dont simulate anything but just compare two already existing sets of data
    {
        MonteCarlo MC(param);
        MC.import_process_data();
    }
    else if(param.get_MCiter()>1) //do a MC simulation
    {
        MonteCarlo MC(param);
        MC.run();
    }
    else //compute and print details of one sample path
    {
        ChemicalSystem* cs;
        Solver* sol;
        param.init_Solver_ChemicalSystem(sol,cs);
        
        if(!sol->run_once_save_path())
            cout<<"ERROR: simulation aborted."<<endl;
        
        delete sol;
        delete cs;
    }
}

