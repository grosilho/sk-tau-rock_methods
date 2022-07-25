#include "headers.h"
#include "MonteCarlo.h"
#include "ProblemsList.h"
#include "SSA.h"

int main(int argc, char** argv) 
{
    Parameters param;
    
    // default values
    param.filename = "sol"; //output file name
    param.solver_name = "str"; //name of integration method
    param.prob_num = 2; //problem/example to solve
    param.MCiter = 1e2; //number of Monte Carlo iterations
    param.tau_max = 0.01; //tau/dt. Can be reduced during integration due to stability restrictions, for instance in standard tau-leap method
    param.nout = 1e4;   //number of output points. Used only when one path is computed, not for Monte Carlo simulations
    param.post_process = true; //postprocess the solution or not
    param.n_bins =30; //number of bins used in the approximation of density function
    param.refsol = "";//name of a reference solutions, if exists
    param.Newton_tol = 1e-2; //tolerance for Newton algorithm in implicit methods
    param.s_add = 1; //in stabilized methods add s_add stages to the estimated value. Used for safety when stiffness increases within one time step
    param.s = 0; //when s>0 then the number of stages is fixed to that value. It is not choses according to current stiffness
    param.damping = 0.05; //the damping parameter for sk-tau-rock. In tau-rock and rev-tau-rock the damping is chosen according to a formula depending on s. If param.s>0 the damping parameter is still chosen accordin to the same formula depending on s. But is param.damping>0 then it is fixed to param.damping, no formula is used.
    
    param.read_command_line(argc, argv);
    param.init();
    
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

