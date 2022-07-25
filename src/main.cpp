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

#include "headers.h"
#include "MonteCarlo.h"
#include "ProblemsList.h"
#include "SSA.h"

int main(int argc, char** argv) 
{
    cout<<"SK-tau-ROCK methods.  Copyright (C) 2022 Giacomo Rosilho de Souza\n"
           "This program comes with ABSOLUTELY NO WARRANTY; for details type run it with the --help option.\n"
           "This is free software, and you are welcome to redistribute it under certain conditions.\n"
           "You should have received a copy of the GNU General Public License along with this program.\n"
           "If not, see <https://www.gnu.org/licenses/>.\n"<<endl;
    
    Parameters param;
    
    // default values
    param.filename = "sol"; //output file name
    param.solver_name = "str"; //name of integration method
    param.prob_num = 4; //problem/example to solve
    param.MCiter = 1e4; //number of Monte Carlo iterations
    param.tau_max = 0.05; //tau/dt. Can be reduced during integration due to stability restrictions, for instance in standard tau-leap method
    param.nout = 1e4;   //number of output points. Used only when one path is computed, not for Monte Carlo simulations
    param.post_process = true; //postprocess the solution or not
    param.n_bins =30; //number of bins used in the approximation of density function
    param.refsol = "";//name of a reference solutions, if exists
    param.Newton_tol = 1e-2; //tolerance for Newton algorithm in implicit methods
    param.s_add = 1; //in stabilized methods add s_add stages to the estimated value. Used for safety when stiffness increases within one time step
    param.s = 0; //when s>0 then the number of stages is fixed to that value. It is not choses according to current stiffness
    param.damping = 0.05; //the damping parameter for sk-tau-rock. In tau-rock and rev-tau-rock the damping is chosen according to a formula depending on s. If param.s>0 the damping parameter is still chosen accordin to the same formula depending on s. But is param.damping>0 then it is fixed to param.damping, no formula is used.
    
    if(param.read_command_line(argc, argv))
        return 0;
    
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
    
    return 0;
}

