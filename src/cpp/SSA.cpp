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

#include "SSA.h"

SSA::SSA(Parameters param_, ChemicalSystem* cs_)
:Solver(param_, cs_), unif_dist(0.0,1.0)
{
    
}

SSA::~SSA()
{
}

bool SSA::run()
{
    nsteps=0;
    t = 0.;
    X = cs->get_X0();
    
    cs->eval_propensity_fun(X,a);
    a0 = a.sum();
    tau = log(1./unif_dist(generator))/a0;
    t = t+tau;
    
    while(t<T)
    {
        r = unif_dist(generator);
        psum_a=a(0);
        j=0;
        while(psum_a<r*a0)
            psum_a += a(++j);
                
        X += Nu.col(j);
        
        nsteps++;
               
        cs->eval_propensity_fun(X,a);
        a0 = a.sum();
        tau = log(1./unif_dist(generator))/a0;
        t = t+tau;
    }
    
    tau_mean = T/nsteps;
    s_mean = 0.;
    rho_or_Nit_mean = 0.;
    
    return true;
}

bool SSA::run_once_save_path()
{
    remove(string(param.get_filename()+string(".bin")).c_str());
    
    nsteps=0;
    t = 0.;
    X = cs->get_X0();
    save_solution(param.get_filename());
    
    unsigned int out_count=1;
    const Real out_dt = T/param.get_nout();
    
    cs->eval_propensity_fun(X,a);
    a0 = a.sum();
    tau = log(1./unif_dist(generator))/a0;
    t = t+tau;
       
    while(t<T)
    {
        r = unif_dist(generator);
        psum_a=a(0);
        j=0;
        while(psum_a<r*a0)
            psum_a += a(++j);
                
        X += Nu.col(j);
        
        nsteps++;
        
        print_step_info();
              
        if(t>=out_count*out_dt)
        {
            save_solution(param.get_filename());
            out_count++;
        }
        
        cs->eval_propensity_fun(X,a);
        a0 = a.sum();
        if(a0!=0.)
        {
            tau = log(1./unif_dist(generator))/a0;
            t = t+tau;
        }
        else
        {
            cout<<"All reactions have zero probability. End of simulation."<<endl;
            tau = T-t;
            t=T;
            break;
        }
    }
    
    
    tau_mean = T/nsteps;
    s_mean=0.;
    rho_or_Nit_mean=0.;
    
    return true;
}

void SSA::print_step_info()
{
    unsigned int w=8;
    cout.setf(ios::scientific, ios::floatfield);
    cout.precision(4); 
            
    cout<<"Time step: "<<left<<setw(3)<<nsteps<<", ";
    cout<<"t = "<<t<<", ";
    cout<<"\u03C4 = "<<tau<<", ";
    cout<<"|X| = "<<X.lpNorm<Eigen::Infinity>()<<".";
    cout<<endl<<flush;
}