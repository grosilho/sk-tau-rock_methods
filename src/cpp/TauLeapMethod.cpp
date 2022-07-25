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

#include "TauLeapMethod.h"

TauLeapMethod::TauLeapMethod(Parameters param_, ChemicalSystem* cs_)
:Solver(param_, cs_), K(VectorXd::Zero(cs_->get_num_reactions())), tau_max(param_.get_tau_max())
{
    eigmax = 0.;
    rho_or_Newton_iter=0;
    n_tmp_vec = 0;
}

TauLeapMethod::~TauLeapMethod()
{
}

void TauLeapMethod::init_eigvec_tmp_vec(bool init_eigvec)
{
    if(init_eigvec)
        eigvec = VectorXd::Zero(param.get_num_species());

    tmp_vec.resize(n_tmp_vec,VectorXd::Zero(param.get_num_species()));
}

void TauLeapMethod::sample_K(const VectorXd& means)
{    
    for(unsigned int i=0;i<K.size();i++)
    {
        poisson_dist.param(poisson_distribution<int>::param_type(means(i)));
        K(i)= poisson_dist(generator);
    }
}

void TauLeapMethod::compute_Jacobian(MatrixXd& J, VectorXd& Xeval, 
                    VectorXd& Xtmp, VectorXd atmp[2])
{
    Real eps = 1e-8;
    Xtmp = Xeval;
    Real dx;
            
    cs->eval_propensity_fun(Xeval,atmp[1]);
    
    for(unsigned int i=0;i<Xeval.size();i++)
    {
        dx = eps*(1.+abs(Xeval(i)));
        
        Xtmp(i) = Xeval(i)+dx;
        cs->eval_propensity_fun(Xtmp,atmp[0]);
       
        J.col(i) = Nu*(atmp[0]-atmp[1])/dx;
        Xtmp(i) = Xeval(i);
    }
}

bool TauLeapMethod::run()
{
    nsteps=0;
    last = false;
    
    t = 0.;
    X = cs->get_X0();
    
    VectorXd::Index minI;
    s_mean=0.;
    damping_mean=0.;
    rho_or_Nit_mean=0.;
    Real normX, normXold;
    Real increase_tol = 1e4;
    normX = X.lpNorm<Eigen::Infinity>();
    
    while(!last)
    {
        normXold = normX;
        tau = tau_max; //try to use the maximal tau allowed by the user
        
        if(t+1.01*tau>T)
        {
            tau = T-t;
            last = true;
        }
       
        if(!this->ensure_stability()) //checks if stability conditions are met, adapts tau or stages if needed
            return false;
        
        this->step();
        
        s_mean += s;
        rho_or_Nit_mean += rho_or_Newton_iter;
        
        if(isNan(X))
        {
            cout<<"ERROR: solution is nan."<<endl;
            return false;
        }
        
        normX = X.lpNorm<Eigen::Infinity>();
        if(normXold*increase_tol<normX)
        {
            cout<<"ERROR: solution norm has increased too much."<<endl;
            return false;
        }
        
//        while(X.minCoeff(&minI)<0)
//            X(minI)=0;
        
        nsteps++;
    }
    
    tau_mean = T/nsteps;
    s_mean /= nsteps;
    rho_or_Nit_mean /= nsteps;
    damping_mean /= nsteps;
    
    if(param.get_post_proc())
        this->post_process();
    
    return true;

}

bool TauLeapMethod::run_once_save_path()
{
    remove(string(param.get_filename()+string(".bin")).c_str());
    
    nsteps=0;
    last = false;
    
    t = 0.;
    X = cs->get_X0();
    
    save_solution(param.get_filename());
    
    VectorXd::Index minI;
    tau_mean = 0.;
    s_mean =0.;
    rho_or_Nit_mean=0.;
    unsigned int out_count=1;
    const Real out_dt = T/param.get_nout();
    
    while(t<T)
    {
        tau = tau_max;
        
        if(t+1.01*tau>T)
        {
            tau = T-t;
            last = true;
        }
       
        if(!this->ensure_stability())
            return false;
                        
        this->step();
        
        s_mean += s;
        rho_or_Nit_mean += rho_or_Newton_iter;
                
        nsteps++;
        
        if(isNan(X))
        {
            cout<<"ERROR: solution is nan."<<endl;
            return false;
        }
        
//        while(X.minCoeff(&minI)<0)
//            X(minI)=0;
//            X(minI)=abs(X(minI));
        
        print_step_info();
                
        if(t>=out_count*out_dt)
        {
            save_solution(param.get_filename());
            out_count++;
        }
    }
    
    tau_mean = T/nsteps;
    s_mean /= nsteps;
    rho_or_Nit_mean /= nsteps;
    
    return true;
}

VectorXd TauLeapMethod::f(VectorXd& Y)
{
    cs->eval_propensity_fun(Y,a);
    return Nu*a;
}

VectorXd TauLeapMethod::Qdt(VectorXd& Y, Real dt)
{
    cs->eval_propensity_fun(Y.array().abs(),a);
    a = a.array().abs();
    sample_K(a*dt);
    return Nu*(K-a*dt);
}

bool TauLeapMethod::ensure_stability()
{
    return true;
}

void TauLeapMethod::print_step_info()
{
    unsigned int w=8;
    cout.setf(ios::scientific, ios::floatfield);
    cout.precision(4); 
            
    cout<<"Time step: "<<left<<setw(3)<<nsteps<<", ";
    cout<<"t = "<<t<<", ";
    cout<<"\u03C4 = "<<tau<<", ";
    
    cout<<"\u03C1 = "<<eigmax<<", ";    
    cout<<"s = "<<left<<setw(3)<<s<<", ";
    cout<<"\u03C1/N iter = "<<left<<setw(3)<<rho_or_Newton_iter<<". ";
    cout<<"#Reactions: "<<K.sum()<<", ";
    cout<<"|X| = "<<X.lpNorm<Eigen::Infinity>()<<", ";
    cout<<"Neg values: "<<((X.minCoeff()<0) ? "Yes":"No")<<".";
    cout<<endl<<flush;
}

bool TauLeapMethod::rho()
{
    Real eigmaxo,sqrtu,znor,ynor,dzyn,dfzfn;
    
    const int maxiter=1e4;
    const Real safe=1.05;
    const Real tol = param.get_Newton_tol();
    bool ok = true;
    
    sqrtu = sqrt(numeric_limits<Real>::epsilon());

// ------ The initial vectors for the power method are yn --------
//       and yn+c*f(v_n), where vn=f(yn) a perturbation of yn 
//       (if n_steps=1) or a perturbation of the last computed
//       eigenvector (if n_steps!=0). 
    
    VectorXd& fn = tmp_vec[0];
    VectorXd& dz = tmp_vec[1];  
    VectorXd& z = eigvec; 
       
    if(nsteps==0)//take a random vector as first vector of power method
    {
        for(unsigned int i=0;i<z.size();i++)
            z[i] = ((rand()%2)-0.5)*((double)rand())/((double)RAND_MAX);
    }
    // else z=eigvec is the eigenvector computed by the last power method call    

    znor = z.norm();
    ynor = X.norm();    

    // ------ Perturbation.--------
    // Building the vector z so that the difference z-yn is small
    if(ynor!=0.0 && znor!=0.0)
    {
        dzyn = ynor*sqrtu;
        z = (dzyn/znor)*z + X;
    }
    else if(ynor!=0.0)
    {
        dzyn = ynor*sqrtu;
        z = (1.+sqrtu)*X;
    }
    else if(znor!=0.0)
    {
        dzyn = sqrtu;
        z *= (dzyn/znor);
    }
    else
    {
        dzyn = sqrtu*sqrt(z.size());
        z += sqrtu*VectorXd::Ones(z.size());
    }
    //here dzyn=||z-yn|| and z=yn+(small perturbation)
    //dzyn=||z-yn|| will be always true, even with new z in the loop
    
    //Start the power method for nonlinear operator
    cs->eval_propensity_fun(X,a);
    fn = Nu*a;
//    eigmax = 0.0;
    for(rho_or_Newton_iter=1;rho_or_Newton_iter<=maxiter;rho_or_Newton_iter++)
    {
        cs->eval_propensity_fun(z,a);
        dz = Nu*a;

        dz -= fn; //dz is the new perturbation, not normalized yet
        dfzfn= dz.norm();

        eigmaxo = eigmax;
        eigmax = dfzfn/dzyn; //approximation of the Rayleigh quotient (not with dot product but just norms)
        eigmax = safe*eigmax;    

//        cout<<setprecision(16)<<"eigmax: "<<eigmax<<", rel incr: "<<abs(eigmax-eigmaxo)/eigmax<<", tol: "<<tol<<endl;
        
        if(abs(eigmax-eigmaxo)<= eigmax*tol)
        {
            //The last perturbation is stored. It will very likely be a
            // good starting point for the next rho call.
            z = (z-X)/dzyn;               
            break;
        }
        else if(rho_or_Newton_iter==maxiter)
        {
            cout<<"ERROR: Reached max iterations in the spectral radius computation."<<endl;
            ok = false;
            break;
        }
        
        if (dfzfn!=0.0)
            z = (dzyn/dfzfn)*dz + X;
        else
        {
            cout<<"ERROR: Convergence failure in the spectral radius computation."<<endl;
            ok = false;
            break;
        }
    }
     
    return ok;
}

bool TauLeapMethod::isNan(const VectorXd& vec)
{
    for(unsigned int i=0;i<vec.size();i++)
        if(isnan(vec(i)))
            return true;
    
    return false;
}