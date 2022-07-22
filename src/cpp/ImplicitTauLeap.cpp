 #include "ImplicitTauLeap.h"

ImplicitTauLeap::ImplicitTauLeap(Parameters param_, ChemicalSystem* cs_)
:TauLeapMethod(param_, cs_)
{
    unsigned int N = cs->get_num_species();
    unsigned int M = cs->get_num_reactions();
    J.resize(N,N);
    I = MatrixXd::Identity(N,N);
    
    n_tmp_vec = 4;
    
    if(param.get_post_proc())
        init_eigvec_tmp_vec(true); //in this case we use the interlaced tau and thus we need to estimate the spectral radius at the end of the simulation
    else
        init_eigvec_tmp_vec(false);
    
    atmp[0].resize(M);
    atmp[1].resize(M);
    
    s=1;
    
    ET_postprocess = false;
    n_ET_steps = 10;
    size_ET_steps = 0.1;
}

ImplicitTauLeap::~ImplicitTauLeap()
{
}

void ImplicitTauLeap::step()
{
    Real eps = param.get_Newton_tol();
    bool ok;
    
    const unsigned int max_newton_it = 100;
    
    // preparing to use the interlaced implicit tau, i.e. doing a sequence of small steps after the large step
    if(param.get_post_proc() && last) 
    {
        for(unsigned int i=0;i<eigvec.size();i++)
            eigvec[i] = ((rand()%2)-0.5)*((double)rand())/((double)RAND_MAX);
        rho();
        
        if(n_ET_steps*size_ET_steps*2./eigmax>tau)
        {
            tau = 0;
            size_ET_steps = eigmax*tau/n_ET_steps/2.;
        }
        else
            tau = tau - n_ET_steps*size_ET_steps*2./eigmax;
    }
          
    if(tau>0.)
    {
        cs->eval_propensity_fun(X,a);
        sample_K(a*tau);

        tmp_vec[0] = X + Nu*(K-a*tau);
        tmp_vec[1] = X;

        rho_or_Newton_iter=0;
        do
        {
            compute_Jacobian(J,tmp_vec[1],tmp_vec[2],atmp);
            J = I - J*tau;
            cs->eval_propensity_fun(tmp_vec[1],a);
            tmp_vec[2] = tmp_vec[0]-tmp_vec[1]+Nu*a*tau;

            tmp_vec[3] = J.colPivHouseholderQr().solve(tmp_vec[2]);
            tmp_vec[1] += tmp_vec[3];

            ok=true;
            for(unsigned int j=0;j<X.size();j++)
                ok = ok && (abs(tmp_vec[3](j))<(1.+abs(tmp_vec[1](j)))*eps);

            rho_or_Newton_iter++;

        }while(!ok && rho_or_Newton_iter<max_newton_it);

        if(rho_or_Newton_iter==max_newton_it && !ok)
            cout<<"WARNING: Newton algorithm did not converge."<<endl;

        X = tmp_vec[1];
        t += tau;
    }
    
    if(param.get_post_proc() && last) // perform some small steps of explicit or implicit tau
    {
        tau = size_ET_steps*2./eigmax;
        
        // EXPLICIT TAU
        if(ET_postprocess)
            for(unsigned int i=1;i<=n_ET_steps;i++)
            {
                cs->eval_propensity_fun(X,a);
                sample_K(a*tau);

                X+= Nu*K;
                t+= tau;
            }
        else//IMPLICIT TAU
            for(unsigned int i=1;i<=n_ET_steps;i++)
            {
                cs->eval_propensity_fun(X,a);
                sample_K(a*tau);

                tmp_vec[0] = X + Nu*(K-a*tau);
                tmp_vec[1] = X;

                rho_or_Newton_iter=0;
                do
                {
                    compute_Jacobian(J,tmp_vec[1],tmp_vec[2],atmp);
                    J = I - J*tau;
                    cs->eval_propensity_fun(tmp_vec[1],a);
                    tmp_vec[2] = tmp_vec[0]-tmp_vec[1]+Nu*a*tau;

                    tmp_vec[3] = J.colPivHouseholderQr().solve(tmp_vec[2]);
                    tmp_vec[1] += tmp_vec[3];

                    ok=true;
                    for(unsigned int j=0;j<X.size();j++)
                        ok = ok && (abs(tmp_vec[3](j))<(1.+abs(tmp_vec[1](j)))*eps);

                    rho_or_Newton_iter++;

                }while(!ok && rho_or_Newton_iter<max_newton_it);

                if(rho_or_Newton_iter==max_newton_it && !ok)
                    cout<<"WARNING: Newton algorithm did not converge."<<endl;

                X = tmp_vec[1];
                t += tau;
            }
    }
}
