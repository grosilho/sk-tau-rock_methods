#include "TrapezoidalTauLeap.h"

TrapezoidalTauLeap::TrapezoidalTauLeap(Parameters param_, ChemicalSystem* cs_)
:TauLeapMethod(param_, cs_)
{
    unsigned int N = cs->get_num_species();
    unsigned int M = cs->get_num_reactions();
    J.resize(N,N);
    I = MatrixXd::Identity(N,N);
    
    n_tmp_vec = 4;
    init_eigvec_tmp_vec(false);
            
    atmp[0].resize(M);
    atmp[1].resize(M);
    
    s=1;
}

TrapezoidalTauLeap::~TrapezoidalTauLeap()
{
}

void TrapezoidalTauLeap::step()
{

    Real eps = param.get_Newton_tol();
    bool ok;
    
    const unsigned int max_newton_it = 100;
       
    cs->eval_propensity_fun(X,a);
    sample_K(a*tau);

    tmp_vec[0] = X + Nu*(K-a*tau*0.5);
    tmp_vec[1] = X;

//        cout<<"tn = "<<t<<endl;
//        cout<<"Xn : "<<X<<endl;
    
    rho_or_Newton_iter=0;
    do
    {
        compute_Jacobian(J,tmp_vec[1],tmp_vec[2],atmp);
        J = I - J*tau*0.5;
        cs->eval_propensity_fun(tmp_vec[1],a);
        tmp_vec[2] = tmp_vec[0]-tmp_vec[1]+Nu*a*tau*0.5;

        tmp_vec[3] = J.colPivHouseholderQr().solve(tmp_vec[2]);

//            cout<<"dX = "<<Xtmp[3]<<endl;
//            cout<<"Xnpu_kmu = "<<Xtmp[1]<<endl;

        tmp_vec[1] += tmp_vec[3];

        ok=true;
        for(unsigned int j=0;j<X.size();j++)
            ok = ok && (abs(tmp_vec[3](j))<(1.+abs(tmp_vec[1](j)))*eps);
        
        rho_or_Newton_iter++;
//            cout<<"ok: "<<ok<<endl;

    }while(!ok && rho_or_Newton_iter<max_newton_it);

    X = tmp_vec[1];
    t += tau;

    if(rho_or_Newton_iter==max_newton_it && !ok)
        cout<<"WARNING: Newton algorithm did not converge."<<endl;

//        cout<<"tnpu = "<<t<<endl;
//        cout<<"Xnpu : "<<X<<endl<<endl;

    
}
