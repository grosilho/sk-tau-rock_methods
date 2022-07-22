#include "tau-ROCK.h"
#include "ChebyshevMethods.h"

tauROCK::tauROCK(Parameters param_, ChemicalSystem* cs_)
:StabilizedTauLeapMethod(param_, cs_)
{

}

tauROCK::~tauROCK()
{
}

void tauROCK::step()
{
    if(s_old!=s)
    {
        ChebyshevMethods::CoefficientsRKC1(mu,nu,kappa,s,damping);
        coeff_recomputation++;
    }
    
    VectorXd& Kjm2 = tmp_vec[0];
    VectorXd& Kjm1 = tmp_vec[1];
    VectorXd& Kj = tmp_vec[2];
    
    Kjm1 = X; //K_0
    Kj   = X + mu[0]*tau*f(Kjm1); // K_1
                                  
    for(unsigned int j=2; j<=s; j++)
    {
        Kjm2 = Kjm1;
        Kjm1 = Kj;
        
        Kj = nu[j-1]*Kjm1 + kappa[j-1]*Kjm2 + mu[j-1]*tau*f(Kjm1);
    }
        
    X = Kj + Qdt(Kjm1,tau);
    t += tau;
}

bool tauROCK::ensure_stability()
{
    unsigned int safe_add = 0;
    
    if(!rho())
        return false;
    
    s_old = s;
    
    // S-ROCK choice for s
    s=1;
    damping = 0.05;
    Real w0,w1;
    bool unstable=true;
    do
    {
        w0 = 1.+damping/s/s;
        w1 = ChebyshevMethods::T(w0,s)/ChebyshevMethods::dT(w0,s);
        
        if(2.*w0/w1>tau*eigmax)
            unstable=false;
        else
        {
            s++;
            damping = 3.3*(log(s)-log(2.))+0.35;
        }
        
    }while(unstable);
    
    
    // standard choice for s
    //    damping = 0.05;
//    beta = 2.-4.*damping/3.;
    // s = ceil(sqrt(tau*eigmax/beta))+ safe_add;

    // Manual choice
//    damping = 3500;
//    s =800;
    
    damping_mean += damping;
    
    return true;
}