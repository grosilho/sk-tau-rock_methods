#include "R-tau-ROCK.h"
#include "ChebyshevMethods.h"

RtauROCK::RtauROCK(Parameters param_, ChemicalSystem* cs_)
:StabilizedTauLeapMethod(param_, cs_)
{
    alpha = 1;
}

RtauROCK::~RtauROCK()
{
}

void RtauROCK::step()
{
    if(s_old!=s)
    {
        ChebyshevMethods::CoefficientsRKC1(mu,nu,kappa,s,damping);
        coeff_recomputation++;
    }
    
    VectorXd& Kjm2 = tmp_vec[0];
    VectorXd& Kjm1 = tmp_vec[1];
    VectorXd& Kj = tmp_vec[2];
    
    Kjm1 = X + Qdt(X,tau);      // K_0
    Kj = Kjm1 + mu[0]*tau*f(Kjm1); // K_1
    
    for(unsigned int j=2; j<=s; j++)
    {
        Kjm2 = Kjm1;
        Kjm1 = Kj;
        
        Kj = nu[j-1]*Kjm1 + kappa[j-1]*Kjm2 + mu[j-1]*tau*f(Kjm1);
    }

    X = Kj;
    t += tau;
}

bool RtauROCK::ensure_stability()
{
    unsigned int safe_add = 0;//4
    Real ell, w0, w1;
    
    if(!rho())
        return false;
    
    s_old = s;
    
    // tentative damping parameters
    damping = 0.05;
    beta = 2.-4.*damping/3.;
    s = ceil(sqrt(tau*eigmax/beta));
    
    s--;
    do
    {
        s++;
        beta = pow(8.*tau*eigmax/alpha,1./2./s);
        damping = (1.-beta)*(1.-beta)*s*s/2./beta;
        w0 = 1+damping/s/s;
        w1 = ChebyshevMethods::T(w0,s)/ChebyshevMethods::dT(w0,s);
        ell = (1.+w0)/w1;
    }while(tau*eigmax>ell);
    
    s+=safe_add;
    beta = pow(8.*tau*eigmax/alpha,1./2./s);
    damping = (1.-beta)*(1.-beta)*s*s/2./beta;
    
    
//    damping = 3500;
//    s = 800;
    
    damping_mean += damping;
    
    return true;
}