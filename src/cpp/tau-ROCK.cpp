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
    if(!rho())
        return false;
    
    s_old = s;
    
    if(param.s>0) //fix s to a user defined value
    {
        s = param.s;
        if(param.damping>=0.)//fix damping to user defined value
            damping=param.damping;
        else//use formula for damping depending on s
            damping = 3.3*(log(s)-log(2.))+0.35;
        
    }
    else// automatically chose s and damping
    {
        // S-ROCK choice for s
        s=1;
        damping = 0.05;
        while(ChebyshevMethods::StabBoundaryRKC1(s,damping)<tau*eigmax)
        {        
            s++;
            damping = 3.3*(log(s)-log(2.))+0.35;
        }
        s += param.s_add;
    }
    
    damping_mean += damping;
    
    return true;
}