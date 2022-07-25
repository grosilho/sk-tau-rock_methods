#include "SK-tau-ROCK.h"
#include "ChebyshevMethods.h"

SKtauROCK::SKtauROCK(Parameters param_, ChemicalSystem* cs_)
:StabilizedTauLeapMethod(param_, cs_)
{
    damping = param.damping;
    beta = 2.-4.*damping/3.; //used only for small dampings
}

SKtauROCK::~SKtauROCK()
{
}

void SKtauROCK::step()
{
    if(s_old!=s)
    {
        ChebyshevMethods::CoefficientsRKC1(mu,nu,kappa,s,damping);
        coeff_recomputation++;
    }
    
    VectorXd& Kjm2 = tmp_vec[0];
    VectorXd& Kjm1 = tmp_vec[1];
    VectorXd& Kj = tmp_vec[2];
    
    alpha = 0.5*sqrt(mu[0]);
    
    Kjm2 = Qdt(X,tau);
    Kjm1 = X + nu[0]*Kjm2; // the argument of f for s=1
    
    Kj = X + mu[0]*tau*f(Kjm1) + kappa[0]*Kjm2; //  K_1
    Kjm1 = X;                                   //  K_0
    
    for(unsigned int j=2; j<=s; j++)
    {
        Kjm2 = Kjm1;
        Kjm1 = Kj;
        
        Kj = nu[j-1]*Kjm1 + kappa[j-1]*Kjm2 + mu[j-1]*tau*f(Kjm1);
    }

    X = Kj;
    t += tau;
}

bool SKtauROCK::ensure_stability()
{    
    if(!rho())
        return false;
    
    s_old = s;
    
    if(param.s>0)// fix s to a user defined value
    {
        s=param.s;
        return true;
    }
    else// chose s according to formula  
    {
        if(damping<=0.1)
            s = ceil(sqrt(tau*eigmax/beta))+ param.s_add;
        else//for large dampings the above formula is unvalid
        {
            s = ceil(sqrt(tau*eigmax/2.));
            while(ChebyshevMethods::StabBoundaryRKC1(s,damping)<tau*eigmax)
                s++;
            s += param.s_add;
        }
    }
    damping_mean += damping;
    
    return true;
}

void SKtauROCK::post_process()
{    
    // Old postprocessor parameter
    //Real alpha = (1./2./s);
    //alpha *= sqrt(1.+(2./3.)*damping*(1.-1./s/s)); //damping correction 
            
    X += alpha*Qdt(X,tau); 
}