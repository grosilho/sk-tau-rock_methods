#include "TauLeap.h"

TauLeap::TauLeap(Parameters param_, ChemicalSystem* cs_)
:TauLeapMethod(param_, cs_)
{
    s=1;
    
    n_tmp_vec = 2;
    init_eigvec_tmp_vec(true);
}

TauLeap::~TauLeap()
{
}


void TauLeap::step()
{
    cs->eval_propensity_fun(X,a);
    sample_K(a*tau);

    X+= Nu*K;
    t+= tau;
}

bool TauLeap::ensure_stability()
{
    if(!rho())
        return false;
    
    if(tau>0.9*2./eigmax)
    {
        tau = 0.9*2./eigmax;
        last = false;
    }
    
    return true;
}