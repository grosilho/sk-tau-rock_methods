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