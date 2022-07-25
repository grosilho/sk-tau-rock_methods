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

#include "StabilizedTauLeapMethod.h"

StabilizedTauLeapMethod::StabilizedTauLeapMethod(Parameters param_, ChemicalSystem* cs_)
:TauLeapMethod(param_, cs_)
{
    s = 0;
    s_old = 0;
    coeff_recomputation=0;
    
    n_tmp_vec = 3;
    init_eigvec_tmp_vec(true);
}

StabilizedTauLeapMethod::~StabilizedTauLeapMethod()
{
}

void StabilizedTauLeapMethod::print_step_info()
{
    unsigned int w=8;
    cout.setf(ios::scientific, ios::floatfield);
    cout.precision(4); 
            
    cout<<"Time step: "<<left<<setw(3)<<nsteps<<", ";
    cout<<"t = "<<t<<", ";
    cout<<"\u03C4 = "<<tau<<", ";
    
    cout<<"\u03C1 = "<<eigmax<<", ";    
    cout<<"s = "<<left<<setw(3)<<s<<", ";
    cout<<"\u03B5 = "<<left<<setw(3)<<damping<<", ";
    cout<<"\u03C1 iter = "<<left<<setw(3)<<rho_or_Newton_iter<<". ";
    cout<<"#Reactions: "<<K.sum()<<", ";
    cout<<"|X| = "<<X.lpNorm<Eigen::Infinity>()<<", ";
    cout<<"Neg values: "<<((X.minCoeff()<0) ? "Yes":"No")<<".";
    cout<<endl<<flush;
}