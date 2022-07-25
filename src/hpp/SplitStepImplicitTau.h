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

#ifndef SPLITSTEPIMPLICITTAU_H
#define SPLITSTEPIMPLICITTAU_H

#include "TauLeapMethod.h"

// This scheme comes from the reference below, but here it is not multi level monte carlo.
// It is similar to Implicit-tau but it first does an implicit Euler step (only drift, no noise) 
// obtaining a solution Y and then to Y it adds the zero-mean noise evaluated at Y.
// Hammouda, C. B., Moraes, A., and Tempone, R. (2017). Multilevel hybrid split-step implicit tau-leap. Numer. Algorithms , 74(2):527-560
// We added as well a post-processing procedure, the same as implicit tau.
class SplitStepImplicitTau: public TauLeapMethod
{
public:
    SplitStepImplicitTau(Parameters param_, ChemicalSystem* cs_);
    ~SplitStepImplicitTau();
         
    
protected:
    void step();
    
    MatrixXd J, I;
    VectorXd atmp[2];
    
    // Options for interlaced tau
    unsigned int n_ET_steps; // number of explicit tau steps in the interlaced implicit tau method
    Real size_EE_steps;      // relative size of the steps, compared to maximal step size 2/rho
};

#endif /* SPLITSTEPIMPLICITTAU_H */

