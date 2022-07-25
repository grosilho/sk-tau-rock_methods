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

#ifndef IMPLICITTAULEAP_H
#define IMPLICITTAULEAP_H

#include "TauLeapMethod.h"

// If the post processing option is enabled then this scheme is the interlaced implicit tau, where 
// a few small steps resolving the fast variables are performed at the end of the simulation.
// Rathinam, M., Petzold, L. R., Cao, Y., and Gillespie, D. T. (2003). Stiffness in stochastic chemically reacting systems: The implicit tau-leaping method. Journal of Chemical Physics, 119(24):12784-12794.
class ImplicitTauLeap: public TauLeapMethod
{
public:
    ImplicitTauLeap(Parameters param_, ChemicalSystem* cs_);
    ~ImplicitTauLeap();
         
    
protected:
    void step();
    
    MatrixXd J, I;
    VectorXd atmp[2];
    
    // Options for interlaced tau
    bool ET_postprocess;     // if true then postprocessing is done with explicit tau, else with implicit tau
    unsigned int n_ET_steps; // number of postprocessing steps
    Real size_ET_steps;      // relative size of the steps, compared to 2/rho
};

#endif /* IMPLICITTAULEAP_H */

