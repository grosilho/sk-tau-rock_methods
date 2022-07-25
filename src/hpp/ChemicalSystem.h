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

#ifndef CHEMICALSYSTEM_H
#define CHEMICALSYSTEM_H

#include "headers.h"

class ChemicalSystem
{
public:
    ChemicalSystem(const int N_, const int M_);
    virtual ~ChemicalSystem();
    
    virtual void eval_propensity_fun(const VectorXd& X, VectorXd& a) =0;
    Real get_final_time();
    VectorXd get_X0();
    const MatrixXd& get_nu();
    int get_num_species();
    int get_num_reactions();
    string get_problem_name();
    
protected:
    const int N;             // number of species
    const int M;             // number of reactions
    Real T;         // final time
    VectorXd X0;          // initial values of state-vector
    MatrixXd nu;    // stoichiometric matrix
    VectorXd c;     // rate constants
    
    string problem_name;
};

#endif /* CHEMICALSYSTEM_H */

