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

#include "headers.h"
#include "ChemicalSystem.h"

ChemicalSystem::ChemicalSystem(const int N_, const int M_)
:N(N_), M(M_), X0(N_), nu(N_,M_), c(M_)
{
    
}

ChemicalSystem::~ChemicalSystem()
{
}

Real ChemicalSystem::get_final_time()
{
    return T;
}

VectorXd ChemicalSystem::get_X0()
{
    return X0;
}

const MatrixXd& ChemicalSystem::get_nu()
{
    return nu;
}

int ChemicalSystem::get_num_species()
{
    return N;
}

int ChemicalSystem::get_num_reactions()
{
    return M;
}

string ChemicalSystem::get_problem_name()
{
    return problem_name;
}