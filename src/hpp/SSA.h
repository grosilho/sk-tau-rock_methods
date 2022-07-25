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

#ifndef SSA_H
#define SSA_H

#include "Solver.h"

class SSA: public Solver
{
public:
    SSA(Parameters param_, ChemicalSystem* cs_);
    ~SSA();
        
    bool run();
    bool run_once_save_path();
    
protected:
    void print_step_info();
    
    int j;         // next reaction index

    Real psum_a; // temporary var for the partial sum of prop fun
    Real r;     //temporary variable holding a sampling
    uniform_real_distribution<Real> unif_dist;
};

#endif /* SSA_H */

