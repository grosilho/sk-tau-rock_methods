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

#ifndef CHEBYSHEVMETHODS_H
#define CHEBYSHEVMETHODS_H

#include "headers.h"

namespace ChebyshevMethods
{    
    void CoefficientsRKC1(vector<Real>& mu, vector<Real>& nu, vector<Real>& kappa, 
                          unsigned int s, Real eps=0.05);
    Real StabBoundaryRKC1(unsigned int s, Real eps=0.05);
    
    Real T(Real x, unsigned int s);
    Real dT(Real x, unsigned int s);
    Real ddT(Real x, unsigned int s);
    Real U(Real x, unsigned int s);
    Real dU(Real x, unsigned int s);
};

#endif /* CHEBYSHEVMETHODS_H */

