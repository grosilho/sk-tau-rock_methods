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

#ifndef R_TAU_ROCK_H
#define R_TAU_ROCK_H

#include "StabilizedTauLeapMethod.h"

class RtauROCK: public StabilizedTauLeapMethod
{
public:
    RtauROCK(Parameters param_, ChemicalSystem* cs_);
    ~RtauROCK();
     

protected:
    void step();
    bool ensure_stability();
  
    Real alpha;
};

#endif /* R_TAU_ROCK_H */

