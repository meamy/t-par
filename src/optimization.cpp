/*--------------------------------------------------------------------
  Tpar - T-gate optimization for quantum circuits
  Copyright (C) 2013  Matthew Amy and The University of Waterloo,
  Institute for Quantum Computing, Quantum Circuits Group

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

Author: Matthew Amy
---------------------------------------------------------------------*/

#include "optimization.h"

void remove_x(int n, vector<exponent> & phase) {
  for (vector<exponent>::iterator it = phase.begin(); it != phase.end(); it++) {
    if (it->second.test(n)) {
      it->second.reset(n);
      it->first = (it->first + 7) % 8;
      phase.push_back(make_pair(it->first, xor_func(n + 1)));
    }
  }
}

