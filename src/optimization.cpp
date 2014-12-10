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

void add_vector(int n, vector<exponent> & A, const vector<exponent> & B) {
  vector<exponent>::const_iterator it;
  for (it = B.begin(); it != B.end(); it++) insert_phase(it->first, it->second, A);
}

//void simple_opt(int n, vector<exponent> & A) {

void remove_x(character & circ) {
  int n = circ.n + circ.h;
  int i, ind;
  list<Hadamard>::iterator it;

  for (i = 0; i < circ.phase_expts.size(); i++) {
    if (circ.phase_expts[i].second.test(n)) {
      xor_func tmp = circ.phase_expts[i].second;
      tmp.reset(n);
      insert_phase(circ.phase_expts[i].first, xor_func(n + 1, 0), circ.phase_expts);
      ind = insert_phase((circ.phase_expts[i].first*7) % 8, tmp, circ.phase_expts);
      for (it = circ.hadamards.begin(); it != circ.hadamards.end(); it++) {
	if (it->in.find(i) != it->in.end()) it->in.insert(ind);
      }
      circ.phase_expts[i].first = 0;
    }
  }
}

