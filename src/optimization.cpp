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

int max_weight(const vector<exponent> & A) {
  int cur_max = 3;
  int cur_ind = -1;
  int i;
  
  for (i = 0; i < A.size(); i++) {
    if ((A[i].first % 2 == 1) && (A[i].second.count() > cur_max)) cur_ind = i;
  }

  return cur_ind;
}

void add_id(vector<exponent> & A, xor_func f, int cur, int n, int val) {
  for (; cur < n; cur++) {
    if (f.test(cur)) {
      add_id(A, f, cur + 1, n, val);
      f.reset(cur);
      insert_phase(val, f, A);
      if (f.count() > 1) add_id(A, f, cur + 1, n, val);
    }
  }
}

void gaussian_opt(character & circ) {
  int cur_size, cur_val;
  int index = max_weight(circ.phase_expts);
  list<Hadamard>::iterator it;
  int i;

  while (index != -1) {
    cur_size = circ.phase_expts.size();
    cur_val  = 8 - circ.phase_expts[index].first;
    circ.phase_expts[index].first = 0;
    add_id(circ.phase_expts, circ.phase_expts[index].second, 0, circ.n + circ.h, cur_val);
    for (it = circ.hadamards.begin(); it != circ.hadamards.end(); it++) {
      if (it->in.find(index) != it->in.end()) {
	for (i = cur_size; i < circ.phase_expts.size(); i++) it->in.insert(i);
      }
    }
    index = max_weight(circ.phase_expts);
  }

  int total = 0;
  for (i = 0; i < circ.phase_expts.size(); i++) {
    if (circ.phase_expts[i].first % 2 == 1) total++;
  }
  cout << "# T-count: " << total << "\n" << flush;
    
}

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

