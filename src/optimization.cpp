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

#include <algorithm>

// m is number of variables
// n = 2^m - 1, length of words of RM(m, r)* code

typedef set<xor_func> word;

void print_word(int m, word & A) {
  for (int i = 1; i < (1 << m); i++) {
    if (A.find(xor_func(m, i)) == A.end()) cout << "0";
    else cout << "1";
  }
  cout << "\n";
}

word quotient(word & A, xor_func q) {
  word res;
  xor_func p;

  for (word::iterator it = A.begin(); it != A.end(); it++) {
    p = *it - q;
    if (res.find(p) == res.end()) res.insert(p);
    else res.erase(p);
  }

  return res;
}

word eval_monomial(int m, xor_func q) {
  word res;
  xor_func x;
  int i;

  for (i = 1; i < (1 << m); i++) {
    x = xor_func(m, i);
    if ((x & q) == q) res.insert(x);
  }

  return res;
}

word add_words(word & A, word & B) {
  word res;
  set_symmetric_difference(A.begin(), A.end(), B.begin(), B.end(), inserter(res, res.begin()));
  return res;
}

void decode(int m, word & y) {
  string bitmask;
  xor_func q;
  word z;
  word code;

  // For each degree <= m-4...
  for (int k = m - 4; k >= 0; k--) {
    code.clear();
    bitmask = string(k, '1');
    bitmask.resize(m, '0');
    // For each monomial (i.e. each string of length m with k bits)
    do {
      // Create a bitset from bitmask
      q = xor_func(bitmask);
      // Now quotient the word y by those bits
      z = quotient(y, q);
      // z.size() gives the number of 1 votes
      // If the 1 votes form a majority (> 2^(m - k)/2) then assign q 1
      if (z.size() > (1 << (m - k - 1))) {
	// That is, compute the word of evaluations of q as a monomial
	// Then add it to the result
	word tmp = eval_monomial(m, q);
	code = add_words(code, tmp);
      }
    } while (prev_permutation(bitmask.begin(), bitmask.end()));
    // Add the degree k codeword to y
    y = add_words(y, code);
  }
}

void minT(int m, const vector<exponent> & A) {
  word y;

  for (vector<exponent>::const_iterator it = A.begin(); it != A.end(); it++) {
    if (it->first % 2 == 1) {
      xor_func tmp = it->second;
      tmp.resize(m);
      y.insert(tmp);
    }
  }
  //  print_word(m, y);
  cout << "Original T-count: " << y.size() << "\n";

  decode(m, y);
  //  print_word(m, y);
  cout << "Optimized T-count: " << y.size() << "\n";

}

#include <time.h>

void test() {
  srand (time(NULL));
  cout << "Starting..." << "\n" << flush;
  int m = 10;
  xor_func mon(m, 1 << 4);
  word w = eval_monomial(m, mon);
  print_word(m, w);
  for (int i = 0; i < 4; i++) {
    int s = rand() % ((1 << m) - 1) + 1;
    xor_func tmp(m, s);
    if (w.find(tmp) == w.end()) w.insert(tmp);
    else w.erase(tmp);
  }
  print_word(m, w);
  decode(m, w);
  print_word(m, w);
}
