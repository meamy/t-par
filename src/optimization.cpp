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

// ------------------------ Majority logic decoding

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

pair<int, word> reduce_vals(int m, word & A) {
  int mp;
  word Ap;
  xor_func p(m, 0);
  xor_func tmp;
  word::iterator it;

  for (it = A.begin(); it != A.end(); it++) p |= *it;
  mp = p.count();
  for (it = A.begin(); it != A.end(); it++) {
    int j = 0;
    tmp = xor_func(mp, 0);
    for (int i = 0; i < m; i++) {
      if (p.test(i)) tmp[j++] = (*it)[i];
    }
    Ap.insert(tmp);
  }
  return make_pair(mp, Ap); 
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

// --------------------------------------- Recursive decoding

#include <list>

typedef list<pair<xor_func, float> > sparse_vec;

void print_vec(const sparse_vec & A) {
  sparse_vec::const_iterator Ai;
  
  for (Ai = A.begin(); Ai != A.end(); Ai++) {
    cout << "(" << (Ai->first).to_ulong() << ", " << Ai->second << "), ";
  }

  cout << "\n";
}

void increment(xor_func & x) {
  for (int loop = 0; loop < x.size(); loop++) {
    x.flip(loop);
    if (x.test(loop)) break;
  }
}

sparse_vec from_word(int m, word & w) {
  sparse_vec res;
  word::const_iterator wi;
  xor_func tmp = xor_func(m, 0);

  res.push_back(make_pair(tmp, 1));

  for (wi = w.begin(); wi != w.end(); wi++) {
    if (*wi == tmp && tmp != xor_func(m, 0)) {
      res.pop_back();
    } else {
      tmp = *wi;
      res.push_back(make_pair(tmp, -1));
    }
    increment(tmp);
    if (tmp != xor_func(m, 0)) res.push_back(make_pair(tmp, 1));
  }

  return res;
}

sparse_vec from_word_mask(int m, word & w) {
  sparse_vec res;
  word::const_iterator wi;
  xor_func tmp = xor_func(m, 0);

  res.push_back(make_pair(tmp, 0));

  for (wi = w.begin(); wi != w.end(); wi++) {
    if (*wi == tmp && tmp != xor_func(m, 0)) {
      res.pop_back();
    } else {
      tmp = *wi;
      res.push_back(make_pair(tmp, 1));
    }
    increment(tmp);
    if (tmp != xor_func(m, 0)) res.push_back(make_pair(tmp, 0));
  }

  return res;
}

word to_word(int m, sparse_vec & A) {
  word res;
  int current;
  sparse_vec::iterator Ai = A.begin();

  for (xor_func tmp = xor_func(m, 0); ; ) {
    if (tmp == Ai->first) {
      if (Ai->second == -1) current = 1;
      else if (Ai->second == 1) current = 0;
      else assert(false);
      Ai++;
    }
    if (current == 1) res.insert(tmp);
    increment(tmp);
    if (tmp == xor_func(m, 0)) break;
  }

  return res;
}

void reset_high_order(int m, sparse_vec & A) {
  sparse_vec::iterator Ai;

  for (Ai = A.begin(); Ai != A.end(); Ai++) (Ai->first).reset(m-1);
}

void set_high_order(int m, sparse_vec & A) {
  sparse_vec::iterator Ai;

  for (Ai = A.begin(); Ai != A.end(); Ai++) (Ai->first).set(m-1);
}

sparse_vec split_vec(int n, int m, sparse_vec & A) {
  sparse_vec res;
  sparse_vec::iterator Ai;
  xor_func test(n, 0);
  test.set(m-1);

  for (Ai = A.begin(); Ai != A.end(); Ai++) {
    if (Ai->first == test) {
      res.splice(res.begin(), A, Ai, A.end());
      reset_high_order(m, res);
      return res;
    } else if (Ai->first.test(m-1)) {
      res.splice(res.begin(), A, Ai, A.end());
      break;
    }
  }
  reset_high_order(m, res);

  assert(not(A.empty()));
  if (res.empty() || (res.begin())->first != xor_func(n, 0)) {
    Ai = A.end();
    Ai--;
    res.push_front(make_pair(xor_func(n, 0), Ai->second));
  }
  return res;
}

sparse_vec vec_vec_mul(const sparse_vec & A, const sparse_vec & B) {
  sparse_vec res;
  sparse_vec::const_iterator Ai, Bi;

  for (Ai = A.begin(), Bi = B.begin(); Ai != A.end() && Bi != B.end(); Ai++, Bi++) {
    if (Ai->first == Bi-> first) {
      res.push_back(make_pair(Ai->first, Ai->second * Bi->second));
    } else if (Ai->first < Bi->first) {
      Bi--;
      res.push_back(make_pair(Ai->first, Ai->second * Bi->second));
    } else {
      Ai--;
      res.push_back(make_pair(Bi->first, Ai->second * Bi->second));
    }
  }

  Bi--;
  for ( ; Ai != A.end(); Ai++) {
    res.push_back(make_pair(Ai->first, Ai->second * Bi->second));
  }
  Bi++;
  Ai--;
  for ( ; Bi != B.end(); Bi++) {
    res.push_back(make_pair(Bi->first, Ai->second * Bi->second));
  }
  Ai++;

  return res;
}

sparse_vec vec_vec_add(const sparse_vec & A, const sparse_vec & B) {
  sparse_vec res;
  sparse_vec::const_iterator Ai, Bi;

  for (Ai = A.begin(), Bi = B.begin(); Ai != A.end() && Bi != B.end(); Ai++, Bi++) {
    if (Ai->first == Bi-> first) {
      res.push_back(make_pair(Ai->first, Ai->second + Bi->second));
    } else if (Ai->first < Bi->first) {
      Bi--;
      res.push_back(make_pair(Ai->first, Ai->second + Bi->second));
    } else {
      Ai--;
      res.push_back(make_pair(Bi->first, Ai->second + Bi->second));
    }
  }

  Bi--;
  for ( ; Ai != A.end(); Ai++) {
    res.push_back(make_pair(Ai->first, Ai->second + Bi->second));
  }
  Bi++;
  Ai--;
  for ( ; Bi != B.end(); Bi++) {
    res.push_back(make_pair(Bi->first, Ai->second + Bi->second));
  }
  Ai++;

  return res;
}

void scaler_vec_mul(float c, sparse_vec & A) {
  sparse_vec::iterator Ai;

  for (Ai = A.begin(); Ai != A.end(); Ai++) {
    Ai->second *= c;
  }
}

void vec_reduce(sparse_vec & A) {
  sparse_vec::iterator Ai = A.begin();
  float current = Ai->second;

  Ai++;
  for(; Ai != A.end();) {
    if (Ai->second == current) Ai = A.erase(Ai);
    else current = (Ai++)->second;
  }
}

long length(int m, sparse_vec & A, sparse_vec::const_iterator Ai) {
  sparse_vec::const_iterator Bi = Ai;
  Bi++;

  if (Bi != A.end()) {
    return (Bi->first).to_ulong() - (Ai->first).to_ulong();
  } else {
    return (1 << m) - (Ai->first).to_ulong();
  }
}

sparse_vec recursive_decode(int n, int m, int r, sparse_vec A) {
  sparse_vec res;
  sparse_vec::const_iterator Ai;

  if (r == 0) {                                   // {m, 0} decoding
    float tot = 0;
    for (Ai = A.begin(); Ai != A.end(); Ai++) {
      tot += length(m, A, Ai)*Ai->second;
    }
    if (tot >= 0) res.push_back(make_pair(xor_func(n, 0), 1));
    else res.push_back(make_pair(xor_func(n, 0), -1));
  } else if (r == m) {                            // {m, m} decoding
    for (Ai = A.begin(); Ai != A.end(); Ai++) {
      if (Ai->second >= 0) res.push_back(make_pair(Ai->first, 1));
      else res.push_back(make_pair(Ai->first, -1));
    }
  } else {
    sparse_vec B = split_vec(n, m, A);
    sparse_vec tmp1 = vec_vec_mul(A, B);
    vec_reduce(tmp1);
    sparse_vec v = recursive_decode(n, m-1, r-1, tmp1);
    tmp1 = vec_vec_mul(B, v);
    sparse_vec tmp2 = vec_vec_add(A, tmp1);
    vec_reduce(tmp2);
    scaler_vec_mul(0.5, tmp2);
    sparse_vec u = recursive_decode(n, m-1, r, tmp2);
    sparse_vec uv = vec_vec_mul(u, v);
    set_high_order(m, uv);
    res.splice(res.end(), u);
    res.splice(res.end(), uv);
    vec_reduce(res);
  }

  return res;
}

sparse_vec recursive_decode_H(int n, int m, int r, sparse_vec A, sparse_vec mask) {
  sparse_vec res;
  sparse_vec::const_iterator Ai;

  if (r == 0) {                                   // {m, 0} decoding
    for (Ai = mask.begin(); Ai != mask.end(); Ai++) {
      if (Ai->second == 0) {
        res.push_back(make_pair(xor_func(n, 0), 1));
        return res;
      }
    }

    float tot = 0;
    for (Ai = A.begin(); Ai != A.end(); Ai++) {
      tot += length(m, A, Ai)*Ai->second;
    }
    if (tot >= 0) res.push_back(make_pair(xor_func(n, 0), 1));
    else res.push_back(make_pair(xor_func(n, 0), -1));
  } else if (r == m) {                            // {m, m} decoding
    sparse_vec tmp = vec_vec_mul(A, mask);

    for (Ai = tmp.begin(); Ai != tmp.end(); Ai++) {
      if (Ai->second >= 0) res.push_back(make_pair(Ai->first, 1));
      else res.push_back(make_pair(Ai->first, -1));
    }
  } else {
    sparse_vec mask2 = split_vec(n, m, mask);
    sparse_vec mask3 = vec_vec_mul(mask, mask2);
    
    sparse_vec B = split_vec(n, m, A);
    sparse_vec tmp1 = vec_vec_mul(A, B);
    vec_reduce(tmp1);
    sparse_vec v = recursive_decode_H(n, m-1, r-1, tmp1, mask3);
    tmp1 = vec_vec_mul(B, v);
    sparse_vec tmp2 = vec_vec_add(A, tmp1);
    vec_reduce(tmp2);
    scaler_vec_mul(0.5, tmp2);
    sparse_vec u = recursive_decode_H(n, m-1, r, tmp2, mask3);
    sparse_vec uv = vec_vec_mul(u, v);
    set_high_order(m, uv);
    res.splice(res.end(), u);
    res.splice(res.end(), uv);
    vec_reduce(res);
  }

  return res;
}

void decode_rec(int m, word & y/*, word & avoid*/) {
  sparse_vec B = from_word(m, y);
  sparse_vec C = recursive_decode(m, m, m-4, B);
  word z = to_word(m, C);
  word e = add_words(y, z);
  y = e;
}

void decode_rec_H(int m, word & y, word & mask) {
  sparse_vec A = from_word_mask(m, mask);
  print_vec(A);
  sparse_vec B = from_word(m, y);
  print_vec(B);
  sparse_vec C = recursive_decode_H(m, m, m-4, B, A);
  word z = to_word(m, C);
  word e = add_words(y, z);
  y = e;
}

void minT_rec(int m, const vector<exponent> & A) {
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

//  cout << "Original number of variables: " << m << "\n";
//  pair<int, word> yp = reduce_vals(m, y);
//  cout << "New number of variables: " << yp.first << "\n";
  decode_rec(m, y);
  //decode_rec(m, y);
  //  print_word(m, y);
  cout << "Optimized T-count: " << y.size() << "\n";

}

void add_all_combs(int n, int m, xor_func * wires, word & w) {
  for (int i = 1; i < (1 << n); i++) {
    int x = i;
    xor_func tmp(m+1, 0);
    for (int j = 0; j < n; j++) {
      if (x % 2 == 1) tmp ^= wires[j];
      x /= 2;
    }
    tmp.resize(m);
    w.insert(tmp);
  }
}
        
word find_all(character & c) {
  word ret;

  for (list<Hadamard>::iterator ti = c.hadamards.begin(); ti != c.hadamards.end(); ti++) {
    add_all_combs(c.n + c.m, c.n + c.h, ti->wires, ret);
  }
  add_all_combs(c.n + c.m, c.n + c.h, c.outputs, ret);

  return ret;
}

void test_rec_decode(int m, character & c) {
  word y;

  for (vector<exponent>::const_iterator it = c.phase_expts.begin(); it != c.phase_expts.end(); it++) {
    if (it->first % 2 == 1) {
      xor_func tmp = it->second;
      tmp.resize(m);
      y.insert(tmp);
    }
  }
  //  print_word(m, y);
  cout << "Original T-count: " << y.size() << "\n";
  c.print();
  cout << "Free values: " << "\n";
  word w = find_all(c);

//  cout << "Original number of variables: " << m << "\n";
//  pair<int, word> yp = reduce_vals(m, y);
//  cout << "New number of variables: " << yp.first << "\n";
//  decode_rec(m, y);
  decode_rec_H(m, y, w);
  //decode_rec(m, y);
  //  print_word(m, y);
  cout << "Optimized T-count: " << y.size() << "\n";
  c.phase_expts.clear();
  for (word::iterator it = y.begin(); it != y.end(); it++) {
    it->resize(m + 1);
    c.phase_expts.push_back(make_pair(1, *it));
  }
  c.print();

  cout << "Testing realizability...";

  xor_func wires[c.n + c.m + 1];

  bool flag = false;
  for (word::iterator it = y.begin(); it != y.end(); it++) {
    flag = false;
    for (list<Hadamard>::iterator ti = c.hadamards.begin(); !flag && (ti != c.hadamards.end()); ti++) {
      int rk = compute_rank(c.n + c.m, m, ti->wires);
      for (int i = 0; i < c.n + c.m; i++) {
        wires[i] = ti->wires[i];
        wires[i].resize(m);
      }
      wires[c.n+c.m] = *it;
      wires[c.n+c.m].resize(m);
      int nrk = compute_rank(c.n + c.m + 1, m, wires);
      if (rk == nrk) flag = true;
    }
    int rk = compute_rank(c.n + c.m, m, c.outputs);
    for (int i = 0; i < c.n + c.m; i++) {
      wires[i] = c.outputs[i];
      wires[i].resize(m);
    }
    wires[c.n+c.m] = *it;
    wires[c.n+c.m].resize(m);
    int nrk = compute_rank(c.n + c.m + 1, m, wires);
    if (rk == nrk) flag = true;
    if (!flag) break;
  }
  if (flag) cout << "Yes!" << "\n";
  else cout << "No..." << "\n";

}

void test_rec() {
  srand (time(NULL));
  cout << "Starting..." << "\n" << flush;
  int m = 5;
  xor_func mon(m, 1 << 4);
  word w = eval_monomial(m, mon);
  print_word(m, w);
  for (int i = 0; i < 3; i++) {
    int s = rand() % ((1 << m) - 1) + 1;
    xor_func tmp(m, s);
    if (w.find(tmp) == w.end()) w.insert(tmp);
    else w.erase(tmp);
  }
  print_word(m, w);
  decode_rec(m, w);
  print_word(m, w);
}
