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

#include "util.h"
#include <map>
#include <cmath>

bool disp_log = false;
synth_type synth_method = PMH;

void print_wires(const vector<xor_func>& wires, int num, int dim) {
  int i, j;
  for (i = 0; i < num; i++) {
    for (j = 0; j < dim; j++) {
      if (wires[i].test(j)) cout << "1";
      else                  cout << "0";
    }
    cout << "\n";
  }
}

// Commands for making certain circuits
gatelist xor_com(int a, int b, const vector<string>& names) {
  list<string> tmp_list;
  gatelist ret;

  tmp_list.push_back(names[a]);
  tmp_list.push_back(names[b]);
  ret.push_back(make_pair("tof", tmp_list));

  return ret;
}

gatelist swap_com(int a, int b, const vector<string>& names) {
  list<string> tmp_list1, tmp_list2;
  gatelist ret;

  tmp_list1.push_back(names[a]);
  tmp_list1.push_back(names[b]);
  tmp_list2.push_back(names[b]);
  tmp_list2.push_back(names[a]);
  ret.push_back(make_pair("tof", tmp_list1));
  ret.push_back(make_pair("tof", tmp_list2));
  ret.push_back(make_pair("tof", tmp_list1));

  return ret;
}

gatelist x_com(int a, const vector<string>& names) {
  gatelist ret;
  list<string> tmp_list;

  tmp_list.push_back(names[a]);
  ret.push_back(make_pair("tof", tmp_list));

  return ret;
}

gatelist om_com(int a, const vector<string>& names) {
  gatelist ret;
  list<string> tmp_list;

  tmp_list.push_back(names[a]);
  ret.push_back(make_pair("H", tmp_list));
  ret.push_back(make_pair("P", tmp_list));
  ret.push_back(make_pair("H", tmp_list));
  ret.push_back(make_pair("P", tmp_list));
  ret.push_back(make_pair("H", tmp_list));
  ret.push_back(make_pair("P", tmp_list));

  return ret;
}

gatelist i_com(int a, const vector<string>& names) {
  gatelist ret;
  list<string> tmp_list;

  tmp_list.push_back(names[a]);
  ret.push_back(make_pair("tof", tmp_list));
  ret.push_back(make_pair("Z", tmp_list));
  ret.push_back(make_pair("Y", tmp_list));

  return ret;
}

// Make triangular to determine the rank (destructive)
int compute_rank_dest(int m, int n, vector<xor_func>& tmp) {
  int i, j;
  int ret = 0;

  // Make triangular
  for (i = 0; i < n; i++) {
    bool flg = false;
    for (j = ret; j < m; j++) {
      if (tmp[j].test(i)) {
        // If we haven't yet seen a vector with bit i set...
        if (!flg) {
          // If it wasn't the first vector we tried, swap to the front
          if (j != ret) swap(tmp[ret], tmp[j]);
          flg = true;
        } else {
          tmp[j] ^= tmp[ret];
        }
      }
    }
    if (flg) ret++;
  }

  return ret;
}

int compute_rank(int m, int n, const vector<xor_func>& bits) {
  int ret;
  // Make a copy of the bitset
  vector<xor_func> tmp(m);
  for(int i = 0; i < m; i++) {
    tmp[i] = bits[i];
  }
  ret = compute_rank_dest(m, n, tmp);
  return ret;
}

int compute_rank(int n, const vector<exponent> & expnts, const set<int> & lst) {
  int ret;
  int m = lst.size();

  vector<xor_func> tmp(m);
  for (int i = 0; i < m; i++) {
    tmp[i] = expnts[i].second;
  }
  ret = compute_rank_dest(m, n, tmp);
  return ret;
}

// Check linear independence of one vector wrt a matrix (destructive)
bool is_indep_dest(int n, const vector<xor_func>& bits, xor_func & a) {
  map<int, int> pivots;
  // Find all pivot columns
  for (int i = 0, j = 0; i < n && j < bits.size();) {
    if (bits[j].test(i)) {
      pivots[i] = j;
      i++;
      j++;
    } else {
      j++;
    }
  }
  
  for (int i = 0; i < n; i++) {
    if (a.test(i)) {
      map<int, int>::iterator it = pivots.find(i);
      if (it == pivots.end()) return true;
      else a ^= bits[(*it).second];
    }
  }

  return false;
}

bool is_indep(int n, const vector<xor_func>& bits, const xor_func & a) {
  xor_func tmp = a;
  return is_indep_dest(n, bits, tmp);
}

// Make echelon form
gatelist to_upper_echelon(int m, int n, vector<xor_func>& bits, vector<xor_func>* mat, const vector<string>& names) {
  gatelist acc;
  int i, j;
  int rank = 0;
  for (j = 0; j < m; j++) {
    if (bits[j].test(n)) {
      bits[j].reset(n);
      if (mat == NULL) acc.splice(acc.end(), x_com(j, names));
      else             (*mat)[j].set(m);
    }
  }

  // Make triangular
  for (i = 0; i < n; i++) {
    bool flg = false;
    for (j = rank; j < m; j++) {
      if (bits[j].test(i)) {
        // If we haven't yet seen a vector with bit i set...
        if (!flg) {
          // If it wasn't the first vector we tried, swap to the front
          if (j != rank) {
            swap(bits[rank], bits[j]);
            if (mat == NULL) acc.splice(acc.end(), swap_com(rank, j, names));
            else             swap((*mat)[rank], (*mat)[j]);
          }
          flg = true;
        } else {
          bits[j] ^= bits[rank];
          if (mat == NULL) acc.splice(acc.end(), xor_com(rank, j, names));
          else             (*mat)[j] ^= (*mat)[rank];
        }
      }
    }
    if (flg) rank++;
  }

  return acc;
}

gatelist to_lower_echelon(int m, int n, vector<xor_func>& bits, vector<xor_func>* mat, const vector<string>& names) {
  gatelist acc;
  int i, j;

  for (i = n-1; i > 0; i--) {
    for (j = i - 1; j >= 0; j--) {
      if (bits[j].test(i)) {
        bits[j] ^= bits[i];
        if (mat == NULL) acc.splice(acc.end(), xor_com(i, j, names));
        else              (*mat)[j] ^= (*mat)[i];
      }
    }
  }
  return acc;
}

// Expects two matrices in echelon form, the second being a subset of the
//   rowspace of the first. It then morphs the second matrix into the first
gatelist fix_basis(
    int m,
    int n,
    int k,
    const vector<xor_func>& fst,
    vector<xor_func>& snd,
    vector<xor_func>* mat,
    const vector<string>& names )
  {
  gatelist acc;
  int j = 0;
  bool flg = false;
  map<int, int> pivots;  // mapping from columns to rows that have that column as pivot
  for (int i = 0; i < n; i++) pivots[i] = -1;

  // First pass makes sure tmp has the same pivots as fst
  for (int i = 0; i < m; i++) {
    // Find the next pivot
    while (j < n && !fst[i].test(j)) j++;
    if (j < n) {
      pivots[j] = i;
      flg = false;
      for (int h = i; !flg && h < k; h++) {
        // We found a vector with the same pivot
        if (snd[h].test(j)) {
          flg = true;
          if (h != i) {
            swap(snd[h], snd[i]);
            if (mat == NULL) acc.splice(acc.end(), swap_com(h, i, names));
            else             swap((*mat)[h], (*mat)[i]);
          }
        }
      }
      // There was no vector with the same pivot
      if (!flg) {
        if (k >= m) {
          cout << "FATAL ERROR: second space not a subspace\n" << flush;
          exit(1);
        }
        snd[k] = fst[i];
        if (k != i) {
          swap(snd[k], snd[i]);
          if (mat == NULL) {
            acc.splice(acc.end(), swap_com(k, i, names));
          } else {
            swap((*mat)[k], (*mat)[i]);
          }
        }
        k++;
      }
    }
  }

  // Second pass makes each row of tmp equal to that row of fst
  for (int i = 0; i < m; i++) {
    for (int j = i +1; j < n; j++) {
      if (fst[i][j] != snd[i][j]) {
        if (pivots[j] == -1) {
          cout << "FATAL ERROR: cannot fix basis\n" << flush;
          exit(1);
        } else {
          snd[i] ^= snd[pivots[j]];
          if (mat == NULL) acc.splice(acc.end(), xor_com(pivots[j], i, names));
          else             (*mat)[i] ^= (*mat)[pivots[j]];
        }
      }
    }
    if (!(snd[i] == fst[i])) {
      cout << "FATAL ERROR: basis differs\n" << flush;
      exit(1);
    }
  }

  return acc;
}

// A := B^{-1} A
void compose(int num, vector<xor_func>& A, const vector<xor_func>& B) {
  auto tmp = vector<xor_func>(num);
  for (int i = 0; i < num; i++) {
    tmp[i] = B[i];
  }
  to_upper_echelon(num, num, tmp, &A, vector<string>());
  to_lower_echelon(num, num, tmp, &A, vector<string>());
}

//------------------------- CNOT synthesis methods

// Gaussian elimination based CNOT synthesis
gatelist gauss_CNOT_synth(int n, int m, vector<xor_func>& bits, const vector<string>& names) {
  gatelist lst;

  for (int j = 0; j < n; j++) {
    if (bits[j].test(n)) {
      bits[j].reset(n);
      lst.splice(lst.begin(), x_com(j, names));
    }
  }

  // Make triangular
  for (int i = 0; i < n; i++) {
    bool flg = false;
    for (int j = i; j < n + m; j++) {
      if (bits[j].test(i)) {
        // If we haven't yet seen a vector with bit i set...
        if (!flg) {
          // If it wasn't the first vector we tried, swap to the front
          if (j != i) {
            swap(bits[i], bits[j]);
            lst.splice(lst.begin(), swap_com(i, j, names));
          }
          flg = true;
        } else {
          bits[j] ^= bits[i];
          lst.splice(lst.begin(), xor_com(i, j, names));
        }
      }
    }
    if (!flg) {
      cout << "ERROR: not full rank\n";
      exit(1);
    }
  }

  //Finish the job
  for (int i = n-1; i > 0; i--) {
    for (int j = i - 1; j >= 0; j--) {
      if (bits[j].test(i)) {
        bits[j] ^= bits[i];
        lst.splice(lst.begin(), xor_com(i, j, names));
      }
    }
  }

  return lst;
}

// Patel/Markov/Hayes CNOT synthesis
gatelist Lwr_CNOT_synth(int n, int m, vector<xor_func>& bits, const vector<string>& names, bool rev) {
  gatelist acc;
  int sec, tmp, row, col, i;
  vector<int> patt(1<<m);

  for (sec = 0; sec < ceil(n / m); sec++) {

    for (i = 0; i < (1<<m); i++) {
      patt[i] = -1;
    }
    for (row = sec*m; row < n; row++) {
      tmp = 0;
      for (i = 0; i < m; i++) {
        if (bits[row].test(sec*m + i)) tmp += (1 << i);
      }
      if (patt[tmp] == -1) {
        patt[tmp] = row;
      } else if (tmp != 0) {
        bits[row] ^= bits[patt[tmp]];
        if (rev) acc.splice(acc.begin(), xor_com(row, patt[tmp], names));
        else acc.splice(acc.end(), xor_com(patt[tmp], row, names));
      }
    }

    for (col = sec*m; col < (sec+1)*m; col++) {
      for (row=col + 1; row < n; row++) {
        if (bits[row].test(col)) {
          if (not(bits[col].test(col))) {
            bits[col] ^= bits[row];
            bits[row] ^= bits[col];
            bits[col] ^= bits[row];
            if (rev) {
              acc.splice(acc.begin(), xor_com(col, row, names));
              acc.splice(acc.begin(), xor_com(row, col, names));
              acc.splice(acc.begin(), xor_com(col, row, names));
            } else {
              acc.splice(acc.end(), xor_com(row, col, names));
              acc.splice(acc.end(), xor_com(col, row, names));
              acc.splice(acc.end(), xor_com(row, col, names));
            }
          } else {
            bits[row] ^= bits[col];
            if (rev) acc.splice(acc.begin(), xor_com(row, col, names));
            else acc.splice(acc.end(), xor_com(col, row, names));
          }
        }
      }
    }
  }

  return acc;
}

gatelist CNOT_synth(int n, vector<xor_func>& bits, const vector<string>& names) {
  gatelist acc;
  int i, j, m = (int)(log((double)n) / (log(2) * 2));
  // When m <= 1, PMH is just Gaussian elimination, so default to it
  if (m <= 1) return gauss_CNOT_synth(n, 0, bits, names);

  for (j = 0; j < n; j++) {
    if (bits[j].test(n)) {
      bits[j].reset(n);
      acc.splice(acc.begin(), x_com(j, names));
    }
  }

  acc.splice(acc.end(), Lwr_CNOT_synth(n, m, bits, names, false));
  for (i = 0; i < n; i++) {
    for (j = i + 1; j < n; j++) {
      bits[j][i] = bits[i][j];
      bits[i].reset(j);
    }
  }
  acc.splice(acc.end(), Lwr_CNOT_synth(n, m, bits, names, true));
  acc.reverse();

  return acc;
}

gatelist global_phase_synth(int n, int phase, const vector<string>& names) {
  gatelist acc;
  int qubit = 0;

  if (phase % 2 == 1) {
    acc.splice(acc.end(), om_com(qubit, names));
    qubit = (qubit + 1) % n;
  }
  for (int i = phase / 2; i > 0; i--) {
    acc.splice(acc.end(), i_com(qubit, names));
    qubit = (qubit + 1) % n;
  }

  return acc;
}

// Construct a circuit for a given partition
gatelist construct_circuit(
    const vector<exponent> & phase,
    const partitioning & part,
    vector<xor_func>& in,
    const vector<xor_func>& out,
    int num,
    int dim,
    const vector<string>& names) {
  gatelist ret, tmp, rev;
  auto bits = vector<xor_func>(num);
  auto pre = vector<xor_func>(num);
  auto post = vector<xor_func>(num);
  set<int>::iterator ti;
  int i;
  bool flg = true;

  for (int i = 0; i < num; i++) {
    flg &= (in[i] == out[i]);
  }
  for (int i = 0; i < num; i++) {
    flg &= (in[i] == out[i]);
    if (synth_method != AD_HOC) {
      pre[i] = xor_func(num + 1, 0);
      post[i] = xor_func(num + 1, 0);
      pre[i].set(i);
      post[i].set(i);
    }
  }
  if (flg && (part.size() == 0)) return ret;

  // Reduce in to echelon form to decide on a basis
  if (synth_method == AD_HOC) {
    ret.splice(ret.end(), to_upper_echelon(num, dim, in, NULL, names));
  } else {
    to_upper_echelon(num, dim, in, &pre, vector<string>());
  }

  // For each partition... Compute *it, apply T gates, uncompute
  for (partitioning::const_iterator it = part.begin(); it != part.end(); it++) {
    for (ti = it->begin(), i = 0; i < num; i++) {
      if (i < it->size()) {
        bits[i] = phase[*ti].second;
        ti++;
      } else {
        bits[i] = xor_func(dim + 1, 0);
      }
    }

    // prepare the bits
    if (synth_method == AD_HOC) {
      tmp = to_upper_echelon(it->size(), dim, bits, NULL, names);
      tmp.splice(tmp.end(), fix_basis(num, dim, it->size(), in, bits, NULL, names));
      rev = tmp;
      rev.reverse();
      ret.splice(ret.end(), rev);
    } else {
      to_upper_echelon(it->size(), dim, bits, &post, vector<string>());
      fix_basis(num, dim, it->size(), in, bits, &post, vector<string>());
      compose(num, pre, post);
      if (synth_method == GAUSS) ret.splice(ret.end(), gauss_CNOT_synth(num, 0, pre, names));
      else if (synth_method == PMH) ret.splice(ret.end(), CNOT_synth(num, pre, names));
    }

    // apply the T gates
    list<string> tmp_lst;
    for (ti = it->begin(), i = 0; ti != it->end(); ti++, i++) {
      tmp_lst.clear();
      tmp_lst.push_back(names[i]);
      if (phase[*ti].first <= 4) {
        if (phase[*ti].first / 4 == 1) ret.push_back(make_pair("Z", tmp_lst));
        if (phase[*ti].first / 2 == 1) ret.push_back(make_pair("P", tmp_lst));
        if (phase[*ti].first % 2 == 1) ret.push_back(make_pair("T", tmp_lst));
      } else {
        if (phase[*ti].first == 5 || phase[*ti].first == 6) ret.push_back(make_pair("P*", tmp_lst));
        if (phase[*ti].first % 2 == 1) ret.push_back(make_pair("T*", tmp_lst));
      }
    }

    // unprepare the bits
    if (synth_method == AD_HOC) ret.splice(ret.end(), tmp);
    else {
      pre = std::move(post);
      post = vector<xor_func>(num);
      // re-initialize
      for (i = 0; i < num; i++) {
        post[i] = xor_func(num + 1, 0);
        post[i].set(i);
      }
    }
  }

  // Reduce out to the basis of in
  for (i = 0; i < num; i++) {
    bits[i] = out[i];
  }
  if (synth_method == AD_HOC) {
    tmp = to_upper_echelon(num, dim, bits, NULL, names);
    tmp.splice(tmp.end(), fix_basis(num, dim, num, in, bits, NULL, names));
    tmp.reverse();
    ret.splice(ret.end(), tmp);
  } else {
    to_upper_echelon(num, dim, bits, &post, vector<string>());
    fix_basis(num, dim, num, in, bits, &post, vector<string>());
    compose(num, pre, post);
    if (synth_method == GAUSS) ret.splice(ret.end(), gauss_CNOT_synth(num, 0, pre, names));
    else if (synth_method == PMH) ret.splice(ret.end(), CNOT_synth(num, pre, names));
  }
  return ret;
}

// Matroid oracle
bool ind_oracle::operator()(const vector<exponent> & expnts, const set<int> & lst) const {
  if (lst.size() > num) return false;
  if (lst.size() == 1 || (num - lst.size()) >= dim) return true;

  set<int>::const_iterator it;
  int i, j, rank = 0;
  auto tmp = vector<xor_func>(lst.size());

  for (i = 0, it = lst.begin(); it != lst.end(); it++, i++) {
    tmp[i] = expnts[*it].second;
  }

  for (i = 0; i < length; i++) {
    bool flg = false;
    for (j = rank; j < lst.size(); j++) {
      if (tmp[j].test(i)) {
        // If we haven't yet seen a vector with bit i set...
        if (!flg) {
          // If it wasn't the first vector we tried, swap to the front
          if (j != rank) swap(tmp[rank], tmp[j]);
          flg = true;
        } else {
          tmp[j] ^= tmp[rank];
        }
      }
    }
    if (flg) rank++;
  }

  return (num - lst.size()) >= (dim - rank);
}

// Shortcut to find a linearly dependent element faster
int ind_oracle::retrieve_lin_dep(const vector<exponent> & expnts, const set<int> & lst) const {
  set<int>::const_iterator it;
  int i, j, rank = 0, tmpr;
  map<int, int> mp;
  auto tmp = vector<xor_func>(lst.size());

  for (i = 0, it = lst.begin(); it != lst.end(); it++, i++) {
    tmp[i] = expnts[*it].second;
    mp[i] = *it;
  }

  for (j = 0; j < lst.size(); j++) {
    if (tmp[j].test(length)) tmp[j].reset(length);
  }

  for (i = 0; i < length; i++) {
    bool flg = false;
    for (j = rank; j < lst.size(); j++) {
      if (tmp[j].test(i)) {
        // If we haven't yet seen a vector with bit i set...
        if (!flg) {
          // If it wasn't the first vector we tried, swap to the front
          if (j != rank) {
            swap(tmp[rank], tmp[j]);
            tmpr = mp[rank];
            mp[rank] = mp[j];
            mp[j] = tmpr;
          }
          flg = true;
        } else {
          tmp[j] ^= tmp[rank];
          if (tmp[j].none()) return mp[j];
        }
      }
    }
    if (flg) rank++;
  }

  assert((num - lst.size()) >= (dim - rank));
  return -1;
}
