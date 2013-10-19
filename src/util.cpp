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

void print_wires(const xor_func * wires, int num, int dim) {
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
gatelist xor_com(int a, int b, const string * names) {
	list<string> tmp_list;
  gatelist ret;

	tmp_list.push_back(names[a]);
	tmp_list.push_back(names[b]);
	ret.push_back(make_pair("tof", tmp_list));

  return ret;
}

gatelist swap_com(int a, int b, const string * names) {
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

gatelist x_com(int a, const string * names) {
  gatelist ret;
	list<string> tmp_list;

  tmp_list.push_back(names[a]);
  ret.push_back(make_pair("tof", tmp_list));

  return ret;
}

// Make triangular to determine the rank
int compute_rank(int m, int n, const xor_func * bits) {
	int k, i, j;
	int ret = 0;
	bool flg;

	// Make a copy of the bitset
	xor_func * tmp = new xor_func[m];
	for(int i = 0; i < m; i++) {
		tmp[i] = bits[i];
	}

	// Make triangular
	for (i = 0; i < n; i++) {
		flg = false;
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

	delete [] tmp;
	return ret;
}

// Make echelon form
gatelist to_echelon(int m, int n, xor_func * bits, const string * names) {
  gatelist acc;
	int k, i, j;
	int rank = 0;
	bool flg;

  for (j = 0; j < m; j++) {
    if (bits[j].test(n)) {
      acc.splice(acc.end(), x_com(j, names));
      bits[j].reset(n);
    }
  }

	// Make triangular
	for (i = 0; i < n; i++) {
		flg = false;
		for (j = rank; j < m; j++) {
			if (bits[j].test(i)) {
				// If we haven't yet seen a vector with bit i set...
				if (!flg) {
					// If it wasn't the first vector we tried, swap to the front
					if (j != rank) {
            swap(bits[rank], bits[j]);
            acc.splice(acc.end(), swap_com(rank, j, names));
          }
					flg = true;
				} else {
					bits[j] ^= bits[rank];
          acc.splice(acc.end(), xor_com(rank, j, names));
				}
			}
		}
		if (flg) rank++;
	}

  return acc;
}


// Expects two matrices in echelon form, the second being a subset of the 
//   rowspace of the first. It then morphs the second matrix into the first
gatelist fix_basis(int m, int n, int k, const xor_func * fst, const xor_func * snd, const string * names) {
  list<pair<string, list<string> > > acc;
  int j = 0;
  bool flg;
  map<int, int> pivots;  // mapping from columns to rows that have that column as pivot
  for (int i = 0; i < n; i++) pivots[i] = -1;

	// Make a copy of the bitset
	xor_func * tmp = new xor_func[m];
	for(int i = 0; i < m; i++) {
    if (i < k) tmp[i] = snd[i];
    else       tmp[i] = xor_func(n + 1, 0);
	}

  // First pass makes sure tmp has the same pivots as fst
  for (int i = 0; i < m; i++) {
    // Find the next pivot
    while (j < n && !fst[i].test(j)) j++;
    if (j < n) {
      pivots[j] = i;
      flg = false;
      for (int h = i; !flg && h < k; h++) {
        // We found a vector with the same pivot
        if (tmp[h].test(j)) {
          flg = true;
          if (h != i) {
            swap(tmp[h], tmp[i]);
            acc.splice(acc.end(), swap_com(h, i, names));
          }
        }
      }
      // There was no vector with the same pivot
      if (!flg) {
        if (k >= m) {
          cout << "FATAL ERROR: second space not a subspace\n" << flush;
          exit(1);
        }
        tmp[k] = fst[i];
        if (k != i) {
          swap(tmp[k], tmp[i]);
          acc.splice(acc.end(), swap_com(k, i, names));
        }
        k++;
      }
    }
  }

  // Second pass makes each row of tmp equal to that row of fst
  for (int i = 0; i < m; i++) {
    for (int j = i +1; j < n; j++) {
      if (fst[i][j] != tmp[i][j]) {
        if (pivots[j] == -1) {
          cout << "FATAL ERROR: cannot fix basis\n" << flush;
          exit(1);
        } else {
          tmp[i] ^= tmp[pivots[j]];
          acc.splice(acc.end(), xor_com(pivots[j], i, names));
        }
      }
    }
    if (!(tmp[i] == fst[i])) {
      cout << "FATAL ERROR: basis differs\n" << flush;
      exit(1);
    }
  }

  delete [] tmp;
  return acc;
}

// Gaussian elimination based CNOT synthesis
gatelist gauss_CNOT_synth(int n, int m, xor_func * bits, const string * names) {
	int i, j, k;
	gatelist lst;
	list<string> tmp_list1, tmp_list2;

	// Make triangular
	for (i = 0; i < n; i++) {
		bool flg = false;
		for (j = i; j < n + m; j++) {
			if (bits[j].test(i)) {
				// If we haven't yet seen a vector with bit i set...
				if (!flg) {
					// If it wasn't the first vector we tried, swap to the front
					if (j != i) {
						swap(bits[i], bits[j]);
            lst.splice(lst.end(), swap_com(i, j, names));
					}
					flg = true;
				} else {
					bits[j] ^= bits[i];
          lst.splice(lst.end(), xor_com(i, j, names));
				}
			}
		}
		if (!flg) {
			cout << "ERROR: not full rank\n";
			exit(1);
		}
	}

	//Finish the job
	for (i = n-1; i > 0; i--) {
		for (j = i - 1; j >= 0; j--) {
			if (bits[j].test(i)) {
				bits[j] ^= bits[i];
        lst.splice(lst.end(), xor_com(i, j, names));
			}
		}
	}

	return lst;
}

// Patel/Markov/Hayes CNOT synthesis
gatelist Lwr_CNOT_synth(int n, int m, xor_func * bits, const string * names, bool rev) {
  gatelist acc;
  int sec, tmp, row, col, i;
  int patt[1<<m];

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

gatelist CNOT_synth(int n, xor_func * bits, const string * names) { 
  gatelist acc, tmp;
  int i, j, m = (int)(log((double)n) / (log(2) * 2));

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

// Construct a circuit for a given partition
gatelist construct_circuit(const vector<exponent> & phase, 
                           const partitioning & part, 
                           xor_func * in, 
                           const xor_func * out,
                           int num,
                           int dim,
                           const string * names) {
  gatelist ret, tmp, rev;
  xor_func * bits;
  set<int>::iterator ti;
  int i;
  bool flg = true;

  for (i = 0; i < num; i++) flg &= (in[i] == out[i]);
  if (flg && (part.size() == 0)) return ret;

  // Reduce in to echelon form to decide on a basis
  ret.splice(ret.end(), to_echelon(num, dim, in, names));

  // For each partition... Compute *it, apply T gates, uncompute
  for (partitioning::const_iterator it = part.begin(); it != part.end(); it++) {
    bits = new xor_func[it->size()];
    for (ti = it->begin(), i = 0; ti != it->end(); ti++, i++) {
      bits[i] = phase[*ti].second;
    }
    
    // prepare the bits
    tmp = to_echelon(it->size(), dim, bits, names);
    tmp.splice(tmp.end(), fix_basis(num, dim, it->size(), in, bits, names));
    rev = tmp;
    rev.reverse();
    ret.splice(ret.end(), rev);

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
    ret.splice(ret.end(), tmp);
    delete [] bits;
  }

  // Reduce out to the basis of in
  bits = new xor_func[num];
  for (i = 0; i < num; i++) {
    bits[i] = out[i];
  }
  tmp = to_echelon(num, dim, bits, names);
  tmp.splice(tmp.end(), fix_basis(num, dim, num, in, bits, names));
  tmp.reverse();
  ret.splice(ret.end(), tmp);
  delete [] bits;

  return ret;
}


void compute_function(xor_func * bits, int num, const gatelist & pre, const gatelist & post, const string * names) { 
  int i;
  gatelist::const_iterator it;
  map<string, int> name_map;

  // initialize the name mapping
  for (i = 0; i < num; i++) {
    name_map[names[i]] = i;
  }

  // initialize bits
  for (i = 0; i < num; i++) {
    bits[i] = xor_func(num + 1, 0);
    bits[i].set(i);
  }

  // iterate through pre and post
  for (it = pre.begin(); it != pre.end(); it++) {
    assert(it->first == "tof");
    if (it->second.size() == 1) {
      bits[name_map[*(it->second.begin())]].flip(num);
    } else if (it->second.size() == 2) {
      bits[name_map[*(++(it->second.begin()))]] ^= bits[name_map[*(it->second.begin())]];
    } else assert(false);
  }
  for (it = post.begin(); it != post.end(); it++) {
    assert(it->first == "tof");
    if (it->second.size() == 1) {
      bits[name_map[*(it->second.begin())]].flip(num);
    } else if (it->second.size() == 2) {
      bits[name_map[*(++(it->second.begin()))]] ^= bits[name_map[*(it->second.begin())]];
    } else assert(false);
  }
}

gatelist construct_circuit_efficient(const vector<exponent> & phase, 
                                     const partitioning & part, 
                                     xor_func * in, 
                                     const xor_func * out,
                                     int num,
                                     int dim,
                                     const string * names) {
  gatelist ret, tmp, rev, pre, post;
  xor_func * bits;
  xor_func * bits2;
  set<int>::iterator ti;
  int i;
  bool flg = true;
  bits2 = new xor_func[num];

  for (i = 0; i < num; i++) flg &= (in[i] == out[i]);
  if (flg && (part.size() == 0)) return ret;

  // Reduce in to echelon form to decide on a basis
  pre = to_echelon(num, dim, in, names);

  // For each partition... Compute *it, apply T gates, uncompute
  for (partitioning::const_iterator it = part.begin(); it != part.end(); it++) {

    bits = new xor_func[it->size()];
    for (ti = it->begin(), i = 0; ti != it->end(); ti++, i++) {
      bits[i] = phase[*ti].second;
    }
    
    // prepare the bits
    tmp = to_echelon(it->size(), dim, bits, names);
    tmp.splice(tmp.end(), fix_basis(num, dim, it->size(), in, bits, names));
    post = tmp;
    post.reverse();

    compute_function(bits2, num, pre, post, names);
    //ret.splice(ret.end(), gauss_CNOT_synth(num, 0, bits2, names));
    ret.splice(ret.end(), CNOT_synth(num, bits2, names));

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
    pre = tmp;

    delete [] bits;
  }

  // Reduce out to the basis of in
  bits = new xor_func[num];
  for (i = 0; i < num; i++) {
    bits[i] = out[i];
  }
  post = to_echelon(num, dim, bits, names);
  post.splice(tmp.end(), fix_basis(num, dim, num, in, bits, names));
  post.reverse();
  compute_function(bits2, num, pre, post, names);
  //ret.splice(ret.end(), gauss_CNOT_synth(num, 0, bits2, names));
  ret.splice(ret.end(), CNOT_synth(num, bits2, names));
  delete [] bits;
  delete [] bits2;

  return ret;
}

bool ind_oracle::operator()(const vector<exponent> & expnts, const set<int> & lst) const {
  if (lst.size() > num) return false;
  if (lst.size() == 1 || (num - lst.size()) >= dim) return true;

  set<int>::const_iterator it;
  int i, j, rank = 0;
  bool flg;
  xor_func * tmp = new xor_func[lst.size()];

  for (i = 0, it = lst.begin(); it != lst.end(); it++, i++) {
    tmp[i] = expnts[*it].second;
  }

  for (i = 0; i < length; i++) {
    flg = false;
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

  delete[] tmp;
  return (num - lst.size()) >= (dim - rank);
}
