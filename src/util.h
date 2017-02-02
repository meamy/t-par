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

#include <vector>
#include <boost/dynamic_bitset.hpp>
#include "partition.h"

typedef boost::dynamic_bitset<>            xor_func;
typedef pair<char, xor_func >              exponent;
typedef list<pair<string, list<string> > > gatelist;

enum synth_type { AD_HOC, GAUSS, PMH };

extern bool disp_log;
extern synth_type synth_method;

class ind_oracle {
  private: 
    int num;
    int dim;
    int length;
  public:
    ind_oracle() { num = 0; dim = 0; length = 0; }
    ind_oracle(int numin, int dimin, int lengthin) { num = numin; dim = dimin; length = lengthin; }

    void set_dim(int newdim) { dim = newdim; }
    int retrieve_lin_dep(const vector<exponent> & expnts, const set<int> & lst) const;

    bool operator()(const vector<exponent> & expnts, const set<int> & lst) const;
};

void print_wires(const vector<xor_func>& wires, int num, int dim);
int compute_rank_dest(int m, int n, vector<xor_func>& bits);
int compute_rank(int m, int n, const vector<xor_func>& bits);
int compute_rank(int n, const vector<exponent> & expnts, const set<int> & lst);
bool is_indep(int n, const vector<xor_func>& bits, const xor_func & a);

gatelist global_phase_synth(int n, int phase, const vector<string>& names);

gatelist construct_circuit(const vector<exponent> & phase, 
    const partitioning & part, 
    vector<xor_func>& in,
    const vector<xor_func>& out,
    int num,
    int dim,
    const vector<string>& names);
