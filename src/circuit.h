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

#include <string>
#include <list>
#include <iostream>
#include <set>
#include <map>
#include "matroid.h"
#include "util.h"

using namespace std;

// Recognized gates are T, T*, P, P*, Z, Z*, Z a b c, tof a b, tof a, X, H

// Internal representation of a .qc circuit circuit
struct dotqc {
  int n;                   // number of unknown inputs
  int m;                   // number of known inputs (initialized to |0>)
  list<string> names;      // names of qubits
  map<string, bool> zero;  // mapping from qubits to 0 (non-zero) or 1 (zero)
  gatelist circ;           // Circuit

  void input(istream& in);
  void output(ostream& out);
  void print() {output(cout);}
  void clear() {n = 0; m = 0; names.clear(); zero.clear(); circ.clear();}
  void append(pair<string, list<string> > gate);
  void remove_swaps();
  int count_depth();
  int count_t_depth();
  void print_stats();
  void remove_ids();
};

// ------------------------- Hadamard version
struct Hadamard {
  int qubit;        // Which qubit this hadamard is applied to
  int prep;         // Which "value" this hadamard prepares

  set<int> in;      // exponent terms that must be prepared before the hadamard
  vector<xor_func> wires; // state of the wires when this hadamard is applied
};

// Characteristic of a circuit
struct character {
  int n;                        // number of unknown inputs
  int m;                        // number of zero-initialized ancilla qubits
  int h;                        // number of hadamards
  vector<string>     names;     // names of qubits
  vector<bool>       zero;      // Which qubits start as 0
  map<int, int>      val_map;   // which value corresponds to which qubit
  vector<exponent> phase_expts; // a list of exponents of \omega in the mapping
  vector<xor_func> outputs;   // the xors computed into each qubit
  // TODO: make this a dependency graph instead
  list<Hadamard>   hadamards;   // a list of the hadamards in the order we saw them

  void output(ostream& out);
  void print() {output(cout);}
  void parse_circuit(dotqc & input);
  void add_ancillae(int num);
  void remove_x();
  dotqc synthesize();
  dotqc synthesize_unbounded();
};

// ------------------------- {CNOT, T} version

enum circuit_type { CNOTT, OTHER, UNKNOWN };

struct metacircuit {
  int n;                        // number of unknown inputs
  int m;                        // number of known inputs (initialized to |0>)
  list<string> names;           // names of qubits
  map<string, bool> zero;       // mapping from qubits to 0 (non-zero) or 1 (zero)
  list<pair<circuit_type, dotqc> > circuit_list;     // A list of subcircuits

  void partition_dotqc(dotqc & input);
  void output(ostream& out);
  void print() {output(cout);}
  void optimize();
  dotqc to_dotqc();
};
