#include <string>
#include <list>
#include <iostream>
#include <set>
#include <map>
#include <boost/dynamic_bitset.hpp>
#include "matroid.h"

using namespace std;

// Recognized gates are T, T*, P, P*, Z, Z*, and tof a b

typedef boost::dynamic_bitset<>              xor_func;
typedef pair<char, xor_func > exponent;

struct Hadamard {
  int qubit;        // Which qubit this hadamard is applied to
  int prep;         // Which "value" this hadamard prepares

  set<int> in;      // exponent terms that must be prepared before the hadamard
  xor_func * wires; // state of the wires when this hadamard is applied
};

// Internal representation of a .qc circuit and a {CNOT, T} circuit
struct dotqc {
	int n;                        // number of unknown inputs
	int m;                        // number of known inputs (initialized to |0>)
	list<string> names;           // names of qubits
	map<string, bool> zero;       // mapping from qubits to 0 (non-zero) or 1 (zero)
	list<pair<string, list<string> > > circ; // Circuit

	void input(istream& in);
	void output(ostream& out);
	void print() {output(cout);}
	void clear() {n = 0; m = 0; names.clear(); zero.clear(); circ.clear();}
  void remove_swaps();
  void print_stats();
};

struct character {
	int n;                        // number of unknown inputs
	int m;                        // number of zero-initialized ancilla qubits
  int h;                        // number of hadamards
	string           * names;     // names of qubits
  bool             * zero;      // Which qubits start as 0
  map<int, int>      val_map;   // which value corresponds to which qubit
	vector<exponent> phase_expts; // a list of exponents of \omega in the mapping
	xor_func         * outputs;   // the xors computed into each qubit
  // TODO: make this a dependency graph instead
  list<Hadamard>   hadamards;   // a list of the hadamards in the order we saw them

	void output(ostream& out);
	void print() {output(cout);}
  void parse_circuit(dotqc & input);
	dotqc synthesize();
};

