#include <string>
#include <list>
#include <iostream>
#include <set>
#include <map>
#include <boost/dynamic_bitset.hpp>

using namespace std;

// Recognized gates are T, T*, P, P*, Z, Z*, and tof a b

typedef boost::dynamic_bitset<>              xor_func;
typedef pair<char, xor_func > exponent;

struct Hadamard {
  string qubit;                 // Which qubit this hadamard is applied to
  int    prep;                  // Which "value" this hadamard prepares

  list<exponent> in; // exponents that must be prepared before the hadamard
};

// Internal representation of a .qc circuit and a {CNOT, T} circuit
struct dotqc {
	int n;                        // number of unknown inputs
	int m;                        // number of known inputs (initialized to |0>)
	set<string> names;            // names of qubits
	map<string, bool> zero;       // mapping from qubits to 0 (non-zero) or 1 (zero)
	list<pair<string, list<string> > > circ; // Circuit

	void input(istream& in);
	void output(ostream& out);
	void print() {output(cout);}
	void clear() {n = 0; m = 0; names.clear(); zero.clear(); circ.clear();}
	void append(pair<string, list<string> > gate);
};

struct character {
	int n;                        // number of unknown inputs
	int m;                        // number of zero-initialized ancilla qubits
  int h;                        // number of hadamards
	string           * names;     // names of qubits
	vector<exponent> phase_expts; // a list of exponents of \omega in the mapping
	xor_func         * outputs;   // the xors computed into each qubit
  // TODO: make this a dependency graph instead
  list<Hadamard>   hadamards;   // a list of the hadamards in the order we saw them

	void output(ostream& out);
	void print() {output(cout);}
  void parse_circuit(dotqc & input);
	//dotqc synthesize();
};

