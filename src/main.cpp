#include "circuit.h"
#include <cstdio>
#include <boost/dynamic_bitset.hpp>

typedef boost::dynamic_bitset<>              xor_func;
typedef pair<char, xor_func > exponent;

class ind_oracle {
	private: 
		int n;
		int m;
	public:
    ind_oracle() { n = 0; m = 0; }
		ind_oracle(int nin, int min) { n = nin; m = min; }

		bool operator()(const vector<exponent> & expnts, const set<int> lst) {
			if (lst.size() > n) return false;
			if (lst.size() == 1 || (n - lst.size()) >= m) return true;

			set<int>::const_iterator it;
			int i, j, rank = 0;
      bool flg;
			xor_func * tmp = new xor_func[lst.size()];

			for (i = 0, it = lst.begin(); it != lst.end(); it++, i++) {
				tmp[i] = expnts[*it].second;
			}

      for (i = 0; i < n; i++) {
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
			return (n - lst.size()) >= (m - rank);
		}
};

int main() {
	dotqc circuit;
	circuit.input(cin);
  circuit.print();
  character c;
  c.parse_circuit(circuit);
  cout << "Done parsing\n" << flush;
  c.print();
/*
  matroid<exponent, ind_oracle> mat(c.phase_expts, ind_oracle(c.n + c.h + c.m, c.n + c.h));
  partitioning part = mat.partition_matroid();
  cout << part << "\n" << "\n" << flush;

  for (list<Hadamard>::iterator it = c.hadamards.begin(); it != c.hadamards.end(); it++) {
    partitioning tmp = freeze_partitions(part, it->in);
    cout << "frozen: " << tmp << "\nremaining: "<< part << "\n" << "\n" << flush;
  }
  */

  dotqc synth = c.synthesize();

	return 0;
}
