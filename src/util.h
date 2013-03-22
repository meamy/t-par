#include <vector>
#include <boost/dynamic_bitset.hpp>
#include "partition.h"

typedef boost::dynamic_bitset<>            xor_func;
typedef pair<char, xor_func >              exponent;
typedef list<pair<string, list<string> > > gatelist;

class ind_oracle {
	private: 
		int num;
		int dim;
    int length;
	public:
    ind_oracle() { num = 0; dim = 0; length = 0; }
		ind_oracle(int numin, int dimin, int lengthin) { num = numin; dim = dimin; length = lengthin; }

    void set_dim(int newdim) { dim = newdim; }

		bool operator()(const vector<exponent> & expnts, const set<int> & lst) const;
};

int compute_rank(int m, int n, const xor_func * bits);

gatelist construct_circuit(const vector<exponent> & phase, 
                           const partitioning & part, 
                           xor_func * in, 
                           const xor_func * out,
                           int num,
                           int dim,
                           const string * names);
