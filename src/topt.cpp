#include "circuit.h"

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

// Given a number of qubits, target rank, an incomplete basis,
//	add vectors to the basis to make it full rank
void make_full_rank(int n, int m, vector<exponent> & basis) {
	int i, j, rank = 0, k = basis.size();
	bool flg;

	// Make a copy of the bitset
	xor_func * tmp = new xor_func[k];
	for(i = 0; i < k; i++) {
		tmp[i] = basis[i].second;
	}

	// Make triangular
	for (i = 0; i < n; i++) {
		flg = false;
		for (j = rank; j < k; j++) {
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

	j = 0;
	for (i = 0; i < n; i++) {
		if (j < k && tmp[j].test(i)) j++;
		else {
			xor_func ins(n, 0);
			ins.set(i);
			basis.push_back(make_pair(0, ins));
		}
	}

	assert(basis.size() <= n + m);
	while (basis.size() < n + m) {
		basis.push_back(make_pair(0, xor_func(n, 0)));
	}
	assert(basis.size() == n + m);

	delete [] tmp;
}

// Synthesize a circuit computing "bits" by gaussian elimination
list<pair<string, list<string> > > gauss_circuit(int n, int m, xor_func * bits, string * names) {
	int i, j, k;
	list<pair<string, list<string> > > lst;
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
						tmp_list1.clear();
						tmp_list2.clear();
						tmp_list1.push_back(names[i]);
						tmp_list1.push_back(names[j]);
						tmp_list2.push_back(names[j]);
						tmp_list2.push_back(names[i]);
					  lst.push_back(make_pair("tof", tmp_list1));
					  lst.push_back(make_pair("tof", tmp_list2));
					  lst.push_back(make_pair("tof", tmp_list1));
					}
					flg = true;
				} else {
					bits[j] ^= bits[i];
					tmp_list1.clear();
					tmp_list1.push_back(names[i]);
					tmp_list1.push_back(names[j]);
					lst.push_back(make_pair("tof", tmp_list1));
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
				tmp_list1.clear();
				tmp_list1.push_back(names[i]);
				tmp_list1.push_back(names[j]);
				lst.push_back(make_pair("tof", tmp_list1));
			}
		}
	}

	return lst;
}

list<pair<string, list<string> > > 
compute_CNOT_network(int n, int m, const vector<exponent> & expnts, string * names) {
	int i, j;

	// Make a copy of the bitset
	xor_func * tmp = new xor_func[n + m];
	for (i = 0; i < n + m; i++) {
		tmp[i] = expnts[i].second;
	}

	// Uncompute cnot network
	list<pair<string, list<string> > > uncompute = gauss_circuit(n, m, tmp, names);
	// Compute cnot network
	list<pair<string, list<string> > > compute = uncompute;
	compute.reverse();

	list<string> tmp_lst;
	// Set the T gates
	for (i = 0; i < n+m; i++) {
		tmp_lst.clear();
		tmp_lst.push_back(names[i]);
		if (expnts[i].first <= 4) {
			if (expnts[i].first / 4 == 1) compute.push_back(make_pair("Z", tmp_lst));
			if (expnts[i].first / 2 == 1) compute.push_back(make_pair("P", tmp_lst));
			if (expnts[i].first % 2 == 1) compute.push_back(make_pair("T", tmp_lst));
		} else {
			if (expnts[i].first == 5 || expnts[i].first == 6) compute.push_back(make_pair("P", tmp_lst));
			if (expnts[i].first % 2 == 1) compute.push_back(make_pair("T*", tmp_lst));
		}
	}

	compute.splice(compute.end(), uncompute);
	delete[] tmp;
	return compute;
}

list<pair<string, list<string> > > 
compute_output_func(int n, int m, const xor_func * outputs, string * names) {
	int i, j;

	// Make a copy of the bitset
	xor_func * tmp = new xor_func[n + m];
	for (i = 0; i < n + m; i++) {
		tmp[i] = outputs[i];
	}
	// Uncompute cnot network
	list<pair<string, list<string> > > uncompute = gauss_circuit(n, m, tmp, names);

	delete[] tmp;
	return uncompute;
}

/*
class ind_oracle {
	private: 
		int dim;
		int mx;
	public:
		ind_oracle(int d, int m) {dim = d; mx = m;}
		bool operator()(const vector<exponent> & expnts, const set<int> lst) {
			if (lst.size() > mx) return false;
			if (lst.size() == 1) return true;

			set<int>::const_iterator it;
			int i, rank;
			xor_func * tmp = new xor_func[lst.size()];

//			cerr << "Testing set {\n";
			for (i = 0, it = lst.begin(); it != lst.end(); it++, i++) {
//				cerr << *it << ", ";
				tmp[i] = expnts[*it].second;
			}
//			cerr << "}\n" << flush;

			rank = compute_rank(lst.size(), dim, tmp);
			delete[] tmp;

			return (mx - lst.size()) >= (dim - rank);
		}
};
*/
