#include "topt.cpp"

//----------------------------------------- DOTQC stuff

void ignore_white(istream& in) {
	while (in.peek() == ' ') in.ignore();
}

void dotqc::input(istream& in) {
	int i, j;
	string buf, tmp;
	list<string> namelist;
  n = 0;

	// Inputs
	while (buf != ".v") in >> buf;
	ignore_white(in);
	while(in.peek() != '\n') {
		in >> buf;
		names.insert(buf);
		zero[buf] = 1;
		ignore_white(in);
	}

	// Primary inputs
	while (buf != ".i") in >> buf;
	ignore_white(in);
	while (in.peek() != '\n') {
		n++;
		in >> buf;
		zero[buf] = 0;
		ignore_white(in);
	}

	m = names.size() - n;

	// Circuit
	while (buf != "BEGIN") in >> buf;
	in >> tmp;
	while (tmp != "END") {
		namelist.clear();
		// Build up a list of the applied qubits
		ignore_white(in);
		while (in.peek() != '\n') {
			in >> buf;
			if (names.find(buf) == names.end()) {
				cout << "ERROR: no such qubit \"" << buf << "\"\n" << flush;
				exit(1);
			} else {
				namelist.push_back(buf);
			}
			ignore_white(in);
		}
		circ.push_back(make_pair(tmp, namelist));
		in >> tmp;
	}
}

void dotqc::output(ostream& out) {
	int i;
	set<string>::iterator name_it;
	list<pair<string, list<string> > >::iterator it; 
	list<string>::iterator ti;

	// Inputs
	out << ".v";
	for (name_it = names.begin(); name_it != names.end(); name_it++) {
		out << " " << *name_it;
	}

	// Primary inputs
	out << "\n.i";
	for (name_it = names.begin(); name_it != names.end(); name_it++) {
		if (zero[*name_it] == 0) out << " " << *name_it;
	}

	// Circuit
	out << "\n\nBEGIN\n";
	for (it = circ.begin(); it != circ.end(); it++) {
		out << it->first;
		for (ti = (it->second).begin(); ti != (it->second).end(); ti++) {
			out << " " << *ti;
		}
		out << "\n";
	}
	out << "END\n";
}

// Append a gate to the end of a dotqc structure
void dotqc::append(pair<string, list<string> > gate) {
	circ.push_back(gate);

	for (list<string>::iterator it = gate.second.begin(); it != gate.second.end(); it++) {
		names.insert(*it);
	}
}

// Count the Hadamard gates
int count_h(dotqc & qc) {
  int ret = 0;
	list<pair<string,list<string> > >::iterator it;

  for (it = qc.circ.begin(); it != qc.circ.end(); it++) {
    if (it->first == "H") ret++;
  }

  return ret;
}

//-------------------------------------- End of DOTQC stuff

void character::output(ostream& out) {
	int i, j;
	bool flag;
	vector<exponent>::iterator it;

	out << "U|";
	for (i = 0; i < n; i++) {
		if (i != 0) out << " ";
		out << names[i];
	}
  for (; i < (n + m); i++) {
		if (i != 0) out << " ";
		out << "()";
  }

	out << "> --> w^(";

	// Print the phase exponents
	for (it = phase_expts.begin(); it != phase_expts.end(); it++) {
		if (it != phase_expts.begin()) out << "+";
		out << (int)(it->first) << "*";
		for (i = 0; i < n; i++) {
			if (it->second.test(i)) out << names[i];
		}
    for (; i < n + h; i++) {
			if (it->second.test(i)) out << names[i + m];
    }
	}
	out << ")|";

	// Print the output functions
	for (i = 0; i < (n + m); i++) {
		flag = false;
		out << "(";
    for (j = 0; j < n; j++) {
			if (outputs[i].test(j)) {
				if (flag) out << " ";
				out << names[j];
				flag = true;
      }
    }
		for (; j < n + h; j++) {
			if (outputs[i].test(j)) {
				if (flag) out << " ";
        out << j;
				flag = true;
			}
		}
		out << ")";
	}
	out << ">\n";

  // Print the Hadamards:
  for (list<Hadamard>::iterator ti = hadamards.begin(); ti != hadamards.end(); ti++) {
    out << "H:" << names[ti->qubit] << "-->" << ti->prep << " | ";
    for (set<int>::iterator iti = ti->in.begin(); iti != ti->in.end(); iti++) {
      for (int i = 0; i < n + h; i++) {
        if (phase_expts[*iti].second.test(i)) out << 1;
        else out << 0;
      }
      out << " | ";
    }
    out << "\n";
  }
}

void insert_phase (unsigned char c, xor_func f, vector<exponent> & phases) {
	vector<exponent>::iterator it;
	bool flg = false;
	for (it = phases.begin(); it != phases.end(); it++) {
		if (it->second == f) {
			it->first = (it->first + c) % 8;
			flg = true;
		}
	}
	if (!flg) {
		phases.push_back(make_pair(c, xor_func(f)));
	}
}

// Parse a {CNOT, T} circuit
void character::parse_circuit(dotqc & input) {
	int i, j, a, b, c;
	n = input.n;
	m = input.m;
  h = count_h(input);

  hadamards.clear();
  map<int, int> qubit_map;
	map<string, int> name_map, val_map;
	val_map["T"] = 1;
	val_map["T*"] = 7;
	val_map["P"] = 2;
	val_map["P*"] = 6;
	val_map["Z"] = 4;

  // Name stuff
	names = new string [n + m + h];
  i = 0;
  j = n;
	for (set<string>::iterator it = input.names.begin(); it != input.names.end(); it++) {
		if (input.zero[*it]) {
			name_map[*it] = j;
      val_map[*it]  = j;
			names[j++] = *it;
		} else {
			name_map[*it] = i;
      val_map[*it]  = i;
			names[i++] = *it;
		}
	}

  xor_func * wires = outputs = new xor_func [n + m];

	// Start each wire out with only its input in the sum
	for (j = 0; j < n; j++) {
		wires[j] = xor_func(n + h, 0);
		wires[j].set(j);
	}
	for (; j < (n + m); j++) wires[j] = xor_func (n + h, 0);

	bool flg;

  cout << h << "\n" << flush;
	list<pair<string,list<string> > >::iterator it;
	for (it = input.circ.begin(); it != input.circ.end(); it++) {
		flg = false;
		if (it->first == "tof") {
			wires[name_map[*(++(it->second.begin()))]] ^= wires[name_map[*(it->second.begin())]];
		} else if (it->first == "T" || it->first == "T*" || 
				       it->first == "P" || it->first == "P*" || 
							 (it->first == "Z" && it->second.size() == 1)) {
			a = name_map[*(it->second.begin())];
			insert_phase(val_map[it->first], wires[a], phase_expts);
		} else if (it->first == "Z" && it->second.size() == 3) {
			list<string>::iterator tmp_it = it->second.begin();
			a = name_map[*(tmp_it++)];
			b = name_map[*(tmp_it++)];
			c = name_map[*tmp_it];
			insert_phase(1, wires[a], phase_expts);
			insert_phase(1, wires[b], phase_expts);
			insert_phase(1, wires[c], phase_expts);
			insert_phase(7, wires[a] ^ wires[b], phase_expts);
			insert_phase(7, wires[a] ^ wires[c], phase_expts);
			insert_phase(7, wires[b] ^ wires[c], phase_expts);
			insert_phase(1, wires[a] ^ wires[b] ^ wires[c], phase_expts);
		} else if (it->first == "H") {
      // This WILL confuse you later on you idiot
      //   You zero the "destroyed" qubit, compute the rank, then replace the
      //   value with each of the phase exponents to see if the rank increases
      //   i.e. the system is inconsistent. This is so you don't have to make
      //   a new matrix -- i.e. instead of preparing the new value and computing
      //   rank, then adding each phase exponent and checking the rank you do it
      //   in place
      Hadamard new_h;
      new_h.qubit = name_map[*(it->second.begin())];
      new_h.prep  = i++;
      new_h.wires = new xor_func[n + m];
      for (j = 0; j < n + m; j++) {
        new_h.wires[j] = wires[j];
      }

      // Check previous exponents to see if they're inconsistent
    cout << "HELLO " << new_h.qubit << "\n" << flush;
      wires[new_h.qubit] = xor_func (n + h, 0);
      int rank = compute_rank(n + m, n + h, wires);
    cout << "GOODBYE\n" << flush;
      for (j = 0; j < phase_expts.size(); j++) {
        if (phase_expts[j].first != 0) {
          wires[new_h.qubit] = phase_expts[j].second;
          if (compute_rank(n + m, n + h, wires) > rank) new_h.in.insert(j);
        }
      }

      // Done creating the new hadamard
      hadamards.push_back(new_h);

      // Prepare the new value
      wires[new_h.qubit].reset();
      wires[new_h.qubit].set(new_h.prep);

      // Give this new value a name
      names[new_h.prep + m] = names[new_h.qubit];
      names[new_h.prep + m].append(to_string(new_h.prep));

		} else {
			cout << "ERROR: not a {CNOT, T} circuit\n";
			phase_expts.clear();
			delete[] outputs;
		}
	}
}

//---------------------------- Synthesis

class ind_oracle {
	private: 
		int num;
		int dim;
    int length;
	public:
    ind_oracle() { num = 0; dim = 0; length = 0; }
		ind_oracle(int numin, int dimin, int lengthin) { num = numin; dim = dimin; length = lengthin; }

    void set_dim(int newdim) { dim = newdim; }

		bool operator()(const vector<exponent> & expnts, const set<int> & lst) const {
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
};

dotqc character::synthesize() {
  partitioning floats, frozen; 
	dotqc ret;
  xor_func mask(n + h, 0);      // Tells us what values we have prepared
  xor_func wires[n + m];        // Current state of the wires
  list<int> remaining;          // Which terms we still have to partition
  int dim = n, tmp;
  ind_oracle oracle(n + m, dim, n + h);

  // initialize some stuff
	ret.n = n;
	ret.m = m;
	for (int i = 0; i < n + m; i++) {
		ret.names.insert(names[i]);
    wires[i] = xor_func(n + h, 0);
    if (i < n) {
      ret.zero[names[i]] = 0;
      wires[i].set(i);
      mask.set(i);
    } else {
      ret.zero[names[i]] = 1;
    }
  }

  // initialize the remaining list
  for (int i = 0; i < phase_expts.size(); i++) remaining.push_back(i);

  // create an initial partition
  cerr << "Adding new functions to the partition... " << flush;
  for (list<int>::iterator it = remaining.begin(); it != remaining.end();) {
    xor_func tmp = (~mask) & (phase_expts[*it].second);
    if (tmp.none()) {
      add_to_partition(floats, *it, phase_expts, oracle);
      it = remaining.erase(it);
    } else it++;
  }
  cerr << floats << "\n" << flush;
  cout << compute_rank(n +m,n + h,wires) << "\n";

  for (list<Hadamard>::iterator it = hadamards.begin(); it != hadamards.end(); it++) {
    // 1. freeze partitions that are not disjoint from the hadamard input
    // 2. construct CNOT+T circuit
    // 3. apply the hadamard gate
    // 4. add new functions to the partition

    cerr << "Freezing partitions... " << flush;
    frozen = freeze_partitions(floats, it->in);
    cerr << frozen << "\n" << flush;

    cerr << "Constructing {CNOT, T} subcircuit... " << flush;
    //circ = construct_circuit(frozen, wires, ith->state);
	  for (int i = 0; i < n + m; i++) {
      wires[i] = it->wires[i];
    }
    cerr << "\n" << flush;

    cerr << "Applying Hadamard gate... " << flush;
    //circ.push_back(make_pair("H", list<string>(1, it->qubit)));
    wires[it->qubit] = xor_func(n + h, it->prep);
    mask.set(it->prep);
    cerr << "\n" << flush;

    cerr << "Checking for increase in dimension... " << flush;
    tmp = compute_rank(n + m, n + h, wires);
    if (tmp > dim) {
      cerr << "Increased to " << tmp << "\n" << "Repartitioning... " << flush;
      dim = tmp;
      oracle.set_dim(dim);
      repartition(floats, phase_expts, oracle);
      cerr << floats << flush;
    }
    cerr << "\n" << flush;

    cerr << "Adding new functions to the partition... " << flush;
    for (list<int>::iterator it = remaining.begin(); it != remaining.end();) {
      xor_func tmp = (~mask) & (phase_expts[*it].second);
      if (tmp.none()) {
        add_to_partition(floats, *it, phase_expts, oracle);
        it = remaining.erase(it);
      } else it++;
    }
    cerr << floats << "\n" << flush;
  }

  cerr << "Constructing final {CNOT, T} subcircuit... " << floats << "\n" << flush;
  //circ = construct_circuit(frozen, wires, outputs);

    /*
    // Partition the phase exponents
    cerr << "Partitioning matroid... " << flush;
    part = mat.partition_matroid();
    cerr << part << "\n" << flush;

    // Fill the partitions to have full rank
    cerr << "Filling out missing rank... " << flush;
    for (partitioning::iterator it = part.begin(); it != part.end(); it++) {
      tmp = new vector<exponent>();
      tmp->reserve(n + m);
      for (set<int>::iterator ti = it->begin(); ti != it->end(); ti++) {
        tmp->push_back(phase_expts[*ti]);
      }
      make_full_rank(n, m, *tmp);
      basis.push_back(*tmp);
      delete tmp;
    }
    cerr << "\n" << flush;

    // Construct the new circuit
    cerr << "Constructing circuit... " << flush;
    for (list<vector<exponent> >::iterator it = basis.begin(); it != basis.end(); it++) {
      net = compute_CNOT_network(n, m, *it, names);
      ret.circ.splice(ret.circ.end(), net);
    }
    net = compute_output_func(n, m, outputs, names);
    ret.circ.splice(ret.circ.end(), net);
    cerr << "\n" << flush;
    */

  return ret;
}
