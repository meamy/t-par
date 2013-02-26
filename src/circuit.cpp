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
	for (i = 0; i < (n + m); i++) {
		if (i != 0) out << " ";
		out << names[i];
	}
  /*
	out << "> --> w^(";

	// Print the phase exponents
	for (it = phase_expts.begin(); it != phase_expts.end(); it++) {
		if (it != phase_expts.begin()) out << "+";
		out << (int)(it->first) << "*";
		for (int i = 0; i < n + h + m; i++) {
			if (it->second[i]) out << names[i];
		}
	}
	out << ")|";
  */
  out << "\n";
	for (it = phase_expts.begin(); it != phase_expts.end(); it++) {
    for (int i = 0; i < n + h + m; i++) {
      if (it->second[i]) out << 1;
      else out << 0;
    }
    out << "\n";
  }


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
    out << "H:" << ti->qubit << "-->" << ti->prep << " | ";
    for (list<exponent>::iterator iti = ti->in.begin(); iti != ti->in.end(); iti++) {
      for (int i = 0; i < n + h + m; i++) {
        if (iti->second[i]) out << 1;
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

	names = new string [n + m];
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

	vector<exponent>::iterator ti;
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
      new_h.qubit = *(it->second.begin());
      new_h.prep = i++;

      int tmp = name_map[*(it->second.begin())];
      int rank;

      wires[tmp] = xor_func (n + h, 0);
      rank = compute_rank(n + m, n + h, wires);
      // Check previous exponents to see if they're inconsistent
      for (ti = phase_expts.begin(); ti != phase_expts.end(); ti++) {
        if (ti->first != 0) {
          wires[tmp] = ti->second;
          if (compute_rank(n + m, n + h, wires) > rank) new_h.in.push_back(*ti);
        }
      }
      // Prepare the new value
      wires[tmp].set(new_h.prep);

      hadamards.push_back(new_h);
		} else {
			cout << "ERROR: not a {CNOT, T} circuit\n";
			phase_expts.clear();
			delete[] outputs;
		}
	}
}

  /*
dotqc cnot_t_circuit::synthesize() {
	list<pair<string, list<string> > > net;
	list<vector<exponent> > basis;
	partitioning part;
	vector<exponent> * tmp;
	dotqc ret;

	ret.n = n;
	ret.m = m;
	for (int i = 0; i < n + m; i++) {
		ret.names.insert(names[i]);
		if (i >= n) ret.zero[names[i]] = 1;
		else ret.zero[names[i]] = 0;
	}

	// Partition the phase exponents
	cerr << "Partitioning matroid... " << flush;
	part = partition_matroid<exponent, ind_oracle>(phase_expts, ind_oracle(n, n + m));
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

	return ret;
}
  */
