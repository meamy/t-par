#include "topt.cpp"
#include <algorithm>

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
		names.push_back(buf);
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
			if (find(names.begin(), names.end(), buf) == names.end()) {
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
	list<string>::iterator name_it;
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
		if (i != 0)  out << " ";
    if (zero[i]) out << "()";
    else out << names[i];
	}

	out << "> --> w^(";

	// Print the phase exponents
	for (it = phase_expts.begin(); it != phase_expts.end(); it++) {
		if (it != phase_expts.begin()) out << "+";
		out << (int)(it->first) << "*";
		for (i = 0; i < (n + h); i++) {
			if (it->second.test(i)) out << names[val_map[i]];
		}
	}
	out << ")|";

	// Print the output functions
	for (i = 0; i < (n + m); i++) {
		flag = false;
		out << "(";
    for (j = 0; j < (n + h); j++) {
			if (outputs[i].test(j)) {
				if (flag) out << " ";
				out << names[val_map[j]];
				flag = true;
      }
    }
		out << ")";
	}
	out << ">\n";

  // Print the Hadamards:
  for (list<Hadamard>::iterator ti = hadamards.begin(); ti != hadamards.end(); ti++) {
    out << "H:" << names[ti->qubit] << "-->" << ti->prep << " | ";
    /*
    for (set<int>::iterator iti = ti->in.begin(); iti != ti->in.end(); iti++) {
      for (int i = 0; i < n + h; i++) {
        if (phase_expts[*iti].second.test(i)) out << 1;
        else out << 0;
      }
      out << " | ";
    }
    */
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
// NOTE: a qubit's number is NOT the same as the bit it begins with set
void character::parse_circuit(dotqc & input) {
	int i, j, a, b, c, name_max = 0, val_max = 0;
	n = input.n;
	m = input.m;
  h = count_h(input);

  hadamards.clear();
	map<string, int> name_map, gate_lookup;
	gate_lookup["T"] = 1;
	gate_lookup["T*"] = 7;
	gate_lookup["P"] = 2;
	gate_lookup["P*"] = 6;
	gate_lookup["Z"] = 4;

  // Initialize names and wires
	names = new string [n + m + h];
  zero  = new bool   [n + m];
  xor_func * wires = outputs = new xor_func [n + m];
	for (list<string>::iterator it = input.names.begin(); it != input.names.end(); it++) {
    // name_map maps a name to a wire
		name_map[*it] = name_max;
    // names maps a wire to a name
		names[name_max] = *it;
    // zero mapping
    zero[name_max]  = input.zero[*it];
    // each wire has an initial value j, unless it starts in the 0 state
    wires[name_max] = xor_func(n + h, 0);
    if (!zero[name_max]) {
      wires[name_max].set(val_max);
      val_map[val_max++] = name_max;
    }
    name_max++;
	}

	bool flg;

	list<pair<string,list<string> > >::iterator it;
	for (it = input.circ.begin(); it != input.circ.end(); it++) {
		flg = false;
		if (it->first == "tof") {
			wires[name_map[*(++(it->second.begin()))]] ^= wires[name_map[*(it->second.begin())]];
		} else if (it->first == "T" || it->first == "T*" || 
				       it->first == "P" || it->first == "P*" || 
							 (it->first == "Z" && it->second.size() == 1)) {
			a = name_map[*(it->second.begin())];
			insert_phase(gate_lookup[it->first], wires[a], phase_expts);
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
      new_h.prep  = val_max++;
      new_h.wires = new xor_func[n + m];
      for (i = 0; i < n + m; i++) {
        new_h.wires[i] = wires[i];
      }

      // Check previous exponents to see if they're inconsistent
      wires[new_h.qubit].reset();
      int rank = compute_rank(n + m, n + h, wires);
      for (i = 0; i < phase_expts.size(); i++) {
        if (phase_expts[i].first != 0) {
          wires[new_h.qubit] = phase_expts[i].second;
          if (compute_rank(n + m, n + h, wires) > rank) new_h.in.insert(i);
        }
      }

      // Done creating the new hadamard
      hadamards.push_back(new_h);

      // Prepare the new value
      wires[new_h.qubit].reset();
      wires[new_h.qubit].set(new_h.prep);

      // Give this new value a name
      val_map[new_h.prep] = name_max;
      names[name_max] = names[new_h.qubit];
      names[name_max++].append(to_string(new_h.prep));

		} else {
			cout << "ERROR: not a {CNOT, T} circuit\n";
			phase_expts.clear();
			delete[] outputs;
		}
	}
}

//---------------------------- Synthesis

dotqc character::synthesize() {
  partitioning floats, frozen; 
	dotqc ret;
  xor_func mask(n + h, 0);      // Tells us what values we have prepared
  xor_func wires[n + m];        // Current state of the wires
  list<int> remaining;          // Which terms we still have to partition
  int dim = n, tmp, tdepth = 0;
  ind_oracle oracle(n + m, dim, n + h);
  list<pair<string, list<string> > > circ;

  // initialize some stuff
	ret.n = n;
	ret.m = m;
	for (int i = 0, j = 0; i < n + m; i++) {
		ret.names.push_back(names[i]);
    ret.zero[names[i]] = zero[i];
    wires[i] = xor_func(n + h, 0);
    if (!zero[i]) {
      wires[i].set(j);
      mask.set(j++);
    }
  }

  // initialize the remaining list
  for (int i = 0; i < phase_expts.size(); i++) {
    if (phase_expts[i].first != 0) remaining.push_back(i);
  }

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

  for (list<Hadamard>::iterator it = hadamards.begin(); it != hadamards.end(); it++) {
    // 1. freeze partitions that are not disjoint from the hadamard input
    // 2. construct CNOT+T circuit
    // 3. apply the hadamard gate
    // 4. add new functions to the partition

    cerr << "Freezing partitions... " << flush;
    frozen = freeze_partitions(floats, it->in);
    for (partitioning::iterator it = frozen.begin(); it != frozen.end(); it++) {
      for (set<int>::iterator ti = it->begin(); ti != it->end(); ti++) {
        if ((phase_expts[*ti].first % 2) != 0) {
          tdepth += 1;
          break;
        }
      }
    }
    cerr << frozen << "\n" << flush;

    cerr << "Constructing {CNOT, T} subcircuit... " << flush;
    if (!frozen.empty()) {
      ret.circ.splice(ret.circ.end(), construct_circuit(*this, frozen, wires, it->wires));
    }
	  for (int i = 0; i < n + m; i++) {
      wires[i] = it->wires[i];
    }
    cerr << "\n" << flush;

    cerr << "Applying Hadamard gate... " << flush;
    ret.circ.push_back(make_pair("H", list<string>(1, names[it->qubit])));
    wires[it->qubit].reset();
    wires[it->qubit].set(it->prep);
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
  tdepth += floats.size();
  cerr << "T-depth: " << tdepth << "\n" << flush;
  return ret;
}

//----------- This will likely need to be finished to get reasonable CNOT networks
void dotqc::remove_swaps() {
  list<pair<string, list<string> > >::iterator it, tt, ttt;
  list<string>::iterator iti;
  string q1, q2;
  bool flg;
  int i;

  for (it = circ.begin(), i = 0; i < circ.size() - 3;) {
    flg = false;
    if (it->first == "tof" && it->second.size() == 2) {
      iti = it->second.begin();
      q1 = *(iti++);
      q2 = *iti;

      tt = it;
      tt++;
      ttt = tt;
      ttt++;
      if (tt->first == "tof" && tt->second.size() == 2 && ttt->first == "tof" && ttt->second.size() == 2) {
        flg = true;
        iti = tt->second.begin();
        flg &= *(iti++) == q2;
        flg &= *(iti) == q1;
        iti = ttt->second.begin();
        flg &= *(iti++) == q1;
        flg &= *(iti) == q2;

        if (flg) {
          circ.erase(it);
          circ.erase(tt);
          it = ttt = circ.erase(ttt);
          while (ttt != circ.end()) {
            for (iti = ttt->second.begin(); iti != ttt->second.end(); iti++) {
              if (*iti == q1) *iti = q2;
              else if (*iti == q2) *iti = q1;
            }
            ttt++;
          }
        }
      }
    }
    if (!flg) {
      i++;
      it++;
    }
  }
}
