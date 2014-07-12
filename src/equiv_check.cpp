#include <z3++.h>
#include "equiv_check.h"

using namespace z3;

typedef boost::dynamic_bitset<> bitvec;

void subvec(bitvec & x, const bitvec & y, int length) {
  for (int i = 0; i < length; i++) x.set(i, y.test(i));
}

// Generate the partial sums on input bits
vector<expr> gen_partial_sums(character & ch, context & c, solver & s, map<bitvec, int> & mp, expr & inputs) {
  vector<exponent>::iterator it;
  list<Hadamard>::iterator ti;
  vector<expr> ret;
  bitvec tmp(ch.n);
  int j = 0;

  for (it = ch.phase_expts.begin(); it != ch.phase_expts.end(); it++) {
    // Copy first n bits
    subvec(tmp, it->second, ch.n);
    // If it hasn't been seen yet
    if (mp.find(tmp) == mp.end()) {
      expr sum = c.bv_val(0, 1);
      // Generate new variable
      string str = "partial_sum" + std::to_string(j);
      ret.push_back(c.bv_const(str.c_str(), 1));
      // Generate sum
      for (int i = 0; i < ch.n; i++) {
        if (tmp.test(i)) 
          sum = sum ^ to_expr(c, Z3_mk_extract(c, i, i, inputs));
      }
      s.add(ret[j] == sum);
      mp[tmp] = j++;
    }
  }
  // For the Hadamards
  for (ti = ch.hadamards.begin(); ti != ch.hadamards.end(); ti++) {
    // Copy first n bits
    subvec(tmp, ti->wires[ti->qubit], ch.n);
    // If it hasn't been seen yet
    if (mp.find(tmp) == mp.end()) {
      expr sum = c.bv_val(0, 1);
      // Generate new variable
      string str = "partial_sum" + std::to_string(j);
      ret.push_back(c.bv_const(str.c_str(), 1));
      // Generate sum
      for (int i = 0; i < ch.n; i++) {
        if (ti->wires[ti->qubit].test(i)) 
          sum = sum ^ to_expr(c, Z3_mk_extract(c, i, i, inputs));
      }
      s.add(ret[j] == sum);
      mp[tmp] = j++;
    }
  }
  // For the outputs
  for (int l = 0; l < ch.n + ch.m; l++) {
    // Copy first n bits
    subvec(tmp, ch.outputs[l], ch.n);
    // If it hasn't been seen yet
    if (mp.find(tmp) == mp.end()) {
      expr sum = c.bv_val(0, 1);
      // Generate new variable
      string str = "partial_sum" + std::to_string(j);
      ret.push_back(c.bv_const(str.c_str(), 1));
      // Generate sum
      for (int i = 0; i < ch.n; i++) {
        if (ch.outputs[j].test(i)) 
          sum = sum ^ to_expr(c, Z3_mk_extract(c, i, i, inputs));
      }
      s.add(ret[j] == sum);
      mp[tmp] = j++;
    }
  }
  return ret;
}

// Evaluate the xor x on inputs and path
expr xor_sum(character & ch, context & c, const xor_func & x, expr & inputs, bitvec & path) {
  int i;
  bool s = x.test(ch.n + ch.h);
  expr sum = c.bv_val(0, 1);
  for (i = 0; i < ch.n + ch.h; i++) {
    if (x.test(i)) {
      if (i < ch.n) sum = sum ^ to_expr(c, Z3_mk_extract(c, i, i, inputs));
      else if (path.test(i - ch.n)) s = !s;
    }
  }
  
  return s?sum ^ c.bv_val(1,1):sum;
}

// Check whether inputs and path satisfy outputs = inputs
expr assert_output_eq(character & ch, context & c, expr & inputs, bitvec path) {
  int i;
  expr ret = c.bv_val(1, 1);
  for (i = 0; i < ch.n + ch.m; i++) {
    if (i < ch.n) 
      ret = ret & ~(to_expr(c, Z3_mk_extract(c, i, i, inputs)) 
                    ^ xor_sum(ch, c, ch.outputs[i], inputs, path));
    else 
      ret = ret & ~(xor_sum(ch, c, ch.outputs[i], inputs, path));
  }
  return ret;
}

// Evaluate the polynomial p on inputs and path
expr gen_sums(character & ch, context & c, expr & inputs, bitvec path) {
  vector<exponent>::iterator it;
  list<Hadamard>::iterator ti;
  expr ret = c.bv_val(0, 3);
  for (it = ch.phase_expts.begin(); it != ch.phase_expts.end(); it++) {
    // add up the bits in *it
    ret = ret + ite(xor_sum(ch, c, it->second, inputs, path) == c.bv_val(1, 1),
                    c.bv_val((unsigned int)it->first, 3),
                    c.bv_val(0, 3));
  }
  
  for (ti = ch.hadamards.begin(); ti != ch.hadamards.end(); ti++) {
    if (path.test(ti->prep - ch.n))
      ret = ret + ite(xor_sum(ch, c, ti->wires[ti->qubit], inputs, path) == c.bv_val(1, 1),
                      c.bv_val(4, 3),
                      c.bv_val(0, 3));
  } 
  return ret;
}

expr gen_sums(character & ch, context & c, vector<expr> & partials, bitvec path, map<bitvec, int> & mp) {
  map<bitvec, int> q;
  bitvec tmp(ch.n + 1);
  vector<exponent>::iterator it;
  list<Hadamard>::iterator ti;
  map<bitvec, int>::iterator iti;
  expr ret = c.bv_val(0, 3);

  for (it = ch.phase_expts.begin(); it != ch.phase_expts.end(); it++) {
    subvec(tmp, it->second, ch.n);
    tmp[ch.n] = it->second.test(ch.n + ch.h);
    for (int i = ch.n; i < ch.n + ch.h; i++) {
      if (it->second.test(i) && path.test(i - ch.n)) tmp.flip(ch.n);
    }
    q[tmp] = (q[tmp] + it->first) % 8;
  }
  for (ti = ch.hadamards.begin(); ti != ch.hadamards.end(); ti++) {
    if (path.test(ti->prep - ch.n)) {
      subvec(tmp, ti->wires[ti->qubit], ch.n);
      tmp[ch.n] = ti->wires[ti->qubit].test(ch.n + ch.h);
      for (int i = ch.n; i < ch.n + ch.h; i++) {
        if (ti->wires[ti->qubit].test(i) && path.test(i - ch.n)) tmp.flip(ch.n);
      }
      q[tmp] = (q[tmp] + 4) % 8;
    }
  }

  tmp = bitvec(ch.n);
  for (iti = q.begin(); iti != q.end(); iti++) {
    subvec(tmp, iti->first, ch.n);
    ret = ret + ite(partials[mp[tmp]] ==
                      (iti->first.test(ch.n)?c.bv_val(0, 1):c.bv_val(1, 1)),
                    c.bv_val((unsigned int)iti->second, 3),
                    c.bv_val(0, 3));
  }

  return ret;
}

// Count the number of path paths to output = input with phase omega^val
expr count(character & ch, context & c, expr & inputs, expr *path_evals[], func_decl & f) {
  expr ret = c.int_val(0);
  for (unsigned long i = 0; i < (1 << ch.h); i++) {
    ret = ite((assert_output_eq(ch, c, inputs, bitvec(ch.h, i)) == c.bv_val(1, 1)),
              ret + f(*(path_evals[i])),
              ret);
  }
  return ret;
}

void unquantified(character & ch) {
  config cfg;
  context c(cfg);
  tactic t = tactic(c, "simplify") &
             tactic(c, "solve-eqs") &
             tactic(c, "bit-blast") &
             tactic(c, "smt");
  solver s = t.mk_solver();
  map<bitvec, int> mp;

  expr inputs = c.bv_const("input", ch.n);
  expr *path_evals[1 << ch.h];
  cout << ch.h << " Hadamards\n";
  if (ch.h <= 62) {
    for (unsigned long i = 0; i < (1 << ch.h); i++) {
      for (int j = 0; j < ch.h; j++) {
        if (i == (1 << j)) cout << j;
      }
      string str = "e" + std::to_string(i);
      path_evals[i] = new expr(c.bv_const(str.c_str(), 3));
      s.add(*(path_evals[i]) == gen_sums(ch, c, inputs, bitvec(ch.h, i)));
    } 
  } else {
    cout << "ERROR: Cannot handle more than 62 Hadamards";
    exit(1);
  }
  cout << "\nDone phase 1\n";

  func_decl f = z3::function("count_fun", c.bv_sort(3), c.int_sort());
  if (ch.h % 2 == 1) {
    // check #w^1 - #w^3 - #w^5 + #w^7 = 2^{h+1}
    s.add(f(c.bv_val(0, 3)) == 0);
    s.add(f(c.bv_val(1, 3)) == 1);
    s.add(f(c.bv_val(2, 3)) == 0);
    s.add(f(c.bv_val(3, 3)) == -1);
    s.add(f(c.bv_val(4, 3)) == 0);
    s.add(f(c.bv_val(5, 3)) == -1);
    s.add(f(c.bv_val(6, 3)) == 0);
    s.add(f(c.bv_val(7, 3)) == 1);
    s.add(count(ch, c, inputs, path_evals, f)
          != c.int_val(1 << ((ch.h+1) / 2)));
  } else {
    // check #w^0 - #w^4 = 2^h
    s.add(f(c.bv_val(0, 3)) == 1);
    s.add(f(c.bv_val(1, 3)) == 0);
    s.add(f(c.bv_val(2, 3)) == 0);
    s.add(f(c.bv_val(3, 3)) == 0);
    s.add(f(c.bv_val(4, 3)) == -1);
    s.add(f(c.bv_val(5, 3)) == 0);
    s.add(f(c.bv_val(6, 3)) == 0);
    s.add(f(c.bv_val(7, 3)) == 0);
    s.add(count(ch, c, inputs, path_evals, f)
          != c.int_val(1 << (ch.h / 2)));
  }

  check_result chk = s.check();
  if (chk == sat) {
     std::cout << "Circuits are not equivalent:\n";
     model m = s.get_model();
     std::cout << "  Counterexample is |" << m.eval(inputs) << ">\n";
  } else if (chk == unsat) {
     std::cout << "Circuits are equivalent\n";
  } else {
     std::cout << "Problem not solved\n";
  }
}

void unquantified_sum(character & ch) {
  config cfg;
  context c(cfg);
  tactic t = tactic(c, "simplify") &
             tactic(c, "solve-eqs") &
     //        tactic(c, "bit-blast") &
             tactic(c, "smt");
  solver s = t.mk_solver();
  map<bitvec, int> mp;

  expr inputs = c.bv_const("input", ch.n);
  vector<expr> partial_sums = gen_partial_sums(ch, c, s, mp, inputs);

  expr *path_evals[1 << ch.h];
  cout << ch.h << " Hadamards\n";
  if (ch.h <= 62) {
    for (unsigned long i = 0; i < (1 << ch.h); i++) {
      for (int j = 0; j < ch.h; j++) {
        if (i == (1 << j)) cout << j;
      }
      string str = "e" + std::to_string(i);
      path_evals[i] = new expr(c.bv_const(str.c_str(), 3));
      s.add(*(path_evals[i]) == gen_sums(ch, c, partial_sums, bitvec(ch.h, i), mp));
    } 
  } else {
    cout << "ERROR: Cannot handle more than 62 Hadamards";
    exit(1);
  }
  cout << "\nDone phase 1\n";

  func_decl f = z3::function("count_fun", c.bv_sort(3), c.int_sort());
  if (ch.h % 2 == 1) {
    // check #w^1 - #w^3 - #w^5 + #w^7 = 2^{h+1}
    s.add(f(c.bv_val(0, 3)) == 0);
    s.add(f(c.bv_val(1, 3)) == 1);
    s.add(f(c.bv_val(2, 3)) == 0);
    s.add(f(c.bv_val(3, 3)) == -1);
    s.add(f(c.bv_val(4, 3)) == 0);
    s.add(f(c.bv_val(5, 3)) == -1);
    s.add(f(c.bv_val(6, 3)) == 0);
    s.add(f(c.bv_val(7, 3)) == 1);
    s.add(count(ch, c, inputs, path_evals, f)
          != c.int_val(1 << ((ch.h+1) / 2)));
  } else {
    // check #w^0 - #w^4 = 2^h
    s.add(f(c.bv_val(0, 3)) == 1);
    s.add(f(c.bv_val(1, 3)) == 0);
    s.add(f(c.bv_val(2, 3)) == 0);
    s.add(f(c.bv_val(3, 3)) == 0);
    s.add(f(c.bv_val(4, 3)) == -1);
    s.add(f(c.bv_val(5, 3)) == 0);
    s.add(f(c.bv_val(6, 3)) == 0);
    s.add(f(c.bv_val(7, 3)) == 0);
    s.add(count(ch, c, inputs, path_evals, f)
          != c.int_val(1 << (ch.h / 2)));
  }

  check_result chk = s.check();
  if (chk == sat) {
     std::cout << "Circuits are not equivalent:\n";
     model m = s.get_model();
     std::cout << "  Counterexample is |" << m.eval(inputs) << ">\n";
  } else if (chk == unsat) {
     std::cout << "Circuits are equivalent\n";
  } else {
     std::cout << "Problem not solved\n";
  }
}

// ----------------------------- Quantified versions
expr xor_sum(character & ch, context & c, const xor_func & x, expr & inputs, expr & path) {
  int i;
  expr sum = c.bv_val(0, 1);
  for (i = 0; i < ch.n + ch.h; i++) {
    if (x.test(i)) {
      if (i < ch.n) sum = sum ^ to_expr(c, Z3_mk_extract(c, i, i, inputs));
      else          sum = sum ^ to_expr(c, Z3_mk_extract(c, i - ch.n, i - ch.n, path));
    }
  }
  
  return x.test(i)?sum ^ c.bv_val(1,1):sum;
}

expr gen_sums(character & ch, context & c, expr & inputs, expr & path) {
  vector<exponent>::iterator it;
  list<Hadamard>::iterator ti;
  expr ret = c.bv_val(0, 3);

  // Phase gates
  for (it = ch.phase_expts.begin(); it != ch.phase_expts.end(); it++) {
    ret = ret + ite(xor_sum(ch, c, it->second, inputs, path) == c.bv_val(1, 1), 
                    c.bv_val((unsigned int)it->first, 3), 
                    c.bv_val(0, 3));
  }
 
  // Hadamard gates
  for (ti = ch.hadamards.begin(); ti != ch.hadamards.end(); ti++) {
    ret = ret + ite((to_expr(c, Z3_mk_extract(c, ti->prep - ch.n, ti->prep - ch.n, path))
                      & xor_sum(ch, c, ti->wires[ti->qubit], inputs, path)) == c.bv_val(1, 1),
                    c.bv_val(4, 3), c.bv_val(0, 3));
  } 
  return ret;
}

expr counter(character & ch, context & c, expr & inputs,
                 expr & path, expr & count, expr & countp, func_decl & f) {
  expr tmp = c.bv_const("_sum", 3);
  return countp == count + f(gen_sums(ch, c, inputs, path));
}

expr assert_output_eq(character & ch, context & c, expr & inputs, expr & path) {
  int i;
  expr ret = c.bv_val(1, 1);
  for (i = 0; i < ch.n + ch.m; i++) {
    if (i < ch.n) 
      ret = ret & ~(to_expr(c, Z3_mk_extract(c, i, i, inputs)) 
                    ^ xor_sum(ch, c, ch.outputs[i], inputs, path));
    else 
      ret = ret & ~(xor_sum(ch, c, ch.outputs[i], inputs, path));
  }
  return ret;
}

expr trans(character & ch, context & c, expr & inputs,
                 expr hvals, expr count, expr countp, func_decl & f) { 
  return ite(assert_output_eq(ch, c, inputs, hvals) == c.bv_val(1, 1),
             counter(ch, c, inputs, hvals, count, countp, f),
             countp == count);
}

expr trans_sqr(character & ch, context & c, expr & inputs, 
                 int i, expr hvals, expr count, expr countp, func_decl & f) {
  if (i == 0) return trans(ch, c, inputs, hvals, count, countp, f);

  std::string names[7] = { "m-hv" + std::to_string(i),
                           "m-cv" + std::to_string(i),
                           "alpha" + std::to_string(i),
                           "x-hv" + std::to_string(i),
                           "x-cv" + std::to_string(i),
                           "y-hv" + std::to_string(i),
                           "y-cv" + std::to_string(i) };
  expr mcount = c.int_const(names[1].c_str());
  expr alpha  = c.bool_const(names[2].c_str());
  expr xhvals = c.bv_const(names[3].c_str(),ch.h);
  expr xcount = c.int_const(names[4].c_str());
  expr ycount = c.int_const(names[6].c_str());

  return exists(mcount, forall(alpha, exists(xhvals, xcount, ycount,
           trans_sqr(ch, c, inputs, i - 1, xhvals, xcount, ycount, f) &&
           implies(alpha, (xhvals == hvals) && (xcount == count) && (ycount == mcount)) &&
           implies(!alpha, (xhvals == hvals + (1 << (i - 1))) && (xcount == mcount) && (ycount == countp)))));
}

void quantified(character & ch) {
  config cfg;
  cfg.set(":mbqi", true);
  context c(cfg);
  tactic t = tactic(c, "simplify") &
             tactic(c, "solve-eqs") &
             tactic(c, "smt");
  solver s = t.mk_solver();

  expr inputs = c.bv_const("input", ch.n);
  expr path0 = c.bv_val(0, ch.h);
  expr count0 = c.int_val(0);
  expr countn = c.int_const("phase");

  func_decl f = z3::function("count_fun", c.bv_sort(3), c.int_sort());
  s.add(trans_sqr(ch, c, inputs, ch.h, path0, count0, countn, f));
  if (ch.h % 2 == 1) {
    // check #w^1 - #w^3 - #w^5 + #w^7 = 2^{h+1}
    s.add(f(c.bv_val(0, 3)) == 0);
    s.add(f(c.bv_val(1, 3)) == 1);
    s.add(f(c.bv_val(2, 3)) == 0);
    s.add(f(c.bv_val(3, 3)) == -1);
    s.add(f(c.bv_val(4, 3)) == 0);
    s.add(f(c.bv_val(5, 3)) == -1);
    s.add(f(c.bv_val(6, 3)) == 0);
    s.add(f(c.bv_val(7, 3)) == 1);
    s.add(countn != (1 << ((ch.h+1) / 2)));
  } else {
    // check #w^0 - #w^4 = 2^h
    s.add(f(c.bv_val(0, 3)) == 1);
    s.add(f(c.bv_val(1, 3)) == 0);
    s.add(f(c.bv_val(2, 3)) == 0);
    s.add(f(c.bv_val(3, 3)) == 0);
    s.add(f(c.bv_val(4, 3)) == -1);
    s.add(f(c.bv_val(5, 3)) == 0);
    s.add(f(c.bv_val(6, 3)) == 0);
    s.add(f(c.bv_val(7, 3)) == 0);
    s.add(countn != (1 << (ch.h / 2)));
  }

  check_result chk = s.check();
  if (chk == sat) {
     std::cout << "Circuits are not equivalent:\n";
     model m = s.get_model();
     std::cout << "  Model: " << m << "\n";
     std::cout << "  Counterexample is |" << m.eval(inputs) << ">\n";
  } else if (chk == unsat) {
     std::cout << "Circuits are equivalent\n";
  } else {
     std::cout << "Problem not solved\n";
  }
}

