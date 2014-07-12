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

#include <cstdio>
#include <iomanip>
#include <fstream>
#include "equiv_check.h"

int main(int argc, char *argv[]) {
  struct timespec start, end;
  dotqc circuit, tmp;
  ifstream in;
  character c;

  for (int i = 1; i < argc; i++) {
    in.open(argv[i], ifstream::in);
    if (i == 1) circuit.input(in);
    else {
      tmp.clear();
      tmp.input(in);
      for (gatelist::iterator it = tmp.circ.begin(); it != tmp.circ.end(); it++) {
        circuit.append(*it);
      }
    }
    in.close();
  }

  circuit.remove_swaps();
  circuit.remove_ids();
//  circuit.print();

  if (disp_log) cerr << "Parsing circuit...\n" << flush;
  c.parse_circuit(circuit);
  if (disp_log) cerr << "Checking circuit...\n" << flush;
  clock_gettime(CLOCK_MONOTONIC, &start);
  unquantified_sum(c);
//  unquantified(c);
  clock_gettime(CLOCK_MONOTONIC, &end);

  cout << fixed << setprecision(3);
  cout << "# Time: " << (end.tv_sec + (double)end.tv_nsec/1000000000) 
                      - (start.tv_sec + (double)start.tv_nsec/1000000000) << " s\n";

  return 0;
}
