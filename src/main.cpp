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

#include "circuit.h"
#include <cstdio>
#include <iomanip>

int main(int argc, char *argv[]) {
  struct timespec start, end;
	dotqc circuit, synth;
  bool full_character = true;
  bool post_process = true;
  // Quick and dirty solution, don't judge me
  for (int i = 0; i < argc; i++) if ((string)argv[i] == "-no-hadamard") full_character = false;
                                 else if ((string)argv[i] == "-no-post-process") post_process = false;

	circuit.input(cin);

  if (full_character) {
    character c;
    c.parse_circuit(circuit);
    c.output(cerr);
    //	clock_gettime(CLOCK_MONOTONIC, &start);
    synth = c.synthesize();
    //	clock_gettime(CLOCK_MONOTONIC, &end);
  } else {
    metacircuit meta;
    meta.partition_dotqc(circuit);
    //	clock_gettime(CLOCK_MONOTONIC, &start);
    meta.optimize();
    //	clock_gettime(CLOCK_MONOTONIC, &end);
    synth = meta.to_dotqc();
  }
  //	cout << fixed << setprecision(3);
  //	cout << "# Time: " << (end.tv_sec + (double)end.tv_nsec/1000000000) - (start.tv_sec + (double)start.tv_nsec/1000000000) << " s\n";

  if (post_process) {
    synth.remove_swaps();
    synth.remove_ids();
  }
  synth.print_stats();
  synth.print();

	return 0;
}
