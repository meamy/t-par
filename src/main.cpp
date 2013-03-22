#include "circuit.h"
#include <cstdio>
#include <iomanip>

int main(int argc, char *argv[]) {
  struct timespec start, end;
	dotqc circuit, synth;
  bool full_character = true;
  // Quick and dirty solution, don't judge me
  for (int i = 0; i < argc; i++) if ((string)argv[i] == "-no-hadamard") full_character = false;

	circuit.input(cin);

  if (full_character) {
    character c;
    c.parse_circuit(circuit);
    c.output(cerr);
    //	clock_gettime(CLOCK_MONOTONIC, &start);
    synth = c.synthesize();
    //	clock_gettime(CLOCK_MONOTONIC, &end);
    //	cout << fixed << setprecision(3);
    //	cout << "Time: " << (end.tv_sec + (double)end.tv_nsec/1000000000) - (start.tv_sec + (double)start.tv_nsec/1000000000) << " s\n";
  } else {
    metacircuit meta;
    meta.partition_dotqc(circuit);
    meta.optimize();
    synth = meta.to_dotqc();
  }

  synth.remove_swaps();
  synth.remove_ids();
  synth.print_stats();
  synth.print();

	return 0;
}
