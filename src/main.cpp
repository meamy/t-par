#include "circuit.h"
#include <cstdio>
#include <iomanip>

int main() {
  struct timespec start, end;
	dotqc circuit;
  character c;

	circuit.input(cin);
  c.parse_circuit(circuit);
  c.output(cerr);

 //	clock_gettime(CLOCK_MONOTONIC, &start);
  dotqc synth = c.synthesize();
//	clock_gettime(CLOCK_MONOTONIC, &end);
//	cout << fixed << setprecision(3);
//	cout << "Time: " << (end.tv_sec + (double)end.tv_nsec/1000000000) - (start.tv_sec + (double)start.tv_nsec/1000000000) << " s\n";
  synth.remove_swaps();
//	synth.print_stats();
  synth.print();

	return 0;
}
