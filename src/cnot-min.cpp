#include "cnot-min.h"

int hamming(xor_func & a, xor_func & b) { return (a ^ b).count(); }

pair<int,int>
min_hamming(list<int> & candidates,
            const vector<exponent> & expnts, 
            const xor_func * bits,
            int n) {
  int min_index = 0;
  int min_bit = 0;
  int min_hamming = -1;
  int tmp;
  list<int>::iterator it;

  for (it = candidates.begin(); it != candidates.end(); it++) {
    for (int i = 0; i < n; i++) {
      tmp = hamming(expnts[*it], bits[i]);
      if (hamming == -1 || tmp < hamming) {
        min_index = *it;
        min_bit = i;
        min_hamming = tmp;
      }
    }
  }

  return make_pair(min_index, min_bit);
}

int min_count(list<int> & candidates, const vector<exponent> & expnts) {
  int min_index = 0;
  int min_count = -1;
  int tmp;
  list<int>::iterator it;

  for (it = candidates.begin(); it != candidates.end(); it++) {
    tmp = expnts[*it].count();
    if (min_count == -1 || expnts[*it].count < min_count) {
      min_index = *it;
      min_count = tmp;
    }
  }

  return min_index;
}
