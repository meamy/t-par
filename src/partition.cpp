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

#include "partition.h"

template<typename T>
bool is_disjoint(const set<T> & A, const set<T> & B) {
  typename set<T>::iterator itA = A.begin(), itB = B.begin();

  while (itA != A.end() && itB != B.end()) {
    if (*itA == *itB) {
      return false;
    }
    else if (*itA < *itB) itA++;
    else itB++;
  }

  return true;
}

ostream& operator<<(ostream& output, const partitioning& part) {
  partitioning::const_iterator Si;
  set<int>::const_iterator yi;

  for (Si = part.begin(); Si != part.end(); Si++) {
    output << "{";
    for (yi = Si->begin(); yi != Si->end(); yi++) {
      output << *yi << ",";
    }
    output << "}";
  }

  return output;
}

// Take a partition and a set of ints, and return all partitions that are not
//   disjoint with the set, also removing them from the partition
partitioning freeze_partitions(partitioning & part, set<int> & st) {
  partitioning ret;
  partitioning::iterator it, tmp;

  for (it = part.begin(); it != part.end();) {
    if (!is_disjoint(*it, st)) {
      tmp = it;
      it++;
      ret.splice(ret.begin(), part, tmp);
    } else {
      it++;
    }
  }
  return ret;
}

int num_elts(partitioning & part) {
  int tot = 0;
  partitioning::iterator it;
  for (it = part.begin(); it != part.end(); it++) {
    tot += it->size();
  }
  return tot;
}

partitioning create(set<int> & st) {
  return partitioning(1, st);
}
