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
