#include <list>
#include <set>
#include <iostream>

using namespace std;

typedef list<set<int> > partitioning;
typedef list<pair <int, partitioning::iterator> >::iterator path_iterator;

ostream& operator<<(ostream& output, const partitioning& part);
partitioning freeze_partitions(partitioning & part, set<int> & st);
