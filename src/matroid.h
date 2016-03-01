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

#include <vector>
#include <deque>
#include "partition.h"

#include <assert.h>

#ifndef MATROID
#define MATROID

using namespace std;

struct path {
  list<pair <int, partitioning::iterator> > lst;

  path() { }
  path(int i, partitioning::iterator ref) { lst.push_front(make_pair(i, ref)); }
  path(const path & p) { lst = p.lst; }
  path(int i, partitioning::iterator ref, path & p) { 
    lst = p.lst;
    lst.push_front(make_pair(i, ref));
  }

  pair<int, partitioning::iterator> head() { return lst.front(); }
  int                               head_elem() { return lst.front().first; }
  partitioning::iterator            head_part() { return lst.front().second; }

  path_iterator begin() { return lst.begin(); }
  path_iterator   end() { return lst.end(); }

  void insert(int i, partitioning::iterator ref) { lst.push_front(make_pair(i, ref)); }
};

//-------------------------------------- Matroids

// Implements a matroid partitioning algorithm
template <class T, typename oracle_type>
void add_to_partition(partitioning & ret, int i, const vector<T> & elts, const oracle_type & oracle) {
  partitioning::iterator Si;
  set<int>::iterator yi, zi;

  // The node q contains a queue of paths and an iterator to each node's location.
  //	Each path's first element is the element we grow more paths from.
  //	If x->y is in the path, then we can replace x with y.
  deque<path> node_q;
  path t;
  path_iterator p;
  vector<bool> marked(elts.size());
  int tmp;
  bool flag = false;

  // Reset everything
  node_q.clear();
  for (int j = 0; j <= elts.size(); j++) {
    marked[j] = false;
  }

  // Insert element to be partitioned
  node_q.push_back(path(i, ret.end()));
  marked[i] = true;

  // BFS loop
  while (!node_q.empty() && !flag) {
    // The head of the path is what we're currently considering
    t = node_q.front();
    node_q.pop_front();

    for (Si = ret.begin(); Si != ret.end() && !flag; Si++) {
      if (Si != t.head_part()) {
        // Add the head to Si. If Si is independent, leave it, otherwise we'll have to remove it
        Si->insert(t.head_elem());

        if (oracle(elts, *Si)) {
          // We have the shortest path to a partition, so make the changes:
          //	For each x->y in the path, remove x from its partition and add y
          for (p = t.begin(); p != --(t.end()); ) {
            Si = p->second;
            (Si)->erase(p->first);
            (Si)->insert((++p)->first);
          }
          flag = true;
        } else {
          // For each element of Si, if removing it makes an independent set, add it to the queue
          for (yi = Si->begin(); yi != Si->end(); yi++) {
            if (!marked[*yi]) {
              // Generate an iterator to the position before yi
              zi = yi;
              if (zi != Si->begin()) zi--;
              // Take yi out
              tmp = *yi;
              Si->erase(yi);
              if (oracle(elts, *Si)) {
                // Put yi back in
                yi = Si->insert(Si->begin(), tmp);
                // Add yi to the queue
                node_q.push_back(path(*yi, Si, t));
                marked[*yi] = true;
              } else {
                yi = Si->insert(Si->begin(), tmp);
              }
            }
          }
          // Remove CURRENT from Si
          Si->erase(t.head_elem());
        }
      }
    }
  }

  // We were unsuccessful trying to edit the current partitions
  if (!flag) {
    set<int> newset;
    newset.insert(i);
    ret.push_front(newset);
  }

}

// Partition the matroid
template <class T, typename oracle_type>
partitioning partition_matroid(const vector<T> & elts, const oracle_type & oracle) {
  partitioning ret;

  // For each element of the matroid
  for (int i = 0; i < elts.size(); i++) {
    add_to_partition(ret, i, elts, oracle);
  }
  return ret;
}

// Repartition according to a new oracle
template <class T, typename oracle_type>
void repartition(partitioning & part, const vector<T> & elts, const oracle_type & oracle) {
  int tmp;
  partitioning::iterator Si;
  set<int>::iterator yi;

  list<int> acc;

  for (Si = part.begin(); Si != part.end(); Si++) {
    tmp = oracle.retrieve_lin_dep(elts, *Si);
    if (tmp != -1) {
      Si->erase(tmp);
      acc.push_back(tmp);
    }
    //assert(oracle(elts, *Si));
  }

  for (list<int>::iterator it = acc.begin(); it != acc.end(); it++) {
    add_to_partition(part, *it, elts, oracle);
  }
}

#endif
