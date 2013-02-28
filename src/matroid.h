#include <vector>
#include <deque>
#include "partition.h"

#ifndef MATROID
#define MATROID

using namespace std;

//-------------------------------------- Matroids
// Implements a matroid partitioning algorithm
template <class T, typename oracle_type>
class matroid {
    vector<T> elts;
    oracle_type oracle;

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

  public:
    matroid() { }
    matroid(vector<T> & S, oracle_type I) { elts = S; oracle = I; }

    void add_to_partition(partitioning & ret, int i) {
      partitioning::iterator Si;
      set<int>::iterator yi;

      // The node q contains a queue of paths and an iterator to each node's location.
      //	Each path's first element is the element we grow more paths from.
      //	If x->y is in the path, then we can replace x with y.
      deque<path> node_q;
      path t;
      path_iterator p;
      bool marked[elts.size()], flag;
      int tmp;
      set<int> * newset;

      // Reset everything
      node_q.clear();
      for (int j = 0; j <= i; j++) {
        marked[j] = false;
      }
      flag = false;

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
        newset = new set<int>;
        newset->insert(i);
        ret.push_front(*newset);
      }

    }

    // Partition the matroid
    partitioning partition_matroid() {
      partitioning ret;

      // For each element of the matroid
      for (int i = 0; i < elts.size(); i++) {
        this->add_to_partition(ret, i);
      }
      return ret;
    }
};

#endif
