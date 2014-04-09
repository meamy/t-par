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

#include <list>
#include <set>
#include <iostream>

using namespace std;

typedef list<set<int> > partitioning;
typedef list<pair <int, partitioning::iterator> >::iterator path_iterator;

ostream& operator<<(ostream& output, const partitioning& part);
partitioning freeze_partitions(partitioning & part, set<int> & st);

int num_elts(partitioning & part);
partitioning create(set<int> & st);
